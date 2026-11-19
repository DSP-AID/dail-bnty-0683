"""Frame loading and mid-exposure timing for the three Seestar FITS layouts.

Layouts (method v1, mirrors the dialect table in seestar-photometry docs, "Seestar data
formats", https://seestar-photometry.readthedocs.io/en/latest/data-format.html):
  A "cube"  : native Seestar on-board stack, primary HDU shape (3, ny, nx) = R,G,B planes.
  B "mef"   : CrowdSky multi-extension; planes resolved by EXTNAME; FOOTPRINT excluded.
  C "bayer" : raw single-exposure Bayer mosaic, tagged BAYERPAT.

EXPTIME means different things per dialect: native EXPTIME = per-sub (TOTALEXP = total);
CrowdSky EXPTIME = total (NIMAGES = sub count). Mid-exposure time per frame:
  bayer : (DATE-OBS + DATE-EXP)/2, fallback DATE-OBS + EXPTIME/2
  cube  : DATE-OBS + TOTALEXP/2          (short by ~1-2 min of wall-clock mid; inter-sub
                                          overhead is unknowable from native headers)
  mef   : OB-START + OBSTOTAL/2          (true span when present)
DATE-OBS is taken as UT. Evidence for this dataset: the device-local file-name clocks in
SPORIGIN run exactly +2h00m ahead of DATE-OBS in two independent files, i.e. the header
time was converted to UTC by the ALP firmware (FITS standard also requires DATE-OBS=UT).
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.time import Time, TimeDelta

IMAGE_EXTENSIONS = ("FOOTPRINT", "MASK", "WEIGHT")


@dataclass
class Frame:
    path: str
    layout: str                 # cube | mef | bayer
    data: np.ndarray            # (3, ny, nx) float for cube/mef; None for bayer
    bayer: np.ndarray | None    # raw mosaic float for bayer; else None
    bayerpat: str | None
    exptime_sub: float          # seconds per sub-exposure
    total_exptime: float        # on-sky total, seconds
    n_exp: int
    t_mid: Time                 # mid-exposure, UTC scale
    time_source: str            # which rule produced t_mid
    header: fits.Header         # primary header (dialect preserved)

    @property
    def shape(self):
        return self.bayer.shape if self.layout == "bayer" else self.data.shape[1:]

    @property
    def unit_id(self) -> str:
        return str(self.header.get("TELESCOP", "?"))


def _hget(h: fits.Header, key, default=None):
    return h[key] if key in h else default


def _time(card) -> Time | None:
    if card is None:
        return None
    return Time(str(card).strip(), scale="utc")  # ISO-8601; FITS-standard UT


def load_frame(path) -> Frame:
    path = str(path)
    with fits.open(path) as hdul:
        ph = hdul[0].header
        names = [hdu.name for hdu in hdul]
        if "GREEN" in names:  # layout B: CrowdSky MEF, planes by EXTNAME
            planes = {}
            for ext in ("RED", "GREEN", "BLUE"):
                planes[ext] = np.asarray(hdul[ext].data, dtype=np.float64)
            data = np.stack([planes["RED"], planes["GREEN"], planes["BLUE"]])
            n_exp = int(_hget(ph, "NIMAGES", 1))
            exptime_sub = float(_hget(ph, "EXPTIME", 0.0)) / max(n_exp, 1)
            total = float(_hget(ph, "EXPTIME", 0.0))
            ob_start, obstotal = _time(_hget(ph, "OB-START")), _hget(ph, "OBSTOTAL")
            if ob_start is not None and obstotal:
                t_mid = ob_start + TimeDelta(float(obstotal) / 2.0, format="sec")
                src = "ob-start+obstotal/2"
            else:  # shortcut: no span cards; upgrade when OB-END appears in an export
                t_mid = _time(_hget(ph, "DATE-OBS")) + TimeDelta(total / 2.0, format="sec")
                src = "date-obs+exptime/2"
            return Frame(path, "mef", data, None, _hget(ph, "BAYERPAT"), exptime_sub,
                         total, n_exp, t_mid, src, ph.copy())

        if ph.get("NAXIS") == 3:  # layout A: native cube (3, ny, nx) channel-first
            arr = np.asarray(hdul[0].data, dtype=np.float64)
            if arr.shape[0] != 3 and arr.shape[-1] == 3:  # channel-last guard
                arr = np.moveaxis(arr, -1, 0)
            n_exp = int(_hget(ph, "STACKCNT", 1))
            exptime_sub = float(_hget(ph, "EXPTIME", _hget(ph, "EXPOSURE", 0.0)))
            total = float(_hget(ph, "TOTALEXP", exptime_sub * n_exp))
            t_mid = _time(_hget(ph, "DATE-OBS")) + TimeDelta(total / 2.0, format="sec")
            return Frame(path, "cube", arr, None, _hget(ph, "BAYERPAT"), exptime_sub,
                         total, n_exp, t_mid, "date-obs+totalexp/2", ph.copy())

        if str(_hget(ph, "IMAGETYP", "")).strip() == "Light" and ph.get("NAXIS") == 2:
            # layout C: raw Bayer sub
            mos = np.asarray(hdul[0].data, dtype=np.float64)
            exptime_sub = float(_hget(ph, "EXPTIME", _hget(ph, "EXPOSURE", 0.0)))
            t_obs, t_exp = _time(_hget(ph, "DATE-OBS")), _time(_hget(ph, "DATE-EXP"))
            if t_exp is not None:
                t_mid, src = t_obs + (t_exp - t_obs) * 0.5, "date-obs+date-exp/2"
            else:
                t_mid, src = t_obs + TimeDelta(exptime_sub / 2.0, format="sec"), "date-obs+exptime/2"
            return Frame(path, "bayer", None, mos, _hget(ph, "BAYERPAT"), exptime_sub,
                         exptime_sub, 1, t_mid, src, ph.copy())

    raise ValueError(f"{path}: unrecognised Seestar FITS layout")


def g_plane(fr: Frame) -> np.ndarray:
    """Green (TG) plane of a frame: cube/mef plane index 1; bayer G sub-lattice filled to
    full resolution.

    For BAYERPAT GRBG the two green sub-lattices sit at (row+col) even (verified for this
    dataset in tests/test_bayer.py against the channel medians documented upstream). The
    missing (row+col)-odd positions are filled by averaging their four orthogonal (N/E/S/W)
    neighbours - which all sit on the G lattice - giving a smooth full-resolution G image so
    the same aperture code serves stacks and raw subs. R/B planes need different
    interpolation (not needed here: the report uses TG only; see WORKFLOW.md).
    """
    if fr.layout == "bayer":
        mos = fr.bayer
        out = mos.copy()
        p = np.pad(mos, 1, mode="edge")
        orth = 0.25 * (p[:-2, 1:-1] + p[2:, 1:-1] + p[1:-1, :-2] + p[1:-1, 2:])
        miss = np.zeros(mos.shape, dtype=bool)
        miss[0::2, 1::2] = True
        miss[1::2, 0::2] = True
        out[miss] = orth[miss]
        return out
    return fr.data[1]
