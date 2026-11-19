"""Per-frame forced aperture photometry + ensemble calibration (method v1).

For each frame: measure target, comparison stars and check star on the green (TG) plane at
catalogue positions (forced photometry, so the ensemble is identical across frames), then
form the ensemble zero point ZP = clipped mean of (V_ref - m_inst) over valid comparisons
and report V_target = m_inst(target) + ZP. This is the upstream lightcurve recipe
(seestar-photometry docs, "Building a light curve", each comparison referenced to its own
catalogue magnitude; 3-sigma clipping x2 refits as in the zeropoint doc). A colour-term
refinement exists (see WORKFLOW.md) but is not applied here: the synthetic-Johnson V
reference colours span too little spread for it to help at this precision.

Reference magnitudes are the Gaia-derived synthetic Johnson V (v_jkc_mag, Jester-Andrews /
Cutri et al.) from the bundled catalogue - the FILTERCODE "TG" convention.

Reported per-frame error: MERR = sqrt( (rms/sqrt(n_comp))^2 + sigma_floor^2 ), sigma_floor
= 0.020 mag, the differential precision floor measured for this instrument class on a bright
target (~20 mmag, seestar-photometry docs front page). rms alone would understate systematic
error; photon noise is negligible for V<10 targets here.
"""
from __future__ import annotations

import numpy as np

from . import aperture as ap

R_AP, R_IN, R_OUT = 7.0, 14.0, 22.0   # px at 2.378"/px; ~1.7 FWHM, annulus outside PSF
SIGMA_FLOOR = 0.020


def measure_at(img: np.ndarray, x: float, y: float) -> dict | None:
    c = ap.centroid(img, x, y)
    if c is None:
        return None
    f = ap.aperture_flux(img, c[0], c[1], R_AP, R_IN, R_OUT)
    if f is None or f["flux"] <= 0:
        return None
    f["m_inst"] = -2.5 * np.log10(f["flux"])
    f["x"], f["y"] = c
    return f


def ensemble_zp(m_inst, v_ref) -> dict | None:
    """ZP = clipped mean(V_ref - m_inst); 3-sigma clip, up to 2 refits."""
    d = np.asarray(v_ref, dtype=np.float64) - np.asarray(m_inst, dtype=np.float64)
    d = d[np.isfinite(d)]
    if len(d) < 4:
        return None
    for _ in range(3):
        med = float(np.median(d))
        mad = 1.4826 * float(np.median(np.abs(d - med))) + 1e-9
        keep = np.abs(d - med) <= 3 * mad
        if keep.all() or keep.sum() < 4:
            break
        d = d[keep]
    return dict(zp=float(np.mean(d)), rms=float(np.std(d)), n=int(len(d)))


def reduce_frame(img: np.ndarray, comps_xy, comp_v, target_xy, check_xy=None):
    """Returns per-frame photometry dict or None if the frame is unusable."""
    sky, sig = ap.sky_stats(img)
    comps = [measure_at(img, x, y) for x, y in np.atleast_2d(comps_xy)]
    ok = [i for i, c in enumerate(comps) if c is not None]
    tgt = measure_at(img, target_xy[0], target_xy[1])
    chk = measure_at(img, check_xy[0], check_xy[1]) if check_xy is not None else None
    if tgt is None or len(ok) < 4:
        return None
    fit = ensemble_zp([comps[i]["m_inst"] for i in ok], np.asarray(comp_v)[ok])
    if fit is None:
        return None

    v_t = tgt["m_inst"] + fit["zp"]
    merr = float(np.hypot(fit["rms"] / np.sqrt(fit["n"]), SIGMA_FLOOR))
    out = dict(sky=sky, sig=sig, n_comp=fit["n"], zp=fit["zp"], fit_rms=fit["rms"],
               rejected=len(ok) - fit["n"],
               target=dict(m_inst=tgt["m_inst"], v=v_t, flux=tgt["flux"], peak=tgt["peak"],
                           x=tgt["x"], y=tgt["y"], merr=merr),
               comps=comps)
    if chk is not None:
        out["check"] = dict(m_inst=chk["m_inst"], flux=chk["flux"], peak=chk["peak"],
                            x=chk["x"], y=chk["y"], v=chk["m_inst"] + fit["zp"])
    return out
