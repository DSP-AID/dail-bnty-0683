"""End-to-end workflow: bundled Seestar dataset -> AAVSO EFF report + summary (method v1).

Pipeline (see WORKFLOW.md for rationale and citations):
  1. verify/extract the upstream example-data-v1.tar.gz against data/manifest.json;
  2. build the comparison/check ensemble from the bundled Gaia synthetic-V catalogue with
     sequence.build_sequence (upstream cuts);
  3. reduce every frame: stacks + CrowdSky MEF via their solved WCS, raw Bayer subs via the
     pattern plate-solver against the stack sidecar (subs carry no WCS, and alt-az pointing
     + field rotation displace stars by ~100 px / ~10 deg between sessions - measured);
  4. forced aperture photometry + ensemble zero point (phot.reduce_frame);
  5. screen frames on saturation (target peak > 90 % of the 65520 ADU full scale of this
     12-bit firmware export - clipping level measured on these files), write one EFF row per
     surviving frame, and validate the report against the v1.2 spec.

Output: aavso_report_mwcam.txt + summary.json (per-frame diagnostics + golden cross-checks
against the numbers published in the upstream docs).
"""
from __future__ import annotations

import hashlib
import json
import tarfile
from pathlib import Path

import numpy as np
from astropy.table import Table
from astropy.wcs import WCS

from . import aavso_eff, phot, plate, sequence
from .frames import load_frame, g_plane

FULL_SCALE = 65520.0          # 12-bit values << 4: this firmware's clipping level (measured)
SAT_FRAC = 0.90               # peaks above 90 % of full scale may be nonlinear -> reject
STARID = "000-BFQ-649"        # MW Cam AUID (data/vsx_mwcam.json, VSX API)
GOLDEN = [                    # upstream lightcurve-doc table rows (same frames, same method)
    {"jd": 2461163.45220, "v": 9.2518, "file": "stack_c17_15min.fits.gz"},
    {"jd": 2461163.48730, "v": 9.3189, "file": "stack_c17_30min.fits.gz"},
]
JD_TOL, V_TOL = 0.0005, 0.03

FRAME_ORDER = ["raw_sub_1.fits.gz", "raw_sub_2.fits.gz", "raw_sub_3.fits.gz",
               "raw_sub_4.fits.gz", "raw_sub_5.fits.gz", "stack_c17_15min.fits.gz",
               "stack_c17_30min.fits.gz", "stack_saturated.fits.gz",
               "crowdsky_mef.fits.gz"]


def sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def ensure_data(data_dir) -> Path:
    """Check the release tarball against manifest.json; extract missing members."""
    data = Path(data_dir)
    man = json.load(open(data / "manifest.json"))
    tar_path = data / man["tarball"]
    got = sha256(tar_path)
    if got != man["tarball_sha256"]:
        raise ValueError(f"{man['tarball']}: sha256 {got} != manifest {man['tarball_sha256']}")
    with tarfile.open(tar_path, "r:gz") as tf:
        for m in tf.getmembers():
            out = data / m.name
            if not out.exists():
                with tf.extractfile(m) as src, open(out, "wb") as dst:
                    dst.write(src.read())
    for name, want in man["members_sha256"].items():
        got = sha256(data / name)
        if got != want:
            raise ValueError(f"extracted {name}: sha256 mismatch vs manifest")
    return data


def _frame_wcs(fr, data_dir):
    """Solved WCS for a frame: sidecar .wcs beside it, else the frame header (CrowdSky MEF
    carries TAN-SIP in the primary). Returns (None, reason) for native stacks/subs."""
    side = Path(str(fr.path).replace(".fits.gz", ".wcs"))
    if side.exists():
        return plate.wcs_from_sidecar(side), "sidecar-wcs"
    try:
        w = WCS(fr.header)
        if w.has_celestial:
            return w, "header-wcs"
    except Exception:  # shortcut: astropy can raise on odd SIP; fall back to pattern solve
        pass
    return None, "no-wcs"


def _bright_cat_xy(wcs, cat, shape):
    """Catalogue stars projected into the frame box (all magnitudes), brightest first."""
    px = plate.project(wcs, cat["ra"], cat["dec"])
    ny, nx = shape
    inside = (px[:, 0] > 0) & (px[:, 0] < nx) & (px[:, 1] > 0) & (px[:, 1] < ny)
    order = np.argsort(np.asarray(cat["g"], dtype=float)[inside])[:40]
    return px[inside][order]


def run(data_dir, out_dir=None):
    data = ensure_data(data_dir)
    out = Path(out_dir) if out_dir else Path(data_dir).resolve().parent / "report"
    out.mkdir(parents=True, exist_ok=True)

    vsx = json.load(open(data / "vsx_mwcam.json"))
    cat = sequence.load_reference(Table.read(data / "gaia_mwcam.ecsv", format="ascii.ecsv"))
    seq = sequence.build_sequence(cat, vsx["ra"], vsx["dec"])
    comp_idx, chk_idx = np.asarray(seq["comps"]), seq["check"]
    if len(comp_idx) < 8:
        raise ValueError(f"only {len(comp_idx)} comparison stars pass the cuts")

    ref_wcs = plate.wcs_from_sidecar(data / "stack_c17_15min.wcs")
    frames_out, rows = [], []

    for name in FRAME_ORDER:
        fr = load_frame(data / name)
        img = g_plane(fr)
        ny, nx = img.shape
        center = np.array([nx / 2.0 - 0.5, ny / 2.0 - 0.5])

        wcs, solve_note = _frame_wcs(fr, data)
        if wcs is not None:
            project = lambda ra, dec: plate.project(wcs, ra, dec)
        else:
            sky, sig = phot.ap.sky_stats(img)
            det = np.array([[x, y] for x, y, _ in
                            phot.ap.detect_bright(img, sky, sig, n=25, nsigma=6.0)])
            ref_xy = _bright_cat_xy(ref_wcs, cat, (ny, nx))
            sol = plate.solve_by_pattern(det, ref_xy, center) if len(det) >= 6 else None
            if sol is None:
                frames_out.append(dict(file=name, layout=fr.layout, used=False,
                                       reason=f"plate solve failed ({len(det)} detections)"))
                continue
            theta, shift, nmatch = sol
            project = lambda ra, dec: plate.project_pattern(ref_wcs, ra, dec,
                                                            theta, shift, center)
            solve_note = f"solve n={nmatch} rot={np.degrees(theta):.2f}deg"

        comps_xy = project(cat["ra"][comp_idx], cat["dec"][comp_idx])
        tgt_xy = project([vsx["ra"]], [vsx["dec"]])[0]
        chk_xy = (project([cat["ra"][chk_idx]], [cat["dec"][chk_idx]])[0]
                  if chk_idx is not None else None)
        r = phot.reduce_frame(img, comps_xy, np.asarray(cat["v"])[comp_idx], tgt_xy, chk_xy)
        if r is None:
            frames_out.append(dict(file=name, layout=fr.layout, used=False,
                                   reason="photometry failed / <4 usable comps"))
            continue

        peak = r["target"]["peak"]
        # saturation screen over everything we photometer (a clipped comp/check biases ZP):
        # stack_saturated passes on the target alone but its brighter check star is clipped.
        all_peaks = [peak] + [c["peak"] for c in r["comps"] if c is not None]
        if "check" in r:
            all_peaks.append(r["check"]["peak"])
        peak_max = max(all_peaks)
        sat = peak_max > SAT_FRAC * FULL_SCALE
        ck_dv = (r["check"]["v"] - cat["v"][chk_idx]) if chk_idx is not None and "check" in r else None
        # check-star consistency gate (upstream lightcurve doc: a variable/blend/clipped star
        # "stands out" against the check): reject when the measured check deviates from its
        # sequence magnitude by more than 3x the frame's own ensemble scatter (floor 0.10 so
        # single-exposure photon noise never rejects good frames).
        ck_bad = (ck_dv is not None and "check" in r
                  and abs(ck_dv) > max(0.10, 3 * r["fit_rms"]))
        thin = r["n_comp"] < 0.5 * len(comp_idx)   # ensemble degraded (e.g. clipping regime)
        entry = dict(file=name, layout=fr.layout, jd=float(fr.t_mid.jd),
                     v=round(r["target"]["v"], 4), merr=round(r["target"]["merr"], 4),
                     n_comp=r["n_comp"], rejected=r["rejected"], zp=round(r["zp"], 4),
                     fit_rms=round(r["fit_rms"], 4), peak=round(peak_max, 1),
                     sky=round(r["sky"], 1), ck_dv=None if ck_dv is None else round(ck_dv, 4),
                     solve=solve_note, used=not (sat or ck_bad or thin))
        reasons = []
        if sat:
            reasons.append(f"ensemble peak {peak_max:.0f} > {SAT_FRAC * FULL_SCALE:.0f} ADU "
                           "(clipping regime)")
        if ck_bad:
            reasons.append(f"check star dV={ck_dv:+.3f} > 3x ensemble rms "
                           f"({r['fit_rms']:.3f})")
        if thin:
            reasons.append(f"only {r['n_comp']}/{len(comp_idx)} comparison stars measurable")
        if reasons:
            entry["reason"] = "; ".join(reasons)
        else:
            rows.append(dict(
                starid=STARID, jd=float(fr.t_mid.jd), mag=r["target"]["v"],
                merr=r["target"]["merr"], filt="TG", trans="NO", mtype="STD",
                cname="ENSEMBLE", cmag=None,
                kname=f"G{cat['source_id'][chk_idx]}" if chk_idx is not None else None,
                kmag=float(cat["v"][chk_idx]) if chk_idx is not None else None,
                amass="na", group=None, chart=fr.unit_id[:20],
                notes=(f"{fr.layout}; {solve_note}; comps={r['n_comp']}; zp={r['zp']:.3f}; "
                       f"rms={r['fit_rms']:.4f}; ap={phot.R_AP:.0f}px; exp={fr.total_exptime:.0f}s"
                       + (f"; ck_dV={ck_dv:+.3f}" if ck_dv is not None else ""))))
        frames_out.append(entry)

    rows.sort(key=lambda o: o["jd"])
    from . import __version__ as _ver
    eff = aavso_eff.write_eff(rows, obscode="TEST", software=f"seestar_eff {_ver} (TG ensemble)",
                              obstype="CCD", delim=",")
    spec_errors = aavso_eff.validate_eff(eff)
    (out / "aavso_report_mwcam.txt").write_text(eff)

    golden = []
    for g in GOLDEN:
        m = [e for e in frames_out if e["file"] == g["file"]]
        ok = bool(m) and m[0].get("used") and abs(m[0]["v"] - g["v"]) <= V_TOL \
            and abs(m[0]["jd"] - g["jd"]) <= JD_TOL
        golden.append(dict(file=g["file"], ref_v=g["v"], ref_jd=g["jd"],
                           got_v=m[0]["v"] if m else None,
                           got_jd=m[0]["jd"] if m else None,
                           dV=None if not m else round(m[0]["v"] - g["v"], 4),
                           pass_=ok))

    summary = dict(target=vsx, n_comp=int(len(comp_idx)),
                   check_star=None if chk_idx is None else
                   dict(kname=f"G{cat['source_id'][chk_idx]}",
                        ra=float(cat["ra"][chk_idx]), dec=float(cat["dec"][chk_idx]),
                        kmag=float(cat["v"][chk_idx])),
                   rows=len(rows), spec_errors=spec_errors, frames=frames_out,
                   golden_checks=golden)
    (out / "summary.json").write_text(json.dumps(summary, indent=1))
    return summary
