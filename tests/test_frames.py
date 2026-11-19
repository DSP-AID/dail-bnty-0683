"""Frame layout detection + mid-exposure timing on the real shipped frames.

Golden anchors are the numbers published in the upstream docs (lightcurve use-case table)
and the exact header arithmetic - no invented values.
"""
import json
from pathlib import Path

from astropy.time import Time

from seestar_eff.frames import load_frame, g_plane

DATA = Path(__file__).resolve().parents[1] / "data"


def test_bayer_sub_timing_is_mid_exposure_utc():
    fr = load_frame(DATA / "raw_sub_1.fits.gz")
    assert fr.layout == "bayer" and fr.bayerpat == "GRBG"
    # DATE-OBS 21:48:45.018210 + (DATE-EXP - DATE-OBS)/2 = 21:48:55.322069 (UT)
    want = Time("2026-05-02 21:48:55.322070", scale="utc")
    assert fr.t_mid.scale == "utc"
    assert abs(fr.t_mid.jd - want.jd) < 1e-7          # < 9 ms
    assert fr.total_exptime == 20.0 and fr.n_exp == 1


def test_cube_stack_timing_and_golden_jd():
    fr = load_frame(DATA / "stack_c17_15min.fits.gz")
    assert fr.layout == "cube" and fr.data.shape == (3, 1000, 1000)
    assert fr.total_exptime == 760.0 and fr.n_exp == 38
    # upstream lightcurve doc row: JD 2461163.45220 for this frame (DATE-OBS + TOTALEXP/2)
    assert abs(fr.t_mid.jd - 2461163.45220) <= 0.0005


def test_mef_timing_uses_ob_start_span():
    from astropy.time import TimeDelta
    fr = load_frame(DATA / "crowdsky_mef.fits.gz")
    assert fr.layout == "mef"
    ob = Time(str(fr.header["OB-START"]), scale="utc")
    want = ob + TimeDelta(float(fr.header["OBSTOTAL"]) / 2.0, format="sec")
    assert abs(fr.t_mid.jd - want.jd) < 1e-9
    assert fr.time_source == "ob-start+obstotal/2"


def test_g_plane_shapes_and_sky_consistency():
    sub = load_frame(DATA / "raw_sub_1.fits.gz")
    stk = load_frame(DATA / "stack_c17_15min.fits.gz")
    g_s, g_k = g_plane(sub), g_plane(stk)
    assert g_s.shape == sub.bayer.shape and g_k.shape == (1000, 1000)
    # same sky band through two independent layouts. Tolerance 60 ADU: the sub was shot
    # 56 min before the stack session and the sky genuinely rose ~43 ADU in that hour
    # (measured medians: 966 vs 1009), so a tight bound would test the weather, not us.
    import numpy as np
    assert abs(np.median(g_s) - np.median(g_k)) < 60.0


def test_cube_plane_order_is_rgb_by_colour():
    """data[1] must be green: cross-plane magnitude differences vs Gaia BP-RP order the
    planes R > G > B in reddening (slopes -0.25 / -0.29 / -0.55 mag per BP-RP, measured)."""
    import numpy as np
    from astropy.table import Table
    from seestar_eff import plate, sequence, aperture as ap
    fr = load_frame(DATA / "stack_c17_15min.fits.gz")
    cat = sequence.load_reference(Table.read(DATA / "gaia_mwcam.ecsv", format="ascii.ecsv"))
    wcs = plate.wcs_from_sidecar(DATA / "stack_c17_15min.wcs")
    px = plate.project(wcs, cat["ra"], cat["dec"])
    sel = ((px[:, 0] > 40) & (px[:, 0] < 960) & (px[:, 1] > 40) & (px[:, 1] < 960)
           & np.isfinite(cat["bp_rp"]) & (np.asarray(cat["g"]) < 15))
    x, y, bprp = px[sel][:40, 0], px[sel][:40, 1], np.asarray(cat["bp_rp"])[sel][:40]
    mags = []
    for pi in range(3):
        row = []
        for xi, yi in zip(x, y):
            c = ap.centroid(fr.data[pi], xi, yi)
            f = ap.aperture_flux(fr.data[pi], c[0], c[1], 6.0, 12.0, 20.0)
            row.append(np.nan if f is None else -2.5 * np.log10(max(f["flux"], 1)))
        mags.append(np.array(row))
    def sl(a, b):
        d = mags[a] - mags[b]
        k = np.isfinite(d) & np.isfinite(bprp)
        return float(np.polyfit(bprp[k], d[k], 1)[0])
    assert sl(0, 1) < -0.1 and sl(1, 2) < -0.1 and sl(0, 2) < sl(0, 1) < 0
