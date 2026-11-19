"""Aperture photometry on analytic Gaussians + ensemble zero-point statistics."""
import numpy as np

from seestar_eff import aperture as ap
from seestar_eff import phot


def gauss(shape=(60, 60), flux=50000.0, sigma=2.0, xc=30.4, yc=29.7, sky=100.0):
    y, x = np.mgrid[0:shape[0], 0:shape[1]]
    return sky + flux / (2 * np.pi * sigma ** 2) * np.exp(
        -((x - xc) ** 2 + (y - yc) ** 2) / (2 * sigma ** 2))


def test_centroid_and_flux_recovery():
    img = gauss()
    m = phot.measure_at(img, 30.0, 30.0)
    assert m is not None
    assert abs(m["x"] - 30.4) < 0.15 and abs(m["y"] - 29.7) < 0.15
    enclosed = 1.0 - np.exp(-(phot.R_AP ** 2) / (2 * 2.0 ** 2))   # analytic r<7 fraction
    assert abs(m["flux"] / (50000.0 * enclosed) - 1.0) < 0.015    # within 1.5 %


def test_detect_and_measure_synthetic_star():
    img = gauss(shape=(100, 100), flux=8000.0, xc=37.5, yc=62.2)
    sky, sig = ap.sky_stats(img)
    pk = ap.detect_bright(img, sky, sig, n=3, nsigma=6.0)
    assert pk and abs(pk[0][0] - 37.5) <= 1.5 and abs(pk[0][1] - 62.2) <= 1.5


def test_ensemble_zp_exact_and_clips_outlier():
    m = np.array([10.0, 10.4, 10.9, 11.3, 12.0, 9.8, 10.6])
    v = m + 2.30
    v[2] += 1.0                                   # one catastrophically bad comp
    fit = phot.ensemble_zp(m, v)
    assert abs(fit["zp"] - 2.30) < 0.02
    assert fit["n"] == 6 and fit["rms"] < 0.01


def test_reduce_frame_end_to_end_synthetic():
    """Injected target + 4 comps on a 100x100 grid, spacing >> annulus: full reduce path."""
    img = np.full((100, 100), 100.0)
    pos = [(50.5, 50.5), (30.5, 30.5), (30.5, 69.5), (69.5, 30.5), (69.5, 69.5)]
    fluxes = [50000.0, 20000.0, 30000.0, 25000.0, 40000.0]   # target first
    for (x, y), f in zip(pos, fluxes):
        img = img + gauss(shape=img.shape, flux=f, sigma=2.0, xc=x, yc=y, sky=0.0)
    enclosed = 1.0 - np.exp(-(phot.R_AP ** 2) / 8.0)
    vref = np.array([-2.5 * np.log10(f * enclosed) for f in fluxes[1:]])
    res = phot.reduce_frame(img, np.array(pos[1:]), vref, pos[0])
    assert res is not None and res["n_comp"] >= 4
    want = -2.5 * np.log10(fluxes[0] * enclosed)
    assert abs(res["target"]["v"] - want) < 0.01   # calibrated target within 10 mmag
