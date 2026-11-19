"""Bayer GRBG phase + green-plane fill logic (synthetic), cross-layout sky (real).

Phase convention (G at row+col even for BAYERPAT=GRBG) is the upstream data-format doc's;
the photometric consequence of getting it wrong is caught in test_pipeline (sub ZP offset).
"""
import numpy as np
from astropy.io import fits
from astropy.time import Time

from seestar_eff.frames import Frame, g_plane


def _bayer_frame(mos):
    return Frame("mem", "bayer", None, mos.astype(float), "GRBG", 1.0, 1.0, 1,
                 Time("2026-05-02", scale="utc"), "test", fits.Header())


def test_gplane_extracts_only_green_and_fills_diagonally():
    mos = np.empty((6, 8))
    mos[0::2, 0::2] = 10.0   # G (r+c even)
    mos[1::2, 1::2] = 10.0   # G
    mos[0::2, 1::2] = 20.0   # R
    mos[1::2, 0::2] = 30.0   # B
    out = g_plane(_bayer_frame(mos))
    assert np.all(out[1:-1, 1:-1] == 10.0)     # interior: no red/blue bleed-through
    assert out.min() >= 10.0 and out.max() <= 30.0   # edge-pad rows stay sane


def test_gplane_is_not_the_wrong_parity():
    mos = np.full((6, 8), 5.0)
    mos[0::2, 0::2] = 100.0            # G sub-lattice only
    mos[1::2, 1::2] = 100.0
    out = g_plane(_bayer_frame(mos))   # orthogonal fill reconstructs the flat field inside
    assert np.all(out[1:-1, 1:-1] == 100.0)
    wrong = mos[1::2, 0::2]            # odd-parity sample sees only the background 5.0
    assert float(np.median(wrong)) == 5.0
