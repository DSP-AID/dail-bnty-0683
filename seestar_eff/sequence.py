"""Comparison / check star sequence construction (method v1).

Cuts follow the versioned recipe in the seestar-photometry docs ("A differential light
curve" + "Calibrating a frame"):
  * reference = Gaia DR3 synthetic Johnson V (v_jkc_mag), homogeneous all-sky;
  * magnitude window V in [10, 14] (below: saturation regime for this unit; above: noise);
  * catalogue-flagged variables removed (phot_variable_flag == 'VARIABLE');
  * isolation: no other catalogue source within 20 arcsec;
  * proximity to target <= 20 arcmin (shared airmass + shared field-rotation aperture loss);
  * ensemble size capped at 30 (10-30 is the sweet spot; more adds noise faster than it averages).
Check star: isolated, non-variable, V closest to 10.8, excluded from the comparison ensemble.
"""
from __future__ import annotations

import numpy as np


def load_reference(table) -> dict:
    """gaia_mwcam.ecsv astropy Table -> dict of arrays for stars with usable synthetic V."""
    ok = np.isfinite(np.asarray(table["v_jkc_mag"], dtype=float))
    return {
        "source_id": np.asarray(table["source_id"])[ok],
        "ra": np.asarray(table["ra"], dtype=float)[ok],
        "dec": np.asarray(table["dec"], dtype=float)[ok],
        "g": np.asarray(table["phot_g_mean_mag"], dtype=float)[ok],
        "bp_rp": np.asarray(table["bp_rp"], dtype=float)[ok],
        "v": np.asarray(table["v_jkc_mag"], dtype=float)[ok],
        "varflag": np.asarray(table["phot_variable_flag"].astype(str))[ok],
    }


def _separation_deg(ra1, dec1, ra2, dec2):
    """Vectorised angular separation (deg), small-angle-safe at poles via 3-D dot."""
    def xyz(ra, dec):
        ra = np.radians(np.asarray(ra, dtype=float))
        dec = np.radians(np.asarray(dec, dtype=float))
        return np.stack([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra),
                         np.sin(dec)], axis=-1)
    a, b = xyz(ra1, dec1), xyz(ra2, dec2)
    d = np.clip((a * b).sum(axis=-1), -1, 1)
    return np.degrees(np.arccos(d))


def build_sequence(ref: dict, target_ra: float, target_dec: float,
                   v_lo=10.0, v_hi=14.0, iso_arcsec=20.0, max_sep_arcmin=20.0,
                   n_max=30) -> dict:
    """Return {'comps': idx array, 'check': idx or None} into ref arrays."""
    v, varflag = ref["v"], ref["varflag"]
    ra, dec = ref["ra"], ref["dec"]

    nonvar = varflag != "VARIABLE"
    inwin = (v >= v_lo) & (v <= v_hi)
    sep_t = _separation_deg(ra, dec, target_ra, target_dec) * 60.0  # arcmin
    near = sep_t <= max_sep_arcmin

    # isolation: nearest catalogue neighbour (any magnitude in the reference table)
    iso = np.ones(len(ra), dtype=bool)
    for i in np.where(inwin & near)[0]:
        d = _separation_deg(ra, dec, ra[i], dec[i]) * 3600.0
        if np.sort(d)[1] < iso_arcsec:  # [0] is self
            iso[i] = False

    cand = np.where(nonvar & inwin & near & iso)[0]
    if len(cand) == 0:
        return {"comps": cand, "check": None}
    # brightest first (closest to the bright target), deterministic by (V, source_id)
    order = cand[np.lexsort((ref["source_id"][cand], v[cand]))]
    comps = order[:n_max].tolist()

    # check star: closest to V=10.8 among isolated non-variable (may sit outside comps window)
    chk_cand = np.where(nonvar & iso & near & (v >= v_lo - 1.0) & (v <= v_hi + 1.0))[0]
    chk_cand = [i for i in chk_cand if i not in comps]
    check = None
    if len(chk_cand):
        check = int(chk_cand[np.argmin(np.abs(v[chk_cand] - 10.8))])
        comps = [i for i in comps if i != check]
    return {"comps": np.array(comps, dtype=int), "check": check}
