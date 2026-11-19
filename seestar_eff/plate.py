"""Plate solution: catalogue stars -> pixel positions, per frame.

Stacks/CrowdSky carry a usable solved WCS (sidecar .wcs / MEF primary incl. SIP); the raw
Bayer subs carry none, so they are registered by an angle+shift pattern vote against the
stack sidecar solution (alt-az field rotation between frames is real: ~0.31 deg per 100 s
on this dataset, several deg over a 15-min bin).

Method v1: brute-force vote over rotation grid x candidate translations from the brightest
detected/reference pairs, then least-squares similarity refinement on matched pairs.
"""
from __future__ import annotations

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS


def wcs_from_sidecar(sidecar_path) -> WCS:
    return WCS(fits.getheader(str(sidecar_path)))


def project(wcs: WCS, ra, dec):
    x, y = wcs.all_world2pix(np.asarray(ra), np.asarray(dec), 0)
    return np.column_stack([np.atleast_1d(x), np.atleast_1d(y)])


def _similarity_lsq(det_pts, ref_pts):
    """LSQ angle+shift (no scale) from matched pairs; returns (theta, dx) in the caller's
    convention det ~ R(theta) @ ref + dx. atan2 gives the det->ref angle, hence the flip."""
    d = det_pts - det_pts.mean(axis=0)
    r = ref_pts - ref_pts.mean(axis=0)
    num = float((d[:, 0] * r[:, 1] - d[:, 1] * r[:, 0]).sum())
    den = float((d[:, 0] * r[:, 0] + d[:, 1] * r[:, 1]).sum())
    theta = -np.arctan2(num, den)
    c, s = np.cos(theta), np.sin(theta)
    R = np.array([[c, -s], [s, c]])
    dx = det_pts.mean(axis=0) - R @ ref_pts.mean(axis=0)
    return theta, dx


def solve_by_pattern(det_xy: np.ndarray, ref_xy: np.ndarray, center: np.ndarray,
                     angles_deg=np.arange(-179.75, 180.0, 0.5), tol: float = 4.0):
    """Return (theta_rad, shift, n_matched) aligning ref_xy (sky->pixel via sidecar WCS)
    onto det_xy (brightest detections first). The angle grid spans the full circle because
    the parallactic rotation rate near dec +81 is extreme (~70-90 deg/h for this station -
    enough to exceed 60 deg across the dataset's ~1 h session gap; measured: a +-25 deg
    grid leaves the Bayer subs unsolvable).
    """
    det = np.asarray(det_xy, dtype=float)
    ref = np.asarray(ref_xy, dtype=float)
    top_det = det[: min(len(det), 12)]
    top_ref = ref[: min(len(ref), 12)]
    best = (0, None, None)
    for ang in angles_deg:
        th = np.radians(ang)
        c, s = np.cos(th), np.sin(th)
        R = np.array([[c, -s], [s, c]])
        rot = (ref - center) @ R.T + center
        # candidate shifts from brightest pairs, evaluated in one broadcast
        cand = (top_det[:, None, :] - rot[None, :, :]).reshape(-1, 2)
        d2 = ((det[:, None, None, :] - (rot[None, None, :, :] + cand[None, :, None, :]))
              ** 2).sum(axis=3)                      # (ndet, nshift, nref)
        counts = (d2.min(axis=2) <= tol * tol).sum(axis=0)   # matches per shift
        k = int(counts.argmax())
        if counts[k] > best[0]:
            best = (int(counts[k]), th, cand[k].copy())
        if best[0] >= max(6, min(len(det), len(ref))):
            break
    n, th, shift = best
    if th is None or n < 5:
        return None
    th0, shift0 = th, shift

    def _count(t, sh):
        c, s = np.cos(t), np.sin(t)
        R = np.array([[c, -s], [s, c]])
        rot = (ref - center) @ R.T + center + sh
        return int((((det[:, None, :] - rot[None, :, :]) ** 2).sum(axis=2)
                    .min(axis=1) <= tol * tol).sum())

    # similarity refinement on matched pairs (iterate)
    for _ in range(3):
        c, s = np.cos(th), np.sin(th)
        R = np.array([[c, -s], [s, c]])
        rot = (ref - center) @ R.T + center + shift
        d2 = ((det[:, None, :] - rot[None, :, :]) ** 2).sum(axis=2)
        j = d2.argmin(axis=1)
        ok = d2.min(axis=1) <= 9.0
        if ok.sum() < 4:
            break
        th2, sh2 = _similarity_lsq(det[ok], ref[j[ok]])
        # rebuild: rot' = (ref-center)R'^T + center + shift' s.t. det ~ rot'
        c2, s2 = np.cos(th2), np.sin(th2)
        R2 = np.array([[c2, -s2], [s2, c2]])
        shift2 = det[ok].mean(axis=0) - ((ref[j[ok]] - center) @ R2.T + center).mean(axis=0)
        th, shift = th2, shift2
    # refinement can diverge on a lucky low-n vote: keep whichever transform matches more
    # (on a tie prefer the refined one - it is the sub-pixel solution)
    if _count(th0, shift0) > _count(th, shift):
        th, shift = th0, shift0
    return th, shift, _count(th, shift)


def project_pattern(wcs: WCS, ra, dec, theta: float, shift, center: np.ndarray):
    """Apply solved pattern transform to catalogue positions -> pixel coords."""
    px = project(wcs, ra, dec)
    c, s = np.cos(theta), np.sin(theta)
    R = np.array([[c, -s], [s, c]])
    return (px - center) @ R.T + center + shift
