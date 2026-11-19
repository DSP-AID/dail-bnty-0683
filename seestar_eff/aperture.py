"""Circular aperture photometry with annulus sky (pure numpy, method v1).

Centroid: iterated flux-weighted first moments inside a box (oversampled PSF ~4 px FWHM
makes moments stable); aperture = pixels whose centre is within r_ap of the centroid;
sky = sigma-clipped median of annulus r_in..r_out scaled by aperture pixel count.
"""
from __future__ import annotations

import numpy as np


def sky_stats(img: np.ndarray, box=101) -> tuple[float, float]:
    """Sigma-clipped median / MAD-sigma of the frame background (central `box` region)."""
    ny, nx = img.shape
    cy, cx = ny // 2, nx // 2
    y0, y1 = max(cy - box // 2, 0), min(cy + box // 2, ny)
    x0, x1 = max(cx - box // 2, 0), min(cx + box // 2, nx)
    v = img[y0:y1, x0:x1].ravel().astype(np.float64)
    for _ in range(5):
        med = np.median(v)
        mad = np.median(np.abs(v - med)) * 1.4826 + 1e-9
        keep = np.abs(v - med) < 3 * mad
        if keep.sum() < 100:
            break
        v = v[keep]
    return float(med), float(mad)


def centroid(img: np.ndarray, x0: float, y0: float, half: int = 6, iters: int = 3):
    ny, nx = img.shape
    x, y = float(x0), float(y0)
    for _ in range(iters):
        xi, yi = int(round(x)), int(round(y))
        xs = np.arange(max(xi - half, 0), min(xi + half + 1, nx))
        ys = np.arange(max(yi - half, 0), min(yi + half + 1, ny))
        if len(xs) < 3 or len(ys) < 3:
            return None
        sub = img[np.ix_(ys, xs)]
        w = np.clip(sub - np.median(sub), 0, None)
        tot = w.sum()
        if tot <= 0:
            return None
        x = float((w * xs[None, :]).sum() / tot)
        y = float((w * ys[:, None]).sum() / tot)
    return x, y


def aperture_flux(img: np.ndarray, x: float, y: float, r_ap: float,
                  r_in: float, r_out: float) -> dict | None:
    """Sky-subtracted aperture sum. Returns None if flux <= 0 or geometry off-frame."""
    ny, nx = img.shape
    xi, yi = int(round(x)), int(round(y))
    pad = int(np.ceil(r_out)) + 1
    x0, x1 = xi - pad, xi + pad + 1
    y0, y1 = yi - pad, yi + pad + 1
    if x0 < 0 or y0 < 0 or x1 > nx or y1 > ny:
        return None
    sub = img[y0:y1, x0:x1]
    yy, xx = np.mgrid[y0:y1, x0:x1]
    rr = np.hypot(xx - x, yy - y)
    ap = rr <= r_ap
    an = (rr >= r_in) & (rr <= r_out)
    if ap.sum() == 0 or an.sum() < 20:
        return None
    skyv = sub[an]
    med = np.median(skyv)
    mad = np.median(np.abs(skyv - med)) * 1.4826 + 1e-9
    good = np.abs(skyv - med) < 3 * mad
    if good.sum() < 10:
        return None
    sky = float(np.mean(skyv[good]))
    flux = float(sub[ap].sum() - sky * ap.sum())
    peak = float(sub[yi - y0, xi - x0])
    return dict(flux=flux, sky=sky, npix=int(ap.sum()), peak=peak)


def curve_of_growth(img: np.ndarray, x: float, y: float, radii, r_in: float, r_out: float):
    out = []
    for r in radii:
        a = aperture_flux(img, x, y, r, r_in, r_out)
        out.append(None if a is None else a["flux"])
    return out


def detect_bright(img: np.ndarray, sky: float, sig: float, n: int = 40,
                  min_sep: int = 25, edge: int = 30, nsigma: float = 8.0):
    """Top-n local maxima above nsigma*sigma (crude but adequate for plate solving)."""
    ny, nx = img.shape
    thresh = sky + nsigma * sig
    # local max via 3x3 maximum filter
    p = np.pad(img, 1, mode="edge")
    mx = np.maximum.reduce([p[:-2, :-2], p[:-2, 1:-1], p[:-2, 2:],
                            p[1:-1, :-2], p[1:-1, 1:-1], p[1:-1, 2:],
                            p[2:, :-2], p[2:, 1:-1], p[2:, 2:]])
    peaks = np.argwhere((img >= mx) & (img > thresh))
    vals = img[peaks[:, 0], peaks[:, 1]]
    order = np.argsort(-vals)
    picked: list[tuple[int, int, float]] = []
    for i in order:
        py, px = int(peaks[i, 0]), int(peaks[i, 1])
        if px < edge or py < edge or px >= nx - edge or py >= ny - edge:
            continue
        if any((px - qx) ** 2 + (py - qy) ** 2 < min_sep ** 2 for qx, qy, _ in picked):
            continue
        picked.append((px, py, float(vals[i])))
        if len(picked) >= n:
            break
    return picked
