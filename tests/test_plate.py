"""Pattern plate-solver: synthetic full-circle rotations (the near-pole regime)."""
import numpy as np

from seestar_eff import plate

C = np.array([500.0, 500.0])


def _scene(theta_deg, shift, seed=7, n_ref=40, keep=30):
    rng = np.random.default_rng(seed)
    ref = rng.uniform(50, 950, size=(n_ref, 2))
    th = np.radians(theta_deg)
    R = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
    ideal = (ref - C) @ R.T + C + shift
    idx = np.argsort(-ideal[:, 0])[:keep]          # "brightest" subset, shuffled order
    det = ideal[idx] + rng.normal(0, 0.2, size=(keep, 2))
    return ref, det


def test_recovers_rotation_beyond_old_25deg_grid():
    """+150 deg: the parallactic rate near dec +81 can sweep this far in an hour."""
    for theta in (-150.0, 16.75, 90.0):
        ref, det = _scene(theta, np.array([-80.0, 55.0]))
        sol = plate.solve_by_pattern(det, ref, C)
        assert sol is not None, theta
        th, sh, n = sol
        assert abs(np.angle(np.exp(1j * (th - np.radians(theta))))) < np.radians(0.3), theta
        assert np.max(np.abs(sh - [-80.0, 55.0])) < 0.5, theta
        assert n >= 28, theta


def test_project_pattern_roundtrip():
    ref, det = _scene(16.75, np.array([-19.0, -72.0]))
    th, sh, n = plate.solve_by_pattern(det, ref, C)
    c, s = np.cos(th), np.sin(th)
    R = np.array([[c, -s], [s, c]])
    got = (ref - C) @ R.T + C + sh
    d2 = ((det[:, None, :] - got[None, :, :]) ** 2).sum(axis=2)
    assert float(np.median(np.sqrt(d2.min(axis=1)))) < 0.5
