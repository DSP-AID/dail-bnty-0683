"""Integration: full pipeline on the shipped dataset (the bounty's acceptance criteria).

Golden values are the two light-curve rows published in the upstream docs for exactly these
frames; tolerances |dV| <= 0.03 mag / |dJD| <= 0.0005 d were fixed before running (they are
also enforced by scripts/verify_report.py).
"""
import json

import pytest

from seestar_eff import pipeline


@pytest.fixture(scope="module")
def result(tmp_path_factory):
    from pathlib import Path
    return pipeline.run(Path(__file__).resolve().parents[1] / "data",
                        tmp_path_factory.mktemp("report"))


def test_report_valid_and_row_count(result):
    assert result["spec_errors"] == []
    assert result["rows"] == 8                       # 9 frames in, saturated one rejected


def test_golden_upstream_lightcurve_rows(result):
    for g in result["golden_checks"]:
        assert g["pass_"], (g["file"], g["dV"])


def test_saturated_frame_rejected_with_reason(result):
    bad = [f for f in result["frames"] if f["file"] == "stack_saturated.fits.gz"]
    assert bad and not bad[0]["used"]
    assert "check star" in bad[0]["reason"]


def test_target_magnitudes_in_physical_range(result):
    # MW Cam: VSX magmax 9.15 / magmin 9.26 -> rows must sit inside (amplitude + noise)
    for f in result["frames"]:
        if f.get("used"):
            assert 9.0 < f["v"] < 9.6, f["file"]


def test_bayer_and_stack_zeropoints_agree(result):
    """Green-lattice Bayer photometry and cube-plane photometry land on the same zero
    point (stacks average their subs): proves GRBG phase + pattern registration."""
    by_file = {f["file"]: f for f in result["frames"]}
    ref_zp = by_file["stack_c17_15min.fits.gz"]["zp"]
    for i in range(1, 6):
        zp = by_file[f"raw_sub_{i}.fits.gz"]["zp"]
        assert abs(zp - ref_zp) < 0.15, (i, zp, ref_zp)


def test_summary_frames_have_diagnostics(result):
    assert all("jd" in f or "reason" in f for f in result["frames"])
    assert result["check_star"]["kname"].startswith("G")
