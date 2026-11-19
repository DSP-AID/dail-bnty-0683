#!/usr/bin/env python3
"""Verification gate: re-derive the report from the shipped data and fail (exit 1) unless

  * every shipped file matches manifest.json sha256,
  * the EFF report re-written from scratch is byte-identical to the committed one,
  * the committed report passes the EFF v1.2 validator (0 errors),
  * both golden cross-checks against the upstream published light curve pass
    (|dV| <= 0.03 mag, |dJD| <= 0.0005 d),
  * the ensemble has >= 8 comparison stars and a check star, and every reported frame
    passes the saturation screen (peak < 90 % of full scale).

Reproducibility: run this on any clean machine with only numpy+astropy installed; it needs
no network access and no API keys. Exit code 0 == all checks passed.

Usage: python scripts/verify_report.py [data_dir] [report_dir]
"""
import json
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # repo root (package dir)
from seestar_eff import pipeline  # noqa: E402

here = Path(__file__).resolve().parents[1]
data = Path(sys.argv[1]) if len(sys.argv) > 1 else here / "data"
report = Path(sys.argv[2]) if len(sys.argv) > 2 else here / "report"

fails = []


def check(ok, msg):
    print(("PASS " if ok else "FAIL ") + msg)
    if not ok:
        fails.append(msg)


# 1. manifest hashes (tarball + every member + aux files)
man = json.load(open(data / "manifest.json"))
check(pipeline.sha256(data / man["tarball"]) == man["tarball_sha256"], "tarball sha256")
pipeline.ensure_data(data)  # ships as tarball only; re-extract members if absent
for name, want in man["members_sha256"].items():
    check(pipeline.sha256(data / name) == want, f"member sha256 {name}")
for name, info in man.get("aux_files", {}).items():
    check(pipeline.sha256(data / name) == info["sha256"], f"aux sha256 {name}")

# 2. fresh re-run into a temp dir -> byte-identical report
with tempfile.TemporaryDirectory() as td:
    fresh = pipeline.run(data, Path(td))
    committed = (report / "aavso_report_mwcam.txt").read_text()
    check((Path(td) / "aavso_report_mwcam.txt").read_text() == committed,
          "report byte-identical on re-run")
    fresh_sum = json.load(open(Path(td) / "summary.json"))

# 3. validator clean
errs = pipeline.aavso_eff.validate_eff(committed)
check(errs == [], f"EFF validator: {errs if errs else 'no errors'}")
n_rows = sum(1 for ln in committed.splitlines() if ln and not ln.startswith("#"))
check(n_rows == fresh_sum["rows"], f"row count {n_rows} == summary rows {fresh_sum['rows']}")

# 4. golden values from the upstream published light curve
for g in fresh_sum["golden_checks"]:
    check(g["pass_"], f"golden {g['file']}: dV={g['dV']} (ref V={g['ref_v']} JD={g['ref_jd']})")

# 5. sanity gates
check(fresh_sum["n_comp"] >= 8, f"comparison ensemble n={fresh_sum['n_comp']} >= 8")
check(fresh_sum["check_star"] is not None, "check star present")
for f in fresh_sum["frames"]:
    if f.get("used"):
        check(f["peak"] < pipeline.SAT_FRAC * pipeline.FULL_SCALE,
              f"{f['file']}: peak {f['peak']:.0f} ADU below saturation screen")
rejected = [f for f in fresh_sum["frames"] if not f.get("used")]
print(f"info: {len(rejected)} frame(s) rejected: "
      + ", ".join(f"{f['file']} ({f.get('reason', '')})" for f in rejected) or "none")

print(f"\n{len(fails)} failure(s)")
sys.exit(1 if fails else 0)
