#!/usr/bin/env python3
"""Run the full example workflow and print a per-frame summary table.

Usage: python scripts/run_example.py [data_dir] [out_dir]
Defaults: src-relative data/ and report/ directories.
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # repo root (package dir)
from seestar_eff import pipeline  # noqa: E402

here = Path(__file__).resolve().parents[1]
data = Path(sys.argv[1]) if len(sys.argv) > 1 else here / "data"
out = Path(sys.argv[2]) if len(sys.argv) > 2 else here / "report"

s = pipeline.run(data, out)
print(f"target: {s['target']['name']} ({s['target']['auid']})  comps={s['n_comp']}  "
      f"check={s['check_star']['kname'] if s['check_star'] else 'none'}")
print(f"{'file':24s} {'layout':6s} {'JD_UTC':>13s} {'V':>7s} {'merr':>6s} {'n':>3s} "
      f"{'zp':>7s} {'rms':>6s} {'peak':>7s} {'ck_dV':>7s} used")
for f in sorted(s["frames"], key=lambda e: e.get("jd", 9e9)):
    print(f"{f['file']:24s} {f['layout']:6s} {f.get('jd', 0):13.5f} {f.get('v', 0):7.4f} "
          f"{f.get('merr', 0):6.4f} {f.get('n_comp', 0):3d} {f.get('zp', 0):7.3f} "
          f"{f.get('fit_rms', 0):6.4f} {f.get('peak', 0):7.0f} "
          f"{(f['ck_dv'] if f.get('ck_dv') is not None else 0):7.3f} {f['used']}"
          + ("" if f["used"] else f"  <- {f.get('reason', '')}"))
print(f"\nrows={s['rows']} spec_errors={s['spec_errors']}")
for g in s["golden_checks"]:
    print("golden", g["file"], f"dV={g['dV']}", "PASS" if g["pass_"] else "FAIL")
print(json.dumps({"written": str(out / "aavso_report_mwcam.txt"),
                  "spec_errors": s["spec_errors"],
                  "golden_pass": all(g["pass_"] for g in s["golden_checks"])}))
