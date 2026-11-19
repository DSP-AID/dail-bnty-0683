# seestar_eff — Seestar smart-telescope frames to AAVSO extended-format photometry

**Method version: v1** (every method constant and cut is pinned in code and listed in
`WORKFLOW.md`; the same inputs always produce byte-identical outputs).

A complete, tested workflow that takes **raw Seestar S50 FITS files** (Bayer-mosaic single
exposures, on-board cube stacks, CrowdSky MEF exports — three different dialects) all the way
to a **validated AAVSO Extended File Format report** for a real variable star, MW Cam
(AUID `000-BFQ-649`, δ-Scuti type, VSX range 9.15–9.26).

It reproduces the two light-curve magnitudes published in the upstream documentation to
−0.021 and −0.022 mag (tolerance ±0.03), on real public data, with no API keys and no network
access at verification time.

## Quickstart

```bash
python -m venv venv && . venv/bin/activate         # or any env with the two deps below
pip install -r requirements.txt                    # numpy + astropy (+ pytest for tests)

python scripts/run_example.py                      # full workflow -> report/
python scripts/verify_report.py                    # the acceptance gate (exit 0 == pass)
python -m pytest tests -q                          # 26 tests, ~2 min incl. integration
```

`verify_report.py` is the contract: it re-hashes every shipped data file against
`data/manifest.json`, re-runs the whole pipeline from scratch in a temp directory, and fails
unless the freshly written report is **byte-identical** to the committed one, the EFF
validator reports zero errors, and both published golden magnitudes match within tolerance.
It was verified to behave identically on Linux (Python 3.12, numpy 2.5.5, astropy 8.0.3) and
Windows (numpy 2.5.3, astropy 8.0.2).

## What you get

| Path | Contents |
|---|---|
| `report/aavso_report_mwcam.txt` | the worked example: 8 EFF rows (5 single 20 s subs, 2 stacks, 1 CrowdSky session) for MW Cam, FILTERCODE TG, JD UTC |
| `report/summary.json` | per-frame diagnostics (ZP, ensemble rms, check-star deviation, plate-solve stats) + golden cross-check results |
| `seestar_eff/` | the library: frames (3 layouts + timing), Bayer handling, plate solving, aperture photometry, sequence cuts, EFF writer/validator, pipeline |
| `data/` | the public example dataset (tarball + manifest, see below) and target metadata |
| `tests/`, `scripts/` | unit + integration tests; run and verification entry points |
| `WORKFLOW.md` | the method: every constant, cut, timing rule and quality gate with its source |
| `SOURCES.md` | every external reference with retrieval dates |

## The data (real, public, checksum-locked)

`data/example-data-v1.tar.gz` is the example dataset published by astronomyk for the
seestar-photometry project (release `example-data-v1`, 18 319 999 B, sha256
`0c6e38cd…c813`, matching the upstream `example_data/SHA256` file). It contains five raw 20 s
Bayer subs, three cube stacks (including one deliberately overexposed), one CrowdSky MEF
session with its STAR-TAB, solved WCS sidecars, and a Gaia catalogue cutout carrying
synthetic Johnson V/B/Rc (J-K-C) magnitudes. `data/manifest.json` pins the sha256 of the
tarball **and of every member**; the pipeline verifies them and extracts members on demand
(so shipping only the tarball is enough — a clean checkout rebuilds everything). Target
metadata comes from the AAVSO VSX API snapshot `data/vsx_mwcam.json`.

License: GPL-3.0-or-later (inherited from the upstream repository that produced both the
dataset and the calibration recipe this workflow mirrors); see `LICENSE` and `SOURCES.md`.

## Headline results (all reproducible via `scripts/verify_report.py`)

* Stack photometry vs the published light curve: **dV = −0.021 / −0.022 mag** at matching JD
  (within ±0.00001 d), tolerance ±0.03 — the offsets are dominated by comparison-ensemble size
  and are quantified in `WORKFLOW.md` §6.
* Per-stack ensemble scatter **0.019–0.025 mag** over 26–29 stars — at the ~20–23 mmag
  precision floor the upstream project measured for this hardware on bright targets.
* Bayer-sub photometry lands on the same zero point as cube photometry (ΔZP ≤ 0.13 mag
  across five single exposures) — independent proof of the GRBG phase convention and of the
  pattern plate-solver, which must absorb **+16.4…+16.8°** of alt-az field rotation between
  sessions at this declination (+81.5°).
* The deliberately overexposed frame is caught by the check-star consistency gate
  (dV = +0.166 > 3× its ensemble rms) and excluded from the report, with the reason recorded
  in `summary.json`.

## Known limits (read before submitting real observations)

* Native cube-stack headers under-report their true mid-exposure time by ~1–2 min
  (inter-sub overhead is not recorded); see `WORKFLOW.md` §3.
* Only TG rows are produced (green channel reduced with V-band comps, per the EFF spec); B/Rc
  from this hardware carries 0.1–0.35 mag
  per-star calibration residuals upstream measured, and no I-band counterpart exists.
* `#OBSCODE=TEST` is a placeholder — replace with your AAVSO observer code before upload.
