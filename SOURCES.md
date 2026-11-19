# SOURCES.md — external references, with retrieval dates

All numbers used in `WORKFLOW.md` trace to these sources or to *[measured]* results inside
this repo (tests, `scripts/`, `report/summary.json`). No API key is needed to verify anything
shipped here.

## Dataset and calibration recipe (the upstream project this workflow mirrors)

* Repository (GPL-3.0-or-later): https://github.com/astronomyk/seestar_photometry
  * Dataset asset: https://github.com/astronomyk/seestar_photometry/releases/download/example-data-v1/example-data-v1.tar.gz
    retrieved 2026-11-19, sha256 `0c6e38cd5f8c68e6179c5eda4ccecdd5d77df61f89897c67fb573b750cc9c813`
    (18 319 999 bytes) — matches the repository's `example_data/SHA256` file.
* Documentation (markdown originals; retrieved 2026-11-19):
  * data formats (dialect table, GRBG sub-lattice medians 965/986/1103):
    https://seestar-photometry.readthedocs.io/en/latest/_sources/data-format.md.txt
  * stacking (field rotation 0.31 deg per 100 s; green sky RMS 41.9 -> 17.6 ADU on 5 subs;
    similarity-transform registration):
    https://seestar-photometry.readthedocs.io/en/latest/_sources/usecases/stacking.md.txt
  * light curve (golden rows 9.2518 @ JD 2461163.45220, 9.3189 @ 2461163.48730, 9.3041 @
    2461218.45136; ensemble-ZP recipe; comparison cuts proximity/brightness/colour; 10-30
    comp sweet spot; ~20 mmag headline precision, 23 mmag tuned floor; 0.95 enclosed-flux
    aperture rationale; P = 0.1294 d for MW Cam):
    https://seestar-photometry.readthedocs.io/en/latest/_sources/usecases/lightcurve.md.txt
  * zero point (ZP 23.644 +- 0.003, k = +0.036 at colour0 1.137, rms 0.016 over 59 stars on
    stack_c17_15min; deep stack ZP 23.656; V in [10,14] window rationale; rms grading
    <0.06 photometric-grade; ~0.03 mag systematic floor; Gaia synthetic J-K-C reference
    choice): https://seestar-photometry.readthedocs.io/en/latest/_sources/usecases/zeropoint.md.txt
  * astrometry & Gaia (B-plane per-star calibration residual 104-347 mmag -> why only TG is
    reported; 9 arcmin constant S50 pointing offset; scale 2.37-2.39 arcsec/px):
    https://seestar-photometry.readthedocs.io/en/latest/_sources/astrometry-and-gaia.md.txt

## AAVSO

* Extended File Format v1.2 specification (latest update 2021-12-02) — field list, limits,
  DELIM rules, FILTERCODE table incl. TG ("Green Filter (or Tri-color green) ... commonly the
  'green-channel' in a DSLR or color CCD camera ... These observations use V-band comp star
  magnitudes."). Original https://www.aavso.org/aavso-extended-file-format returned HTTP 403
  to automated fetches at retrieval time (2026-11-19); text verified via Internet Archive
  snapshot https://web.archive.org/web/2024id_/https://www.aavso.org/aavso-extended-file-format.
* VSX target metadata for MW Cam (AUID 000-BFQ-649, position, V range 9.15-9.26, DSCT):
  https://apps.aavso.org/v2/api/stars/search/?name=MW%20CAM retrieved 2026-11-19
  (snapshot kept as `data/vsx_mwcam.json`; AAVSO data licence CC-BY-SA 4.0). Note: the legacy
  VSX web API returned 404 at retrieval time; this anonymous v2 endpoint worked.

## License text

* GNU GPL-3.0-only canonical text:
  https://raw.githubusercontent.com/spdx/license-list-data/main/text/GPL-3.0-only.txt
  retrieved 2026-11-19 (kept as `LICENSE`). The repository and the bundled dataset are
  distributed under GPL-3.0-or-later per the upstream repository's own licence.

## Environment notes (verification-time, not API keys)

* Verification ran on: Linux (Python 3.12, numpy 2.5.5, astropy 8.0.3, pytest 9.2.0) and
  Windows 10 (Python 3.13, numpy 2.5.3, astropy 8.0.2) — byte-identical reports on both.
* During development the VizieR/CDS and Gaia archive TAP services were intermittently down
  (2026 migration) and the AAVSO old VSX API returned 404; none are needed by this repo —
  all catalogue data used is bundled in the dataset tarball.
