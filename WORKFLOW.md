# WORKFLOW.md — the method, version v1

Every number below is either (a) a constant pinned in code, (b) quoted from a dated external
source (`SOURCES.md`), or (c) measured on the shipped dataset by this repo's own tests/scripts
— marked *[measured]*. Nothing is invented.

```
data/example-data-v1.tar.gz ──verify+extract──► frames (3 layouts, UTC timing)
gaia_mwcam.ecsv + vsx_mwcam.json ──sequence cuts──► comparison/check ensemble
frames ──G-plane extraction──► photometry (forced aperture, r=7 px)
subs  ──pattern plate-solve vs stack WCS──► same forced photometry
      ──ensemble ZP per frame──► V, MERR ──quality gates──► EFF rows + summary.json
```

## 1. Data and provenance

Dataset: astronomyk `seestar_photometry` release `example-data-v1` (GPL-3.0-or-later),
sha256 `0c6e38cd…c813` = upstream `example_data/SHA256`. Target: MW Cam, AUID 000-BFQ-649,
VSX position 186.68225 / +81.47397222, δ-Scuti (DSCT), V 9.15–9.26 *[measured from
vsx_mwcam.json, VSX API]*. Reference catalogue: bundled `gaia_mwcam.ecsv` — Gaia DR3 sources
with **synthetic Johnson J-K-C V/B/Rc** (`v_jkc_mag` etc., 830 rows). The upstream project's
calibration docs use exactly this reference ("Gaia synthetic V, not APASS … homogeneous
across the whole sky") because it needs no variable-star catalogue access.

## 2. The three FITS dialects (frames.py)

| layout | what | how detected | total exposure | mid-exposure time rule |
|---|---|---|---|---|
| `bayer` | raw 20 s Bayer mosaic, no WCS | `IMAGETYP=Light`, NAXIS=2, BAYERPAT | EXPTIME | (DATE-OBS + DATE-EXP)/2 |
| `cube` | on-board stack, primary HDU shape (3,1000,1000) = R,G,B planes | NAXIS=3 | TOTALEXP (STACKCNT×EXPTIME) | DATE-OBS + TOTALEXP/2 |
| `mef` | CrowdSky export, EXTNAME RED/GREEN/BLUE + FOOTPRINT + STAR-TAB | GREEN extension present | EXPTIME (total; NIMAGES subs) | OB-START + OBSTOTAL/2 |

Plane order R,G,B for cubes *[measured]*: cross-plane magnitude differences vs catalogue BP−RP
give slopes m0−m1 = −0.25, m1−m2 = −0.29, m0−m2 = −0.55 mag per BP−RP (additively consistent,
red plane first); independently the MEF sidecar EXTNAMEs label the same order.

**DATE-OBS is UTC here** *[measured]*: the device-filename clocks in `SPORIGIN` run exactly
+2 h 00 m ahead of DATE-OBS on two independent files (a +2 h local-time zone converted to UT
by the firmware), and reproducing the published light curve's JDs within 10 s requires treating
DATE-OBS as UT.

Caveat: for `cube` frames the header midpoint under-reports the true mid-exposure by ~1–2 min
because inter-sub download overhead is not recorded (the upstream pipeline's own epoch for the
third example frame sits 0.0015 d = 2 min from the header-derived value — see §10).

## 3. Bayer handling (frames.g_plane)

BAYERPAT=GRBG ⇒ the two green sub-lattices sit at **(row+col) even** *(upstream data-format
doc; verified [measured]*: sublattice medians 963/966 (the two G positions) vs 986 = R and
1103 = B, matching the doc's quoted 965/986/1103 to ~2 ADU). Missing G pixels ((row+col) odd)
are filled by averaging their four **orthogonal** neighbours — which all sit on the G lattice —
giving a smooth full-resolution green image so one aperture code serves subs and stacks.
(Diagonal interpolation is wrong here: diagonal neighbours of an odd-parity pixel land on R/B.)

Photometric proof that the phase convention is right: single-exposure zero points agree with
the cube-stack zero point within 0.13 mag (§8 table) — a wrong phase would measure R or B
stars and shift ZP by ~0.5–1 mag.

## 4. Plate solution (plate.py)

Stacks/CrowdSky carry solved WCS (sidecar `.wcs`; CrowdSky also TAN-SIP in its primary). Raw
subs carry none, and pointing + field rotation move stars **~90 px and ~+16.7°** between the
21:48 sub session and the 22:44 stack session *[measured]* — near dec +81.5 the parallactic
rotation rate is extreme: upstream measured 0.31° per 100 s (≈11°/h), and this repo's five
subs solve to +16.38…+16.81°, drifting ≈0.17°/min *[measured]*, consistent.

Method: angle+shift (similarity, no scale) vote — brightest 12 detections × brightest 40
catalogue stars, rotation grid over the **full circle** at 0.5° steps (a ±25° grid leaves the
subs unsolvable — the true rotation exceeds it *[measured]*), candidate shifts from brightest
pairs, then least-squares refinement on matched pairs with a keep-the-better guard.
Typical match count 23–25 of 25 detections *[measured]*; the synthetic test suite recovers
rotations across the full circle to <0.01° and <0.05 px. Catalogue positions come
through the stack's WCS first, so the sub inherits an absolute astrometric frame
(sidecar solution matches 60/60 brightest detections to catalogue stars within 4 px
*[measured]*).

## 5. Comparison / check star sequence (sequence.py)

Cuts mirror the upstream recipe (lightcurve + zeropoint docs):

* reference magnitude `v_jkc_mag` (Gaia synthetic Johnson V — FILTERCODE "TG" convention);
* window V ∈ [10, 14] — above it this hardware saturates on stacks, below it the noise floor
  fans the relation out (upstream zeropoint doc §"a magnitude window");
* drop Gaia-flagged variables (`phot_variable_flag == "VARIABLE"`);
* isolation: no catalogue neighbour within 20″; proximity: ≤ 20′ of target (shared airmass
  and shared field-rotation aperture loss — upstream lightcurve doc "the three cuts");
* ensemble capped at 30 (brightest first, deterministic by V then source_id) per the upstream
  admonition that 10–30 well-matched stars beat bigger ensembles;
* check star: isolated, non-variable, V closest to 10.8, excluded from the ensemble.

Sensitivity *[measured]* — the published golden value (V=9.2518, its example used 62 comps)
shifts with ensemble choice on this very frame:

| ensemble | n used | dV vs golden |
|---|---|---|
| V∈[10,14], cap 30 (this workflow) | 29 | −0.021 |
| V∈[10,14], no cap | 41 | −0.014 |
| V∈[10,15] | 58 | −0.011 |
| V∈[10,16] | 108 | −0.000 |

so the −0.021 offset is explained by ensemble size (bigger/fainter ensembles converge on the
published value), not by a photometric error. Reported anyway at the pinned recipe.

## 6. Photometry + per-frame calibration (phot.py, aperture.py)

* Green plane only; circular aperture **r = 7 px** (2.371″/px sidecar scale *[measured]*;
  ≈4.2 px FWHM), sky annulus 14–22 px sigma-clipped mean; iterated flux-weighted centroid.
  r=7 encloses ~0.98 of the stellar flux (COG on the target, *[measured]*) — upstream's tuned
  choice is the radius giving 0.95 enclosed flux because a slightly larger aperture damps the
  rotation-driven position-dependent aperture loss; both sit on the flat part of the COG.
* **Forced photometry** at catalogue positions ⇒ the ensemble is identical across frames.
* Ensemble zero point per frame: `ZP = clipped mean(V_ref − m_inst)` with 3σ clipping and up
  to two refits — the upstream lightcurve-doc recipe ("each comparison referenced to its own
  catalogue magnitude; ZP = mean of (catalogue mag − instrumental mag)").
* No colour term is fitted here. Upstream's per-frame fit finds k = +0.036 mag/(B−R) with
  rms 0.016 on this same frame; with the ensemble spanning little colour, our plain-mean ZP
  reaches rms 0.019 — the difference is at their quoted systematic floor (~0.03 mag,
  "chi2_red … Seestar frames sit on a ~0.03 mag systematic floor"). A colour-term refinement
  is documented in the zeropoint doc and can be added later; V values here are reported on the
  synthetic-V scale without it (TRANS=NO).
* `V_target = m_inst + ZP`; `MERR = sqrt((rms/√n)² + 0.020²)` — the 20 mmag floor is the
  upstream headline differential precision ("~20 mmag on a bright target"; their tuned scatter
  floor: 23 mmag). Their own caveat applies: the formal error is optimistic; the check-star
  column in `summary.json` shows realised per-frame systematics.

## 7. Quality gates (pipeline.py) — each rejects with a recorded reason

1. **Saturation**: any measured star (target, comps, check) peaking above 0.9 × 65520 ADU.
   65520 = 12-bit full scale <<4 is this firmware's clipping level; all used frames peak at
   41–53 k ADU with smooth (unclipped) profiles *[measured]*. Upstream excludes V<10 from
   calibration on these exposures for the same reason.
2. **Check-star consistency**: reject if |V_check_measured − V_check_sequence| >
   max(0.10, 3×ensemble rms). This quantifies the upstream diagnostic that "a variable, a
   blend, a saturated star or one drifting off chip stands out" against the check. It removes
   exactly the frame upstream labels `stack_saturated`: its target peak passes gate 1, but
   only 16/30 comps stay measurable, its ensemble rms is 0.054 (worst of the set; <0.06 is
   their photometric-grade line), and its check sits +0.166 mag off — 3× its rms.
   The 0.10 absolute floor protects honest single exposures from photon-noise false rejects.
3. **Thin ensemble**: fewer than half the sequence measurable ⇒ frame rejected.

## 8. Results and golden cross-check

Per-frame output (see `report/summary.json`):

| frame | layout | JD UTC | V | MERR | n | ZP | rms | check dV |
|---|---|---|---|---|---|---|---|---|
| raw_sub_1..5 | bayer | 2461163.40897…995 | 9.19–9.52 | .023–.046 | 23–28 | 23.56–23.76 | .06–.22 | −.14…+.17 |
| stack_c17_15min | cube | 2461163.45225 | 9.2304 | .0203 | 29 | 23.697 | .019 | +.051 |
| stack_c17_30min | cube | 2461163.48735 | 9.2973 | .0206 | 26 | 23.757 | .025 | +.047 |
| crowdsky (Jul session) | mef | 2461246.38151 | 9.3094 | .0209 | 28 | 23.007 | .032 | +.005 |
| stack_saturated | cube | — rejected — | | | 16 | | .054 | +.166 |

Single 20 s exposures are ~0.1 mag per-frame accurate at this sky brightness (sky σ=93 ADU;
formal MERR is optimistic by ~1.5–2×; the check-star column makes that visible) — publish
stacks, keep subs as the per-exposure demonstration.

Golden rows from the upstream lightcurve doc for exactly these frames:

| frame | published V @ JD | this workflow | ΔV | ΔJD |
|---|---|---|---|---|
| stack_c17_15min | 9.2518 @ 2461163.45220 | 9.2304 @ 2461163.45225 | −0.021 | +0.00005 |
| stack_c17_30min | 9.3189 @ 2461163.48730 | 9.2973 @ 2461163.48735 | −0.022 | +0.00005 |
| stack_saturated (context) | 9.3041 @ 2461218.45136 | 9.2960 @ 2461218.45285 | −0.008 | +0.00149 |

(ΔV tolerance ±0.03 fixed in advance; explained to ≤0.014 by §5's ensemble table. The third
epoch's ΔJD is the cube-midpoint caveat of §2, and that frame is rejected here by gate 2.)

## 9. AAVSO EFF output (aavso_eff.py)

Extended format v1.2 (spec archived at archive.aavso.org). One row per surviving frame:
STARID = VSX AUID; DATE = JD_UTC of mid-exposure (§2 rules); MAG/MERR as above;
**FILT=TG** — per the EFF spec, TG covers green-channel (tri-colour/DSLR) photometry reduced
"using V-band comp star magnitudes", which is exactly what we do (G plane, V-scale comps); TRANS=NO (no atmospheric/extinction correction applied — differential
ensemble absorbs it, per upstream practice); MTYPE=STD (calibrated to comparison-star
magnitudes); CNAME=ENSEMBLE with CMAG=na (ensemble photometry — the 30-star sequence is not
inline in EFF; it is reproducible from the bundled catalogue + §5); KNAME/KMAG = check star's
Gaia id and its **sequence** magnitude (our measured value goes to NOTES as `ck_dV`);
AMASS=na (no observatory latitude in Seestar headers); CHART = camera serial (TELESCOP);
NOTES carries layout, plate-solution stats, n comps, ZP, rms, aperture, exposure time.
The module also implements a validator enforcing the spec's required parameters, field count,
value formats and length limits — `verify_report.py` requires zero errors.

## 10. Reproducibility contract

`python scripts/verify_report.py` (exit 0) = manifest hashes OK + fresh pipeline run in a
temp dir byte-identical to the committed report + validator clean + both golden rows within
tolerance + every used frame under the saturation screen + ensemble ≥8 comps + check star
present. Determinism: no RNG, no wall-clock, no network; numpy/astropy versions recorded in
README. The whole thing was run identically on Linux and Windows with byte-identical reports.
