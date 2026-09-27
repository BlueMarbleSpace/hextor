# FILLET code-comparison review: HEXTOR findings and roadmap

> **Status (2026-09-27, later the same day):** executed on branch `fillet-review`; the outcome is in `notes/fillet_refile_2026.md`. Decisions 1–4 of §5 were taken as: corrected numerics on for the FILLET files with the code defaults unchanged, `diffadj` off so D = 0.5 exactly, the current branch chain kept as the base, Benchmark 1 tuned to 288 K alone; the ice line stays at 263.15 K pending v1.2 (decision 5: filed without waiting, as AVALON did).

*2026-09-27. Source: Rory's audit "Where the FILLET models differ, and why"
(25 Sep 2026, `/models/avalon/docs/FilletCodeComparison.pdf`), §4.2, §5, §6,
§9. The audit read HEXTOR at `f12c8d9` (the tip of `master`); every finding
below was re-checked against `HEAD` of `ethane-table` (12 commits ahead of
`master`) and, where it is a number, recomputed independently. AVALON's
reply and v1.2 release are the pattern to follow (`docs/response_to_rory.txt`
in the AVALON repo).*

No development has been done. This note records what was verified, what
should change, and in which order, under the constraint that HEXTOR is
published and in use: **every published configuration must keep producing
bit-identical output unless a namelist switch says otherwise.**

---

## 1. Verdict

- Every HEXTOR finding in the audit checks out against the source. The three
  HIGH items (zenith argument, nearest-level CO2 lookup, CO2 coordinate) were
  fixed in April–July 2026 and are in the regenerated outputs under `fillet/`
  (commit `3574830`, release 4.2.0). **The FILLET archive still holds the June
  2024 files (`852fd03`)**, so the re-file is the main deliverable.
- What remains splits into three classes:
  - **A. Diagnostics, output, scripts, table loader** — no physics. Fix
    outright. Includes one item the audit rated LOW that is more serious than
    it looks: the committed Makefiles still cannot build (§2, item 8).
  - **B. Initialisation and convergence** — fix outright (ghost cells) or add
    tighter, switchable tests (convergence, seasonal clock).
  - **C. Numerics and physics** (non-conservative diffusion operator, the
    sea-ice jump at 263.15 K) — meaningful, and each moves every published
    result. Implement behind switches defaulting to the published behaviour,
    evaluate on a regression set, then decide the defaults with a version
    bump. In the tidally locked configuration the affected belts are the
    substellar and antistellar points, i.e. the THAI/SAMOSA calibration
    targets, so adopting them means re-calibrating.
- Separately, the **FILLET configuration** needs decisions that are not code:
  untuned Benchmarks 2/3, D handling, Experiment 2a grid, Benchmark 1
  re-tune, table choice, ice-line threshold (§4, Phase 4).

---

## 2. Findings, verified

Line numbers are for `HEAD` of `ethane-table`; the audit's numbers refer to
`f12c8d9`.

### 2.1 Defects in the source

| # | Finding (audit severity) | Status at HEAD | What was verified | Class |
|---|---|---|---|---|
| 1 | Albedo table queried at `mu*180/pi` (HIGH) | fixed `91240ab` | `driver.f:888` uses `acos(mu)*180/pi` | done |
| 2 | OLR from nearest CO2 level; albedo along one diagonal (HIGH) | fixed `f6acf6f` | full multilinear `bracket`/`axfrac` in `radiation.f90` | done |
| 3 | CO2 coordinate stored as pCO2/(X+pCO2) (HIGH) | fixed `45627dd` | The audit measured X = 44/28 (mass-weighted mixing ratio), not π/2 as the commit message says; the two differ by 0.04 %, so the applied repair is correct in effect. Only the commit message / notes are wrong. | doc |
| 4 | Diffusion evaluated at x-midpoints while areas use true edges; pole ghost cell of width 0.0038 (MED) | present, `driver.f:526, 559, 690–696` | Reproduced the stencil offline. For T = 288 − 40 x² and D = 0.55 the belt-area sum of D·t2prime is **+0.78 W m⁻²** (a finite-volume operator on the true edges gives 10⁻¹⁴); the polar belt's divergence is **1.70×** the flux-form value (the audit says 37 %; the ratio depends on the profile), interior belts are within 0.5–3 %. Cause: the three-point stencil is fed a spacing ratio of 8 at the pole (dx = 0.0038 ghost against 0.030 belt). | C |
| 5 | Weak convergence test; N–S asymmetric ice edges (MED) | present, `driver.f:1946` | Halt is a single global-OLR criterion, `fluxcnvg` = 0.1 W m⁻² by default; nothing tests per-belt drift or hemispheric symmetry. In the regenerated outputs the N–S ice-edge asymmetry exceeds 0.5° in exp4 cases 65 (6.3°) and 64 (2.4°), exp1 case 101 (2.7°), exp3 case 92 (1.3°), all near S ≈ 1.05–1.07 / 288 K, i.e. at the cap transition. | B |
| 6 | Pole ghost cells start at 273 K regardless of `tempinit` (MED) | present, `driver.f:146, 412, 1531` | `data temp/20*273./`; the `tempinit` loop covers k = 1..nbelts; the ghost cells are only copied from their neighbours at the end of the first step. On a 233 K start the first step diffuses from a 273 K ghost across a 0.0038-wide cell. | B |
| 7 | Sea-ice fraction reaches 0.63 at 263.15 K then jumps to 1 (MED) | present, `driver.f:902–908` | fice = 1 − exp((T − 273.15)/10) reaches 0.632 at `icetemp`, then 1. With the FILLET heat capacities the belt heat capacity drops 1.18×10⁸ → 1.0×10⁷ (12×) and the surface albedo steps by 0.11 at the very temperature that defines the ice line. | C |
| 8 | Makefiles build `m_ffhash.f90`, deleted in `f6acf6f` (LOW) | **still true in the committed tree** | `git show HEAD:model/Makefile` and `HEAD:model/radiation/Makefile` list `m_ffhash.o`. The working copies were fixed locally but both are flagged *assume-unchanged* (`git ls-files -v` → `h`), so `git status` hides the edit and a fresh clone cannot build. | A, urgent |
| 9 | Run scripts do not check for errors (LOW) | present, `fillet/*/fillet_*.sh` | csh, `./driver > /dev/null`, then `tail -n 1` of the previous `fillet_global.out`; a failed case silently repeats the previous line under a new case number. | A |
| 10 | Table duplicates a CO2 level; 380 albedo cells unset; −1 sentinels interpolated (LOW) | present in both legacy tables (Sun and 2600 K) | fco2 = 2.4989×10⁻⁴ appears at indices 7 and 8. `read_table_v1` maps albedo rows with `minloc`, which always picks index 7, so index 8's 380 cells are never written (uninitialised memory, not −1). `bracket()` returns lo = 8 for **249.9 ≤ ppm < 274.9**, so any lookup in that range reads them. The 1300 albedo and 65 OLR sentinels (fco2 ≥ 0.41, T ≤ 230 K) pass the init check because OLR is negated (+1 mW m⁻²) and albedo is only tested for NaN. No FILLET case touches either (the Exp 4 grid brackets the gap at 222.3 and 281.2 ppm). | A |
| 11 | Ice-line array overflows past five crossings (LOW) | present, `driver.f:107, 1757–1762` | `iceline(0:5)`, `nedge` unbounded; nedge = 3 or 4 (cap plus belt) falls into the single-crossing branch and is misreported. Scanned all 974 regenerated per-case latitude files: none has more than two crossings, so latent for FILLET. | A |

### 2.2 Departures from the protocol

| # | Finding | Status | Verified | Class |
|---|---|---|---|---|
| 12 | Filed D is `d0`, used D is rescaled | present, `driver.f:1723, 1891` | `d = d0·pg0·(avemol0/avemol)²·(hcp/hcp0)·(rot0/rot)²` = 1.101·d0 for 1 bar N2 at 280 ppm; every FILLET namelist has `diffadj = .true.`, so Ben2/3 and all experiments ran at 0.5505 and Ben1 at 0.418 while filing 0.50/0.38. In Exp 4 the true D also falls to 0.485 at 10 % CO2. | A + config |
| 13 | Ben1 tuning carried into untuned runs | present | `cloudir = −6.3` in Ben2/3 and every experiment; Ben1 itself now has `+3.0` (and `d0 = 0.38`, Fresnel ocean, Earth geography). Headers still say "consistent with Ben1". | config |
| 14 | Exp 2a on the Exp 1a grid | present, `fillet/exp2a/fillet_exp2a.sh:16–18` | 0.875–1.1 au (190 cases) instead of 0.8–0.975 au (150 cases); the same in the archived script. | A |
| 15 | Exp 4 clamped below the table floor | table property | Legacy Sun table: 100 ppm floor, 0.91 ceiling, 190–370 K. Regenerated warm branch glaciates at ≤ 139 ppm and is flat (225.14 K) below 100 ppm. | doc/config |
| 16 | ATOA is an unweighted time mean including polar night | present, `driver.f:1739` | `znalbsum += atoa(k)` every step. During polar night `mu = 0` → zenith 90°, where the table gives 0.73 for a 0.6 surface; the regenerated Ben2 polar ATOA is exactly 0.73 against Asurf 0.60. The global mean (`albsum`) *is* insolation-weighted. | A |
| 17 | Exp 3/4 on the Ben2 base, 365-d year, fixed-start branches, undeclared | present | Exp 3/4 namelists: Ben2 parameters, `a` fixed, every case from 233 or 288 K. None stated in headers. | A (headers) |
| 18 | Print precision | present, `driver.f:1892` | `Inst` as `f4.2` cannot resolve the 0.0125 grid; `XCO2` as `f8.1` cannot resolve 1.26 ppm; exp4 lat headers use `{:.1F}`. | A |
| 19 | No experiment latitude files in the archive | files exist locally | 974 per-case files under `fillet/exp*/`; the archive README's layout does not ask for them (only `ben*/case_0/lat_output.dat`). | keep in repo |
| 20 | Provenance | — | Archived = `852fd03`; regenerated = `3574830` (4.2.0), not filed. HEAD differs from `3574830` by `4299f50` (`pg0` is now total dry pressure: Exp 4 at 10 % CO2 sees D lower by 10 % through `pg0`), plus moistdiff/CH4/RH tables, all off by default. | Phase 0 |

### 2.3 Found here, not in the audit

| # | Finding | Verified | Class |
|---|---|---|---|
| 21 | **Odd, non-integer step count per orbit.** With `a`, `msun` and `dt = 4.32×10⁴ s` the orbit is 730.62 half-day steps; the year boundary (`t ≥ 2π/w`, `driver.f:1596`) fires after 731 steps, so each model year is 365.5 d, overshoots the equinox by 0.19 d, and then restarts at equinox. The hemispheres are sampled at different phases and unequal counts — the mechanism AVALON traced its 0.4 K asymmetries to (365 vs 366 steps). Candidate cause of item 5. | B |
| 22 | Ben1 uses `solarcon = 1360` and `a = 1.49597892e13`; the others 1361 and 1.495978707e13. Inst = 0.9993, printed as 1.00. | trivial, align at re-tune |
| 23 | `fillet/ben1/lat_header.txt` first line is truncated ("…Haywor>"). | A |
| 24 | Item 4 in `do_longitudinal` mode: belts 1 and 18 are the antistellar and substellar points, whose temperatures are the THAI calibration targets. A conservative operator changes them most, so (d0, cloudir) for THAI/SAMOSA and (fcloud, cloudir) for Earth must be re-derived if it becomes the default. | C |
| 25 | `icetemp` (263.15 K, hard-coded at `driver.f:283`) is both the ice-fraction jump and the ice-line diagnostic threshold. Filing a 273.15/271.15 K ice line (audit rec. 1) needs a separate diagnostic threshold, or it silently turns the ice ramp into a 0 °C step. | A |

Not re-verified: the audit's Exp 1/2 state-change counts and its 37 % polar figure (profile-dependent; the mechanism and the +0.7 W m⁻² are confirmed).

---

## 3. Principles for the changes

1. **Bit-reproducibility of published configurations.** Physics or numerics changes go behind namelist switches whose defaults reproduce the published behaviour, as `diffadj_rot` and `moistdiff` already do. Establish a regression set and run it before and after every change with the switches off, requiring bit-identical output: FILLET Ben1–3; `namelists/input.nml.earth.ch4`; `namelists/input.nml.thai.hab1.calibrated`; SAMOSA Cases 4 and 11.
2. **Diagnostics are added, not redefined.** `zonal.out` albedo is read by ten plotting scripts; add an insolation-weighted column rather than change the existing one.
3. **Versioning.** Fixes plus new switches are a minor release (4.3.0). Changing a default (conservative operator, continuous ice ramp) is 5.0, shipped together with re-calibrated THAI/SAMOSA/Earth constants.
4. **Order.** Decide the FILLET configuration → regenerate → compare against `852fd03` and `3574830` → tag + Zenodo → PR → reply to Rory. AVALON's `docs/before_after_fillet.py` and `response_to_rory.txt` are the templates.

---

## 4. Roadmap

### Phase 0 — Baseline and provenance (½ day)

- **0.1** Clear the assume-unchanged flags and commit Makefiles that build (`m_ffhash.o` removed; `FC`/`WDIR` left as the documented placeholders). Verify a fresh clone builds. Fixes item 8.
- **0.2** Re-run Ben1–3 at HEAD with the current `fillet/` namelists and diff against the `3574830` files; record the drift from `4299f50` (expected ≲ 0.01 K at 280 ppm, up to 10 % in D at the Exp 4 top).
- **0.3** Branch plan: the FILLET work should sit on `master` after `ch4-lookup-table`/`ethane-table` are merged (their features are off by default and the release will carry them). `tools/run_samosa.py` has an uncommitted change on `ethane-table` to resolve first.

### Phase 1 — Output, scripts, loader, initialisation (no physics; 1–2 days)

- **1.1** FILLET global file: write the effective `d`, `Inst` with four decimals, `XCO2` with four; header states the D rule. Items 12, 18.
- **1.2** Per-belt insolation-weighted planetary albedo (1 − ASR/S̄; `znasrsum` exists, add `znssum`) as a new column in `zonal.out`/`fillet.out`; FILLET `ATOA` uses it and the header says so. Item 16.
- **1.3** Convergence record per case: orbits run, final |ΔOLR|, max N–S asymmetry max_k |T_k − T_{N+1−k}|; written to the run log and to the lat-file header. Audit rec. 7.
- **1.4** Ice-line diagnostic: array bounded by `nbelts`; nedge > 2 classified (cap coexisting with belt → report the cap, log it); new namelist `icelinetemp` (diagnostic only, default 263.15) separate from `icetemp`. Items 11, 25.
- **1.5** Run scripts rewritten as one resumable runner (bash with `set -euo pipefail` or Python in `tools/`, after `run_samosa.py`): exit-code check on `driver`, non-finite guard, corrected Exp 2a grid (0.8–0.975 au, 15 × 10), and direct output in the archive layout (`ben*/case_0/lat_output.dat` + `global_output.dat`; `exp1…exp2a/global_output.dat`; `exp3_cold`, `exp3_warm`, `exp4_cold`, `exp4_warm`; tarball). Headers state per experiment: `cloudir`, D, table and its ranges, branch construction, 365-d year, Ben2 base, ice-line definition, heat capacities; fix the truncated Ben1 line; `solarcon = 1361` and one `a`. Items 9, 13, 14, 17, 22, 23.
- **1.6** Legacy-table loader: collapse the duplicated CO2 level on load (bit-identical outside 249.9–274.9 ppm, where it replaces uninitialised memory) and halt with a message on any lookup that would interpolate a sentinel cell; note it in `notes/` and correct the π/2 attribution. Items 3, 10.
- **1.7** Initialise `temp(0)` and `temp(nbelts+1)` from the initial profile in every `resfile` branch. Changes only the first step; regression set confirms equilibria are unchanged. Item 6.

### Phase 2 — Convergence and the seasonal clock (2–3 days; switchable)

- **2.1** Namelist `nstepyr` (even integer; 0 = current behaviour): `dt = P/nstepyr` computed from the orbit, year boundary by step count, so Exp 1a/2a stay exact when `a` changes. FILLET: 730. Item 21.
- **2.2** Convergence: `fluxcnvg = 1e-3` for FILLET; optional per-belt criterion (max_k |ΔT_k| over an orbit); symmetry diagnostic printed at halt. Item 5.
- **2.3** Test on the asymmetric cases (exp4 63–66, exp1 101/111, exp3 91–93): tighter convergence alone, then with even steps. Decide the FILLET settings from what removes the asymmetry.

### Phase 3 — Numerics and physics under switches (1–2 weeks with evaluation)

- **3.1** Flux-conservative diffusion, `diffcons` (default `.false.`): fluxes (1 − x²)∂T/∂x at the true belt edges sin(lat ± 5°), divergence over the true belt width (the same widths as `area(k)`), zero flux at the poles; same path for `moistdiff`. Validate Σ area·D·t2prime = 0 to round-off; compare Ben1–3, THAI Hab1 (substellar/antistellar belts), SAMOSA Cases 4/11/16, Earth. Expected: ≈ −0.3 K globally (0.78 W m⁻² / 2.35 W m⁻² K⁻¹) plus ice feedback, largest at the poles. Items 4, 24.
- **3.2** Continuous sea-ice ramp, `iceramp` (default = published): fice = (1 − exp((T − 273.15)/10)) / (1 − e⁻¹), reaching 1 at 263.15 K without a jump; heat capacity follows. Raises fice by up to 58 % inside the ramp, so calibrations move; evaluate on Exp 3/4 widths, THAI ice line, Earth. Item 7.
- **3.3** Decision after 3.1/3.2: file FILLET with the switches on (the archive is being replaced anyway and both are numerical defects), keep the code defaults published until a 5.0 with re-calibrated constants. Not for the re-file: a v2-table FILLET configuration (12 zenith nodes, CO2 to 1 ppm) — the v2 tables are clear-sky with saturating OLR, an untuned Earth on them runs to 367 K, and the cloud parameterisation they need is exactly what the untuned experiments forbid.

### Phase 4 — FILLET configuration decisions (Jacob / Rory)

- **4.1** Untuned Ben2/3 and all experiments: `cloudir = 0`, `fcloud = 0`, and `diffadj = .false.` so D = 0.5 exactly (recommended; alternatively keep the rescaling and file `d`). HEXTOR then becomes the clear-sky member of the ensemble (Ben2 ≈ 275 K against 272 K now and 287 K archived); say so in the headers.
- **4.2** Ben1 re-tune of (d0, cloudir) to 288.0 K on a clean configuration (1361 W m⁻², one `a`, D reported). If v1.2 adds an albedo target (rec. 5), add `cloudalb`/`fcloud` with a two-knob fit as in `tools/calibrate_earth.py`.
- **4.3** Table: stay on the legacy 1-bar Sun table for continuity; declare its ranges (100 ppm floor → Exp 4 clamped below it; 190–370 K, extrapolated outside) per rec. 10.
- **4.4** Ice line: keep 263.15 K until v1.2 defines it, with `icelinetemp` ready to file 273.15 K (one temperature per belt, as AVALON argued). Optionally record the secondary definition in the header.
- **4.5** Exp 3/4 base (Ben2) and fixed-start branches: keep and declare; align with whatever Rory tells AVALON on the Ben1-vs-Ben2 question.

### Phase 5 — Regenerate, compare, release, re-file (2–3 days)

- **5.1** Full regeneration: Ben1–3 plus 934 experiment cases (190 + 190 + 190 + 150 + 114 + 100), as a low-priority Slurm array like the table builds.
- **5.2** Before/after figure and drift table (archived `852fd03` → regenerated 4.2.0 → new): state census, Exp 3/4 thresholds, N–S asymmetry, branch-order check (warm start never colder than cold start).
- **5.3** Release 4.3.0: tag, Zenodo (concept DOI 10.5281/zenodo.20074008), `CITATION.cff`, README, `notes/fillet_refile_2026.md`, CLAUDE.md.
- **5.4** PR replacing `Results/hextor` in `projectcuisines/fillet` (layout as in its `Results/README.md`), and the reply to Rory in the form of AVALON's.

---

## 5. Decisions needed before Phase 3–5

1. Adopt the conservative operator and continuous ramp for the re-file after evaluation, or file with published numerics and only the Phase 1–2 fixes.
2. D for FILLET: `diffadj` off (protocol value exactly) or rescaled and reported.
3. Merge order: land `ethane-table` on `master` first, then branch.
4. Ben1 targets: 288 K alone, or 288 K plus albedo 0.29.
5. Whether to wait for Protocol v1.2 (ice-line definition, Exp 3/4 base) before regenerating. AVALON did not wait; it filed and asked.
