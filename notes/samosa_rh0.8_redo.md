# Redoing SAMOSA on the RH 0.8 M-dwarf table

Started 2026-09-15.  The previous HEXTOR SAMOSA submission
(`samosa/samosa_summary.csv`, `global_output_HEXTOR.dat`, committed 2026-09-14)
was computed on `radiation_N2_CO2_3000K_p.h5`, whose columns are saturated
(RH = 1), with dry temperature diffusion and the constants
`(d0, cloudir) = (3.33, -45.0)` refitted at THAI Hab1 against the RH = 1
2600 K table.

The new table `radiation_N2_CO2_3000K_p_rh0.8.h5` (built 2026-09-15, see
`subsaturated_rce_tables.md`) is the same grid at RH = 0.8.  Two decisions
were taken before any of the runs below (Jacob, 2026-09-15):

1. **Recalibrate rather than transplant.**  Part of the -45 W/m^2 compensates
   for the RH = 1 table's overstated clear-sky greenhouse, so carrying it onto
   RH = 0.8 radiative transfer would double-count.  The THAI Hab1 calibration
   needs its own RH = 0.8 table, around a 2600 K star, which did not exist and
   was built for this (`radiation_N2_CO2_2600K_p_rh0.8.h5`, 70 columns, 1.3 h).
2. **Turn on moist static energy diffusion** (`&ebm::moistdiff`,
   `rhmoist = 0.8`), which is what the RH < 1 tables were built to pair with.

## What changed in the tooling

- `tools/moist_stability.py` (new).  Mirrors `driver.f`'s `moistprops`,
  `esatw` and `qmoist` so the harness can compute
  `beta = d(h/cp)/dT = 1 + (L/cp) dq/dT` and divide the explicit-diffusion
  timestep by it.  beta is 2.5 at 290 K, 3.5 at 300 K and 13 at 343 K for
  1 bar of N2, so this is not a small correction.
- `tools/run_samosa.py` and `tools/calibrate_thai.py` take `--moistdiff` /
  `--rhmoist`, and set `dt` per run by a **fixed point**: run at a guessed
  substellar temperature, then redo at the beta the run actually produced.
  A first pass that goes non-finite is retried once at 420 K, the tightest
  in-range beta, so an unstable timestep is not mistaken for a runaway -- in
  the output the two are indistinguishable.
- `namelists/input.nml.samosa` and `namelists/input.nml.thai.hab1` gained
  `moistdiff`, `rhmoist` and (for Hab1) an explicit `fluxcnvg`.
- The summary CSV records `table`, `rh_table`, `moistdiff`, `rhmoist`,
  `beta`, `t_hot`, `dt_passes` and `dt_capped`.

### Two defects found and fixed on the way

**The drift metric mislabelled converged runs.**  `classify()` called a case
`drifting` on the change in global mean over its last ten years.  These cases
approach equilibrium exponentially and the flux halt fires at ~20 years, so a
ten-year window still spans most of the transient.  Case 11 is settled to
0.001 K/yr by year 19 and yet scored -0.51 K/decade, which labelled it
`drifting` and would have dropped it from the submission
(`samosa_submit.py` reports only `equilibrium` cases).  Verified against runs
with the halt disabled: case 9 and case 10 sit at the flux-halted value to
0.001 K out to 2000 years, so the **halt is sound and only the label was
wrong**.  Replaced by `dT_rate10`, ten times the median of the last three
annual increments, which equals the ten-year change for anything drifting
linearly and goes to zero once the exponential has died.  Separation is now
0.02 K/decade for equilibria against 41+ for runaways, against 0.13 vs 0.6
before.

**`--years` never bounded a run.**  With `seasons = .true.` the halt at `tend`
is commented out (`driver.f:593-600`), so reaching `tend` only stops the
orbital clock; the integration continues to flux convergence or to the
hard-coded `niter = 5000`.  A case run with the halt disabled reached 5000
years with `--years 300`.  This changes nothing for SAMOSA -- these cases are
synchronous at zero obliquity and eccentricity, so a frozen orbital clock
leaves the insolation constant -- but it is a live trap for any seasonal
configuration, and `--years` is now documented for what it does.  Not fixed
in `driver.f`: that halt has been commented out since before this work and
re-enabling it would change published behaviour.

## Decomposition (published constants, so the two changes can be separated)

`--sequence all16 --transport constant --d0-ref 3.33 --cloudir -45`, warm
start, global mean surface temperature [K]:

| case | S | p [bar] | RH 1, dry | RH 0.8, dry | RH 1, moist | RH 0.8, moist |
|---|---|---|---|---|---|---|
| 1 | 500 | 0.70 | 174.01 | 174.01 | 174.01 | 174.01 |
| 4 | 1200 | 2.34 | 308.34 | 300.93 | 294.38 | 289.77 |
| 8 | 800 | 6.16 | 220.17 | 218.58 | 220.05 | 218.53 |
| 9 | 1100 | 0.70 | 278.17 | 273.85 | 268.97 | 266.95 |
| 10 | 400 | 4.83 | 153.52 | 153.52 | 153.52 | 153.52 |
| 11 | 900 | 0.10 | 232.00 | 231.11 | 224.94 | 224.53 |
| 14 | 900 | 1.44 | 241.47 | 240.09 | 240.87 | 239.77 |
| 15 | 600 | 0.43 | 189.40 | 189.31 | 189.40 | 189.31 |
| 16 | 1400 | 10.0 | 465.15 | 405.37 | 396.08 | 375.47 |

The other seven cases (2, 3, 5, 6, 7, 12, 13) run away in every configuration.
**Which cases have a climate solution is unchanged by either change** -- 9 of
16 in all four.  Both changes cool, and they are roughly additive.  The
largest single effect is on case 16, the 10 bar case, which falls from 465 K
(above the table's 420 K ceiling, so the old value was extrapolated) to 375 K,
inside the tabulated range.

Note the RH 1 / dry column reproduces the committed
`samosa/samosa_summary.csv` to 0.000 K in every case, so none of the harness
changes above affect a dry run.

## Calibration

`tools/calibrate_thai.py`, THAI Hab1 (TRAPPIST-1 e, 2600 K), targets
240.9 / 194.7 / 291.6 K (global / antistellar / substellar, 4-GCM ensemble).

| table | diffusion | d0 | cloudir | T_glob | T_anti | T_sub | contrast |
|---|---|---|---|---|---|---|---|
| RH 1 | dry | 3.33 | -45.00 | 240.41 | 203.53 | 299.62 | 96.08 |
| RH 0.8 | dry | 3.26 | -42.57 | 240.90 | 203.69 | 300.67 | 96.97 |
| RH 1 | moist 0.8 | 2.50 | -46.98 | 240.90 | 197.43 | 294.03 | 96.59 |
| **RH 0.8** | **moist 0.8** | **2.46** | **-45.09** | **240.90** | **197.31** | **294.19** | **96.88** |
| target | | | | 240.9 | 194.7 | 291.6 | 96.9 |

Two things to take from this.

**`d0` moves, `cloudir` barely does, and that is a coincidence of two
cancelling effects.**  Moist diffusion supplies part of the transport `d0` was
carrying, so `d0` falls 3.33 -> 2.46.  Going to RH 0.8 raises the clear-sky
OLR, so less of an OLR correction is needed: `cloudir` -45.0 -> -42.6 at fixed
dry diffusion.  Moist diffusion pushes the other way, -45.0 -> -47.0 at fixed
RH.  Together they nearly cancel and `cloudir` lands at -45.09, within 0.1 of
the published value.  **That near-equality is not a reason to have skipped the
recalibration** -- it is the outcome of it, and neither half of it is small.

**The calibrated configuration now reproduces the profile, not only the two
numbers it was fitted to.**  Against the THAI 4-GCM ensemble the antistellar
and substellar means are 2.6 K too warm, where the published dry/RH 1 fit was
8.8 and 8.0 K too warm.  The day-night contrast error falls from 0.8 K to
0.02 K.  Nothing in the procedure targets the two extremes separately -- it
fits the global mean with `cloudir` and the contrast with `d0` -- so this is
an independent check, and moist transport is what moved it.

The `dt` used over the moist sweep ran 820-870 s, from beta = 9.5.  Note that
is the **conservative first guess** (`T_HOT_GUESS = 330 K`) and not a
refinement: the fixed point only tightens `dt` when a run comes out hotter than
the guess, never loosens it when the run is colder, and these runs peak near
294 K.  Verified dt-independent: the calibrated point gives 240.904 / 240.908 K
at dt = 1350 / 844 / 400 / 200 / 100 s, a spread of 0.004 K.

## Production run

    python tools/run_samosa.py --sequence all16 --init both --transport constant \
        --d0-ref 2.46 --cloudir -45.09 --moistdiff --rhmoist 0.8 \
        --table ./radiation/radiation_N2_CO2_3000K_p_rh0.8.h5 --outdir samosa

32 runs, 0 failed.  Results in `samosa/samosa_summary.csv`, the protocol table
in `samosa/global_output_HEXTOR.dat`, the figure in `samosa/samosa_hextor.png`.

| case | S | p [bar] | T_glob | T_min | T_max | TOA imb. | drift [K/dec] | dt [s] | beta |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 500 | 0.70 | 173.10 | 148.46 | 210.53 | +0.63 | -0.02 | 1350 | 4.5 |
| 4 | 1200 | 2.34 | 292.22 | 244.33 | 334.69 | +2.21 | -0.01 | 1350 | 5.9 |
| 8 | 800 | 6.16 | 224.40 | 184.26 | 287.50 | +1.27 | -0.01 | 1350 | 1.4 |
| 9 | 1100 | 0.70 | 267.86 | 224.33 | 306.81 | +1.98 | -0.01 | 1338 | 6.0 |
| 10 | 400 | 4.83 | 152.42 | 131.72 | 183.61 | +0.51 | -0.01 | 1350 | 1.5 |
| 11 | 900 | 0.10 | 225.52 | 188.82 | 265.15 | +1.40 | -0.02 | 415 | 19.4 |
| 14 | 900 | 1.44 | 241.46 | 197.22 | 297.56 | +1.49 | -0.02 | 1350 | 2.7 |
| 15 | 600 | 0.43 | 188.62 | 160.11 | 232.25 | +0.75 | -0.01 | 1234 | 6.5 |
| 16 | 1400 | 10.0 | 376.07 | 354.67 | 399.96 | +2.45 | +0.04 | 952 | 8.4 |

Cases 2, 3, 5, 6, 7, 12 and 13 run away, and are omitted from the submission
under the protocol's allowance for incipient runaway.  **The set of cases with
a climate solution is the same nine as the published submission.**

**No bistability.**  Warm (300 K) and cold (233 K) starts land on the same
state in all nine, agreeing to 0.003 K.  The published run tested only the
warm start.

**Case 11 got dt = 415 s**, the smallest in the set: it has the largest beta
(19.4), because q ~ esat/p_dry and its dry pressure is 0.1 bar, the thinnest in
the protocol.

**How much the reduced `dt` was actually worth: nothing measurable, here.**
Every case was re-run at a quarter of its chosen `dt` (cases 1, 4, 11, 16:
agreement 0.000, 0.001, 0.001, 0.014 K) *and* at the published 1350 s (cases
9, 11, 15, 16: 0.000, 0.002, 0.000, 0.011 K).  So at D = 2.46 these runs are
stable at 1350 s even where the beta-corrected criterion puts the CFL number at
1.3 (case 11).  The criterion is conservative: it charges the full beta against
a coefficient that diffuses h/cp while the prognostic update is on T.  The
beta correction is therefore **insurance rather than necessity at this D** --
it costs a factor of ~3 in steps on one case and buys protection against the
silent integration to NaN that the dry criterion already documents.  It would
bite for real under `--transport perbar` or `diffadj`, where D reaches 33 and
197 and the dry criterion alone already demands a smaller step.

**Case 16 is now inside the table.**  It equilibrates at 376 K, where on the
RH 1 table with the published constants it sat at 465 K -- above the 420 K
ceiling, so that value was extrapolated and should not have been read as a
climate.  The RH 0.8 table's runaway plateau is 307 W/m^2 rather than 292,
which is what buys the headroom.

## Cloud scaling variant

`--cloud-scaling instellation` (cloudir = -45.09 x S/900) was run as well, into
`samosa/samosa_summary_cloudscale_instellation.csv`, warm start.  It gives
**eleven** equilibria rather than nine: cases 5 (293 K) and 12 (336 K) no
longer run away, and case 7 stalls at 418 K on the plateau.  Subtracting the
per-case cloudir from the reported OLR puts cases 5 and 12 at a clear-sky OLR
of ~287 W/m^2, about 20 W/m^2 below the 307 W/m^2 plateau, so unlike the RH 1
result (where case 12 "equilibrated" at 290 W/m^2 against a 292 W/m^2 plateau,
i.e. on it) these two are genuinely constrained.  The primary result above
keeps `constant` scaling, matching the published submission and the definition
the calibration itself uses.

## Files

- `model/radiation/radiation_N2_CO2_2600K_p_rh0.8.h5` (new, 1.5 MB, 70
  columns, 1 h 14 min, 0 failed).  Axes identical to the RH 1 2600 K table
  except for 10-digit truncation of `fco2` (max 3.6e-10 relative).
  `check_table_sanity`: all hard checks OK; runaway plateau 306.8 W/m^2
  against ~292 at RH 1.  Build log `model/radiation/build_2600K_rh0.8.log`.
- `samosa/samosa_summary.csv`, `global_output_HEXTOR.dat`,
  `samosa_hextor.png/.pdf` -- the new submission.
- `samosa/samosa_summary_decomp_*.csv` -- the four decomposition runs.
- `samosa/samosa_summary_cloudscale_instellation.csv` -- the variant.
- The previous (RH 1, dry) submission's summary and figure are in git as of
  commit 358d484.  Its `global_output_HEXTOR.dat` and zonal files, which git
  never tracked, are archived in `/models/data/samosa/hextor/old_rh1_drydiff/`.

## Still open

- The 2600 K and 3000 K RH 0.8 tables still use the fixed `min(200 K, Ts)`
  stratosphere.  That is step 2 of the plan in `subsaturated_rce_tables.md`
  and is untouched here.
- `tools/samosa_sensitivity.py` / `samosa_summarize.py` have not been re-run
  against the new table, so the cloudir-sweep robustness statement in
  CLAUDE.md still refers to the RH 1 results.
- Organic haze is still absent from the tables (ExoRT `calc_opd_mod.F90`), but
  these cases carry only 400 ubar of CO2 and no CH4, so it does not bite here.

## Comparison with the other SAMOSA submissions

`tools/samosa_compare.py` (new) scores any set of `run_samosa.py` output
directories against the other models' submissions:

    python tools/samosa_compare.py "published=OLD" "D const=samosa"

Two references, answering different questions.

- **ExoCAM** (`/models/data/samosa/exocam/`), a 3-D GCM with clouds.  The
  closest thing to a reference answer, and the only one that constrains the
  day-night contrast.  Global means come from its `analysis.py` digest
  (area-weighted); the extrema come from the `TS` maps.
- **ExoColumn** (`/models/data/samosa/exocolumn/`), a 1-D RCE column using the
  **same radiation core** HEXTOR's tables are built from (ExoRT n68equiv).  It
  has no horizontal dimension, so it says nothing about contrast, but for the
  global mean it isolates everything *except* the spectroscopy.

Neither is a calibration target: `(d0, cloudir)` are fitted at THAI Hab1 around
a 2600 K star, so all of this is out of sample.

### Global mean

| | vs ExoCAM (n=9) | vs ExoColumn (n=8) |
|---|---|---|
| published, RH 1, dry | RMSE 45.2 K, bias +9.3 | RMSE 26.1 K, bias −16.7 |
| new, RH 0.8, moist, D constant | **RMSE 24.1 K**, bias −4.1 | RMSE 25.8 K, bias −20.7 |
| new, cloudir × S/900 | **RMSE 15.8 K**, bias −5.6 | **RMSE 16.6 K**, bias −15.7 |

Most of the ExoCAM improvement is Case 16, whose published 465 K was
extrapolated past the table ceiling and was never a temperature.  Excluding it,
RMSE goes 28.9 → 24.6 K.

### Day-night contrast, vs ExoCAM

RMSE 76.0 → 33.1 K (44.8 for the instellation-scaled clouds, which is worse —
it over-cools Case 16 and re-inflates its contrast).  Cases 4 and 9 go from
40–60 K wrong to within 4 K, and `d0` *fell* to 2.46 while that happened, so it
is latent transport doing it and not a stronger diffusion coefficient.

### The two findings that matter

**On the warm cases HEXTOR now reproduces the same-RT 1-D model almost
exactly.**  Against ExoColumn, Case 4 is −1.0 K and Case 9 −1.8 K, where the
published configuration had them at +15.1 and +8.5.  Since ExoColumn shares the
radiative transfer, that is a statement that HEXTOR's transport and surface
treatment are right for those cases, not that its spectroscopy happens to
cancel an error.

**The whole remaining disagreement is the cold, glaciated end, and it is in the
substellar temperature, not the night side.**  Decomposing the contrast error
against ExoCAM:

| case | T_max err | T_min err |
|---|---|---|
| 1 | **−58.4** | +1.0 |
| 15 | **−53.3** | +2.2 |
| 10 | **−46.4** | −34.8 |
| 11 | **−25.9** | −6.4 |
| 16 | +39.3 | +0.9 |

On Cases 1 and 15 the antistellar temperature is right to 2 K and the
substellar point is 53–58 K too cold.  ExoCAM keeps a warm substellar pool at
269–286 K under a fully glaciated surface (ICEFRAC = 1.0, CLDTOT 0.69–0.79);
HEXTOR's clear-sky column at 0.4–0.7 bar and ~210 K has almost no greenhouse,
and the uniform `cloudir` = −45 W/m² costs it a further ~20 K where the
absorbed flux is only ~100–125 W/m².  The instellation-scaled clouds recover
about 15 K of that (Case 1 substellar 225.8 against 210.5) and no more.  This
is a structural limit of a clear-sky table plus a uniform OLR offset, not a
tuning error.

### Energy budget: the right temperature from two large compensating errors

| case | albedo HEXTOR / ExoCAM | OLR HEXTOR / ExoCAM |
|---|---|---|
| 4 | 0.038 / 0.307 | 287.6 / 207.7 |
| 9 | 0.037 / 0.342 | 263.9 / 181.1 |
| 14 | 0.068 / 0.286 | 208.9 / 161.2 |
| 16 | 0.042 / 0.165 | 334.3 / 291.4 |

On the warm cases HEXTOR absorbs ~80 W/m² too much shortwave and emits ~80
W/m² too much longwave.  That is exactly what `cloudir` is: a longwave stand-in
for a missing cloud shortwave albedo.  The surface temperature comes out close
because the two errors are tied together by the calibration, but the fluxes
themselves are not comparable with the GCMs, and the protocol asks for them.

## Pressure-dependent transport (2026-09-16)

Tested at Jacob's request, because the constant-D contrast residual had a
pressure signature: too little contrast below 1 bar (Case 15, 0.43 bar:
−55.5 K), too much at 10 bar (Case 16: +38.4 K).  All runs: RH 0.8 table,
`moistdiff`, `cloudir` −45.09, warm start.

| `--transport` | D | calibration | runs |
|---|---|---|---|
| `constant` | `d0` | 2.46 | `samosa/samosa_summary.csv` |
| `psqrt` (new) | `d0 · √p` | 2.46, unchanged | `samosa_summary_transport_psqrt.csv` |
| `perbar` | `d0 · p · comp` | **2.24** (effective D 2.466 at Hab1) | `samosa_summary_transport_perbar.csv` |

`perbar` was recalibrated with `calibrate_thai.py --transport perbar`, which
needed a fix first: the Hab1 template carried no `diffadj_rot`, so the old
`--diffadj` flag picked up the driver default (`.true.`) and would have
calibrated the full rotation-scaled transport — a factor 40.9 at Hab1, not
1.10.  `calibrate_thai.py` now takes the same `--transport` names as
`run_samosa.py`.  `psqrt` needs no recalibration: the reference point is 1 bar,
where √p = 1.

### Scores

| vs ExoCAM | constant | √p | p |
|---|---|---|---|
| global mean RMSE | 24.1 K (n=9) | **24.0 K** (n=9) | 25.6 K (n=8) |
| contrast RMSE | 33.1 K | **28.3 K** | 33.4 K |
| contrast bias | −8.8 K | −14.1 K | −23.6 K |
| vs ExoColumn, mean RMSE | **25.8 K** | 26.4 K | 28.7 K |

**D ∝ p is a wash, and it costs Case 11.**  It is spectacular at the two
extremes — Case 16 contrast 45.3 → 5.1 K (ExoCAM 6.9), Case 15 72.1 → 128.9 K
(ExoCAM 127.7) — and breaks the middle: Case 8 (6.16 bar) goes from +35 to −49,
Case 10 (4.83 bar) from −12 to −51, Case 4 from +2 to −38.  At 0.1 bar it gives
D = 0.246, the night side collapses to 92 K (ExoCAM 195 K) and **Case 11 loses
its climate solution**, the same failure the dry-diffusion analysis found.
Latent transport does not rescue it.

**D ∝ √p is a modest, genuine improvement.**  Contrast RMSE 33.1 → 28.3 K with
the global mean unchanged and all nine equilibria kept.  Case 16 contrast falls
to 15.4 K and Case 11's global mean lands at 232.9 K against ExoCAM's 234.0.
It is not a clean win — Cases 4 and 10 get worse, and Case 11's contrast
overshoots to 130 K — it redistributes the error more evenly.

**Transport barely touches the global mean** (24.0–25.6 K across all three),
as the dry analysis also found: diffusion moves energy around, it does not
change how much the planet absorbs or emits.

### Why no transport choice can close the gap

The cold-case contrast error is in the substellar temperature, and transport
moves it only a little:

| case | T_max error, constant | √p | p |
|---|---|---|---|
| 1 (0.70 bar) | −58.4 | −53.2 | −47.6 |
| 15 (0.43 bar) | −53.3 | −39.2 | −24.8 |
| 10 (4.83 bar) | −46.4 | −60.4 | −68.1 |

Even D ∝ p, which starves the thin cases of transport and so keeps the most
heat under the star, leaves Case 1's substellar point 48 K too cold.  What
ExoCAM has and HEXTOR lacks there is a greenhouse and cloud warming over the
substellar pool, not a slower drain of heat away from it.

### Decision left open

The primary submission stays at `constant`.  √p is better against ExoCAM, but
it was chosen *after* comparing with ExoCAM, and the exponent is a round
interpolation between the two published options rather than something derived.
Adopting it on this evidence would be tuning to one GCM in the ensemble the
submission is meant to be compared against.  If there is an independent
argument for sub-linear scaling (e.g. from the Koll & Abbot 2015 or Wordsworth
2015 day-night scalings), that would be the basis for switching; fitting α to
ExoCAM would not.

### A harness defect this exposed

`perbar` first ran for two hours and still failed three cases.  Case 16 went
**NaN at year 21** from a timestep that was stable for its equilibrium but not
for its transient (342 K at year 6), then integrated NaN to `niter = 5000`.
The dt refinement should have retried it at a smaller step, but the timeout
returned `rc != 0`, which broke out of the retry loop first — so an instability
was reported as a runaway.  Two fixes in `run_samosa.py`:

- `_integrate()` now polls the annual series and kills a run the moment it goes
  non-finite, returning a distinct code (`RC_NONFINITE`).  The whole `perbar`
  sequence then took seconds, and Case 16 equilibrated at 373.2 K.
- A failed or timed-out pass now feeds the retry instead of aborting it, and
  `classify()` labels the outcome `unstable` or `no result (timeout)` rather
  than folding it into `runaway`.

Audited the production (`constant`) run for the same contamination: no NaN in
any case's series, and every runaway climbs monotonically to a large finite
temperature, so **the submission is unaffected**.

## Constant D or D ∝ √p: held-out test against the whole ensemble (2026-09-16)

√p was chosen after looking at ExoCAM, so its ExoCAM score cannot justify it.
The other resolved GCMs were not used, and are a held-out test.  Scored with the
validated readers in `~/research/samosa/figures/allcases/extract_maxtemp.py`
(each model's own area weighting; every global mean checked against the values
embedded in the figure scripts, all within 1.5 K):

| model | n | mean RMSE const / √p | contrast RMSE const / √p |
|---|---|---|---|
| ExoCAM (chosen against) | 9 | 24.1 / 24.0 | 33.1 / **28.3** |
| ExoPlaSim | 9 | 21.3 / 20.7 | 53.6 / **51.5** |
| ROCKE-3D | 9 | 33.0 / 32.9 | **28.7** / 29.3 |
| Generic PCM | 7 | 29.6 / 30.3 | 33.4 / **31.0** |
| LFRic | 9 | 24.9 / 24.4 | 38.3 / **31.8** |
| PlaHab | 9 | 36.6 / 37.5 | **52.8** / 55.3 |
| **4 held-out GCMs, pooled** | 34 | 27.4 / 27.4 | 40.0 / **37.3** |

- **The global mean does not distinguish them** at all (27.4 / 27.4 K).
- **The contrast gain is real but small, and halves out of sample**: 15 % on
  ExoCAM, 7 % on the held-out GCMs.  Three of the four prefer √p; ROCKE-3D and
  PlaHab prefer constant.  That shrinkage is the signature of partly fitting the
  model the option was picked against.
- **Against the GCM spread it is a tie.**  Each option falls outside the range
  of the resolved GCMs' contrasts on 5 of 9 cases, just different ones: √p
  brings Cases 11 and 16 inside and pushes 8 and 14 out.  Cases 1, 10 and 15 are
  outside under both — the structural cold substellar deficit.

Caveat common to both: HEXTOR's Tmax/Tmin are substellar and antistellar belt
means, the GCMs' are field extrema, so HEXTOR's absolute contrast is biased low.
That affects both options equally and does not change the comparison between
them.

## Submitted (2026-09-16)

**Constant D**, chosen by Jacob on the held-out test above.  Installed in
`/models/data/samosa/hextor/`: `global_output_HEXTOR.dat` and the nine
`zonal_output_HEXTOR_caseNN.dat`, copied byte-for-byte from `samosa/`.  The
previous version is archived byte-identical in `old_rh1_drydiff/`, following the
`old_pco2_400ppm/` convention, and `README_zonal.txt` has a revision entry.

Checked before installing: the README's claim that the sin(theta)-weighted belt
means reproduce the global file's Tglob, Tmax and Tmin holds for all nine cases
(worst 0.002 K).  Its Case 1 meridian example was recomputed with a method that
first reproduced the old numbers exactly (202.9 / 195.9 → 210.5 / 201.3 K).

Two README cautions, both measured with cloudir ±5 W/m²:
- **Case 16**: day side still at the plateau (substellar clear-sky OLR 306 W/m²
  against ~307), 16 K per 10 W/m² against 5–9 K for the other cases.  In the
  previous version −5 W/m² moved it 31 K and +5 W/m² made it a runaway.
- **Case 10**: at cloudir −50 it drops into CO2 condensation (91.7 K), which
  HEXTOR detects but cannot represent.

Not yet updated: `~/research/samosa/figures/allcases/` still embeds the old
HEXTOR global means, and `extract_maxtemp.py` still leaves Case 16 out of
`ACCEPTED`.

