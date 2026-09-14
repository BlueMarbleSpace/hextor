# Adding CH4 to the HEXTOR radiation lookup tables

Scoping note, 2026-09-01; **implemented 2026-09-12**, see the addendum at the
end for what the measurements changed.  Everything below was measured on this
machine, not estimated from first principles, except where marked.

## Status

Implemented as **Route B** (CH4 as a table axis) with CH4 held constant within
a run: `&ebm::fch4`, no CH4 cycle.  The Sun table is
`model/radiation/radiation_N2_CO2_CH4_Sun_p.h5`, built with
`tools/make_radiation_table.py --ch4-axis`.  See `namelists/input.nml.earth.ch4`.

## Short answer

The radiative transfer side needs **no new work at all** — ExoRT and ExoColumn
already do CH4 today. The cost is compute, and how much code you touch depends
entirely on one question:

**Is CH4 a scenario parameter (fixed per experiment) or a prognostic variable
(varying within a run)?**

- **Scenario** → generate one table per CH4 abundance and switch with `radfile`.
  Zero changes to `radiation.f90` or `driver.f`. Roughly a day, most of it
  waiting for the generator.
- **Prognostic** → add CH4 as a sixth table axis. Half a day of code, plus a
  CH4 cycle in `driver.f`, which does not exist in any form and is the larger
  piece of work.

## What already exists

`ExoRT/source/src.n68equiv/kabs.F90` reads CH4 k-coefficients from
`data/kdist/n68ch4/hitran2016/` (`n68_8gpt_ch4_hitran16_...nc`) and also carries
CO2-CH4 collision-induced absorption (`CO2-CH4_cia_n68.nc`). ExoColumn plumbs
`&exocol_composition::ch4_vmr` through `exocol_coldstart` to `aerad_driver`.

So the generator works unmodified. Measured with the existing
`exocol_sweepts200.exe` at 1 bar, 400 ppm CO2, Ts = 280 K, 3000 K blackbody,
mu = 0.5, surface albedo 0.3:

| CH4 vmr | OLR (W/m2) | planetary albedo |
|---------|-----------|------------------|
| 0       | 239.777   | 0.1823           |
| 1e-6    | 237.055   | 0.1780           |
| 1e-4    | 227.493   | 0.1503           |
| 1e-3    | 221.973   | 0.1273           |
| 1e-2    | 217.058   | 0.0963           |

Both effects are large enough to be worth tabulating: 23 W/m2 in OLR and 0.086
in albedo across 0 to 1%. The albedo response matters as much as the
greenhouse one here, because CH4 absorbs in the near-infrared where a 3000 K
star emits most of its flux — this would be a weaker effect under the Sun.

## Route A: separate tables per CH4 abundance (recommended if scenario)

Generate one table per CH4 value and select it with `&radiation::radfile`.

    python tools/make_radiation_table.py --star <SED> \
        --out model/radiation/radiation_N2_CO2_CH4_<ppm>_<star>.h5 ...

`make_radiation_table.py` needs one small change: `ch4_vmr` is currently
hard-wired to 0.0 in the `NML` template, so add a `CH4` config constant and a
`--ch4` override. Nothing else changes anywhere.

Cost is one full table generation per CH4 value: about 3.5 h each at 4 workers
on this machine for the standard grid (15 pressures x 14 CO2 x 31 temperatures
x 12 zenith x 9 albedo).

## Route B: CH4 as a sixth table axis

### radiation.f90

- `olr_table` rank 3 -> 4, `palb_table` rank 5 -> 6, plus a `/ch4` axis vector
  and `nch4`.
- Interpolation 8 -> 16 corners for OLR, 32 -> 64 for the albedo. The corner
  loops are already nested `do` loops with `merge()`, so this is mechanical.
- Interpolate in log10 along the CH4 axis, as for pressure and CO2, but note
  the axis must then reach a nonzero floor: decide whether CH4 = 0 is a level
  (breaking log interpolation) or whether the floor is, say, 1e-9.
- Backward compatibility is free using the pattern already in place for
  pressure: detect the absence of `/ch4`, load with `nch4 = 1`, and the CH4
  argument clamps away. v1 and v2 tables keep working through the same path.

### driver.f

This is the real gap. **HEXTOR has no CH4 concept at all** — no `fch4`, no
partial pressure, nothing (`grep -i ch4 model/driver.f` returns nothing).
Needed:

- a namelist variable in `&ebm`, alongside `fco2` and `fh2`;
- inclusion in the partial-pressure bookkeeping near driver.f:305, which
  currently sets `pco2 = pg0*fco2` and `pn2 = pg0 - pco2`;
- pass-through at both radiation call sites (the `getOLR` and `getPALB` calls);
- a decision on whether CH4 enters `avemol` and `hcp` for the diffusion scaling
  at driver.f:1599 — it should, since CH4 has mw 16 against N2's 28.

If CH4 is to be prognostic, that is a further piece of work analogous to
`do_cs_cycle` or `do_h2_cycle`, and larger than the table itself.

### tools

`make_radiation_table.py` gains an outer axis; `check_table_reader.py`,
`check_table_sanity.py` and `check_table_interpolation.py` each need the extra
dimension in their reference implementations. `merge_radiation_tables.py`
would need to learn which axis it is joining along.

## Cost and size

CH4 is an outer loop like pressure and CO2, so cost scales linearly in the
number of levels. Measured baseline: the current 6510-profile table took about
3.5 h at 4 workers, plus 35 min for the 6-level hot-end extension.

| CH4 levels | profiles | radiation calls | wall time (4 workers) | file size |
|-----------|----------|-----------------|----------------------|-----------|
| 3  | 19530 | 2.1M | ~10 h | 17 MB |
| 4  | 26040 | 2.8M | ~14 h | 23 MB |
| 5  | 32550 | 3.5M | ~18 h | 28 MB |
| 6  | 39060 | 4.2M | ~21 h | 34 MB |

Note the tables barely compress — gzip-6 on the current table gives only 1.1x,
because the data is smooth but stored as doubles. A five-level table at 28 MB
is past what belongs in git; keep it out of the repository and distribute
separately, which CLAUDE.md notes was the original convention for these tables
anyway.

Worker count: this machine saturates at 3-4 workers (6 physical cores, and the
job is memory-bandwidth bound at ~240 MB per process). Do not raise it.

## Cautions carried over from the pressure work

- **Timestep.** If CH4 changes the calibrated diffusion, re-check
  `D*dt/(C*dx^2) < 0.5`. Exceeding it does not raise an error; the model
  integrates quietly to NaN over thousands of years.
- **Re-calibration.** Any table change invalidates the existing `(d0, cloudir)`
  calibrations. Re-run `tools/calibrate_thai.py`, which solves for `cloudir` at
  fixed `D` because the global mean barely depends on `D`.
- **Stale probe.** Run `make clean` in `model/radiation/` before
  `check_table_reader.py`, or it can validate a stale binary against a stale
  checker and report OK while testing nothing.

## Likely motivation

Probably `~/research/ebmEarlyEarth` — Archean Earth is the classic CH4
application, and would want a solar SED (`G2V_SUN_n68.nc`) rather than the
M-dwarf tables built for SAMOSA. If CH4 is a prescribed scenario value there,
Route A is the whole job.

---

# Addendum, 2026-09-12: what implementation changed

Decisions taken: CH4 constant within a run (namelist `fch4`, no cycle); solar
SED only; axis capped at 1e-1; CO2 axis trimmed to pay for it; the CH4-floor
slice used as the acceptance test; table kept in `model/radiation/`.

## The node count was the thing the note got wrong

The note guessed 3 to 6 CH4 levels.  The axis actually needs **13**, and the
binding constraint is the PLANETARY ALBEDO, not the OLR -- CH4 absorbs in the
near-infrared, so the albedo responds as strongly as the greenhouse does.
Measured against a 0.25 dex reference scan over 1e-8 to 1e-1, across all 12
zenith x 9 albedo nodes at several temperatures, under the solar SED:

| spacing | nodes | max dOLR (W/m2) | max d(albedo) |
|---|---|---|---|
| 0.50 dex uniform | 15 | 0.19-0.26 | 0.0017-0.0023 |
| 0.75 dex uniform | 11 | 0.25-0.31 | 0.0040-0.0061 |
| 1.00 dex uniform | 8 | 0.65-0.67 | 0.0059-0.0095 |
| 1 dex below 1e-6, 0.5 dex above | **13** | **0.43-0.50** | **0.0019-0.0023** |

The budgets to beat are ~1 W/m2 (the vertical-resolution uncertainty, per
`check_table_convergence.py`) and 0.003 in albedo (the zenith grid, per
`check_table_interpolation.py`).  Uniform 0.75 dex misses the albedo budget by
2x; the non-uniform spacing meets both on 13 nodes rather than 15, because the
CH4 response is nearly flat below 1e-6 and steep above it.  That spacing is
`CH4_NODES` in `make_radiation_table.py`.

Under a 3000 K blackbody the albedo errors are about twice as large again, so
an M-dwarf CH4 table would need a finer axis than this solar one.

## A nonzero floor replaces CH4 = 0

The axis floors at 1e-8 rather than including zero, so the log10 interpolation
is well posed, and `fch4 = 0` simply clamps to it.  That is radiatively
indistinguishable from CH4-free (1 bar, 400 ppm CO2, 280 K, mu = 0.5,
alpha_s = 0.3):

| ch4_vmr | dOLR (W/m2) | d(albedo) |
|---------|------------|-----------|
| 1e-9    | -0.008     | -0.00001  |
| 1e-8    | -0.076     | -0.00011  |
| 1e-7    | -0.475     | -0.00078  |

## The CO2 axis paid for part of the CH4 axis

Measured with CH4 present, the CO2 axis was trimmed from 14 nodes (0.44 dex)
to 11 (0.57 dex): 0.57 W/m2 and 0.0023 in albedo, against 0.375 and 0.0015
before.  Both inside budget.  12 nodes (0.52 dex) is the middle option at
0.497 and 0.0020.  Since the axis multiplies the whole build, this bought back
about 11 h.

Note the consequence: three log axes now each contribute ~0.002 in albedo
where the zenith axis contributes 0.003.  If that ever needs tightening, the
CO2 axis is the cheapest to put back.

## Cost, as built

2145 columns (15 p x 11 CO2 x 13 CH4), 5.79M radiation calls, ~41 h of wall
time on 4 workers, 45 MB.  The note's 3.5 h per CH4 level was right: the
measured Sun build was 210 columns in 4.0 h, i.e. 0.101 s per radiation call
rather than the 0.049 s the cost estimator used to assume (now corrected).

## Two bugs found on the way

- **Uninitialised annual accumulators** (`driver.f`).  `tempavesum`, `nstep`
  and the rest are zeroed at label 790, at the *end* of each year, so the
  first year accumulated onto whatever was in memory and the year-0 row of
  `tempseries.out` was a function of the variable layout.  Adding the CH4
  variables changed it from 6.295 K to 288.051 K -- the latter is correct,
  which is exactly the problem.  Now zeroed at startup.  Years 1 and up are
  bit-identical either way.
- **Relative `--out` broke the workers** (`make_radiation_table.py`).  Worker
  scratch paths were derived from `--out` without being made absolute, and
  each worker runs ExoColumn with `cwd` set to its own directory, so a
  relative path left the symlinked executable unresolvable.  Every previous
  build had used an absolute path and missed it.

## Cache names are now physical, not positional

Trimming the CO2 axis exposed a hazard: cache files were named by axis index,
so `col_000_005` still existed after the regrid but referred to a different
CO2 mixing ratio -- a resume would have loaded the wrong physics and reported
a full cache hit.  Names now carry the column's values, so a regridded axis
misses and recomputes.  `tools/migrate_cache_names.py` renames an old cache
given the table whose axes it was built on (the 210-column Sun cache was
migrated this way rather than discarded).

## Acceptance test on the finished table (2026-09-14)

Built 2026-09-13, 2145/2145 columns, no failures, 46.8 MB.
`check_table_reader.py` agrees to 4.8e-11; `check_table_sanity.py` is OK, with
35/49500 CH4 non-monotonic pairs, all <= 0.05 W/m2 and in the cold,
high-pressure, CO2-rich columns that also produce the existing CO2 notes.

**Table level.** Pressure, temperature, zenith and albedo axes are identical
to `radiation_N2_CO2_Sun_p.h5`; the CO2 axes share only their end nodes.  At
fCO2 = 0.5 the CH4-floor slice matches the CH4-free table to 0.02 W/m2 and
0.001 in albedo.  At fCO2 = 1e-6 the maximum is 0.77 W/m2 and 0.0055, and it is
physics rather than a defect: zero at 1 bar and below, growing with column
mass to 20 bar, peaking at 250-290 K, i.e. 10 ppb of CH4 in a thick column
with no CO2 to overlap its bands.  The CH4 axis runs 201.0 -> 198.0 -> 194.3
W/m2 through 1e-8, 1e-7, 1e-6 at that point, against 201.8 CH4-free.  So the
floor is "CH4-free" to within budget everywhere except thick, CO2-poor
atmospheres, and marginally over the albedo budget only at extreme corners
(20 bar, 87 deg zenith, white surface).  A lower floor level (1e-10) would
cost one more CH4 level, ~165 columns, if that regime ever matters.

**Model level.** `input.nml.earth.ch4` at fch4 = 0 and S/S0 = 0.90, which
equilibrates (285.95, 285.67, 285.55, 285.52 K): the CH4 table reproduces the
CH4-free v2 table to **0.06 K** at every year.  So calibrations made against
`radiation_N2_CO2_Sun_p.h5` carry over.

**But there is no Earth calibration for the pressure-resolved Sun tables to
carry over.**  `input.nml.earth` (cloudir = 3.0, d0 = 0.58, no clouds) was tuned
on the v1 table.  Against the v2 or v3 Sun table at S/S0 = 1 it equilibrates at
367 K (an earlier run halted at 348 K still warming 0.7 K/yr: the default
`fluxcnvg = 0.1` is too loose there; 1e-4 was used for the numbers here).
The response to instellation is smooth: 0.90 -> 285.5, 0.92 -> 297, 0.94 -> 307,
0.96 -> 318, 0.98 -> 335, 1.00 -> 367 K.

This is the old calibration, not the new table.  Where both tables see the same
surface they agree to 0.003 in albedo.  The v1 table differs in two ways that
had been acting as hidden calibration: its surface-albedo axis starts at 0.2,
so the 0.06 ocean was clamped to 0.2, and its OLR never saturates (341 W/m2 at
320 K, against ~285 for v2).  Global means at equilibrium:

| run | T (K) | global albedo | OLR (W/m2) |
|---|---|---|---|
| v1 table, no clouds | 285.8 | 0.261 | 252 |
| v2 table, no clouds | 366.9 | 0.152 | 289 |
| v2, ocean albedo pinned at 0.2 | 304.0 | 0.217 | 267 |
| v2, S/S0 = 0.90 | 285.5 | 0.204 | 244 |

(Global albedo is the model's own insolation-weighted `planet average albedo`
from model.out.)  v1 gave cloud-free Earth 0.26 because the clamp acted as
partial cloud cover; v2 gives the clear-sky 0.15.  Pinning the ocean at 0.2
recovers 63 of the 81 K; the rest is v1's missing OLR saturation.  The
namelist's own default cloud option (fcloud = 0.5, zenith-dependent cloud
albedo, cloudir = 3) over-corrects into a snowball.

## Cloud calibration (2026-09-14)

`tools/calibrate_earth.py`, targets 288 K and global albedo 0.30, at 280 ppm
CO2 and 730 ppb CH4 on the CH4 table: **fcloud = 0.4427, cloudir = 10.28 W/m2**
(d0 = 0.58 unchanged) -> 287.85 K, albedo 0.300.  Written into
`namelists/input.nml.earth.ch4`.  Pre-industrial CH4 is included on purpose:
calibrating at zero would fold its greenhouse into the clouds and shift the
zero of every CH4 experiment.

- Warm branch for cloudir ~0-22 W/m2 at this fcloud; snowball below S/S0 ~0.95,
  runaway above ~1.05.
- Local sensitivity 1.65 K per W/m2 (~6 K per CO2 doubling), rising toward
  the runaway -- about twice Earth's.
- cloudir ~10 against a real LW cloud effect of ~+25-30 W/m2: the RH = 1
  column most likely over-states the clear-sky greenhouse, and cloudir absorbs it.
- Zonal: equator 303 K, poles 257-262 K, ice lines -77 / +64 deg; the
  pole-equator contrast is somewhat weak with the v1-era d0.
- Same clouds on the CH4-free v2 table: 283.7 K, vs 283.8 K for the CH4 table
  at fch4 = 0 -- so pre-industrial CH4 is worth ~4 K in this model.

CH4 sweep on the calibration (S/S0 = 1, 280 ppm CO2), dT against 730 ppb:
0 -> -4.0 K, 1 ppm -> +0.5, 10 ppm -> +11.1, 100 ppm -> +25.6 K; 1000 ppm and
above run away.  Those runs report "Flux converged" at 430-630 K with no
drift, non-monotonically in CH4, because above the table's 420 K ceiling the
OLR is extrapolated flat -- they are runaways, not climates.  The narrow warm
band (+-5% in S) matters for Archean work: at S/S0 ~0.8 this model needs a
large greenhouse to avoid the snowball and not much more to run away.

A CH4 sweep at S/S0 = 0.90 is monotonic with no steps at the axis nodes
(+4.6 K at 1 ppm, +23 K at 100 ppm, +48 K at 1%), but every case from 10 ppm
up was still warming at halt and the top of the range is heading for the
runaway plateau, so these are not equilibria and not calibrated results.
The implied sensitivity (~1.4 K per W/m2 from 0 to 1 ppm) is high, as
expected for a clear-sky RH = 1 column in a base state near the runaway.

## Still not done

- **Organic haze remains the known gap.**  `calc_opd_mod.F90` carries no haze
  optics -- only a note to hook up CARMA aerosols -- so above CH4/CO2 ~ 0.1 a
  real Archean atmosphere would have an anti-greenhouse the table cannot
  represent.  The axis reaches 1e-1 deliberately, and the table records the
  caveat in its own `ch4_caveat` attribute; `check_table_sanity.py` prints it.
- **SAMOSA and THAI may need revisiting** in light of the cloud treatment: they
  use cloudir alone, and the published THAI calibration went through the v1
  ocean-albedo clamp, which pushes cloudir the same way as the zenith-grid
  error the RNAAS note credits with the -35 -> -45 shift.
- **d0 was not re-fitted** for the pressure-resolved Sun table; the
  pole-equator contrast is somewhat weak.
- **No M-dwarf CH4 table.**  Needs a finer axis than this one (see above).
- **A prognostic CH4 cycle**, if it is ever wanted, needs no table work.
- **`co2cldavesum` is never reset** at label 790 with the other accumulators,
  so `ann_co2cldave` is cumulative across years.  Pre-existing, untouched,
  and unrelated to CH4 -- noted here only because it was found nearby.
