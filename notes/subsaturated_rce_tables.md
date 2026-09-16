# Sub-saturated, RCE-consistent radiation tables — plan

Written 2026-09-15 at the end of an ebmEarlyEarth session. **Status: planned,
not started.** Nothing in HEXTOR or ExoColumn has been changed for this yet.

## The request (Jacob)

Rebuild **all the pressure-resolved lookup tables** with

1. a **sub-saturated troposphere** (RH < 1), and
2. a **stratospheric temperature consistent with ExoColumn's
   radiative-convective equilibrium (RCE) solutions**.

The legacy v1 tables (`radiation_N2_CO2_Sun.h5`, `radiation_N2_CO2_2600K.h5`)
are out of scope.

## What the tables assume now

All four pressure-resolved tables carry the attribute
`water = "variable_ps: total ps = p_dry + esat(Ts); prescribed moist adiabat, RH = 1.00"`
and `t_strato = "min(200.0 K, Ts)"`. Verified in `tools/make_radiation_table.py`.
`RH = 1.0` is a module constant with no CLI override, and it has been 1.0 since
the builder was first committed (`e95bc1c`). `t_strato_for(ts) = min(200, ts)`
is used for every pressure, CO2, CH4 and host star.

| table | build cost (4 workers) |
|---|---|
| `radiation_N2_CO2_Sun_p.h5` (v2) | ~4 h |
| `radiation_N2_CO2_CH4_Sun_p.h5` (v3, 2145 columns) | ~41 h |
| `radiation_N2_CO2_3000K_p.h5` (SAMOSA) | ~4 h |
| `radiation_N2_CO2_2600K_p.h5` | ~4 h |

Why it matters: an RH = 1 column overstates the clear-sky greenhouse. This is
probably why the calibrated `cloudir` is ~10 W/m² against a real LW cloud
effect of 25–30, and why sensitivity is ~2× Earth's (see CLAUDE.md,
"Pre-industrial Earth calibration"). It also biases the ebmEarlyEarth v3 CH4
solutions, whose global means are 45–66 °C, far from where the calibration was
fitted.

## Feasibility tests already done (2026-09-15)

A single column: Sun (G2V_SUN_n68), S/S0 = 0.75 (`msdist = 1.3333`), 1 bar
dry, fCO2 = 0.05, no O2/O3, `coszrs = 0.5`, surface albedo 0.25, Ts = 300 K,
`exocol.exe` (PVER = 70), `h2o_eos = 'nonideal'`, `variable_ps`,
`cold_trap_phase = 'ice'`, `p_top = 1 Pa`.

### 1. Fixed-surface-temperature RCE works

The RCE loop has no fixed-Ts switch, but `&exocol_nml::dz_slab = 1.0e5` holds
Ts at 300.00 K while the atmosphere equilibrates (SBM, `rh_sbm = 0.7`,
`moisture_scheme = 'prognostic'`, `max_model_days = 1500`). Cost: **6–8 min per
column**. The convergence test never fires, because TOA stays out of balance by
design, so a stratospheric-drift criterion is needed. The 1500-day results
below are therefore not certified converged.

Unlike free-Ts RCE, fixed-Ts RCE can in principle cover the whole table
temperature axis. Free-Ts RCE cannot reach runaway or snowball states.

### 2. The RCE stratosphere depends strongly on composition

| CH4 | cold point | mass-weighted T above cold point | RCE tropospheric RH (median) |
|---|---|---|---|
| 0 | 145.5 K at 1.9 hPa | 146.5 K | ~1.0 (0.71–1.00) |
| 0.1 % | 164.4 K at 1.8 hPa | 164.5 K | ~0.9 |
| 1 % | 219.7 K at 100 hPa (tropopause); coldest 197 K at 1.1 hPa | 214.0 K, **inverted** | ~0.7 |

The fixed 200 K is wrong **in both directions**. CO2-rich, O3-free columns
have a much colder, higher stratosphere, and CH4 shortwave absorption makes a
warm inversion. RCE humidity is also not a uniform 0.7. (These RH values come
from a rough Clausius-Clapeyron diagnostic; recheck with ExoColumn's own
`esat`.)

### 3. Flux impact (prescribed flux_only columns vs the RCE columns)

| column | OLR [W/m²] | albedo |
|---|---|---|
| CH4 0: current assumption (RH 1, 200 K) | 232.47 | 0.2288 |
| CH4 0: RH 0.7, 200 K | 241.28 | 0.2338 |
| CH4 0: RH 1, 146.5 K | 222.67 | 0.2288 |
| CH4 0: **RH 0.7, 146.5 K** | **231.46** | 0.2338 |
| CH4 0: full RCE column | 230.99 | 0.2349 |
| CH4 1 %: current assumption | 205.24 | 0.2083 |
| CH4 1 %: RH 0.7, 200 K | 211.70 | 0.2111 |
| CH4 1 %: RH 0.7, 214 K (mass-weighted) | 218.17 | 0.2111 |
| CH4 1 %: RH 0.7, 197 K (coldest) | 210.56 | 0.2111 |
| CH4 1 %: full RCE column | 211.39 | 0.2113 |

Conclusions:

- **Without CH4 the two changes nearly cancel** (+8.8 and −9.8 W/m²). RH 0.7
  with the RCE stratospheric T reproduces the RCE column to 0.5 W/m² OLR and
  0.001 albedo.
- **With CH4 they do not**: the net is ~+6 W/m² more OLR than the current
  table, an anti-greenhouse from stratospheric warming (cf. Haqq-Misra et al.
  2008).
- **For an inverted stratosphere the mass-weighted mean is the wrong
  effective temperature** (+6.8 W/m²). A single isothermal `t_strato` must be
  *fitted* to the RCE fluxes, or the RCE profile shape must be carried.

## RH = 0.8 tables for the moist EBM (option A, started 2026-09-15)

HEXTOR now diffuses moist static energy (`&ebm::moistdiff`, `rhmoist = 0.8`;
Jacob adopted RH = 0.8 after an RH 0.4–0.9 calibration sweep in ebmEarlyEarth).
That settles **D1**: the table RH is the EBM's RH.  The tables were still RH = 1,
so Jacob chose **option A**: build RH = 0.8 tables now with the existing
`min(200 K, Ts)` stratosphere, and do the RCE/zenith stratosphere as a later,
second rebuild.

**D5 resolved: ExoColumn `&exocol_nml::variable_ps_rh`** (uncommitted, in
`src/exocol_config.F90` and `src/exocol_coldstart.F90`; built at PVER = 200 as
`run/exocol_sweepts200_psrh.exe`, the builder's new default executable).
- `variable_ps` alone puts esat(Ts) on top of p_dry, while the column holds
  only rh·qsat. The pressure coordinate then carries (1−RH)·esat(Ts) of extra
  dry gas, up to 7× the real dry gas in a 0.05 bar column at 390 K.
- The naive fix, surface vapour = RH·esat, is wrong. The non-ideal adiabat takes
  Pn = P − esat(T), clamped at 1e-6·P, so it would silently switch to a steam
  lapse rate in hot, thin columns.
- The switch integrates the saturated adiabat as before, then lowers each
  interface by (1−RH)·[esat(T) − esat(t_strato)]:
  - T(Pn) is unchanged;
  - dry mass = p_dry and ps = p_dry + RH·esat(Ts), within 0.03 Pa;
  - the stratosphere and model top stay put;
  - water is w = ε·RH·esat/(p − RH·esat).
- Validation:
  - switch off: bit-identical to `exocol_sweepts200.exe` (all records, RH 1 and
    0.8); switch on at RH 1: identical;
  - ps − p_dry = 0.8·esat exactly at 0.05–20 bar and 250–420 K, and surface
    humidity is unchanged;
  - smooth to 420 K;
  - versus variable_ps alone: ≤ 0.4 W/m² OLR up to 330 K at ≥ 1.5 bar; +1.8 at
    1.5 bar/360 K; +5–8.5 at 1.5 bar/390–420 K; up to +14 W/m² in 0.05–0.1 bar
    hot columns. Albedo changes by ≤ 0.0045.
- Size of the RH effect itself (RH 0.8 vs 1, the old convention): +3–8 W/m² OLR
  at 300–390 K; albedo up to +0.017 over bright surfaces.

**Builder:** `tools/make_radiation_table.py --rh 0.8` writes `variable_ps_rh`
when RH < 1, tags cache names `_rh0.8000` (so RH = 1 caches cannot be reused),
and records `rh`, `water` and `exocolumn_exe` in the HDF5 attributes.

**Pilot result (2026-09-15, `radiation_N2_CO2_Sun_p_rh0.8.h5`, 165 columns, 0
failed; the build took 2 h 56 min):**
- `check_table_sanity`: all hard checks OK. The soft exceptions are fewer than
  in the RH 1 table: OLR vs p 526/3850 against 810/4900; OLR vs CO2 6/3750
  against 26/4875.
- The runaway plateau rises from ~292 to ~308 W/m².
- Moist RH 0.8 Earth calibration at CH4 = 0 (ebmEarlyEarth
  `calibrate_present_earth_ch4.py`):

  | table | d0 | fcloud | cloudir | RMSE | ΔT 2×CO2 | ΔT 4×CO2 | polar amp. |
  |---|---|---|---|---|---|---|---|
  | RH 1 | 0.2938 | 0.4451 | 13.73 | 1.96 K | 6.33 K | 15.0 K | 1.85 |
  | RH 0.8 | 0.2818 | 0.4315 | 19.66 | 1.99 K | 5.05 K | 11.2 K | 1.86 |

  `cloudir` rises 5.9 W/m² toward the observed LW cloud effect, as predicted,
  and sensitivity falls 20–25 %.
- CH4 RH 0.8 build launched 2026-09-15 14:22
  (`model/radiation/build_ch4_sun_rh0.8.log`, ~41 h). **Paused 15:23 at Jacob's
  request** with 56 of 2145 columns cached. Resume by re-running the same
  command (`--ch4-axis --star G2V_SUN_n68.nc --rh 0.8 --out
  model/radiation/radiation_N2_CO2_CH4_Sun_p_rh0.8.h5 --workers 4`).

**M-dwarf (SAMOSA) CH4-free RH 0.8 table, launched 2026-09-15 15:27:**
`radiation_N2_CO2_3000K_p_rh0.8.h5`, log `build_3000K_rh0.8.log`, 210 columns,
~5 h.
- Star: `blackbody_3000K_n68.nc`, per Jacob.
- Axes are exactly those of `radiation_N2_CO2_3000K_p.h5`: 15 p, 14 CO2, 31 T
  from 180 to 620 K, passed via `--pressures/--fco2/--temps`. That table was
  merged from a 180–420 K and a 440–620 K build; one build over all 31 T is
  equivalent.
- Checks before launch:
  - at RH 1 the current builder and `_psrh` executable reproduce the existing
    table exactly (0.05 and 1 bar; 300, 420, 440 and 620 K);
  - at RH 0.8, ps − p_dry = 0.8·esat to within 0.05 bar up to 127 bar (620 K),
    and OLR flattens smoothly onto a ~303 W/m² plateau (RH 1: ~289).
- **Done 2026-09-15 20:23** (4 h 56 min, 210/210 columns, 0 failed). Axes are
  bit-identical to the RH 1 table.
  - `check_table_sanity`: all hard checks OK; soft exceptions match the RH 1
    table (OLR vs CO2 31/6045 against 30; OLR vs p 1265/6076 against 1435).
  - Runaway plateau 307 W/m² against 292, and it starts later in Ts (1 bar:
    350 K against 340; 10 bar: 440 K against 420).
  - RH 1 → 0.8: OLR +1 W/m² at 250 K, +5–8 at 280–300 K, +10–15 above 320 K;
    planetary albedo +0.005 to +0.010 (less near-IR absorption).

**Plan:**
1. Pilot `radiation_N2_CO2_Sun_p_rh0.8.h5` (CH4-free, ~3 h), then
   `check_table_sanity`.  **Done, see above.**
2. Moist Earth calibration at CH4 = 0 on it, against the same calibration on the
   RH 1 CH4-free table. Baseline on RH 1: `d0 = 0.2938, fcloud = 0.4451,
   cloudir = 13.73`, RMSE 1.96 K. Prediction: `cloudir` rises.
3. Then build `radiation_N2_CO2_CH4_Sun_p_rh0.8.h5` (~41 h), recalibrate with
   730 ppb CH4, and re-run the ebmEarlyEarth sweep.

## Prior method: Williams & Kasting (1997) / Haqq-Misra et al. (2016)

Read 2026-09-15 at Jacob's request. Sources: Haqq-Misra et al. 2016, ApJ 827:120,
Appendix (`~/Documents/library/haqqmisra_etal_2016_limitcycles.pdf`), which
cites Williams & Kasting 1997, Icarus 129:254, Appendix B/F
(`williams_kasting1997.pdf`). Those fits came from the Kasting/Kopparapu RC
model at 1 bar N2, not ExoColumn, and were polynomials rather than lookup
tables.

- **Humidity: the old sweeps also used RH = 1.** WK97 says the troposphere was
  "fully saturated with water vapor", and that a realistic RH "would have
  complicated the model without changing any of our basic conclusions". The
  HM16 appendix does not mention humidity, and its RC model (Kopparapu et al.
  2013) assumes a fully saturated troposphere. So there is **no precedent for
  D1**: sub-saturation is new.
- **Stratosphere: isothermal, at the gray skin temperature** (Kasting 1991
  Eq. 1; HM16 Eq. 5; WK97 A11),
  `T_strat = 2^(-1/4) [S (1 - α) / 4σ]^(1/4)`. It is valid only for the global
  mean (z = 60°), so both papers scaled it with zenith angle as
  `T_strat(z) = T_strat(60°) [F_s(z) / F_s(60°)]^(1/4)` (HM16 Eq. 6; WK97 A12),
  where F_s is the absorbed fraction of incident solar flux from the RC model.
  WK97 also fitted `T_strat(pCO2, T)` at z = 60° (A19), a precedent for a
  composition- and Ts-dependent `t_strato`. The reference S used in Eq. 5 is not
  stated.
- **Why they needed the zenith scaling.** Applying A11 with the local α made TOA
  albedo *decrease* with zenith angle (the Caldeira & Kasting 1992 error). Our
  builder does not tie `t_strato` to the local α, and computes palb over
  zenith × albedo on one profile, so that artifact cannot arise here.

### Does the gray skin temperature fix the bias? No.

For the feasibility columns above (S/S0 = 0.75, `coszrs` = 0.5):

| column | gray T_skin | RCE stratosphere |
|---|---|---|
| CH4 0 (α = 0.2288) | 204.1 K | 146.5 K |
| CH4 1 % (α = 0.2083) | 205.4 K | inverted: 214 K mass-weighted, 197 K coldest |
| present Earth, S0, α = 0.30 | 214.1 K | — |

The gray formula knows nothing about composition except through α, so it
reproduces the current fixed 200 K to within ~5 K and misses the ~60 K
CO2/CH4 spread found in the RCE tests. It is the "old table" assumption in
different words.

### What carries over

- **Parameterize an isothermal `t_strato`** by composition and Ts, as WK97 A19
  did, but fit it to fixed-Ts RCE fluxes, not the gray formula. This is D2's
  fitted-isothermal option, and needs no ExoColumn change.
- **Keep a latitude-dependent stratosphere (Jacob, 2026-09-15)**, because the
  instellation reaching each latitude depends on zenith angle. My first draft
  recommended dropping the scaling; that is overruled. Open design questions:
  - *Zenith angle alone does not fix the local instellation.* HEXTOR's
    diurnal-mean insolation is `s(k) = q·(h/π)·mu(k)` (`driver.f:1386`), so day
    length `h/π` and `q` (S/S0, orbit) enter as well as `mu`. Example, solstice
    at δ = 23.5°: the pole gets `mu` = 0.40 and `s` = 0.40 q, while the equator
    gets `mu` = 0.58 and `s` = 0.29 q. Scaling by zenith angle makes the summer
    pole's stratosphere colder than the equator's, but absorbed flux makes it
    warmer. The difference grows with obliquity, which the study sweeps to 90°.
  - *WK97/HM16 Eq. 6 scales by F_s, "the absorbed fraction of incident solar
    flux".* Read literally that is (1 − α(z)), a weak dependence that carries
    no instellation magnitude. Check what Kasting's RC model actually varied
    with z before copying it.
  - *Options.* (a) A zenith axis on `/olr`, with one profile per zenith node:
    12× the profiles. This carries zenith angle but not day length or S.
    (b) A `t_strato` axis on `/olr` and `/palb`, with HEXTOR computing each
    belt's `t_strato` from its absorbed flux `s(k)(1 − α)` and composition,
    using a small RCE-fitted relation. This carries zenith angle, day length
    and S, so it also settles D3. It costs ~N_tstrat × the build; check whether
    ExoColumn's paired `sweep_ts` / `sweep_t_strato` lists accept repeated Ts.
  - *Polar night*: `mu = 0`, `s = 0`, so the skin formula gives 0 K. That
    needs a floor or an LW-only RCE value. The table's smallest `mu` node is
    0.05.
  - *Driver ordering*: `getOLR` (`driver.f:759`) is called before `mu(k)` is
    computed (`driver.f:828–837`), so any local stratosphere needs a reorder.
- **D3 is the same open question they had.** Eq. 5 scales as S^(1/4), so gray
  T_skin at S/S0 = 0.75 is 7 % below S0 (~15 K). The RCE sensitivity to S has
  not been measured yet.
- **Cheap baseline.** The gray T_skin column costs nothing extra, so include it
  in the step-4 flux test as a "2016 method" reference.

## Decisions to make before building

**D1. Tropospheric RH.** Options:
- constant 0.7, equal to ExoColumn's `rh_sbm`, the SBM target (simple, and
  self-consistent with the RCE convection scheme);
- a Manabe–Wetherald profile (RH_s = 0.77, as in CLIMA);
- RH diagnosed from the RCE runs, which test 2 shows varies with composition.

Recommendation: diagnose RH properly from the RCE runs first, then decide
whether a constant is adequate by the flux test in step 4 below.

**D2. Stratosphere representation.** Options:
- an effective isothermal `t_strato(P, fCO2, fCH4, Ts[, star])`, fitted so the
  prescribed column matches RCE OLR and albedo. This needs **no ExoColumn
  change**: `sweep_t_strato` already takes one value per Ts;
- a spliced profile: moist adiabat up to the RCE tropopause, then the RCE T(p)
  above. This needs an ExoColumn change, and has to handle the cold point and
  inversion separately.

Recommendation: try the fitted isothermal value first. Fall back to splicing
only where it misses the budget (likely high CH4 / M-dwarf).

**D3. RCE forcing.** The stratosphere depends on absorbed shortwave, but the
table has no instellation axis. Choose the instellation, `coszrs` and surface
albedo for the RCE runs, e.g. per-table reference S, global-mean zenith.
Test S/S0 = 0.75 vs 1.0: if the sensitivity is small, one reference is enough;
if not, S must be a build-time choice (an early-Earth table at 0.75).

**D4. Stratospheric humidity closure mismatch.** The cold start sets
stratospheric q = `rh_init × qsat(t_strato, p_tropo)`
(`exocol_coldstart.F90` ~l.575–580), so RH < 1 also dries the stratosphere.
The RCE loop instead floors it at `qsat(T_coldpoint)`
(`apply_stratospheric_coldtrap`). For consistency, add an ExoColumn option
applying `rh_init` to the troposphere only, or using the RCE cold-point q.

**D5. Surface vapour pressure with RH < 1.** `variable_ps` sets total
ps = p_dry + **esat(Ts)** regardless of `rh_init` (ocean-saturated surface),
while the column holds only rh·qsat. The pressure coordinate then contains
vapour mass that is not in the column. This is negligible when esat ≪ p_dry,
but **not** at the hot end or at low pressure (esat(330 K) ≈ 0.17 bar;
esat(420 K) ≈ 4.4 bar; p_dry goes down to 0.05 bar). Decide on a convention,
e.g. ps = p_dry + RH_surface·esat, or keep the saturated ocean surface.

**D6. Edges of the grid.**
- Cold end: CO2 condensation (`co2_condense`) exists only in the cold start,
  not in the RCE loop, so cold CO2-rich RCE columns can supersaturate in CO2.
- Hot end: steam-dominated columns up to 420 K; check that fixed-Ts RCE
  behaves there.
- Ts ≤ t_strato gives an isothermal column.

**D7. Cost / subgrid.** The full grid (15 P × 11 CO2 × 13 CH4 × 25 T = 53,625
columns) at ~7 min each is out of reach. Use a coarse RCE subgrid (e.g.
5 P × 4 CO2 × 5 CH4 × 6 Ts ≈ 600 runs, ~15–20 h on 4 workers) and interpolate
the effective `t_strato` (and RH if needed). Validate on held-out nodes.
Starting from a nearby converged state shortens runs. Define a convergence
criterion on stratospheric T drift.

**D8. Vertical resolution.** The RCE tests used PVER = 70; tables are built at
PVER = 200 (`exocol_sweepts200.exe`). Check that the stratospheric T and the
fitted `t_strato` are resolution-robust.

## Implementation steps

1. **ExoColumn** (ExoRT stays read-only): troposphere-only RH / RCE-consistent
   stratospheric humidity (D4); an explicit fixed-Ts RCE switch (or document
   `dz_slab`) plus a stratospheric-drift stop; optional spliced stratosphere
   (D2).
2. **RCE driver** (new, e.g. `tools/make_strato_table.py`): run the fixed-Ts
   RCE subgrid per host star; diagnose tropospheric RH and the cold point;
   fit the effective `t_strato` by matching prescribed and RCE fluxes; write a
   small table with its interpolation error.
3. **`tools/make_radiation_table.py`**:
   - make `RH` configurable;
   - replace `t_strato_for(ts)` with the RCE lookup (it needs P, fCO2, fCH4
     too);
   - record RH, the stratosphere method and the RCE source in the HDF5
     attributes.
   - **Trap:** cache files are keyed only by (p, CO2, CH4), so a build pointed
     at an existing `*.h5.cache` would silently reuse RH = 1 / 200 K columns
     and report a full cache hit. Use a new `--cache` directory, or put an
     RH/stratosphere tag in `cache_name()`.
4. **Acceptance**:
   - `check_table_reader.py` (after `make clean`) and `check_table_sanity.py`;
   - **new test**: prescribed-column vs full-RCE-column fluxes at held-out
     nodes, within ~1 W/m² OLR and 0.003 albedo.
   - Note the runaway plateau (~292 W/m²) will move with RH.
5. **Build order**:
   - pilot `radiation_N2_CO2_Sun_p.h5` (4 h);
   - then `radiation_N2_CO2_CH4_Sun_p.h5` (41 h);
   - then the 3000 K and 2600 K tables (their stratospheres need their own RCE
     subgrid; the M-dwarf SED heats CH4-bearing stratospheres more).
6. **Downstream** (all calibrations are table-specific):
   - HEXTOR: `tools/calibrate_earth.py` → `namelists/input.nml.earth.ch4`,
     THAI (`calibrate_thai.py`), SAMOSA.
   - ebmEarlyEarth: `calibrate_present_earth_ch4.py`, then re-run
     `ch4_sweep.py` into a new output directory (not a resume).
   - Prediction to check: `cloudir` moves toward the observed 25–30 W/m².

## State of ebmEarlyEarth at the halt

The v3 CH4 sweep (`experiments/early_earth_v3/`) was stopped on 2026-09-15.
aqua warm is complete (27,720 rows); earth warm is partial (11,148 rows); the
rest were not started. Its results (reconciliation needs ≥ 0.1 % CH4, but
almost only in the organic-haze regime) are **provisional**: they were
computed on RH = 1 / 200 K tables. **Do not resume that sweep on the old
tables.** Haze optics are a separate follow-up (with Eric Wolf).
