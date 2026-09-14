# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Is

**HEXTOR** (Habitable EBM for eXoplaneT ObseRvations) is a Fortran-based latitudinal energy balance model (EBM) for simulating seasonal and latitudinal temperature evolution on Earth-like exoplanets. It supports varying atmospheric compositions (N2, O2, CO2, H2O), orbital parameters, obliquity, eccentricity, and different host stars.

## Build and Run

Before building, load the Intel compiler environment and set machine-specific paths:

```bash
# 1. Load the Intel oneAPI compiler environment (required — driver links against libimf)
source /opt/intel/oneapi/setvars.sh

# 2. Set FC (Fortran compiler) and WDIR in model/Makefile
#    FC is typically an HDF5-aware wrapper: h5fc (wrapping ifort or gfortran)
#    WDIR must be the absolute path to model/

# 3. Set wdir in runEBM.sh to match your system

# 4. Copy a namelist template and edit parameters
cp namelists/input.nml.earth.aqua.23 input.nml
# Edit input.nml as needed

# 5. Build and run
./runEBM.sh
```

If `source /opt/intel/oneapi/setvars.sh` is not run first, `./driver` will fail silently with a missing `libimf.so` error and `runEBM.sh` (a csh script with no error-halting) will appear to succeed while returning stale output files.

**Bifurcation analysis** (sweeps solar constant across multiple runs):
```bash
./bifurcation.sh
```

**Outputs** are written to `model/out/`:
- `tempseries.out` — annual time series: year, ann_tempave (K), pg0 (bar), pco2 (bar), pco2soil (bar), gammaout, q (W/m²), d
- `zonal.out` — per-belt zonal statistics (final year): lat, Tave, Tmin, dec@Tmin, Tmax, dec@Tmax, albedo, OLR (W/m²), ASR (W/m²)
- `geog.out` — per-belt geography: lat, ocean fraction (focean)

**Plotting**: NCL scripts in `plots/` (e.g., `plotTempSeries.ncl`, `plotBistability.ncl`). A Python 4-panel summary script is at `plots/summary_plot.py` — panels are: surface temperature (annual mean + seasonal envelope), planetary albedo, energy balance (ASR and OLR vs latitude), and land fraction. Saves PNG + EPS. Run as `python plots/summary_plot.py` (reads `model/out/`) or `python plots/summary_plot.py <subdir>` (reads `/models/hextor/experiments/<subdir>/` and saves output there).

`tools/calibrate_earth.py` calibrates pre-industrial Earth clouds (`fcloud` → global albedo 0.30, `cloudir` → 288 K) against a pressure-resolved table; ~1 s per run, resume-safe, `--table/--fch4` selectable.

A radiation module comparison script is at `model/radiation/compare_radiation.py` — run as `python compare_radiation.py out_old.txt out_new.txt [figure.png]`. Prints statistics and saves a 3-panel publication figure (PNG + EPS) showing mean |ΔOLR| and |ΔPALB| over CO₂ × temperature, and |ΔPALB| over zenith × surface albedo.

**THAI calibration**: `calibrate_thai_hab1.py` (repo root) does a 2D sweep over (d0, cloudir) to match THAI Hab1 GCM ensemble mean targets; reads T_global from `tempseries.out` and T_min/T_max from the ZONAL STATISTICS block of `model.out` (not `zonal.out`, which is only written by `runEBM.sh`). `plots/thai_hab1_comparison.py` generates a 3-panel zonal comparison figure vs GCM reference values — run as `python plots/thai_hab1_comparison.py [subdir]`.

**THAI Hab1 instellation sweep**: `thai_instellation_sweep.py` (repo root) sweeps `relsolcon` (S/S₀, 0.60–1.42) × 3 initial temperatures (233/273/300 K) using `namelists/input.nml.thai.hab1.calibrated`. Writes `thai_hab1_instellation.csv`. Resume-safe. `plots/thai_iceline_figure.py` reads that CSV and saves `thai_hab1_iceline.png/.pdf` — ice line longitude from substellar vs S/S₀. Ice line longitude (0° = substellar, 90° = terminator, 180° = antistellar): `90° − icelineN` when the ice line is on the dayside (icelineN < 90°); `90° − icelineS` when icelineN = 90 (no dayside ice, use nightside ice line; icelineS is negative so lon > 90°). Sentinel icelineN=90, icelineS=−90 → 180° (truly ice-free); icelineN=0, icelineS=0 → 0° (ice-ball).

**Building a radiation table** (`tools/`, all using ExoColumn at `/models/ExoColumn`):

```bash
python tools/make_radiation_table.py --dry-run            # grid size + cost estimate
python tools/make_radiation_table.py --star blackbody_3000K_n68.nc \
    --out model/radiation/radiation_N2_CO2_3000K_p.h5 --workers 4
```

Each grid point is an ExoColumn run in `flux_only` + `sweep_mode` on a *prescribed* profile — a moist adiabat from the tabulated surface temperature, not an RCE solution (solving for Ts is the EBM's job). `variable_ps` puts the water vapour on top of the dry pressure, matching SAMOSA's definition of surface pressure. The script is resume-safe (one cache file per pressure × CO2 × CH4 column) and takes `--pressures/--fco2/--ch4/--temps` overrides for reduced tables. **Cache files are named by the column's physical values, not by axis index** — index names silently survive a regrid and would resume the wrong physics while reporting a full cache hit. `tools/migrate_cache_names.py TABLE.h5 --apply` renames a pre-2026-09-12 cache, taking the axes from the table the cache was built for. **Worker count:** this machine is a 6-core Xeon with hyperthreading; throughput saturates near 38 records/s at 3–4 workers and *degrades* beyond (12 workers measured 9 records/s, worse than one process).

Three validation scripts, each answering a question the table's accuracy depends on:
- `tools/check_table_convergence.py` — vertical resolution and model top. At 10 bar / 360 K, OLR is 255.6 / 253.1 / 252.5 / 252.1 W/m² for 70 / 140 / 200 / 400 layers, so the 200-layer build sits within ~1 W/m² of converged.
- `tools/check_table_interpolation.py` — grid spacing of the zenith and surface-albedo axes. The curvature is concentrated at the limb, which is why `MU_NODES` is spaced geometrically rather than uniformly (max albedo error 0.011 → 0.003).
- `tools/check_table_reader.py` — that `radiation.f90` interpolates a table the way its axes say it should, against an independent NumPy implementation (agreement ~1e-11 including clamping and OLR extrapolation). **Run `make clean` in `model/radiation/` first.** `make` only relinks `tabletest` when `radiation.f90` is newer than the object file, so a stale probe can validate against a stale checker and report OK while testing nothing — that happened once here, hiding a 0.4% disagreement above the table ceiling.
- `tools/check_table_sanity.py` — physical checks on a finished table plus an inspection figure. Note the OLR-versus-temperature check applies only *below* the runaway plateau: the table saturates near 292 W/m² (Simpson-Nakajima), above which OLR wobbles ~1 W/m² either way. That plateau is why `n_eff_upper` is floored at zero in `radiation.f90` — half the columns otherwise fit a negative exponent and extrapolate OLR *downward* as the model heats.

**Re-calibrating after a table change** (`tools/calibrate_thai.py`): the published `(d0, cloudir) = (3.10, −35.0)` was fitted against the old 1 bar 2600 K table, so part of that −35 W/m² compensates for that table rather than for clouds. Re-running the published procedure against a BT-Settl 2600 K table built with ExoRT gives `(3.33, −45.0)` and fits the THAI ensemble day-night contrast six times better (0.7 K error against 4.3 K). Consistency with a published calibration means re-running its *procedure*, not transplanting its constants onto different radiative transfer.

**SAMOSA intercomparison**: `namelists/input.nml.samosa` is the protocol template (3000 K blackbody — the M-dwarf table `radiation_N2_CO2_3000K_p.h5`, not the Sun; 15 d synchronous, aquaplanet, CO2 at a fixed 400 µbar partial pressure per the 2024 erratum — not 400 ppm, so `run_samosa.py` sets `pg0 = pN2 + 4e-4` and `fco2 = 4e-4/pg0` per case) and `tools/run_samosa.py` runs the case sequences, each in its own scratch directory:

```bash
python tools/run_samosa.py --sequence all16 --d0-ref 3.10 --init both
```

It writes `samosa/samosa_summary.csv` plus a per-case zonal profile in the substellar-angle coordinate, and labels each outcome `equilibrium` / `drifting` / `runaway` — HEXTOR's own `converged` flag halts on the year-to-year change in global OLR, which is a false positive on the runaway plateau where OLR is insensitive to surface temperature.

**Heat transport** is selected with `--transport`, and the three options are calibrated (`tools/calibrate_thai.py`) to the same effective D ≈ 3.3 at THAI Hab1, so they are indistinguishable there and differ only in how that is carried across the SAMOSA parameter space:

| `--transport` | D | at 0.1 → 10 bar (15 d rotation) |
|---|---|---|
| `constant` | `d0` | 3.33 everywhere |
| `perbar` | `d0 · p · comp` | 0.33 → 33 |
| `diffadj` | `d0 · p · comp · (rot0/rot)²` | 2.0 → 197 |

`diffadj_rot` (new, `&ebm`, default `.true.` = published behaviour) gates just the rotation factor, so the pressure and composition scaling can be kept without it. The rotation term is 37× for TRAPPIST-1e but 225× for a 15 d rotator, so a `d0` calibrated for one rotator moves the transport by that ratio when carried to another. `tools/samosa_sensitivity.py` sweeps all three against a range of `cloudir`.

**Timestep and diffusion stability:** HEXTOR integrates the diffusion term explicitly on 18 belts in x = sin(lat), so it needs `D·dt/(C·dx²)` below about ½. At the published `dt = 1350 s` that is 0.9 for `perbar` at 10 bar and 5.4 for `diffadj`, and exceeding it does not raise an error — the model integrates quietly to NaN over thousands of years. The same case gives NaN at 1350 s, 289.34 K at 400 s and 289.35 K at 135 s. `run_samosa.py` therefore sets `dt` per case from that limit, capped at 1350 s. **Any run with D much above ~3 needs the timestep reduced.**

The broadband surface albedos are the protocol's two-channel ice/snow values weighted by the fraction of a 3000 K blackbody below 0.7 µm (f_vis = 0.083 → ice 0.21, snow 0.50, against 0.40/0.71 under the Sun).

**CH4 tables** are implemented, and `notes/ch4_lookup_table.md` carries the scoping note plus an addendum on what the measurements changed. CH4 is a table axis held constant within a run (`&ebm::fch4`); there is no CH4 cycle, and adding one would need no table work. Build one with:

```bash
python tools/make_radiation_table.py --ch4-axis --star G2V_SUN_n68.nc \
    --out model/radiation/radiation_N2_CO2_CH4_Sun_p.h5 --workers 4
```

Without `--ch4-axis` (or an explicit `--ch4` list) the output is a CH4-free v2 table exactly as before. `namelists/input.nml.earth.ch4` is a worked example. **`fch4` needs a table with a CH4 axis:** `driver.f` halts rather than silently returning CH4-free fluxes, and also halts if `fch4` is above the axis, since the bracket would otherwise clamp it and model a different atmosphere than the namelist asked for.

The CH4 axis needs **13 nodes** (1 dex below 1e-6, 0.5 dex above, floor 1e-8, ceiling 1e-1). The binding constraint is the **planetary albedo, not the OLR** — CH4 absorbs in the near-infrared, so uniform 0.75 dex leaves 0.006 in albedo against the zenith axis's 0.003 budget. Under a 3000 K blackbody the albedo errors are ~2× larger again, so an M-dwarf CH4 table would need a finer axis than the solar one. The CO2 axis was trimmed from 14 to 11 nodes (0.57 dex, 0.57 W/m² and 0.0023 in albedo) to pay for part of the cost.

**Results above 420 K are not climates.** That is the table's temperature ceiling; beyond it OLR is extrapolated with `n_eff_upper` floored at zero, so a runaway can report "Flux converged" with no drift at 430–630 K (and non-monotonically in the forcing). Treat any T > 420 K as a runaway. On the calibrated Earth this happens from ~1000 ppm CH₄ at S/S₀ = 1.

**Organic haze is the known gap:** ExoRT's `calc_opd_mod.F90` carries no haze optics, so above CH4/CO2 ≈ 0.1 the table misses the anti-greenhouse a real Archean atmosphere would have and is biased warm. The table records this in its own `ch4_caveat` attribute.

There is no traditional test suite; correctness is verified by comparing simulation outputs to known results.

## Architecture

### Core Model: `model/driver.f`

The main Fortran program (`energy_balance_climate_model`) implements a 1D latitudinal EBM across 18 climate belts. The per-timestep loop:

1. Computes solar insolation for each latitude/season using orbital parameters (obliquity, eccentricity, argument of perihelion)
2. Calls the radiation module to get OLR and planetary albedo via lookup tables
3. Solves the energy balance equation for each belt, including thermal diffusion between adjacent belts
4. Updates CO2 partial pressure via carbonate-silicate cycle (outgassing + weathering)
5. Checks CO2 condensation conditions
6. Writes output at configured intervals

All physical parameters are read from Fortran namelists at startup — no recompilation needed to change scenarios.

### Radiation Module: `model/radiation/radiation.f90`

Provides two public subroutines used by `driver.f`:
- `getOLR(pdry, fco2, fch4, tg0, olr)` → outgoing longwave radiation
- `getPALB(pdry, fco2, fch4, tg0, zy, surfalb, palb)` → planetary albedo
- `radiation_nch4()` → CH4 levels in the loaded table, or 0 if it has no CH4 axis
- `radiation_ch4_range(lo, hi)` → the CH4 axis bounds

`pdry` is the **dry** surface pressure in bar (pN2 + pCO2 + pCH4 = HEXTOR's `pg0`); `fco2` and `fch4` are the mixing ratios pCO2/pg0 and pCH4/pg0. Both query precomputed HDF5 lookup tables, and `radiation_init(radfile)` auto-detects the format:

- **v3 (pressure- and CH4-resolved)** — as v2 with a `/ch4` axis vector and one extra dimension on both data arrays, so `/olr` is rank 4 and `/palb` rank 6. Interpolation gains the CH4 axis (16 corners for OLR, 64 for albedo), in log10 like pressure and CO2. The axis has a **nonzero floor** (1e-8) rather than a zero level, so the log interpolation is well posed; `fch4 = 0` clamps to it, which is radiatively indistinguishable from CH4-free (0.08 W/m² in OLR, 1e-4 in albedo).
- **v2 (pressure-resolved)** — rank-3 `/olr` and rank-5 `/palb` datasets plus explicit axis vectors `/pressure`, `/fco2`, `/temperature`, `/zenith`, `/surfalb`. Grid dimensions are read from the file, so the table can be regridded without recompiling. Loaded with a single degenerate CH4 level, so v2 and v3 share one code path and `fch4` is clamped away — verified bit-identical to the pre-CH4 module over 4000 sample points spanning the clamping and extrapolation paths, and on a full pre-industrial Earth run.
- **v1 (legacy, 1 bar)** — the original flat row-per-sample datasets on the fixed 92 × 19 (× 4 × 5) grid. Loaded into the same arrays with a single 1 bar pressure level and a single CH4 level, so there is only one downstream code path and the pressure and CH4 arguments are simply clamped away. Existing namelists and tables keep working unchanged; the rewrite was verified bit-identical to the previous module over 10,976 sample points and on a full pre-industrial Earth run.

**OLR units:** every format stores OLR in mW/m² (W/m² × 1000) — `driver.f` divides by 1000. A v2 or v3 table must follow the same convention. Outside the tabulated temperature range OLR is extrapolated by a power law whose exponent is fitted per (pressure, CO2, CH4) column from the table's own boundary gradient.

The array index order is the reverse of the axis order, so the fastest-varying axis is innermost in Fortran storage; the HDF5 Fortran interface maps that onto datasets written in natural axis order.

`model/radiation/tabletest.f90` (`make tabletest`) is a standalone probe: it reads a table and answers `pdry fco2 fch4 tg0 zenith surfalb` queries on stdin, which is how `tools/check_table_reader.py` checks the Fortran interpolation against an independent implementation. The same query file works for v1, v2 and v3, since a table without a CH4 axis clamps `fch4` away.

### Namelist Configuration

All model parameters are controlled through `input.nml` (5 groups):
- `&ebm` — orbital properties, heat capacity, diffusion coefficient, CO2 levels, `fch4` (CH4 mixing ratio of dry air, constant within a run; needs a `radfile` with a CH4 axis), timestep
- `&radiation` — solar constant, albedo values, cloud effects, host star; `radfile` sets the HDF5 lookup table path
- `&co2cycle` — outgassing rate, weathering parameters
- `&h2cycle` — H2 cycling (partially implemented)
- `&stochastic` — stochastic noise amplitude and seed

The `namelists/` directory contains 30+ pre-configured scenarios (Earth aquaplanet at various obliquities, land surfaces, Mars, habitability/limit-cycle runs, exoplanet variants).

### Machine Configuration

`config/machine.sh` is a symlink to a machine-specific shell script (currently `merlin.sh`) that sources the correct compiler environment (Intel oneAPI). An alternative `discover.sh` exists for NASA Discover cluster. To port to a new machine, add a new config script and update the symlink.

## Key Technical Notes

- The radiation HDF5 lookup tables must be present at the path `model/radiation/` points to for the model to run. In practice the small ones (1.5–6 MB) *are* tracked in git; the CH4-resolved table is ~45 MB and is deliberately left untracked. The per-column build caches (`*.h5.cache/`) and build logs are gitignored.
- Compiler flags include `-parallel` (ifort) — the model supports light multi-threading through the HDF5 layer.
- `model/driver` (the compiled binary) and output files in `model/out/` are gitignored.
- `input.nml` at the repo root is gitignored; `model/input.nml` is the copy used at runtime (written by `runEBM.sh`).
- **Zenith angle fix (driver.f:1130):** `getPALB` receives `zendeg` (zenith angle in degrees). The correct expression is `acos(mu(k))*180./pi`; the earlier form `mu(k)*180/pi` passed cos(z)×(180/π) instead, producing incorrect planetary albedo for `radparam=3`.
- **cloudir behavior:** `cloudir` applies globally to all belts (both dayside and nightside). A former nightside correction that undid `cloudir` on the nightside has been removed. Positive `cloudir` reduces OLR everywhere (warms); negative `cloudir` increases OLR everywhere (cools).
- **OLR extrapolation (radiation.f90):** Outside the table range [190K, 370K], OLR is extrapolated using a data-driven power-law exponent `n_eff` fitted from the boundary gradient at init time (`n_eff_upper ≈ 2.35`, `n_eff_lower ≈ 3.77`). Not a simple T⁴ law.
- **Pre-industrial Earth calibration** (`radparam=3`, `igeog=1`): `fco2=2.8e-4`, `d0=0.58`, `cloudir=3.0`, `diffadj=.false.` converges to T≈288 K. See `namelists/input.nml.earth` as the base template (the older `input.nml.earth.pres.23` now lives under the gitignored `namelists/old/`). Note that template as it stands converges to 285.8 K, not 288 K, against the v1 Sun table. Against the pressure-resolved Sun tables (`radiation_N2_CO2_Sun_p.h5`, `radiation_N2_CO2_CH4_Sun_p.h5`) it is **not calibrated at all**: at S/S₀ = 1 it equilibrates at **367 K** (0.90 → 285.5, 0.94 → 307, 0.98 → 335 K; smooth, no bifurcation). The tables are not at fault — where both see the same surface, v1 and v2 agree to 0.003 in albedo. The v1 calibration leaned on two artefacts of the old table: its surface-albedo axis starts at 0.2, so the 0.06 ocean was clamped to 0.2, and its OLR never saturates (341 W/m² at 320 K against ~285 for v2). With no clouds in the namelist, v1 still gave Earth a global (insolation-weighted) albedo of 0.261 — the clamp was acting as partial cloud cover — while v2 gives the clear-sky 0.152. Pinning ocean albedo at 0.2 on v2 recovers 304 K; the rest is the missing saturation. **Calibrated 2026-09-14** with `tools/calibrate_earth.py`, in `namelists/input.nml.earth.ch4` (use it, not `input.nml.earth`, with any pressure-resolved Sun table): cloud albedo on (`cloudalb=.true.`, `fcloud=0.4427`, the zenith-dependent form), `cloudir=10.28`, `d0=0.58`, at 280 ppm CO₂ and 730 ppb CH₄ → 287.85 K, global albedo 0.300. Two observed targets, two knobs, solved in the well-conditioned order (`cloudir` for T at each `fcloud`, then `fcloud` for albedo), refined by a local plane fit. Properties: snowball below S/S₀ ≈ 0.95, runaway above ≈ 1.05; warm branch for `cloudir` ≈ 0–22 W/m²; sensitivity ~1.65 K per W/m² (~6 K per CO₂ doubling, about twice Earth's); `cloudir` ≈ 10 against a real LW cloud effect of ~+25–30 W/m², consistent with the RH = 1 column over-stating the clear-sky greenhouse. The same clouds on the CH₄-free v2 table give 283.7 K, matching the CH₄ table at `fch4=0` (283.8 K): pre-industrial CH₄ is worth 4 K here. **Also:** the default `fluxcnvg = 0.1` halts runs that are still warming ~0.7 K/yr near the plateau (one stopped at 348 K on its way to 367 K); use ~1e-4 and check the drift. (With `diffadj=.true.`, effective D rises to ~0.639 because the model has no explicit O₂ — the N₂-dominated atmosphere is lighter and has higher Cp than the reference, so `cloudir=2.5` was the prior tuning for that case.)
