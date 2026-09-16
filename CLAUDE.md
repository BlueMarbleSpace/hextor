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

**Sub-saturated tables (`--rh`, default 1.0).** Use `--rh 0.8` to match `&ebm::rhmoist` when the EBM diffuses moist static energy.
- With RH < 1 the builder sets ExoColumn's `variable_ps_rh`, so the column carries RH·esat (not esat) on top of `p_dry` and its dry mass equals the table's pressure coordinate. Without it, a hot, thin column holds up to 7× its real dry gas.
- That needs `run/exocol_sweepts200_psrh.exe`, now `DEFAULT_EXE`. With the switch off it is bit-identical to the old `exocol_sweepts200.exe`.
- RH < 1 cache names carry `_rh<RH>`, so an RH = 1 cache is never reused.
- Details and validation are in `notes/subsaturated_rce_tables.md`.

Three validation scripts, each answering a question the table's accuracy depends on:
- `tools/check_table_convergence.py` — vertical resolution and model top. At 10 bar / 360 K, OLR is 255.6 / 253.1 / 252.5 / 252.1 W/m² for 70 / 140 / 200 / 400 layers, so the 200-layer build sits within ~1 W/m² of converged.
- `tools/check_table_interpolation.py` — grid spacing of the zenith and surface-albedo axes. The curvature is concentrated at the limb, which is why `MU_NODES` is spaced geometrically rather than uniformly (max albedo error 0.011 → 0.003).
- `tools/check_table_reader.py` — that `radiation.f90` interpolates a table the way its axes say it should, against an independent NumPy implementation (agreement ~1e-11 including clamping and OLR extrapolation). **Run `make clean` in `model/radiation/` first.** `make` only relinks `tabletest` when `radiation.f90` is newer than the object file, so a stale probe can validate against a stale checker and report OK while testing nothing — that happened once here, hiding a 0.4% disagreement above the table ceiling.
- `tools/check_table_sanity.py` — physical checks on a finished table plus an inspection figure. Note the OLR-versus-temperature check applies only *below* the runaway plateau: the table saturates near 292 W/m² (Simpson-Nakajima), above which OLR wobbles ~1 W/m² either way. That plateau is why `n_eff_upper` is floored at zero in `radiation.f90` — half the columns otherwise fit a negative exponent and extrapolate OLR *downward* as the model heats.

**Re-calibrating after a table change** (`tools/calibrate_thai.py`): the published `(d0, cloudir) = (3.10, −35.0)` was fitted against the old 1 bar 2600 K table, so part of that −35 W/m² compensates for that table rather than for clouds. Re-running the published procedure against a BT-Settl 2600 K table built with ExoRT gives `(3.33, −45.0)` and fits the THAI ensemble day-night contrast six times better (0.7 K error against 4.3 K). Consistency with a published calibration means re-running its *procedure*, not transplanting its constants onto different radiative transfer. `--moistdiff` / `--rhmoist` calibrate with moist static energy diffusion on, and then `dt` is set per run by a fixed point on the substellar beta rather than left at 1350 s.

**SAMOSA intercomparison**: `namelists/input.nml.samosa` is the protocol template (3000 K blackbody — an M-dwarf table, not the Sun; the current submission uses `radiation_N2_CO2_3000K_p_rh0.8.h5`, see "SAMOSA on the RH 0.8 tables" below; 15 d synchronous, aquaplanet, CO2 at a fixed 400 µbar partial pressure per the 2024 erratum — not 400 ppm, so `run_samosa.py` sets `pg0 = pN2 + 4e-4` and `fco2 = 4e-4/pg0` per case) and `tools/run_samosa.py` runs the case sequences, each in its own scratch directory:

```bash
python tools/run_samosa.py --sequence all16 --d0-ref 3.10 --init both
```

It writes `samosa/samosa_summary.csv` plus a per-case zonal profile in the substellar-angle coordinate, and labels each outcome `equilibrium` / `drifting` / `runaway` — HEXTOR's own `converged` flag halts on the year-to-year change in global OLR, which is a false positive on the runaway plateau where OLR is insensitive to surface temperature. **The drift used for that label is `dT_rate10`**, ten times the median of the last three annual increments, *not* the change over the last ten years: these cases approach equilibrium exponentially and the flux halt fires around year 20, so a ten-year window still spans the transient and labelled converged runs `drifting` (Case 11 is settled to 0.001 K/yr by year 19 and scored −0.51 K/decade, which dropped it from the submission). Verified against runs with the halt disabled — the halted value holds to 0.001 K out to 2000 years, so only the label was wrong. Separation is now 0.02 K/decade for equilibria against 41+ for runaways.

**`--years` does not bound a run.** With `seasons = .true.` the halt at `tend` is commented out (`driver.f:593–600`), so reaching `tend` only freezes the orbital clock and the integration continues to flux convergence or the hard-coded `niter = 5000`. Harmless for SAMOSA (synchronous, zero obliquity and eccentricity, so insolation is constant within a year) but a live trap for any seasonal configuration.

**Heat transport** is selected with `--transport`, and the three options are calibrated (`tools/calibrate_thai.py`) to the same effective D ≈ 3.3 at THAI Hab1, so they are indistinguishable there and differ only in how that is carried across the SAMOSA parameter space:

| `--transport` | D | at 0.1 → 10 bar (15 d rotation) |
|---|---|---|
| `constant` | `d0` | 3.33 everywhere |
| `psqrt` | `d0 · √p` (applied by the runner, driver scaling off) | 1.05 → 10.5 |
| `perbar` | `d0 · p · comp` | 0.33 → 33 |
| `diffadj` | `d0 · p · comp · (rot0/rot)²` | 2.0 → 197 |

(The D values in that column are for the old `d0` = 3.33; the current submission uses 2.46.) **Calibrate with the same `--transport` the run will use**: `calibrate_thai.py` now takes the same names. Before that it had only `--diffadj`, and because the Hab1 template carried no `diffadj_rot` it picked up the driver default and calibrated the *full* rotation-scaled transport — a factor 40.9 at Hab1 where `perbar` means 1.10. `psqrt` needs no recalibration, since √p = 1 at the 1 bar reference.

**Measured against ExoCAM** (RH 0.8, moist, `samosa/samosa_summary_transport_*.csv`): global-mean RMSE is 24.1 / 24.0 / 25.6 K for constant / √p / p — transport barely moves the mean. Contrast RMSE is 33.1 / **28.3** / 33.4 K. `perbar` is spectacular at the extremes (Case 16 contrast 45 → 5 K against ExoCAM's 7; Case 15 72 → 129 against 128) but breaks the mid-pressure cases and **destroys Case 11's climate** (night side 92 K at D = 0.25), even with latent transport on. √p keeps all nine equilibria and is a modest real improvement, but was chosen *after* seeing ExoCAM, so adopting it would be tuning to one member of the ensemble — the submission stays at `constant`. No transport choice closes the cold-case gap: Case 1's substellar point stays 48–58 K too cold under all three, because what is missing there is substellar greenhouse and cloud warming, not a slower drain of heat.

`diffadj_rot` (new, `&ebm`, default `.true.` = published behaviour) gates just the rotation factor, so the pressure and composition scaling can be kept without it. The rotation term is 37× for TRAPPIST-1e but 225× for a 15 d rotator, so a `d0` calibrated for one rotator moves the transport by that ratio when carried to another. `tools/samosa_sensitivity.py` sweeps all three against a range of `cloudir`.

**Cloud correction** is scaled the same way with `--cloud-scaling`: `constant` (default) writes `cloudir` unchanged into every case; `instellation` writes `cloudir × S/900`, anchored at the THAI Hab1 instellation so the calibration point is untouched (Cases 11 and 14, at 900 W/m², reproduce bit-identically). `cloudir` stands in for the shortwave cooling of the clouds the clear-sky table lacks, and a fixed W/m² offset is not shortwave-like: −45 W/m² is 45% of the absorbed flux at S = 500 and 15% at 1200, which drove a bias monotonic in S (−48 K against ExoColumn at 400 W/m², +15 K at 1200). Scaled, it is a fixed +0.20 planetary-albedo increment, still uniform over both hemispheres. The summary CSV records `cloudir` (per case), `cloudir_ref` and `cloud_scaling`, and `samosa_submit.py` describes the rule in the output header. Consequence at the warm end: the extra OLR grows with S, so Case 16 (1400 W/m², 10 bar) leaves the runaway plateau (465 → 365 K) and Case 12 (1500 W/m², 2.98 bar) now "equilibrates" at 412 K with clear-sky OLR 290 W/m², i.e. *on* the ~292 W/m² plateau, where the temperature is not well constrained. (Those two numbers are RH 1 results; on the RH 0.8 table the plateau is 307 W/m² and Case 12 clears it — see "SAMOSA on the RH 0.8 tables" below.)

**Timestep and diffusion stability:** HEXTOR integrates the diffusion term explicitly on 18 belts in x = sin(lat), so it needs `D·dt/(C·dx²)` below about ½. At the published `dt = 1350 s` that is 0.9 for `perbar` at 10 bar and 5.4 for `diffadj`, and exceeding it does not raise an error — the model integrates quietly to NaN over thousands of years. The same case gives NaN at 1350 s, 289.34 K at 400 s and 289.35 K at 135 s. `run_samosa.py` therefore sets `dt` per case from that limit, capped at 1350 s. **Any run with D much above ~3 needs the timestep reduced.** With `moistdiff` the limit tightens by a further factor `beta = 1 + (L/cp)·dq/dT` (`tools/moist_stability.py`, mirroring `driver.f`'s own `moistprops`/`esatw`/`qmoist`), which `run_samosa.py` and `calibrate_thai.py` find by a **fixed point**: run at a guessed substellar temperature, then redo at the beta the run actually produced, capped at 420 K. A first pass that goes non-finite is retried once at that cap, because an unstable timestep and a runaway are indistinguishable in the output. beta is largest where the atmosphere is thin (q ~ esat/p_dry): SAMOSA Case 11, at 0.1 bar, gets `dt = 415 s` from beta = 19.4. **The beta correction is conservative** — it charges the full beta against a coefficient that diffuses h/cp while the prognostic update is on T. At the SAMOSA D of 2.46 the cases are in fact stable at 1350 s, agreeing with the reduced-`dt` answer to ≤0.011 K (and to ≤0.014 K at a quarter of the chosen `dt`), so it is insurance rather than necessity there; it matters under `perbar`/`diffadj`, where D reaches 33 and 197.

**A blow-up does not fail the driver.** It produces NaN and keeps integrating NaN to `niter = 5000` — two hours per case — so `run_samosa.py`'s `_integrate()` polls `out/tempseries.out` and kills the run the moment it goes non-finite (`RC_NONFINITE`), and a failed pass now feeds the dt retry instead of aborting it. Without that, `perbar` Case 16 went NaN at year 21 (stable for its equilibrium, not for its 342 K transient), timed out, and was reported as a runaway; with it, the case equilibrates at 373 K in 15 s. `classify()` now labels these `unstable` or `no result (timeout)`, never `runaway`.

**Moist static energy diffusion** (`&ebm::moistdiff`, default `.false.`; `rhmoist`, default 0.8). Diffuses h/cp = T + (L/cp)·q instead of T (Frierson et al. 2007), so latent heat transport follows Clausius–Clapeyron. q is the surface specific humidity at `rhmoist`, over dry pressure `pg0`, with the vapour sitting on top of the dry air as in the tables (`variable_ps`). esat mirrors ExoColumn's steam mode: Wagner & Pruss over liquid, Clausius–Clapeyron over ice. `diffadj` still scales D. With the flag off, output is bit-identical to the pre-change binary. Things to know:
- h/cp responds to T faster than T by β = 1 + (L/cp)·dq/dT: 3.5 at 300 K, 7 at 320 K, 13 at 343 K. **Divide the stable dt by β at the warmest belt.**
- The moist run flickers ±0.001 K on a 3-year cycle, so `fluxcnvg = 1e-4` never fires. Use 1e-3.
- With `fillet`, `out/transport.out` gives the final-orbit annual-mean northward transport and its latent part at the belt edges, in PW for Earth's radius.
- Earth calibration: ebmEarlyEarth `calibrate_present_earth_ch4.py --moist-rh 0.8`.
- `ph2` is now initialised to 0. `radparam = 3` without the H2 cycle never assigned it, and D (via `avemol`/`hcp`) had only read zero by luck of the memory layout; adding variables exposed it as a year-2 blow-up.

**SAMOSA on the RH 0.8 tables (2026-09-15, current submission).** Redone on `radiation_N2_CO2_3000K_p_rh0.8.h5` with `moistdiff = .true., rhmoist = 0.8`. Because the calibration is table-specific, a matching 2600 K RH 0.8 table was built for it (`radiation_N2_CO2_2600K_p_rh0.8.h5`, 70 columns, 1 h 14 min, sanity OK, plateau 306.8 W/m²) and the THAI Hab1 procedure re-run on it:

| table | diffusion | d0 | cloudir | T_anti | T_sub | contrast |
|---|---|---|---|---|---|---|
| RH 1 | dry | 3.33 | −45.00 | 203.53 | 299.62 | 96.08 |
| RH 0.8 | dry | 3.26 | −42.57 | 203.69 | 300.67 | 96.97 |
| RH 1 | moist 0.8 | 2.50 | −46.98 | 197.43 | 294.03 | 96.59 |
| **RH 0.8** | **moist 0.8** | **2.46** | **−45.09** | **197.31** | **294.19** | **96.88** |
| THAI 4-GCM | | | | 194.7 | 291.6 | 96.9 |

`cloudir` landing within 0.1 of the published −45 is a **coincidence of two cancelling effects**, not a reason to have skipped the recalibration: RH 0.8 raises clear-sky OLR so less correction is needed (−45.0 → −42.6), and moist diffusion pushes the other way (−45.0 → −47.0). `d0` genuinely moves, 3.33 → 2.46, because latent transport supplies part of what `d0` was carrying. The calibrated configuration now reproduces the *profile*, not only the two numbers fitted: both extremes are 2.6 K too warm against 8–9 K before, and the contrast error is 0.02 K against 0.8 K — an independent check, since the procedure fits only the global mean and the contrast.

Production run (`--sequence all16 --init both --transport constant --d0-ref 2.46 --cloudir -45.09 --moistdiff`): **the same nine cases equilibrate** as the published RH 1 / dry submission (1, 4, 8, 9, 10, 11, 14, 15, 16); the other seven run away. **No bistability** — warm (300 K) and cold (233 K) starts agree to 0.003 K in all nine. **Case 16 is now inside the table**, equilibrating at 376 K where the RH 1 run sat at 465 K, above the 420 K ceiling and therefore extrapolated. Both changes cool and are roughly additive; at the published constants Case 4 runs 308.3 → 300.9 (table) → 294.4 (moist) → 289.8 K (both). The four decomposition runs are in `samosa/samosa_summary_decomp_*.csv`, and the RH 1 / dry one reproduces the committed result to 0.000 K. With `--cloud-scaling instellation` eleven cases equilibrate instead of nine (`samosa_summary_cloudscale_instellation.csv`); Cases 5 and 12 then sit at a clear-sky OLR of ~287 W/m², ~20 W/m² below the RH 0.8 plateau, so unlike the RH 1 result they are genuinely constrained. Full account in `notes/samosa_rh0.8_redo.md`.

**Comparing against the other SAMOSA submissions**: `tools/samosa_compare.py` scores any set of `run_samosa.py` output directories against the other models, reading `/models/data/samosa/`:

```bash
python tools/samosa_compare.py "published=OLD" "D const=samosa"
```

Two references, answering different questions. **ExoCAM** is a 3-D GCM with clouds — the closest thing to a reference answer, and the only submission that constrains the day-night contrast (global means from its `analysis.py` digest, extrema from the `TS` maps). **ExoColumn** is a 1-D RCE column using the *same* radiation core HEXTOR's tables are built from (ExoRT n68equiv), so it says nothing about contrast but isolates everything *except* the spectroscopy: a HEXTOR–ExoColumn gap is transport, clouds or surface albedo, never the line lists. Neither is a calibration target — `(d0, cloudir)` are fitted at THAI Hab1 around a 2600 K star, so every number is out of sample. Current scores in `notes/samosa_rh0.8_redo.md`; the headline is RMSE 45.2 → 24.1 K against ExoCAM on the global mean and 76.0 → 33.1 K on the contrast, with Cases 4 and 9 landing within 2 K of ExoColumn.

**Where HEXTOR still disagrees, and why it is structural.** The residual is the cold, glaciated end, and it sits in the *substellar* temperature, not the night side: on Cases 1 and 15 the antistellar temperature is right to 2 K while the substellar point is 53–58 K too cold. ExoCAM keeps a warm substellar pool at 269–286 K under a fully glaciated surface (ICEFRAC 1.0, CLDTOT 0.69–0.79); a clear-sky column at 0.4–0.7 bar and ~210 K has almost no greenhouse to do that with, and the uniform `cloudir` costs a further ~20 K where the absorbed flux is only ~100–125 W/m². Related: on the *warm* cases the energy budget is right for the wrong reason — planetary albedo 0.04 against ExoCAM's 0.31, OLR 288 against 208, i.e. ~80 W/m² too much shortwave absorbed and ~80 W/m² too much longwave emitted. That is what `cloudir` is: a longwave stand-in for a missing cloud shortwave albedo. Surface temperature survives it because the calibration ties the two together; the fluxes do not, and the protocol asks for them.

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
