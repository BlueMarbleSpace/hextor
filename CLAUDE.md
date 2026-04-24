# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Is

**HEXTOR** (Habitable EBM for eXoplaneT ObseRvations) is a Fortran-based latitudinal energy balance model (EBM) for simulating seasonal and latitudinal temperature evolution on Earth-like exoplanets. It supports varying atmospheric compositions (N2, O2, CO2, H2O), orbital parameters, obliquity, eccentricity, and different host stars.

## Build and Run

Before building, set machine-specific paths in two files:

```bash
# 1. Set FC (Fortran compiler) and WDIR in model/Makefile
#    FC is typically an HDF5-aware wrapper: h5fc (wrapping ifort or gfortran)
#    WDIR must be the absolute path to model/

# 2. Set wdir in runEBM.sh to match your system

# 3. Copy a namelist template and edit parameters
cp namelists/input.nml.earth.aqua.23 input.nml
# Edit input.nml as needed

# 4. Build and run
./runEBM.sh
```

**Bifurcation analysis** (sweeps solar constant across multiple runs):
```bash
./bifurcation.sh
```

**Outputs** are written to `model/out/`:
- `tempseries.out` — time series: time, Tmax, Tmin, pg0, pco2, co2_condensation, fh2, diffusion
- `zonal.out` — zonal (latitudinal) temperature profiles per timestep
- `geog.out` — geographic temperature profiles

**Plotting** uses NCL scripts in `plots/` (e.g., `plotTempSeries.ncl`, `plotBistability.ncl`).

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
- `getOLR(fco2, tg0, olr)` → outgoing longwave radiation (bilinear interpolation in log-CO2 × T)
- `getPALB(fco2, tg0, zy, surfalb, palb)` → planetary albedo (quadrilinear interpolation)

Both query precomputed HDF5 lookup tables indexed by CO2 fraction, temperature, solar zenith angle, and surface albedo. At startup, `radiation_init(radfile)` reads the HDF5 file at `radfile` into two in-memory arrays (`olr_table(nco2,ntmp)` and `palb_table(nco2,ntmp,nzen,nsab)`) indexed directly by grid position. The HDF5 file path is passed as an argument — configurable via `radfile` in the `&radiation` namelist (default: `./radiation/radiation_N2_CO2_Sun.h5`).

### Namelist Configuration

All model parameters are controlled through `input.nml` (5 groups):
- `&ebm` — orbital properties, heat capacity, diffusion coefficient, CO2 levels, timestep
- `&radiation` — solar constant, albedo values, cloud effects, host star; `radfile` sets the HDF5 lookup table path
- `&co2cycle` — outgassing rate, weathering parameters
- `&h2cycle` — H2 cycling (partially implemented)
- `&stochastic` — stochastic noise amplitude and seed

The `namelists/` directory contains 30+ pre-configured scenarios (Earth aquaplanet at various obliquities, land surfaces, Mars, habitability/limit-cycle runs, exoplanet variants).

### Machine Configuration

`config/machine.sh` is a symlink to a machine-specific shell script (currently `merlin.sh`) that sources the correct compiler environment (Intel oneAPI). An alternative `discover.sh` exists for NASA Discover cluster. To port to a new machine, add a new config script and update the symlink.

## Key Technical Notes

- The radiation HDF5 lookup tables are **not** stored in the repo; they must be present at the path `model/radiation/` points to for the model to run.
- Compiler flags include `-parallel` (ifort) — the model supports light multi-threading through the HDF5 layer.
- `model/driver` (the compiled binary) and output files in `model/out/` are gitignored.
- `input.nml` at the repo root is gitignored; `model/input.nml` is the copy used at runtime (written by `runEBM.sh`).
