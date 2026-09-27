# HEXTOR (Habitable EBM for eXoplaneT ObseRvations)

README file for release 4.3.0

## Building

1. Source the Intel oneAPI environment: `source /opt/intel/oneapi/setvars.sh`
   (unset `SETVARS_COMPLETED` first in a shell that inherited it, or the script
   silently does nothing). The driver links against Intel's `libimf` and will
   not start without it.

2. `cd model && make`. The Makefile names this machine's HDF5-aware Fortran
   wrapper (`/opt/hdf5/bin/h5fc`, an HDF5 built against Intel ifx). On another
   machine override it on the command line, `make FC=h5pfc`; paths are relative
   to `model/`, so nothing else needs editing.

## Running

3. Edit `runEBM.sh` and set `wdir` to this directory.

4. Copy a namelist template from `./namelists/` to `./input.nml` and edit it.
   `radfile` in `&radiation` selects the HDF5 radiation lookup table under
   `model/radiation/`; the legacy 1 bar tables and the pressure-resolved tables
   (with or without a CH4 axis) are detected automatically.

5. The integration runs to the convergence test (`iterhalt = .true.`,
   `fluxcnvg`, optionally `cnvgcycle`) with a hard cap of `niter = 5000` orbits
   set in `model/driver.f`.

6. Run `./runEBM.sh`. Output is in `./model/out/`:
   - `model.out`: zonal statistics of the final orbit, ice lines, geography
   - `tempseries.out`: year, global mean temperature (K), pg0 (bar), pco2 (bar),
     pco2soil (bar), gammaout, instellation (W m^-2), diffusion coefficient
   - `zonal.out`, `geog.out`: extracted from `model.out` by `runEBM.sh`
   - `convergence.out`: orbits to the halt, which test fired (1 year-to-year,
     2-4 a cycle of that many years, 0 the iteration cap), the final
     year-to-year changes, and the largest north-south difference of the
     annual-mean belt temperatures
   - `fillet.out`, `fillet_global.out`, `transport.out` when `fillet = .true.`

## Regression test

`python tests/regression/run_regression.py` runs eight published
configurations and compares every output file byte for byte with the recorded
baseline. Run it after any change to `driver.f` or `radiation.f90`. Physics or
numerics changes go behind namelist switches whose defaults keep every case
identical; `tools/switch_sensitivity.py` measures what a switch does.

## Numerics switches (4.3.0)

All default to the published behaviour: `nstepyr` (time step as a whole
number of steps per orbit), `diffcons` (flux-conservative diffusion on the
true belt edges), `icecont` (sea-ice fraction reaching 1 continuously at
`icetemp`), `cnvgcycle` (accept a periodic multi-year oscillation as
converged) and `icelinetemp` (threshold of the reported ice line, separate
from the ice-ramp threshold `icetemp`). CLAUDE.md describes each with its
measured effect.

## FILLET

`python tools/run_fillet.py --exp all` runs the FILLET benchmarks and
experiments and writes the `projectcuisines/fillet` archive layout under
`fillet/Results/hextor/`, with per-case latitude files and convergence records
beside it. `tools/tune_fillet_ben1.py` tunes Benchmark 1, `tools/fillet_compare.py`
scores two states of the submission, and `plots/fillet_before_after.py` draws
them. `notes/fillet_review_roadmap.md` and `notes/fillet_refile_2026.md`
record the September 2026 code comparison and the re-file.

## Notes

- H2 is not implemented in the lookup tables.
- The legacy 1 bar tables cannot represent the CO2-condensing corner
  (fco2 >= 0.41 with T <= 230 K); a lookup there halts with a message. Use a
  pressure-resolved table for that regime.
