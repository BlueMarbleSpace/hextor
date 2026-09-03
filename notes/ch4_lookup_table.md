# Adding CH4 to the HEXTOR radiation lookup tables

Scoping note, 2026-09-01. Written after the pressure-axis work (merged in
`14e8667`); everything below was measured on this machine, not estimated from
first principles, except where marked.

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
