# Draft reply to Rory (FILLET code comparison, HEXTOR)

*Draft, 2026-09-27. Numbers from `notes/fillet_refile_2026.md`. The pull
request is a placeholder until it exists.*

Rory,

Thanks for the audit. Every HEXTOR finding in §4.2 checked out against the
source, and we have worked through your §9 list. The result is HEXTOR 4.3.0
(https://github.com/BlueMarbleSpace/hextor/releases/tag/4.3.0, Zenodo
https://doi.org/10.5281/zenodo.22999531). The attached figure compares the archived submission with the new
files for every benchmark and experiment.

What changed

- The three table errors you list (zenith argument, nearest-level CO2
  lookup, CO2 coordinate) were fixed in April to July; the archive never
  received those files. Your measurement of the coordinate constant is
  right: it was 44/28, a mass-weighted mixing ratio, not π/2; the repair is
  the same to 0.04 %.
- The committed Makefiles did still build the deleted hash module; the
  working copies had been fixed locally and hidden from git. A fresh clone
  builds again, and the rebuilt driver reproduces the binary behind the
  published results bit for bit on eight configurations that now form a
  regression test.
- The Diff column is the D the run used. The files before this one carried
  d0 while the model ran 1.10 d0; Benchmarks 2/3 and every experiment now
  run at D = 0.5 exactly (the rescaling is off), and Benchmark 1 files its
  effective 0.418.
- Asurf and ATOA are insolation-weighted (1 − ASR/S), stated in every file;
  the unweighted mean counted polar night at full weight where the table is
  brightest, which is why ATOA sat above Asurf at the poles.
- The diffusion operator was not flux-conservative: on a smooth profile it
  added 0.8 W m⁻² of spurious global heating and overstated the polar
  convergence by ~70 %, from the 0.0038-wide ghost cell you found. A flux
  form on the true belt edges is now available and is used for these files
  (declared in every header). It stays off by default in the code because
  our THAI and SAMOSA calibrations were fitted with the old operator and
  move by 1–2 K under the new one; re-deriving them is the next release.
- The ghost cells start from the initial profile rather than 273 K.
- The model year was 731 half-day steps for a 730.4-step orbit, so the two
  hemispheres were sampled at different phases; the files now use 730 equal
  steps per orbit. With that and a convergence tolerance of 0.001 W m⁻²
  (from 0.1), the ice-edge asymmetries you saw (6.3° in Experiment 4) are
  below 0.3°. A run whose annual means cycle with a 2–4-year period (high
  obliquity, glaciated states) halts on the cycle and says so; five of 937
  cases reached the 5000-orbit cap and are flagged.
- Benchmarks 2/3 and the experiments are untuned (cloudir = 0). Benchmark 1
  keeps its Earth configuration and is tuned through cloudir alone, to
  288.0 K.
- Experiment 2a is on 0.8–0.975 au (150 cases); the archived run used the
  1a range.
- Instellation and CO2 are printed to four decimals, the ice-line array can
  hold every crossing (a cap coexisting with a belt is reported as the cap
  and logged), and every header states the table and its ranges, the cloud
  offset, the transport rule, surface values, orbit, initial state,
  convergence rule and numerics switches. The old csh scripts are replaced
  by a runner that checks every case.

Results, archived -> now

- Benchmark 1: 288.03 K, edge 57.2°N -> 288.00 K, edge 56.7°N.
- Benchmark 2: 287.15 K, caps to 70.1° -> 279.16 K, caps to 52.0°.
- Benchmark 3: 288.89 K -> 285.83 K, ice-free.
- Experiment 1 (free / caps / belt / snowball): 113/16/1/60 -> 86/34/3/67;
  30 of 190 states change. Experiment 2: 151/0/5/34 -> 162/3/1/24.
- Experiment 3: glaciation 0.950 (unchanged), deglaciation 1.160 -> 1.125,
  so the bistable width narrows from 0.210 to 0.175 S⊕.
- Experiment 4: the cold branch deglaciates at 1.2 × 10⁴ ppm (was 1.9 ×
  10⁴). The warm branch still does not glaciate: the table's CO2 axis ends
  at 100 ppm and is clamped below it, where the untuned planet keeps caps to
  43.5°. That range is now stated in the file, as your recommendation 10
  asks; only a table with a lower floor would change it.
- Warm starts are never colder than cold starts on any matched grid.

What we kept

- The legacy 1 bar table. The pressure-resolved tables are clear sky with a
  saturating OLR; an untuned Earth on them runs to 367 K, and the cloud
  parameterisation they need is what the untuned experiments forbid.
- The 263.15 K ice line, interpolated. The reported threshold is now a
  separate parameter from the one that ends the sea-ice ramp, so we can
  file 273.15 K (or 271.15 at sea; we have one temperature per belt) as
  soon as v1.2 fixes the convention.
- The sea-ice ramp, including its jump to full cover at 263.15 K. A
  continuous version exists as an option; it cools Benchmark 2 by 4 K and
  we did not want to change the physics under the intercomparison.
- Experiments 3/4 on the Benchmark 2 base, fixed starts (288 K warm, 233 K
  cold, no continuation), both branches on one grid and orbit; 1a/2a scale
  the year as a^1.5.

We will open the pull request replacing Results/hextor once you confirm the
Experiment 3/4 base you settled with the AVALON submission.

Jacob
