# FILLET re-file, September 2026: what changed and why

*Branch `fillet-review`, worktree `/hugespace/models/hextor-fillet`, 2026-09-27.
Follows `notes/fillet_review_roadmap.md`, the analysis of Rory's FILLET code
comparison (25 Sep 2026). Constraint from Jacob: nothing that undermines prior
publications, so every physics or numerics change is a namelist switch whose
default reproduces the published behaviour, verified byte for byte.*

## 1. Summary

- The archive still held the June 2024 files (`852fd03`). This re-file
  regenerates every benchmark and experiment with the fixed code, an untuned
  Benchmark 2/3 configuration, the corrected Experiment 2a grid, a tighter
  convergence test, and two numerics corrections declared in every header
  (flux-conservative diffusion and an even step count per orbit).
- Nothing changes for any published configuration: the eight-case regression
  set (FILLET Ben1–3, pre-industrial Earth with and without CH4, THAI Hab1, two
  SAMOSA moist cases) is byte-identical with the switches off.
- Headline numbers against the archived submission and the 4.2.0 regeneration:

  | | archived (852fd03) | 4.2.0 regeneration | **now** |
  |---|---|---|---|
  | Ben1 Tglob / N edge | 288.03 K / 57.2° | 283.81 / 51.7° | **288.00 / 56.7°** |
  | Ben2 Tglob / N edge | 287.15 / 70.1° | 271.96 / 44.0° | **279.16 / 52.0°** |
  | Ben3 Tglob | 288.89 (ice-free) | 282.77 | **285.83** |
  | Exp 1 census free/caps/belt/snowball | 113/16/1/60 | 85/25/2/78 | **86/34/3/67** (30 of 190 states differ from archived) |
  | Exp 2 census | 151/0/5/34 | 151/0/7/32 | **162/3/1/24** (15 differ) |
  | Exp 3 glaciation / deglaciation / width | 0.950 / 1.160 / 0.210 | 0.99 / 1.16 / 0.170 | **0.950 / 1.125 / 0.175** S⊕ |
  | Exp 4 warm branch at the 100 ppm floor | caps, never glaciates | glaciates at 139 ppm | **caps to 43.5°, 272.0 K, never glaciates** |
  | Exp 4 deglaciation | 1.93e4 ppm | 2.44e4 | **1.21e4** |
  | branch-order violations (warm start colder than cold start) | 0 | 14 (≤ 0.02 K) | **0** |
  | largest N–S ice-edge asymmetry, Exp 1 / 4 | 0.7° / 0.1° | 2.7° / 6.3° | **0.3° / 0.1°** |

## 2. What the audit found, and what was done

| Audit item (§4.2) | Change | Verification |
|---|---|---|
| Zenith argument, nearest-level CO2, CO2 coordinate (HIGH) | already fixed in April–July 2026 | in the 4.2.0 files; the coordinate constant was 44/28, not π/2 (0.04 % apart) |
| Makefiles build the deleted `m_ffhash.f90` (LOW) | committed Makefiles fixed; paths relative; `-parallel` dropped (ifx never honoured it) | fresh build reproduces the published binary bit for bit on all eight regression cases |
| Filed D is `d0`, used D is 1.10 `d0` | FILLET file writes the D used; Inst/XCO2 printed to 4 decimals | `fillet_global.out` only |
| ATOA an unweighted time mean including polar night | per-belt Asurf/ATOA insolation-weighted (1 − ASR/S); stated in the file | polar ATOA 0.73 → 0.64 with Asurf 0.60; `model.out` untouched |
| Weak convergence test, N–S asymmetric ice edges | `out/convergence.out` (orbits, which test fired, final drift, N–S asymmetry); `fluxcnvg` 1e-3 for FILLET; `cnvgcycle` accepts 2–4-year cycles; `nstepyr = 730` | asymmetry results in §4 |
| Ghost cells start at 273 K | seeded from the initial profile | ≤ 0.012 K global, ≤ 0.023 K per belt, only in loosely converged runs |
| Sea-ice jump 0.63 → 1 at 263.15 K | `icecont` option (normalised ramp) | −4.0 K on Ben2; offered, not adopted |
| Diffusion not flux-conservative | `diffcons` option (flux form on true belt edges) | closes the budget; adopted for the FILLET files |
| Ice-line array overflow; cap + belt misreported | array sized to the grid; ≥ 3 crossings classified and logged; `icelinetemp` separates the reported threshold from the ramp | no regenerated case has > 2 crossings |
| Legacy table duplicate CO2 level, −1 sentinels | loader collapses the level; sentinel lookups halt | affects only 249.9–274.9 ppm (uninitialised memory before) and the CO2-condensing corner |
| Run scripts unchecked; Exp 2a on the 1a grid | `tools/run_fillet.py`: exit-code and non-finite checks, resume-safe, archive layout, Exp 2a on 0.8–0.975 au | dry-run case lists; benchmark test |
| Ben1 tuning carried into untuned runs | Ben2/3 and experiments at `cloudir = 0`, `D = 0.5` exactly; Ben1 retuned through `cloudir` alone | §3 |
| Undeclared configuration, headers stale | every header carries table and ranges, cloud offset, transport rule, surface values, orbit, initial state, convergence rule and switches | see any `global_output.dat` |

## 3. Numerics switches: what they do

`tools/switch_sensitivity.py` on the published configurations (Tglob in K,
end belts = poles, or antistellar/substellar in longitudinal mode; ASR − OLR
is the final orbit's global imbalance):

| case | variant | Tglob | ΔTglob | T(1) | T(18) | ASR−OLR |
|---|---|---|---|---|---|---|
| earth | published | 285.832 | | 256.51 | 252.34 | −0.911 |
| earth | diffcons | 284.982 | −0.850 | 252.68 | 249.12 | −0.152 |
| earth | icecont | 284.969 | −0.863 | 254.66 | 251.44 | −0.960 |
| earth_ch4 | published | 287.857 | | 261.70 | 256.82 | −0.762 |
| earth_ch4 | diffcons | 286.831 | −1.026 | 258.21 | 253.74 | −0.060 |
| earth_ch4 | icecont | 286.146 | −1.711 | 259.51 | 254.90 | −0.796 |
| fillet_ben1 | published | 283.832 | | 246.41 | 241.39 | −0.856 |
| fillet_ben1 | diffcons | 282.919 | −0.913 | 242.81 | 236.67 | −0.128 |
| fillet_ben1 | icecont | 281.898 | −1.934 | 241.83 | 240.02 | −0.910 |
| fillet_ben2 | published | 271.993 | | 239.62 | 239.62 | −1.050 |
| fillet_ben2 | diffcons | 270.832 | −1.161 | 236.40 | 236.39 | −0.340 |
| fillet_ben2 | icecont | 268.000 | −3.993 | 235.77 | 235.76 | −1.008 |
| fillet_ben3 | published | 282.807 | | 287.06 | 287.11 | −0.157 |
| fillet_ben3 | diffcons | 282.864 | +0.057 | 287.50 | 287.55 | −0.271 |
| samosa_case04_warm | published | 292.217 | | 244.33 | 334.69 | +2.213 |
| samosa_case04_warm | diffcons | 293.928 | +1.711 | 244.81 | 336.61 | −0.002 |
| samosa_case11_cold | published | 225.518 | | 188.82 | 265.15 | +1.396 |
| samosa_case11_cold | diffcons | 226.397 | +0.879 | 188.58 | 266.51 | −0.001 |
| samosa_case11_cold | icecont | 224.207 | −1.311 | 187.71 | 264.03 | +1.281 |
| thai_hab1 | published | 240.942 | | 204.84 | 300.81 | +1.474 |
| thai_hab1 | diffcons | 242.633 | +1.691 | 204.96 | 306.46 | −0.018 |
| thai_hab1 | icecont | 240.946 | +0.004 | 204.84 | 300.82 | +1.470 |

Reading it:

- **`diffcons`.** The published operator's spurious source is the imbalance
  of a converged run: 0.8–1.0 W m⁻² on Earth-like profiles (a heating), 1.4–2.2
  W m⁻² of the opposite sign in the tidally locked configuration, where the
  steep substellar and antistellar gradients sit in the belts the stencil
  gets wrong. With the flux form the tightly converged THAI and SAMOSA runs
  balance to 0.02 W m⁻²; the Earth-like residuals of 0.1–0.3 W m⁻² are
  drift at the loose 0.1 W m⁻² halt of those namelists. Global means move by
  about −1 K on Earth-like cases and +0.9 to +1.7 K on the tidally locked
  ones, with THAI Hab1's substellar point +5.7 K. That is why it stays off by
  default: (d0, cloudir) for THAI and SAMOSA and (fcloud, cloudir) for Earth
  were fitted with the published operator and must be re-derived before it
  can become the default (a 5.0 change).
- **`icecont`.** Removing the jump raises the ice fraction inside the ramp by
  up to 58 %, a stronger ice–albedo feedback: Benchmark 2 cools by 4.0 K,
  Benchmark 1 by 1.9 K, Earth by 0.9 K. It is a physics change larger than
  the numerical one, and the audit's to-do list did not ask for it, so it is
  offered as an option and not used in the re-file.
- **Ghost cells.** The 273 K DATA value through the first step changed
  converged results by at most 0.012 K globally (Ben2, halted at 0.1 W m⁻²
  while still drifting 0.04 K/yr) and 0.023 K in a belt; tightly converged
  runs were unchanged. Fixed outright.

Benchmarks under each variant, untuned configuration (`cloudir = 0`, D = 0.5,
`fluxcnvg` 1e-3, published clock; Ben1 at the placeholder `cloudir = 3`):

| variant | Ben1 Tglob / NMin | Ben2 Tglob / NMin | Ben3 Tglob |
|---|---|---|---|
| published numerics | 283.88 / 51.7 | 279.90 / 53.1 | 285.77 |
| diffcons | 283.17 / 50.7 | 279.12 / 51.9 | 285.83 |
| icecont | 281.94 / 50.9 | 278.08 / 50.4 | 285.77 |
| both | 281.15 / 50.0 | 277.47 / 49.6 | 285.83 |
| 4.2.0 files (cloudir −6.3, D 0.55, fluxcnvg 0.1) | 283.81 / 51.7 | 271.96 / 44.0 | 282.77 |
| archived (852fd03) | 288.03 / 57.2 | 287.15 / 70.1 | 288.89 |

## 4. Convergence and the seasonal clock

The published time step is 730.4 steps per orbit, so every model year ran
731 steps, overshot the equinox by 0.3 d and sampled the two hemispheres'
seasons at different phases and counts. Experiments 1 and 4 were run twice
(`experiments/fillet_eval/clock_*`), at `fluxcnvg` 1e-3 with the published
clock and with `nstepyr = 730`:

| | published clock | 730 even steps |
|---|---|---|
| Exp 1 N–S asymmetry of annual-mean belt T: median / 90th pct / max (K) | 0.0016 / 0.174 / 0.454 | 0.0000 / 0.034 / 0.443 |
| Exp 1 cases above 0.05 K | 43 | 18 |
| Exp 4 cold: median / max | 0.0007 / 0.072 | 0.0000 / 0.023 |
| Exp 4 warm: median / max | 0.0189 / 0.068 | 0.0001 / 0.087 |
| climate states changed | | 0 of 290 |
| Tglob change, rms / max | | 0.04 / 0.24 K |

Tight convergence alone already removed the large ice-edge asymmetries of
the 4.2.0 files (6.3° in Exp 4 case 65 → ≤ 0.2°). The even clock reduces the
asymmetry in 161 of 190 Exp 1 cases and all 50 Exp 4 cold cases and is
adopted. The residue (up to 0.44 K) sits in the cases that never met the
year-to-year test: at high obliquity and in glaciated states the ice
thresholds switch periodically and the annual means cycle with a period of
two to four years (S = 0.825, obliquity 60°: 220.007, 220.046, 220.021 K,
repeating), which is why 26 of 190 Exp 1 cases ran to the 5000-orbit cap.
`cnvgcycle = 4` halts those on the cycle and records its period and
amplitude.

## 5. Configuration of the re-file

- **Benchmarks 2/3 and all experiments:** protocol Table 4 values; `cloudir = 0`
  (untuned, clear-sky table), `fcloud = 0`; `diffadj` off so D = 0.5 exactly;
  heat capacities land 1e7, water 4e8, ice 1e7 J m⁻² K⁻¹; albedo 0.20 / 0.30 /
  0.60; ocean fraction 0.75 in every belt.
- **Benchmark 1:** published Earth configuration (geography, Fresnel ocean,
  HEXTOR heat capacities, d0 = 0.38 with `diffadj`, filed as the effective
  0.418) at 1361 W m⁻²; `cloudir` tuned to 288.0 K with
  `tools/tune_fillet_ben1.py`: **9.37 W m⁻²** (1.65 K per W m⁻²; northern edge
  56.7°). The 4.2.0 file carried +3.0 and reached 283.8 K.
- **Numerics:** `diffcons = .true.`, `nstepyr = 730`, `icecont` off,
  `cnvgcycle = 4`, `fluxcnvg` 1e-3, cap 5000 orbits.
- **Table:** the legacy 1 bar Sun table, as filed before, with its ranges
  declared (100 ppm floor, so Experiment 4 is clamped below it; 190–370 K with
  power-law extrapolation outside). The pressure-resolved tables are clear-sky
  with saturating OLR, an untuned Earth on them runs to 367 K, and the cloud
  parameterisation they need is what the untuned experiments forbid.
- **Ice line:** 263.15 K crossing, interpolated, as before; `icelinetemp` can
  file 273.15 K when Protocol v1.2 defines the convention.
- **Starts and branches:** warm 288 K, cold 233 K, every case from that state;
  Experiments 3/4 on the Benchmark 2 base, both branches on the same grid and
  orbit. Experiments 1a/2a vary the orbital period with a^1.5.
- **Experiment 2a:** 0.8–0.975 au, 150 cases (the archived and 4.2.0 files used
  the Experiment 1a range).

## 6. Results

`tools/fillet_compare.py archived now`, `4.2.0 now` and
`experiments/fillet_published_operator now`; the figure is
`plots/fillet_before_after.png`.

### Benchmarks

| | archived | 4.2.0 | now |
|---|---|---|---|
| Ben1: Tglob, edges N/S, Diff, OLR | 288.03, 57.2/66.6, 0.38 (filed; 0.418 used), 263.9 | 283.81, 51.7/57.6, 0.38, 247.7 | 288.00, 56.7/61.5, 0.4184, 251.0 |
| Ben2 | 287.15, 70.1, 0.50 (0.551 used), 261.4 | 271.96, 44.0, 0.50, 230.2 | 279.16, 52.0, 0.5000, 240.0 |
| Ben3 | 288.89, ice-free, 265.5 | 282.77, ice-free, 254.0 | 285.83, ice-free, 254.9 |

Benchmark 1 is back on target with a northern edge of 56.7° (the audit's
common-0 °C recomputation put HEXTOR 7–10° equatorward of the others; that
was the archived zenith and CO2 errors). Untuned Benchmark 2 sits 8 K below
the archived file and 7 K above the 4.2.0 one: the archived value carried
the zenith and CO2-coordinate errors, the 4.2.0 value the −6.3 W m⁻² Benchmark
1 tuning and 10 % extra transport. Against the audit's Table 5 (AVALON
299.6, VPLanet/POISE 302.9, OPS 284.5, Shields-Bitz 284.2 K) HEXTOR is now
the coldest member by 5 K, which is what a clear-sky table does: its OLR at
288 K is 260 W m⁻² against 234–240 for the linear fits, whose intercepts
carry Earth's clouds.

### Experiments 1, 1a, 2, 2a

Drift and state changes (matched cases; Exp 2a on the 90 cases the archived
1a-range grid shares with the corrected one):

| exp | vs archived: mean / rms / max ΔT (K), states changed | vs 4.2.0 | census now (free/caps/belt/snow) |
|---|---|---|---|
| exp1 | −2.9 / 10.9 / −49.1, 30 of 190 | +5.9 / 11.1 / +51.0, 19 | 86/34/3/67 |
| exp1a | −3.4 / 11.7 / −48.8, 30 | +5.9 / 11.1 / +50.8, 20 | 88/32/3/67 |
| exp2 | +1.2 / 15.3 / +64.7, 15 | +6.2 / 14.2 / +63.4, 15 | 162/3/1/24 |
| exp2a | +5.7 / 20.2 / +64.6, 13 of 90 | +10.7 / 20.8 / +63.3, 15 | 124/4/1/21 (150 cases) |

Warm starts glaciate more often than in the archive (67 snowballs against
60) and cold starts deglaciate more often (24 against 34): both directions
of the transition are easier, the signature of the ramp and the corrected
polar physics. The mean Exp 1/Exp 2 bistable width (largest warm-start
snowball instellation to smallest cold-start non-snowball, per obliquity,
bounded by the 1.05–1.25 S⊕ overlap) is 0.168 S⊕ against 0.212 archived and
0.163 in the 4.2.0 files. No matched warm start is colder than its cold
start (14 such cases, of ≤ 0.02 K, in the 4.2.0 files were convergence noise).

### Experiments 3 and 4

Experiment 3: glaciation at 0.950 S⊕ as archived, deglaciation 1.125 (was
1.160), width 0.175 S⊕ (was 0.210; AVALON now 0.26, VPLanet/POISE 0.42).
Experiment 4 does not close its loop on the warm side: the table's CO2 axis
ends at 100 ppm and everything below is clamped to it, where the untuned
planet keeps caps to 43.5° at 272 K. The 4.2.0 file glaciated at 139 ppm
only because it carried the −6.3 W m⁻² tuning. The cold branch deglaciates at
1.2 × 10⁴ ppm (archived 1.9 × 10⁴). Both facts are stated in the headers, as
the audit's recommendation 10 asks.

### Convergence and symmetry

Every case halts on the year-to-year OLR test at 1e-3 W m⁻² or on a 2–4-year
cycle, except five of 937 (Exp 1: 1, 1a: 1, 2: 2, 2a: 1) that reached the
5000-orbit cap, flagged in their latitude-file headers. Cases halted on a
cycle: Exp 1 15, 1a 15, 2 2, 2a 1, 3 warm 2, 4 warm 25 (periods 2, 3 and 4;
the amplitude is the recorded dOLR, ≤ 0.1 W m⁻²). The north–south ice-edge
asymmetry is at most 0.3° in Experiment 1 and 0.1° in Experiment 4 (2.7° and
6.3° in the 4.2.0 files); one Experiment 3 warm-branch case at the cap
transition has 1.3°. The annual-mean belt temperatures are symmetric to
0.05 K in all but 19 (Exp 1), 16 (1a), 5 (2), 3 (2a) and 1 (3 warm) cases,
worst 0.48 K, all at high obliquity near the glaciation threshold.

### What the conservative operator did to the set

Companion set with the published operator, everything else identical
(`experiments/fillet_published_operator/`): global means move by −0.1 to
−0.8 K on average (largest −5.4 K, an Exp 1 case at the cap transition), 25
of 937 climate states change (7 in each of Exp 1 and 1a, 3 in Exp 2, 4 in
2a, 2 each in Exp 3 and 4 warm), the Experiment 3 and 4 thresholds are
unchanged, and Benchmark 1 at fixed cloudir would be 0.84 K colder, which is
why it was retuned (9.37 W m⁻² instead of the 8.3 the old operator would
need). The correction is real but small on this configuration; it matters
more where the polar or substellar belts carry steep gradients.

## 7. Not done, and why

- The conservative operator and the continuous ramp stay off by default. Making
  either the default means re-deriving the THAI Hab1 (d0, cloudir), SAMOSA and
  Earth (fcloud, cloudir) calibrations and re-running those submissions: a 5.0
  change, planned, not started.
- The ice-line convention stays at 263.15 K until Protocol v1.2 fixes one.
- Per-case latitude files for the experiments are kept in the repository
  (`fillet/latfiles/`), not filed: the archive README asks for none.
- Tagging, pushing, the Zenodo deposit and the pull request wait for Jacob.
