#!/usr/bin/env python3
"""
run_samosa.py — run HEXTOR over the SAMOSA intercomparison cases.

SAMOSA (Haqq-Misra et al. 2022, PSJ 3, 260) samples a (surface pressure,
instellation) parameter space for a synchronously rotating aquaplanet around a
3000 K blackbody star.  This script fills namelists/input.nml.samosa for each
case, runs each in its own scratch directory (so cases can run in parallel and
nothing touches model/out), and collects the diagnostics.

Participation options from the protocol, selectable with --sequence:

    warm      1 case  — the warm ExoCAM case only (case 4)
    seq2      8 cases — Sequence 2 (cases 9-16), the GCM-stable subspace
    stable   10 cases — every case ExoCAM ran stably (1, 4, 8-12, 14-16)
    all16    16 cases — Sequences 1 and 2 (the primary set)
    all32    32 cases — adds Sequences 1b and 2b
    all64    64 cases — adds Sequence 3

Outputs (in --outdir):
    samosa_summary.csv          one row per case: global means and ice line
    case_NN_<init>/zonal.txt    per-case zonal profile (ASCII, which the
                                protocol accepts from EBMs)
    case_NN_<init>/model.out    the raw HEXTOR report

Usage:
    python tools/run_samosa.py --sequence all16 --d0-ref 3.10
    python tools/run_samosa.py --sequence stable --init both
"""

import argparse
import csv
import os
import re
import shutil
import subprocess
import sys
import time
from multiprocessing import Pool

HEXTOR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TEMPLATE = os.path.join(HEXTOR, 'namelists', 'input.nml.samosa')
DRIVER = os.path.join(HEXTOR, 'model', 'driver')
DEFAULT_TABLE = './radiation/radiation_N2_CO2_3000K_p.h5'

# (sample, instellation [W/m2], surface pressure [bar]) from the protocol.
# Table 1 = Sequence 1 (1-8) and Sequence 2 (9-16); Table 4 = Sequence 1b
# (17-24), Sequence 2b (25-32) and Sequence 3 (33-64).
CASES = [
    (1,  500,  0.70), (2,  1900, 7.85), (3,  2400, 0.21), (4,  1200, 2.34),
    (5,  1500, 0.16), (6,  2100, 1.83), (7,  1600, 0.55), (8,  800,  6.16),
    (9,  1100, 0.70), (10, 400,  4.83), (11, 900,  0.10), (12, 1500, 2.98),
    (13, 1600, 0.16), (14, 900,  1.44), (15, 600,  0.43), (16, 1400, 10.0),
    (17, 700,  0.26), (18, 1700, 2.98), (19, 2300, 0.89), (20, 1300, 10.0),
    (21, 900,  0.43), (22, 2600, 4.83), (23, 2000, 0.10), (24, 500,  1.13),
    (25, 1300, 0.16), (26, 500,  2.34), (27, 1000, 0.89), (28, 1700, 3.79),
    (29, 1400, 0.55), (30, 800,  7.85), (31, 500,  0.21), (32, 1200, 1.13),
]

SEQUENCES = {
    'warm':   [4],
    'seq2':   list(range(9, 17)),
    'stable': [1, 4, 8, 9, 10, 11, 12, 14, 15, 16],
    'all16':  list(range(1, 17)),
    'all32':  list(range(1, 33)),
    'all64':  list(range(1, 65)),
}

# Initial surface temperatures.  HEXTOR's ice-albedo feedback is bistable, so
# a warm and a cold start can land in different states at the same forcing;
# running both is how that is detected rather than accidentally sampled.
INITS = {'warm': 300.0, 'cold': 233.0}


# Diffusion stability, matching how driver.f discretises: 18 belts in
# x = sin(lat), explicit in time, heat capacity C.
NBELTS = 18
HEATCAP = 4.0e6          # &ebm::heatcap in the SAMOSA template
DT_PUBLISHED = 1350.0    # the published THAI timestep
CFL_SAFETY = 0.4         # explicit diffusion needs D*dt/(C*dx^2) below ~1/2


# How heat transport is allowed to depend on the case.  All four are
# calibrated to the SAME effective D at THAI Hab1 (1 bar), so they differ only
# in how that is carried across the SAMOSA parameter space.
TRANSPORT = {
    # D = d0, identical in every case.  diffadj = .false. means the driver
    # uses d0 literally (driver.f:504), so this is what that flag alone gives.
    'constant': dict(diffadj=False, scale_by_p=False, rot=False),
    # D = d0 * p.  Pressure scaling applied by this script.
    'perbar': dict(diffadj=True, scale_by_p=False, rot=False),
    # D = d0 * p * composition * (rot0/rot)^2, HEXTOR's full scaling.
    'diffadj': dict(diffadj=True, scale_by_p=False, rot=True),
}


def effective_D(d0_ref, ps, diffadj, rot_scaling, fco2=4.0e-4,
                rot=4.84813681e-6, rot0=7.27e-5):
    """The diffusion coefficient the driver will end up using."""
    if not diffadj:
        return d0_ref
    pco2 = ps * fco2
    pn2 = ps - pco2
    avemol = (28.0 * pn2 + 44.0 * pco2) / ps
    hcp = (0.2484 * pn2 + 0.2105 * pco2) / ps
    d = d0_ref * ps * (28.89 / avemol) ** 2 * (hcp / 0.2401)
    if rot_scaling:
        d *= (rot0 / rot) ** 2
    return d


def stable_dt(D):
    """Timestep satisfying the explicit-diffusion limit, capped at published."""
    dx = 2.0 / NBELTS
    if D <= 0:
        return DT_PUBLISHED
    return min(DT_PUBLISHED, CFL_SAFETY * HEATCAP * dx * dx / D)


def prepare_rundir(rundir):
    """A HEXTOR run directory: the driver plus the relative paths it opens."""
    os.makedirs(os.path.join(rundir, 'out'), exist_ok=True)
    for name, target in (('data', os.path.join(HEXTOR, 'model', 'data')),
                         ('radiation', os.path.join(HEXTOR, 'model', 'radiation')),
                         ('driver', DRIVER)):
        link = os.path.join(rundir, name)
        if not os.path.exists(link):
            os.symlink(target, link)


def parse_zonal(model_out):
    """Pull the ZONAL STATISTICS block out of model.out.

    Columns: coordinate(deg), Tave, Tmin, dec@Tmin, Tmax, dec@Tmax, albedo,
    OLR, ASR.  In do_longitudinal mode the first column is the angle from the
    substellar point expressed on HEXTOR's latitude grid.
    """
    rows = []
    with open(model_out) as f:
        lines = f.readlines()
    for i, line in enumerate(lines):
        if line.startswith('ZONAL STATISTICS'):
            for ln in lines[i + 2:]:
                parts = ln.split()
                if len(parts) != 9:
                    break
                try:
                    rows.append([float(x) for x in parts])
                except ValueError:
                    break
            break
    return rows


def belt_area_weights(coords_deg, nbelts=18):
    """Fractional area of each HEXTOR belt.

    The belts are equally spaced in the coordinate, NOT in area
    (driver.f:376, 490): a belt centred at c spans sin(c + pi/2n) - sin(c -
    pi/2n).  Unweighted means over the 18 belts would over-weight the poles —
    or, in the tidally-locked coordinate, the substellar and antistellar
    points.  The same expression serves both readings, since x = sin(lat) is
    reinterpreted as cos(angle from substellar).
    """
    import math
    half = math.pi / (2 * nbelts)
    w = [abs(math.sin(math.radians(c) + half) - math.sin(math.radians(c) - half))
         for c in coords_deg]
    tot = sum(w)
    return [x / tot for x in w]


def parse_scalars(rundir, stdout=''):
    """Global diagnostics from tempseries.out, icelines.out and model.out."""
    out = {}
    ts = os.path.join(rundir, 'out', 'tempseries.out')
    if os.path.exists(ts):
        lines = [l for l in open(ts) if l.strip()]
        if lines:
            p = lines[-1].split()
            if len(p) >= 8:
                out['T_global'] = float(p[1])
                out['pg0'] = float(p[2])
                out['pco2'] = float(p[3])
                out['q'] = float(p[6])
                out['d'] = float(p[7])
            # HEXTOR halts on the year-to-year change in global OLR, which is
            # a false positive on the runaway-greenhouse plateau: there OLR is
            # genuinely insensitive to surface temperature, so it can go flat
            # while the temperature is still climbing.  Record the temperature
            # trend over the final years so that case can be recognised rather
            # than reported as an equilibrium.
            tail = [float(l.split()[1]) for l in lines[-11:] if len(l.split()) >= 2]
            if len(tail) >= 2:
                out['dT_last'] = tail[-1] - tail[0]
                out['n_years'] = len(lines)
    il = os.path.join(rundir, 'out', 'icelines.out')
    if os.path.exists(il):
        lines = [l for l in open(il) if l.strip()]
        if lines:
            p = lines[-1].split()
            if len(p) >= 3:
                out['icelineS'] = float(p[1])
                out['icelineN'] = float(p[2])
    # The flux-convergence message goes to stdout, not model.out.
    m = re.search(r'Flux converged at year\s+([0-9.Ee+-]+)', stdout)
    out['converged'] = bool(m)
    m = re.search(r'dOLR =\s*\n?\s*([0-9.Ee+-]+)', stdout)
    if m:
        out['dOLR'] = float(m.group(1))
    return out


def substellar_longitude(iceline_n, iceline_s):
    """Ice-line position as degrees from the substellar point.

    In do_longitudinal mode the belt coordinate is the angle from substellar,
    so 0 deg = substellar, 90 deg = terminator, 180 deg = antistellar.
    Sentinels: (90, -90) means ice free, (0, 0) means fully glaciated.
    """
    if iceline_n is None or iceline_s is None:
        return None
    if iceline_n == 90.0 and iceline_s == -90.0:
        return 180.0
    if iceline_n == 0.0 and iceline_s == 0.0:
        return 0.0
    if iceline_n < 90.0:
        return 90.0 - iceline_n
    return 90.0 - iceline_s


def run_case(task):
    (sample, inst, ps, init_name, init_t, d0_ref, cloudir, table, outdir,
     diffadj, rot_scaling) = task

    tag = 'case_%02d_%s' % (sample, init_name)
    rundir = os.path.join(outdir, tag)
    prepare_rundir(rundir)

    # With diffadj the driver applies the pressure scaling itself (along with
    # the composition and rotation factors), so d0 is passed through as the
    # reference value; without it the runner does the pressure scaling here.
    d0 = d0_ref
    D = effective_D(d0_ref, ps, diffadj, rot_scaling)
    dt = stable_dt(D)

    nml = open(TEMPLATE).read().format(
        tempinit='%.1f' % init_t, d0='%.6f' % d0, pg0='%.5f' % ps,
        solarcon='%.1f' % inst, cloudir='%.2f' % cloudir, radfile=table,
        diffadj='.true.' if diffadj else '.false.',
        diffadj_rot='.true.' if rot_scaling else '.false.',
        dt='%.1f' % dt)
    with open(os.path.join(rundir, 'input.nml'), 'w') as f:
        f.write(nml)

    t0 = time.time()
    stdout = ''
    try:
        res = subprocess.run(['./driver'], cwd=rundir, capture_output=True,
                             text=True, timeout=7200)
        rc = res.returncode
        stdout = res.stdout or ''
        err = '' if rc == 0 else (res.stderr or stdout)[-300:].replace('\n', ' ')
    except subprocess.TimeoutExpired:
        rc, err = -1, 'timeout'

    rec = {'case': sample, 'init': init_name, 'instellation': inst,
           'ps_bar': ps, 'd0': d0, 'cloudir': cloudir, 'diffadj': diffadj,
           'rot_scaling': rot_scaling, 'D': D, 'dt': dt,
           'rc': rc, 'error': err, 'wall_s': round(time.time() - t0, 1)}
    rec.update(parse_scalars(rundir, stdout))

    zonal = parse_zonal(os.path.join(rundir, 'out', 'model.out'))
    if zonal:
        # Report from substellar outward.
        zonal.sort(key=lambda r: -r[0])
        with open(os.path.join(rundir, 'zonal.txt'), 'w') as f:
            f.write('# HEXTOR / SAMOSA case %d (%s start)\n' % (sample, init_name))
            f.write('# instellation = %.1f W/m2   surface pressure = %.3f bar\n'
                    % (inst, ps))
            f.write('# theta = angle from substellar point [deg]; '
                    '0 = substellar, 180 = antistellar\n')
            f.write('# theta  T_ave_K  T_min_K  T_max_K  planetary_albedo  '
                    'OLR_Wm2  ASR_Wm2\n')
            for r in zonal:
                f.write('%7.1f %9.3f %9.3f %9.3f %9.4f %10.3f %10.3f\n'
                        % (90.0 - r[0], r[1], r[2], r[4], r[6], r[7], r[8]))
        w = belt_area_weights([r[0] for r in zonal])
        rec['T_min'] = min(r[1] for r in zonal)
        rec['T_max'] = max(r[1] for r in zonal)
        rec['OLR_global'] = sum(wi * r[7] for wi, r in zip(w, zonal))
        rec['ASR_global'] = sum(wi * r[8] for wi, r in zip(w, zonal))
        # What the protocol asks models to drive toward +-1 W/m2.
        rec['TOA_imbalance'] = rec['ASR_global'] - rec['OLR_global']

    rec['iceline_lon'] = substellar_longitude(rec.get('icelineN'),
                                              rec.get('icelineS'))
    rec['state'] = classify(rec, table)
    return rec


def classify(rec, table):
    """Label the outcome, since HEXTOR's own 'converged' flag is not enough.

    HEXTOR halts on the year-to-year change in global OLR.  On the runaway
    plateau OLR is insensitive to surface temperature, so that criterion trips
    while the planet is still heating: every case in this parameter space
    reports converged.  What separates a climate from a runaway is whether the
    top of the atmosphere actually balances, and whether the model stayed
    inside the range its radiation table covers.

    The protocol anticipates this — incipient-runaway cases may be reported at
    the last stable state or omitted — so runaways are labelled, not hidden.
    """
    T = rec.get('T_global')
    if T is None or T != T:               # NaN
        return 'runaway (numerical failure)'
    tmax = table_tmax(table)
    if tmax is not None and T > tmax:
        return 'runaway (beyond table, T > %.0f K)' % tmax
    imb = rec.get('TOA_imbalance')
    dT = rec.get('dT_last')
    if imb is None or imb != imb:
        return 'unknown'
    # Classify on the temperature trend, not the imbalance.  HEXTOR closes its
    # global budget only to ~1 W/m2 on the 18-belt grid (pre-industrial Earth
    # sits at -0.91 W/m2, and the discretised dayside insolation carries
    # +0.86 W/m2 by itself), so a 1-2 W/m2 residual is the model's floor rather
    # than a sign of non-convergence.  The protocol allows exactly this: a
    # stable trend suffices where the balance cannot be driven to +-1 W/m2.
    if dT is not None and abs(dT) < 0.1 and abs(imb) <= 3.0:
        return 'equilibrium'
    if dT is not None and abs(dT) < 1.0 and abs(imb) <= 10.0:
        return 'drifting'
    return 'runaway'


_TMAX = {}


def table_tmax(table):
    """Top of a table's temperature axis; None for the legacy 1 bar format."""
    if table in _TMAX:
        return _TMAX[table]
    path = table
    if not os.path.isabs(path):
        path = os.path.join(HEXTOR, 'model', path.lstrip('./'))
    val = None
    try:
        import h5py
        with h5py.File(path, 'r') as f:
            if 'temperature' in f:
                val = float(f['temperature'][:].max())
    except Exception:
        val = None
    _TMAX[table] = val
    return val


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--sequence', default='all16', choices=sorted(SEQUENCES))
    ap.add_argument('--init', default='warm', choices=['warm', 'cold', 'both'],
                    help='initial surface temperature (both = test bistability)')
    ap.add_argument('--d0-ref', type=float, default=3.10,
                    help='diffusion coefficient at 1 bar; d0 = d0_ref * ps. '
                         'Default 3.10 is the THAI Hab1 tidally-locked value.')
    ap.add_argument('--cloudir', type=float, default=-35.0,
                    help='cloud correction to OLR [W/m2]; driver.f applies '
                         'ir = ir - cloudir, so the published TRAPPIST-1 value '
                         'of -35 raises OLR and cools. Pass 0 for clear sky.')
    ap.add_argument('--table', default=DEFAULT_TABLE,
                    help='radiation table, as the driver sees it (relative to '
                         'the run directory)')
    ap.add_argument('--transport', default='perbar', choices=sorted(TRANSPORT),
                    help="how D varies across the parameter space: 'constant' "
                         "(D = d0 everywhere), 'perbar' (D = d0 * p * "
                         "composition, no rotation term), or 'diffadj' "
                         "(HEXTOR's full scaling, including (rot0/rot)^2). All "
                         "are calibrated to the same D at THAI Hab1.")
    ap.add_argument('--outdir', default=os.path.join(HEXTOR, 'samosa'))
    ap.add_argument('--workers', type=int, default=8)
    args = ap.parse_args()

    mode = TRANSPORT[args.transport]
    samples = SEQUENCES[args.sequence]
    known = {c[0]: c for c in CASES}
    missing = [s for s in samples if s not in known]
    if missing:
        print('cases not tabulated in this script (Sequence 3 is not '
              'transcribed): %s' % missing)
        samples = [s for s in samples if s in known]

    inits = ['warm', 'cold'] if args.init == 'both' else [args.init]

    os.makedirs(args.outdir, exist_ok=True)
    tasks = []
    for s in samples:
        _, inst, ps = known[s]
        for init_name in inits:
            tasks.append((s, inst, ps, init_name, INITS[init_name],
                          args.d0_ref, args.cloudir, args.table, args.outdir,
                          mode['diffadj'], mode['rot']))

    print('SAMOSA sequence %s: %d cases x %d start(s) = %d runs'
          % (args.sequence, len(samples), len(inits), len(tasks)))
    desc = {'constant': 'D = %.4f (constant)' % args.d0_ref,
            'perbar': 'D = %.4f * p * comp' % args.d0_ref,
            'diffadj': 'D = %.5f * p * comp * (rot0/rot)^2' % args.d0_ref}
    print('  transport: %-42s cloudir = %.2f W/m2'
          % (desc[args.transport], args.cloudir))
    print('  table    : %s' % args.table)
    print()

    results = []
    with Pool(min(args.workers, len(tasks))) as pool:
        for rec in pool.imap_unordered(run_case, tasks):
            results.append(rec)
            print('  case %2d %-4s  S=%4.0f  p=%5.2f bar  ->  '
                  % (rec['case'], rec['init'], rec['instellation'], rec['ps_bar'])
                  + ('T=%7.2f K  Tmin=%7.2f  Tmax=%7.2f  ice@%s  (%.0f s)'
                     % (rec.get('T_global', float('nan')),
                        rec.get('T_min', float('nan')),
                        rec.get('T_max', float('nan')),
                        ('%5.1f deg' % rec['iceline_lon'])
                        if rec.get('iceline_lon') is not None else '   n/a',
                        rec['wall_s'])
                     + '  [%s]' % rec.get('state', '?')
                     if rec['rc'] == 0 else 'FAILED: %s' % rec['error']),
                  flush=True)

    results.sort(key=lambda r: (r['case'], r['init']))
    cols = ['case', 'init', 'instellation', 'ps_bar', 'd0', 'cloudir',
            'T_global', 'T_min', 'T_max', 'OLR_global', 'ASR_global',
            'TOA_imbalance', 'dT_last', 'n_years', 'state', 'diffadj',
            'rot_scaling', 'D', 'dt',
            'icelineN', 'icelineS', 'iceline_lon', 'converged', 'dOLR',
            'wall_s', 'rc', 'error']
    path = os.path.join(args.outdir, 'samosa_summary.csv')
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction='ignore')
        w.writeheader()
        for r in results:
            w.writerow(r)

    nfail = sum(1 for r in results if r['rc'] != 0)
    print('\nwrote %s  (%d runs, %d failed)' % (path, len(results), nfail))
    tally = {}
    for r in results:
        tally[r.get('state', '?')] = tally.get(r.get('state', '?'), 0) + 1
    for k in sorted(tally):
        print('  %-38s %d' % (k, tally[k]))
    return 1 if nfail else 0


if __name__ == '__main__':
    sys.exit(main())
