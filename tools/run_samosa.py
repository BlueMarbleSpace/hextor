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


def parse_scalars(rundir):
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
    il = os.path.join(rundir, 'out', 'icelines.out')
    if os.path.exists(il):
        lines = [l for l in open(il) if l.strip()]
        if lines:
            p = lines[-1].split()
            if len(p) >= 3:
                out['icelineS'] = float(p[1])
                out['icelineN'] = float(p[2])
    mo = os.path.join(rundir, 'out', 'model.out')
    if os.path.exists(mo):
        txt = open(mo).read()
        m = re.search(r'Flux converged at year\s+([0-9.Ee+-]+)', txt)
        out['converged'] = bool(m)
        m = re.search(r'\(dOLR =\s*\n?\s*([0-9.Ee+-]+)', txt)
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
    sample, inst, ps, init_name, init_t, d0_ref, cloudir, table, outdir = task

    tag = 'case_%02d_%s' % (sample, init_name)
    rundir = os.path.join(outdir, tag)
    prepare_rundir(rundir)

    # Heat transport scales with surface pressure; see the template's notes on
    # why HEXTOR's built-in diffadj is not used for a slow synchronous rotator.
    d0 = d0_ref * ps

    nml = open(TEMPLATE).read().format(
        tempinit='%.1f' % init_t, d0='%.5f' % d0, pg0='%.5f' % ps,
        solarcon='%.1f' % inst, cloudir='%.2f' % cloudir, radfile=table)
    with open(os.path.join(rundir, 'input.nml'), 'w') as f:
        f.write(nml)

    t0 = time.time()
    try:
        res = subprocess.run(['./driver'], cwd=rundir, capture_output=True,
                             text=True, timeout=7200)
        rc = res.returncode
        err = '' if rc == 0 else (res.stderr or res.stdout or '')[-300:].replace('\n', ' ')
    except subprocess.TimeoutExpired:
        rc, err = -1, 'timeout'

    rec = {'case': sample, 'init': init_name, 'instellation': inst,
           'ps_bar': ps, 'd0': d0, 'cloudir': cloudir,
           'rc': rc, 'error': err, 'wall_s': round(time.time() - t0, 1)}
    rec.update(parse_scalars(rundir))

    zonal = parse_zonal(os.path.join(rundir, 'out', 'model.out'))
    if zonal:
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
        rec['T_min'] = min(r[1] for r in zonal)
        rec['T_max'] = max(r[1] for r in zonal)
        rec['OLR_global'] = sum(r[7] for r in zonal) / len(zonal)
        rec['ASR_global'] = sum(r[8] for r in zonal) / len(zonal)

    rec['iceline_lon'] = substellar_longitude(rec.get('icelineN'),
                                              rec.get('icelineS'))
    return rec


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--sequence', default='all16', choices=sorted(SEQUENCES))
    ap.add_argument('--init', default='warm', choices=['warm', 'cold', 'both'],
                    help='initial surface temperature (both = test bistability)')
    ap.add_argument('--d0-ref', type=float, default=3.10,
                    help='diffusion coefficient at 1 bar; d0 = d0_ref * ps. '
                         'Default 3.10 is the THAI Hab1 tidally-locked value.')
    ap.add_argument('--cloudir', type=float, default=0.0,
                    help='uniform OLR offset standing in for clouds [W/m2]')
    ap.add_argument('--table', default=DEFAULT_TABLE,
                    help='radiation table, as the driver sees it (relative to '
                         'the run directory)')
    ap.add_argument('--outdir', default=os.path.join(HEXTOR, 'samosa'))
    ap.add_argument('--workers', type=int, default=8)
    args = ap.parse_args()

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
                          args.d0_ref, args.cloudir, args.table, args.outdir))

    print('SAMOSA sequence %s: %d cases x %d start(s) = %d runs'
          % (args.sequence, len(samples), len(inits), len(tasks)))
    print('  d0 = %.3f * ps    cloudir = %.2f W/m2    table = %s'
          % (args.d0_ref, args.cloudir, args.table))
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
                     if rec['rc'] == 0 else 'FAILED: %s' % rec['error']),
                  flush=True)

    results.sort(key=lambda r: (r['case'], r['init']))
    cols = ['case', 'init', 'instellation', 'ps_bar', 'd0', 'cloudir',
            'T_global', 'T_min', 'T_max', 'OLR_global', 'ASR_global',
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
    return 1 if nfail else 0


if __name__ == '__main__':
    sys.exit(main())
