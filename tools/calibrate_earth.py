#!/usr/bin/env python3
"""
calibrate_earth.py — calibrate pre-industrial Earth clouds against a
pressure-resolved (v2/v3) radiation table.

WHY THIS EXISTS
---------------
namelists/input.nml.earth reaches ~286 K on the legacy v1 Sun table with no
clouds at all.  That worked only because of two artefacts of the v1 table: its
surface-albedo axis starts at 0.2, so the 0.06 ocean was clamped to 0.2 (acting
as partial cloud cover, global albedo 0.26), and its OLR never saturates.  The
ExoColumn tables resolve the dark ocean and the runaway OLR limit, so the same
cloud-free namelist gives the clear-sky albedo 0.15 and equilibrates at 367 K.
Clouds have to be put in explicitly.  See notes/ch4_lookup_table.md.

METHOD
------
Two observed targets, two knobs:

  * global (insolation-weighted) albedo  -> 0.30, set mainly by fcloud, the
    cloud fraction of the model's zenith-dependent cloud albedo
    (acloud = -0.078 + 0.65 z, driver.f);
  * global mean surface temperature      -> 288 K, set mainly by cloudir, the
    longwave cloud offset (ir = ir - cloudir).

The knobs are only loosely coupled -- cloudir moves the albedo through the ice
line, fcloud moves the temperature strongly -- so the search is ordered the way
that is well conditioned, as in tools/calibrate_thai.py: for each fcloud,
interpolate the cloudir that gives 288 K; then interpolate, among those
solutions, the fcloud whose albedo is 0.30.  A local sweep around the estimate
refines it, and a final run verifies.

Only warm-branch, converged runs enter the interpolation.  A warm start can
still fall into a snowball at high cloud albedo, and interpolating across that
jump would return a meaningless cloudir.

A physical cross-check, not a target: CERES puts Earth's longwave cloud
radiative effect at roughly +25 to +30 W/m2, so a calibrated cloudir far
outside that range means the cloud albedo form is doing something odd.

    python tools/calibrate_earth.py
    python tools/calibrate_earth.py --table ./radiation/radiation_N2_CO2_Sun_p.h5 --fch4 0
"""

import argparse
import csv
import json
import os
import re
import shutil
import subprocess
import sys
import time
from multiprocessing import Pool

import numpy as np

HEXTOR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TEMPLATE = os.path.join(HEXTOR, 'namelists', 'input.nml.earth.ch4')
DRIVER = os.path.join(HEXTOR, 'model', 'driver')

TARGET_T = 288.0        # K, pre-industrial global mean surface temperature
TARGET_ALB = 0.30       # global Bond albedo (CERES ~0.29)

# Pre-industrial composition.  CH4 ~730 ppb is the 1750 value; using zero here
# would fold the pre-industrial CH4 greenhouse into the cloud parameters.
FCO2_PI = 2.8e-4
FCH4_PI = 7.3e-7

FCLOUD_SWEEP = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
CLOUDIR_SWEEP = [0.0, 10.0, 20.0, 30.0, 40.0, 50.0]

FLUXCNVG = 1.0e-4       # the default 0.1 halts runs still drifting ~0.7 K/yr
SNOWBALL_T = 255.0      # below this a run is treated as off the warm branch
DRIFT_TOL = 0.2         # K over the last 10 years, for "converged"


def edit_namelist(text, values, insert_after=None):
    """Set namelist keys line by line.

    Whole lines are rewritten by key rather than by regex on the value, since a
    value such as './radiation/...' contains '/', which also ends a namelist
    group -- a partial substitution once truncated &radiation silently.
    """
    out, seen = [], set()
    for ln in text.split('\n'):
        m = re.match(r'^\s*([A-Za-z_0-9]+)\s*=', ln)
        if m and m.group(1) in values:
            k = m.group(1)
            tail = ' /' if ln.rstrip().endswith('/') else ','
            out.append('   %-13s = %s%s' % (k, values[k], tail))
            seen.add(k)
        else:
            out.append(ln)
        if insert_after and m and m.group(1) == insert_after[0]:
            k, v = insert_after[1], insert_after[2]
            out.append('   %-13s = %s,' % (k, v))
            seen.add(k)
    missing = set(values) - seen
    if missing:
        raise KeyError('namelist keys not found in template: %s' % sorted(missing))
    return '\n'.join(out)


def run_case(task):
    """Run one (fcloud, cloudir) case in its own scratch directory."""
    fcloud, cloudir, table, fch4, workroot = task
    name = 'f%.4f_c%+08.3f' % (fcloud, cloudir)
    rundir = os.path.join(workroot, name)
    result = os.path.join(rundir, 'result.txt')
    if os.path.exists(result):
        with open(result) as f:
            return json.load(f)

    shutil.rmtree(rundir, ignore_errors=True)
    os.makedirs(os.path.join(rundir, 'out'))
    os.symlink(os.path.join(HEXTOR, 'model', 'data'), os.path.join(rundir, 'data'))
    os.symlink(os.path.join(HEXTOR, 'model', 'radiation'),
               os.path.join(rundir, 'radiation'))
    shutil.copy(DRIVER, os.path.join(rundir, 'driver'))

    template = open(TEMPLATE).read()
    nml = edit_namelist(template, {
        'radfile': "'%s'" % table,
        'fco2': '%.4e' % FCO2_PI,
        'fch4': '%.4e' % fch4,
        'fcloud': '%.4f' % fcloud,
        'cloudalb': '.true.',
        'cloudir': '%.3f' % cloudir,
    }, insert_after=('iterhalt', 'fluxcnvg', '%.1e' % FLUXCNVG))
    with open(os.path.join(rundir, 'input.nml'), 'w') as f:
        f.write(nml)

    t0 = time.time()
    try:
        proc = subprocess.run(['./driver'], cwd=rundir, capture_output=True,
                              text=True, timeout=3600)
        log = proc.stdout + proc.stderr
    except subprocess.TimeoutExpired:
        log = 'TIMEOUT'
    with open(os.path.join(rundir, 'run.log'), 'w') as f:
        f.write(log)

    res = dict(fcloud=fcloud, cloudir=cloudir, T=np.nan, albedo=np.nan,
               olr=np.nan, years=0, drift=np.nan, halt='failed',
               seconds=round(time.time() - t0, 1))
    try:
        mo = open(os.path.join(rundir, 'out', 'model.out')).read()
        res['T'] = float(re.search(r'planet average temperature =\s*([-\d.]+)', mo).group(1))
        res['albedo'] = float(re.search(r'planet average albedo =\s*([-\d.]+)', mo).group(1))
        res['olr'] = float(re.search(r'planet average outgoing infrared =\s*([-\d.]+)', mo).group(1))
        ts = np.loadtxt(os.path.join(rundir, 'out', 'tempseries.out'), ndmin=2)
        res['years'] = int(ts[-1, 0])
        k = max(0, len(ts) - 11)
        res['drift'] = float(ts[-1, 1] - ts[k, 1])
        if 'Flux converged' in log:
            res['halt'] = 'converged'
        elif 'Maximum iterations' in log:
            res['halt'] = 'max_iterations'
        elif 'TIMEOUT' in log:
            res['halt'] = 'timeout'
    except Exception as exc:          # a crashed run is recorded, not fatal
        res['halt'] = 'failed: %s' % exc

    # Cache the result (json keeps NaN for failed runs); the namelist and log
    # stay beside it for inspection.
    with open(result, 'w') as f:
        json.dump({k: (float(v) if isinstance(v, np.floating) else v)
                   for k, v in res.items()}, f)
    return res


def usable(r):
    """Warm-branch, converged run whose numbers can be interpolated."""
    return (r['halt'] == 'converged' and np.isfinite(r['T'])
            and r['T'] > SNOWBALL_T and abs(r['drift']) < DRIFT_TOL)


def solve(results):
    """Interpolate (fcloud, cloudir) hitting both targets, or None."""
    by_f = {}
    for r in results:
        if usable(r):
            by_f.setdefault(r['fcloud'], []).append(r)

    # For each fcloud: the cloudir giving TARGET_T, and the albedo there.
    rows = []
    for f, rs in sorted(by_f.items()):
        rs.sort(key=lambda r: r['cloudir'])
        c = np.array([r['cloudir'] for r in rs])
        T = np.array([r['T'] for r in rs])
        a = np.array([r['albedo'] for r in rs])
        for i in range(len(rs) - 1):
            if (T[i] - TARGET_T) * (T[i + 1] - TARGET_T) <= 0 and T[i] != T[i + 1]:
                w = (TARGET_T - T[i]) / (T[i + 1] - T[i])
                rows.append((f, c[i] + w * (c[i + 1] - c[i]),
                             a[i] + w * (a[i + 1] - a[i])))
                break

    for f, c, a in rows:
        print('    fcloud %.3f : cloudir %7.2f W/m2 gives %.0f K, albedo %.4f'
              % (f, c, TARGET_T, a))

    for i in range(len(rows) - 1):
        (f0, c0, a0), (f1, c1, a1) = rows[i], rows[i + 1]
        if (a0 - TARGET_ALB) * (a1 - TARGET_ALB) <= 0 and a0 != a1:
            w = (TARGET_ALB - a0) / (a1 - a0)
            return f0 + w * (f1 - f0), c0 + w * (c1 - c0)
    return None


def sweep(tasks, workers, label):
    print('%s: %d runs on %d workers' % (label, len(tasks), workers), flush=True)
    t0 = time.time()
    with Pool(workers) as pool:
        results = pool.map(run_case, tasks, chunksize=1)
    print('  done in %.1f min' % ((time.time() - t0) / 60.0), flush=True)
    bad = [r for r in results if not usable(r)]
    for r in bad:
        print('  excluded fcloud %.3f cloudir %6.2f: T %.1f K, drift %+.2f K, %s'
              % (r['fcloud'], r['cloudir'], r['T'], r['drift'], r['halt']))
    return results


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--table', default='./radiation/radiation_N2_CO2_CH4_Sun_p.h5',
                    help='radiation table, relative to model/')
    ap.add_argument('--fch4', type=float, default=FCH4_PI,
                    help='CH4 mixing ratio (default: pre-industrial %.1e)' % FCH4_PI)
    ap.add_argument('--workers', type=int, default=6)
    ap.add_argument('--work', default=os.path.join(HEXTOR, 'calibrate_earth_runs'),
                    help='scratch directory for the runs (resume-safe)')
    ap.add_argument('--csv', default=None, help='results CSV (default: <work>/results.csv)')
    args = ap.parse_args()

    os.makedirs(args.work, exist_ok=True)
    print('table   : %s' % args.table)
    print('fco2    : %.2e   fch4: %.2e' % (FCO2_PI, args.fch4))
    print('targets : T = %.1f K, global albedo = %.3f' % (TARGET_T, TARGET_ALB))

    tasks = [(f, c, args.table, args.fch4, args.work)
             for f in FCLOUD_SWEEP for c in CLOUDIR_SWEEP]
    results = sweep(tasks, args.workers, 'coarse sweep')
    print('  solving on the coarse grid:')
    est = solve(results)
    if est is None:
        print('no bracketing solution on the coarse grid -- widen the sweeps')
        return 1
    print('  coarse estimate: fcloud %.4f, cloudir %.2f' % est)

    # Local refinement.  Interpolating row by row, as the coarse solve does,
    # fails here: the neighbouring fcloud rows need not bracket the target
    # within a small cloudir window, and the estimate then silently stays at
    # the coarse value.  Near the solution both targets are close to linear in
    # both knobs, so fit T and albedo as planes over a local 3 x 3 sweep and
    # solve the 2 x 2 system instead.  Repeat until the verification run lands.
    for it in range(4):
        f0, c0 = est
        local = [(round(f0 + df, 4), round(c0 + dc, 3), args.table, args.fch4,
                  args.work)
                 for df in (-0.02, 0.0, 0.02) for dc in (-2.0, 0.0, 2.0)
                 if 0.0 < f0 + df < 1.0]
        loc = [r for r in sweep(local, args.workers, 'local sweep %d' % (it + 1))
               if usable(r)]
        results += loc
        if len(loc) < 3:
            print('  too few usable local runs to fit; keeping the estimate')
            break
        A = np.array([[1.0, r['fcloud'], r['cloudir']] for r in loc])
        cT = np.linalg.lstsq(A, np.array([r['T'] for r in loc]), rcond=None)[0]
        cA = np.linalg.lstsq(A, np.array([r['albedo'] for r in loc]), rcond=None)[0]
        M = np.array([[cT[1], cT[2]], [cA[1], cA[2]]])
        rhs = np.array([TARGET_T - cT[0], TARGET_ALB - cA[0]])
        est = tuple(np.linalg.solve(M, rhs))
        print('  local fit: dT/dfcloud %.1f K, dT/dcloudir %.3f K per W/m2, '
              'dalb/dfcloud %.3f, dalb/dcloudir %.5f per W/m2'
              % (cT[1], cT[2], cA[1], cA[2]))
        print('  refined estimate %d: fcloud %.4f, cloudir %.3f' % ((it + 1,) + est))
        chk = run_case((round(est[0], 4), round(est[1], 3), args.table,
                        args.fch4, args.work))
        if (usable(chk) and abs(chk['T'] - TARGET_T) < 0.1
                and abs(chk['albedo'] - TARGET_ALB) < 0.001):
            break

    final = run_case((round(est[0], 4), round(est[1], 3), args.table, args.fch4,
                      args.work))
    results.append(final)
    print()
    print('VERIFY  fcloud %.4f  cloudir %.3f  ->  T %.2f K  albedo %.4f  '
          'OLR %.1f W/m2  (%d yr, drift %+.3f K, %s)'
          % (final['fcloud'], final['cloudir'], final['T'], final['albedo'],
             final['olr'], final['years'], final['drift'], final['halt']))
    ok = (usable(final) and abs(final['T'] - TARGET_T) < 0.5
          and abs(final['albedo'] - TARGET_ALB) < 0.005)
    print('CALIBRATION OK' if ok else 'CALIBRATION MISSED TARGETS')

    path = args.csv or os.path.join(args.work, 'results.csv')
    keys = ['fcloud', 'cloudir', 'T', 'albedo', 'olr', 'years', 'drift', 'halt', 'seconds']
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in sorted(results, key=lambda r: (r['fcloud'], r['cloudir'])):
            w.writerow({k: r[k] for k in keys})
    print('wrote %s' % path)
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
