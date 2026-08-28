#!/usr/bin/env python3
"""
calibrate_thai.py — re-derive HEXTOR's THAI Hab1 calibration against a given
radiation table.

This reproduces the published calibration procedure (calibrate_thai_hab1.py,
HEXTOR 4.2.1) with two changes: the radiation table is selectable, and runs
happen in isolated scratch directories in parallel so nothing touches
model/out.  The point is that the published (d0, cloudir) = (3.10, -35.0) was
fitted against the OLD 1 bar 2600 K lookup table, so part of that -35 W/m^2
compensates for that table rather than for clouds.  Re-running the same
procedure against a new table separates the two.

Method, as published:
  * for each cloudir, sweep d0 and interpolate the d0 that reproduces the THAI
    4-GCM ensemble global mean;
  * among those solutions, choose the cloudir whose day-night contrast best
    matches the ensemble contrast.

Targets are the THAI 4-GCM ensemble means for Hab1 (Turbet et al. 2022,
PSJ Part II).

    python tools/calibrate_thai.py --table ./radiation/radiation_N2_CO2_2600K_p.h5
"""

import argparse
import csv
import os
import re
import subprocess
import sys
from multiprocessing import Pool

import numpy as np

HEXTOR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TEMPLATE = os.path.join(HEXTOR, 'namelists', 'input.nml.thai.hab1')
DRIVER = os.path.join(HEXTOR, 'model', 'driver')

# THAI 4-GCM ensemble means for Hab1 [K]
TARGET = 240.9          # global mean
TARGET_MIN = 194.7      # antistellar minimum
TARGET_MAX = 291.6      # substellar maximum
TARGET_CONTRAST = TARGET_MAX - TARGET_MIN

D0_SWEEP = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0]
CLOUDIR_SWEEP = [-50., -45., -40., -35., -30., -25., -20., -15., -10., -5., 0., 5., 10.]


def prepare_rundir(rundir):
    os.makedirs(os.path.join(rundir, 'out'), exist_ok=True)
    for name, target in (('data', os.path.join(HEXTOR, 'model', 'data')),
                         ('radiation', os.path.join(HEXTOR, 'model', 'radiation')),
                         ('driver', DRIVER)):
        link = os.path.join(rundir, name)
        if not os.path.exists(link):
            os.symlink(target, link)


def patch_nml(text, d0, cloudir, table):
    text = re.sub(r'([ \t]+d0[ \t]*=[ \t]*)[\d.eE+-]+', r'\g<1>%s' % d0, text)
    text = re.sub(r'([ \t]*cloudir[ \t]*=[ \t]*)[-\d.eE+]+', r'\g<1>%s' % cloudir, text)
    text = re.sub(r"([ \t]*radfile[ \t]*=[ \t]*)'[^']*'", r"\g<1>'%s'" % table, text)
    return text


def run_one(task):
    d0, cloudir, table, work = task
    rundir = os.path.join(work, 'd%.3f_c%+.1f' % (d0, cloudir))
    prepare_rundir(rundir)
    with open(TEMPLATE) as f:
        text = patch_nml(f.read(), d0, cloudir, table)
    with open(os.path.join(rundir, 'input.nml'), 'w') as f:
        f.write(text)

    nan = (d0, cloudir, float('nan'), float('nan'), float('nan'))
    try:
        res = subprocess.run(['./driver'], cwd=rundir, capture_output=True,
                             text=True, timeout=1800)
        if res.returncode != 0:
            return nan
    except subprocess.TimeoutExpired:
        return nan

    ts = os.path.join(rundir, 'out', 'tempseries.out')
    mo = os.path.join(rundir, 'out', 'model.out')
    if not (os.path.exists(ts) and os.path.exists(mo)):
        return nan
    lines = [l for l in open(ts) if l.strip()]
    if not lines:
        return nan
    T = float(lines[-1].split()[1])

    # Zonal block: first row is the antistellar belt, last is the substellar.
    zon = []
    seen = False
    for line in open(mo):
        if 'ZONAL STATISTICS' in line:
            seen = True
            continue
        if seen:
            if line.strip().startswith('latitude') or not line.strip():
                continue
            try:
                float(line.split()[0])
                zon.append(line)
            except (ValueError, IndexError):
                if zon:
                    break
    if not zon:
        return nan
    return (d0, cloudir, T, float(zon[0].split()[1]), float(zon[-1].split()[1]))


def interpolate_d0(d0s, Ts):
    """The d0 reproducing TARGET on the first ascending crossing (as published)."""
    pairs = [(d, t) for d, t in zip(d0s, Ts) if not np.isnan(t)]
    for (dl, tl), (dh, th) in zip(pairs, pairs[1:]):
        if tl <= TARGET <= th and th > tl:
            return float(np.interp(TARGET, [tl, th], [dl, dh]))
    for (dl, tl), (dh, th) in zip(pairs, pairs[1:]):
        if th <= TARGET <= tl and tl > th:
            return float(np.interp(TARGET, [th, tl], [dh, dl]))
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--table', default='./radiation/radiation_N2_CO2_2600K_p.h5',
                    help='radiation table, as the driver sees it')
    ap.add_argument('--work', default=None, help='scratch directory')
    ap.add_argument('--out', default=None, help='CSV of every sweep point')
    ap.add_argument('--workers', type=int, default=4)
    args = ap.parse_args()

    work = args.work or os.path.join('/tmp/claude-1000/-hugespace-models-hextor/'
                                     '380c2fc7-6424-4334-8480-5cc50da5b2ff/'
                                     'scratchpad', 'thaical',
                                     os.path.basename(args.table))
    os.makedirs(work, exist_ok=True)

    print('THAI Hab1 calibration against %s' % args.table)
    print('  targets: global %.1f K   antistellar %.1f K   substellar %.1f K'
          '   contrast %.1f K' % (TARGET, TARGET_MIN, TARGET_MAX, TARGET_CONTRAST))
    print('  sweep  : %d d0 x %d cloudir = %d runs'
          % (len(D0_SWEEP), len(CLOUDIR_SWEEP), len(D0_SWEEP) * len(CLOUDIR_SWEEP)))
    print()

    tasks = [(d, c, args.table, work) for c in CLOUDIR_SWEEP for d in D0_SWEEP]
    results = []
    with Pool(args.workers) as pool:
        for r in pool.imap_unordered(run_one, tasks):
            results.append(r)
    results.sort(key=lambda r: (r[1], r[0]))

    if args.out:
        with open(args.out, 'w', newline='') as f:
            w = csv.writer(f)
            w.writerow(['d0', 'cloudir', 'T_global', 'T_antistellar', 'T_substellar'])
            w.writerows(results)

    # For each cloudir: the d0 that hits the global-mean target, then verify.
    print('%9s %9s %9s %9s %9s %10s' %
          ('cloudir', 'd0*', 'T_global', 'T_anti', 'T_sub', 'contrast'))
    print('-' * 60)
    curve = []
    verify = []
    for c in CLOUDIR_SWEEP:
        sub = [r for r in results if r[1] == c]
        d0s = [r[0] for r in sub]
        Ts = [r[2] for r in sub]
        d0star = interpolate_d0(d0s, Ts)
        if d0star is None:
            print('%9.1f %9s   target not bracketed (T range %.1f-%.1f K)'
                  % (c, '-', np.nanmin(Ts), np.nanmax(Ts)))
            continue
        verify.append((d0star, c, args.table, work))
        curve.append((c, d0star))

    if verify:
        with Pool(args.workers) as pool:
            vres = list(pool.imap_unordered(run_one, verify))
        vres.sort(key=lambda r: r[1])
        best = None
        for d0, c, T, tmin, tmax in vres:
            contrast = tmax - tmin
            print('%9.1f %9.3f %9.2f %9.2f %9.2f %10.2f'
                  % (c, d0, T, tmin, tmax, contrast))
            if best is None or abs(contrast - TARGET_CONTRAST) < abs(best[5] - TARGET_CONTRAST):
                best = (c, d0, T, tmin, tmax, contrast)
        print('-' * 60)
        print('%9s %9s %9.1f %9.1f %9.1f %10.1f  <- THAI ensemble'
              % ('target', '', TARGET, TARGET_MIN, TARGET_MAX, TARGET_CONTRAST))
        if best:
            print()
            print('best match on day-night contrast:')
            print('   cloudir = %.1f W/m2   d0 = %.3f' % (best[0], best[1]))
            print('   T_global = %.2f K (target %.1f)   contrast = %.2f K (target %.1f)'
                  % (best[2], TARGET, best[5], TARGET_CONTRAST))
    return 0


if __name__ == '__main__':
    sys.exit(main())
