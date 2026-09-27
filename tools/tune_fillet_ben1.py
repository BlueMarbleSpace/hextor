#!/usr/bin/env python3
"""Tune FILLET Benchmark 1's cloudir so that HEXTOR reaches 288.0 K.

Benchmark 1 is the protocol's one tuned case: present-day Earth geography
and orbit, tuned to 288 K.  HEXTOR keeps its published Earth transport
(d0 = 0.38 with diffadj) and adjusts the longwave cloud offset cloudir alone,
one knob for one target, so the tuning is well posed.  Everything else in
the case is whatever tools/run_fillet.py writes, including any numerics
switches passed through --extra-ebm, so the tuning is redone whenever they
change.

    python tools/tune_fillet_ben1.py --extra-ebm 'diffcons = .true.' --extra-ebm 'nstepyr = 730'

Prints the secant iterations and the cloudir to pass to run_fillet.py as
--cloudir-ben1.  Scratch runs go to experiments/fillet_ben1_tune/ (gitignored).
"""
import argparse
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RUNNER = os.path.join(HERE, 'run_fillet.py')


def tglob(cloudir, extra, workdir, driver):
    outdir = os.path.join(workdir, 'cloudir_%+.4f' % cloudir)
    cmd = [sys.executable, RUNNER, '--exp', 'ben1', '--outdir', outdir, '--jobs', '1',
           '--cloudir-ben1', '%.4f' % cloudir, '--driver', driver, '--label', 'Ben1 tuning']
    for e in extra:
        cmd += ['--extra-ebm', e]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)
    rows = [l.split() for l in open(os.path.join(outdir, 'Results', 'hextor', 'ben1', 'global_output.dat'))
            if l.strip() and not l.startswith('#')]
    return float(rows[0][4]), float(rows[0][6])   # Tglob, northern ice edge


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--target', type=float, default=288.0)
    ap.add_argument('--start', type=float, nargs=2, default=[3.0, 8.0],
                    help='two starting cloudir values for the secant iteration')
    ap.add_argument('--tol', type=float, default=0.005, help='|Tglob - target| to stop at (K)')
    ap.add_argument('--extra-ebm', action='append', default=[])
    ap.add_argument('--driver', default=os.path.join(ROOT, 'model', 'driver'))
    ap.add_argument('--workdir', default=os.path.join(ROOT, 'experiments', 'fillet_ben1_tune'))
    args = ap.parse_args()

    c0, c1 = args.start
    t0, _ = tglob(c0, args.extra_ebm, args.workdir, args.driver)
    t1, edge = tglob(c1, args.extra_ebm, args.workdir, args.driver)
    print('cloudir %+8.4f -> Tglob %8.3f K' % (c0, t0))
    print('cloudir %+8.4f -> Tglob %8.3f K' % (c1, t1))
    for it in range(12):
        if abs(t1 - args.target) < args.tol:
            break
        if t1 == t0:
            sys.exit('flat response; widen --start')
        c2 = c1 + (args.target - t1) * (c1 - c0) / (t1 - t0)
        c0, t0 = c1, t1
        c1 = c2
        t1, edge = tglob(c1, args.extra_ebm, args.workdir, args.driver)
        print('cloudir %+8.4f -> Tglob %8.3f K' % (c1, t1))
    print('\nBenchmark 1: cloudir = %.4f W/m2 gives Tglob = %.3f K, northern ice edge %.1f deg'
          % (c1, t1, edge))
    print('sensitivity dT/dcloudir = %.3f K per W/m2' % ((t1 - t0) / (c1 - c0) if c1 != c0 else float('nan')))
    print('pass to run_fillet.py as: --cloudir-ben1 %.4f' % c1)


if __name__ == '__main__':
    main()
