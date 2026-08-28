#!/usr/bin/env python3
"""
samosa_sensitivity.py — how much of HEXTOR's SAMOSA answer depends on the two
parameters that are not fixed by the protocol.

Sweeps the cloud correction (cloudir) and the treatment of heat transport, and
reports how the equilibrium/runaway split and the global mean temperatures
move.  Both parameters are calibrated at a single point — THAI Hab1, 1 bar,
900 W/m2, around a 2600 K star — and then applied across 0.1-10 bar and
400-2400 W/m2 around a 3000 K star, so knowing the spread matters more than
knowing any single number.

Three transport treatments, each calibrated to give the SAME effective
diffusion at THAI conditions, so they are indistinguishable there and differ
only in how that calibration is carried across SAMOSA's parameter space:

  constant  D = d0 in every case.  This is what diffadj = .false. gives on its
            own (driver.f:504 sets d = d0 and nothing else touches it).
  perbar    D = d0 * p * composition.  Transport per bar carries over from the
            THAI calibration; the rotation term is dropped via diffadj_rot.
  diffadj   D = d0 * p * composition * (rot0/rot)^2, HEXTOR's full scaling.
            The rotation factor is 37 for TRAPPIST-1e but 225 for a 15 d
            rotator, so this gives ~6x more transport here than perbar.

Note the timestep: HEXTOR integrates diffusion explicitly, and run_samosa.py
sets dt per case from the stability limit.  At the published 1350 s the
rotation-scaled runs integrate quietly to NaN rather than failing.

    python tools/samosa_sensitivity.py
"""

import argparse
import csv
import os
import subprocess
import sys

HEXTOR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RUNNER = os.path.join(HEXTOR, 'tools', 'run_samosa.py')

# d0 for each transport treatment, all calibrated against the THAI Hab1
# ensemble with the new ExoRT 2600 K table (tools/calibrate_thai.py) so that
# every treatment gives the SAME effective D ~ 3.3 at THAI conditions (1 bar,
# 6.1 d rotation).  They differ only in how that is carried across SAMOSA's
# 0.1-10 bar and 15 d rotation:
#
#   constant  D = 3.33 everywhere                    (no scaling at all)
#   perbar    D = 3.33 * p                           (pressure only)
#   diffadj   D = 19.7 * p at SAMOSA's rotation rate (pressure x rotation)
D0 = {'constant': 3.33, 'perbar': 3.0262, 'diffadj': 0.07971}
MODES = ['constant', 'perbar', 'diffadj']

CLOUDIR = [-30.0, -35.0, -40.0, -45.0, -50.0, -55.0, -60.0]


def run(cloudir, mode, sequence, table, workdir, workers):
    tag = '%s_c%+.0f' % (mode, cloudir)
    outdir = os.path.join(workdir, tag)
    cmd = [sys.executable, RUNNER, '--sequence', sequence, '--init', 'warm',
           '--transport', mode, '--d0-ref', str(D0[mode]),
           '--cloudir', str(cloudir), '--table', table,
           '--outdir', outdir, '--workers', str(workers)]
    subprocess.run(cmd, cwd=HEXTOR, capture_output=True, text=True, timeout=14400)

    rows = []
    csvpath = os.path.join(outdir, 'samosa_summary.csv')
    if os.path.exists(csvpath):
        with open(csvpath) as f:
            rows = list(csv.DictReader(f))
    return tag, rows


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--sequence', default='all16')
    ap.add_argument('--table', default='./radiation/radiation_N2_CO2_3000K_p.h5')
    ap.add_argument('--workdir', default=os.path.join(HEXTOR, 'samosa_sensitivity'))
    ap.add_argument('--workers', type=int, default=4)
    ap.add_argument('--out', default=None)
    args = ap.parse_args()
    os.makedirs(args.workdir, exist_ok=True)

    print('SAMOSA sensitivity: %d cloudir values x %d transport treatments'
          % (len(CLOUDIR), len(MODES)))
    print('  constant : D = %.3f everywhere' % D0['constant'])
    print('  perbar   : D = %.3f * p' % (D0['perbar'] * 1.1004))
    print('  diffadj  : D = %.3f * p  (includes (rot0/rot)^2 = 224.9)'
          % (D0['diffadj'] * 1.1004 * 224.9))
    print()

    all_rows = []
    print('%-16s %5s %5s %6s %8s %8s   %s'
          % ('config', 'eq', 'run', 'nan', 'T_min', 'T_max', 'equilibrium cases'))
    print('-' * 96)
    for mode in MODES:
        for c in CLOUDIR:
            tag, rows = run(c, mode, args.sequence, args.table,
                            args.workdir, args.workers)
            eq, run_, nan = [], 0, 0
            for r in rows:
                st = r.get('state', '')
                if st == 'equilibrium':
                    eq.append((int(r['case']), float(r['T_global'])))
                elif 'numerical' in st:
                    nan += 1
                    run_ += 1
                elif st.startswith('runaway'):
                    run_ += 1
                r['config'] = tag
                all_rows.append(r)
            Ts = [t for _, t in eq]
            print('%-16s %5d %5d %6d %8s %8s   %s'
                  % (tag, len(eq), run_, nan,
                     '%.1f' % min(Ts) if Ts else '-',
                     '%.1f' % max(Ts) if Ts else '-',
                     ' '.join(str(c) for c, _ in sorted(eq))), flush=True)
        print('-' * 96)

    out = args.out or os.path.join(args.workdir, 'sensitivity_summary.csv')
    if all_rows:
        cols = ['config'] + [k for k in all_rows[0] if k != 'config']
        with open(out, 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=cols, extrasaction='ignore')
            w.writeheader()
            w.writerows(all_rows)
        print('wrote %s' % out)

    # Per-case spread across the cloudir sweep, at fixed transport treatment.
    print()
    print('global mean temperature across cloudir = %.0f .. %.0f W/m2:'
          % (CLOUDIR[0], CLOUDIR[-1]))
    print('%-6s %-21s %-21s %-21s'
          % ('case', 'constant D', 'D ~ p', 'D ~ p, rotation'))
    for case in sorted({int(r['case']) for r in all_rows}):
        cells = []
        for pref in MODES:
            Ts = [float(r['T_global']) for r in all_rows
                  if int(r['case']) == case and r['config'].startswith(pref)
                  and r.get('state') == 'equilibrium']
            cells.append('%6.1f-%6.1f K (%d)' % (min(Ts), max(Ts), len(Ts))
                         if Ts else ' runaway throughout')
        print('%-6d %-21s %-21s %-21s' % (case, cells[0], cells[1], cells[2]))
    return 0


if __name__ == '__main__':
    sys.exit(main())
