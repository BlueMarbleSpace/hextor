#!/usr/bin/env python3
"""
samosa_summarize.py — read a completed sensitivity sweep and say which parts of
HEXTOR's SAMOSA answer are robust and which are statements about the two free
parameters.

Reads every per-configuration summary written by tools/samosa_sensitivity.py
and reports:

  * per case and transport treatment, whether it equilibrates and over what
    range of global mean temperature across the cloudir sweep;
  * the cases that equilibrate in every configuration (the robust core), the
    ones that never do, and the ones that flip;
  * the instellation above which no case equilibrates, per treatment, which is
    the runaway threshold the model actually exhibits.

    python tools/samosa_summarize.py [sweep_dir]
"""

import csv
import os
import sys
from collections import defaultdict

MODES = ['constant', 'perbar', 'diffadj']
MODE_LABEL = {'constant': 'D constant',
              'perbar': 'D ~ p',
              'diffadj': 'D ~ p, rotation'}


def co2_fraction(case_dir):
    """Fraction of belts with CO2 condensing, from the run's own output."""
    path = os.path.join(case_dir, 'out', 'co2clouds.out')
    if not os.path.exists(path):
        return 0.0
    lines = [l for l in open(path) if l.strip()]
    if not lines:
        return 0.0
    parts = lines[-1].split()
    if len(parts) < 19:
        return 0.0
    n = 0
    for x in parts[1:19]:
        try:
            n += 1 if float(x) > 0.5 else 0
        except ValueError:
            pass
    return n / 18.0


def relabel(r, case_dir):
    """Re-derive the outcome uniformly across a sweep.

    Configurations in one sweep can have been produced by different versions of
    the classifier, so the label is recomputed here from the stored numbers
    plus the run's own CO2 output, and everything is compared on equal terms.
    """
    if co2_fraction(case_dir) > 0.5:
        return 'CO2 condensing'
    T = fnum(r, 'T_global')
    if T != T:
        return 'runaway (out of range)'
    tmax = 420.0                      # top of the 3000 K table
    if T > tmax:
        return 'runaway'
    imb, dT = fnum(r, 'TOA_imbalance'), fnum(r, 'dT_last')
    if imb != imb or dT != dT:
        return 'unknown'
    if abs(dT) < 0.5 and abs(imb) <= 3.0:
        return 'equilibrium'
    if abs(dT) < 2.0 and abs(imb) <= 10.0:
        return 'drifting'
    return 'runaway'


def load(sweep_dir):
    rows = []
    for name in sorted(os.listdir(sweep_dir)):
        path = os.path.join(sweep_dir, name, 'samosa_summary.csv')
        if not os.path.exists(path):
            continue
        mode = name.split('_c')[0]
        cloudir = float(name.split('_c')[1])
        with open(path) as f:
            for r in csv.DictReader(f):
                r['mode'] = mode
                r['cloudir_val'] = cloudir
                r['state'] = relabel(r, os.path.join(
                    sweep_dir, name, 'case_%02d_%s' % (int(r['case']), r['init'])))
                rows.append(r)
    return rows


def fnum(r, key):
    try:
        return float(r[key])
    except (KeyError, ValueError, TypeError):
        return float('nan')


def main():
    sweep = sys.argv[1] if len(sys.argv) > 1 else 'samosa_sensitivity'
    rows = load(sweep)
    if not rows:
        print('no configuration summaries found in %s' % sweep)
        return 1

    modes = [m for m in MODES if any(r['mode'] == m for r in rows)]
    ncfg = {m: len({r['cloudir_val'] for r in rows if r['mode'] == m}) for m in modes}
    cases = sorted({int(r['case']) for r in rows})
    print('sweep: %s' % sweep)
    print('  %d rows   %s' % (len(rows),
                              '   '.join('%s: %d configs' % (m, ncfg[m]) for m in modes)))
    print()

    # ---- per case x treatment -------------------------------------------
    eq = defaultdict(list)
    noneq = defaultdict(lambda: defaultdict(int))
    for r in rows:
        if r.get('state') == 'equilibrium':
            eq[(int(r['case']), r['mode'])].append(fnum(r, 'T_global'))
        else:
            noneq[(int(r['case']), r['mode'])][r['state']] += 1

    def outcome(c, m):
        """Label a cell by its dominant non-equilibrium outcome, not 'runaway'
        for everything -- a CO2-condensing collapse and a runaway greenhouse
        are opposite ends of the parameter space."""
        d = noneq.get((c, m), {})
        if not d:
            return 'runaway'
        return max(d, key=d.get)

    inst = {int(r['case']): fnum(r, 'instellation') for r in rows}
    pres = {int(r['case']): fnum(r, 'ps_bar') for r in rows}

    print('%-5s %6s %7s   %s' % ('case', 'S', 'p_bar',
                                 '   '.join('%-22s' % MODE_LABEL[m] for m in modes)))
    print('-' * (22 + 25 * len(modes)))
    for c in cases:
        cells = []
        for m in modes:
            Ts = eq.get((c, m), [])
            cells.append('%6.1f-%6.1f K %2d/%d' % (min(Ts), max(Ts), len(Ts), ncfg[m])
                         if Ts else '%-16s0/%d' % (outcome(c, m)[:15], ncfg[m]))
        print('%-5d %6.0f %7.2f   %s'
              % (c, inst[c], pres[c], '   '.join('%-22s' % x for x in cells)))
    print('-' * (22 + 25 * len(modes)))

    # ---- robustness ------------------------------------------------------
    total = sum(ncfg[m] for m in modes)
    always, never, flips = [], [], []
    for c in cases:
        n = sum(len(eq.get((c, m), [])) for m in modes)
        (always if n == total else never if n == 0 else flips).append(c)

    # Regimes other than equilibrium, counted so they are not lumped together.
    print()
    other = defaultdict(int)
    for r in rows:
        if r['state'] != 'equilibrium':
            other[r['state']] += 1
    print('non-equilibrium outcomes across all configurations:')
    for k in sorted(other):
        print('   %-26s %d' % (k, other[k]))

    print()
    print('equilibrates in every configuration : %s'
          % (' '.join(map(str, always)) or 'none'))
    print('never equilibrates                  : %s'
          % (' '.join(map(str, never)) or 'none'))
    print('depends on the free parameters      : %s'
          % (' '.join(map(str, flips)) or 'none'))

    # ---- runaway threshold ----------------------------------------------
    print()
    print('highest instellation reaching equilibrium, and lowest that never does:')
    for m in modes:
        hi = [inst[c] for c in cases if eq.get((c, m))]
        lo = [inst[c] for c in cases if not eq.get((c, m))]
        print('  %-18s equilibrium up to %4.0f W/m2   |   runaway from %4.0f W/m2'
              % (MODE_LABEL[m], max(hi) if hi else float('nan'),
                 min(lo) if lo else float('nan')))

    # ---- sensitivity to cloudir at fixed treatment -----------------------
    print()
    print('temperature sensitivity to cloudir (K per 10 W/m2 of cloud correction):')
    for m in modes:
        slopes = []
        for c in cases:
            pts = sorted((r['cloudir_val'], fnum(r, 'T_global')) for r in rows
                         if int(r['case']) == c and r['mode'] == m
                         and r.get('state') == 'equilibrium')
            if len(pts) >= 2:
                dT = pts[-1][1] - pts[0][1]
                dc = pts[-1][0] - pts[0][0]
                if dc:
                    slopes.append(abs(dT / dc) * 10.0)
        if slopes:
            print('  %-18s %.2f - %.2f K   (median %.2f)'
                  % (MODE_LABEL[m], min(slopes), max(slopes),
                     sorted(slopes)[len(slopes) // 2]))
    return 0


if __name__ == '__main__':
    sys.exit(main())
