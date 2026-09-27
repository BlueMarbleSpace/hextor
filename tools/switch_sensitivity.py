#!/usr/bin/env python3
"""Effect of the 4.3.0 numerics switches on the published configurations.

Runs each namelist in tests/regression/cases/ as published and once per
switch setting, and tabulates the global mean, the two end belts (the poles,
or the antistellar and substellar points in longitudinal mode), the
year-to-year drift at the halt, the north-south asymmetry and the global
energy imbalance ASR - OLR of the final orbit.  A conservative diffusion
operator should bring the imbalance of a converged run to the convergence
tolerance; the published operator leaves a spurious source of order 1 W/m2.

    python tools/switch_sensitivity.py --switch 'diffcons = .true.' \
        --switch 'icecont = .true.' --switch 'diffcons = .true., icecont = .true.'

Scratch runs go to tests/regression/runs/switch_*/ (gitignored).
"""
import argparse
import concurrent.futures
import glob
import math
import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
CASES = os.path.join(ROOT, 'tests', 'regression', 'cases')
RUNS = os.path.join(ROOT, 'tests', 'regression', 'runs')
ENV = os.path.join(ROOT, 'config', 'machine.sh')
NBELTS = 18


def inject(nml_text, switch):
    """Insert '<switch>,' before the line that closes the &ebm group."""
    if not switch:
        return nml_text
    lines = nml_text.split('\n')
    in_ebm = False
    for i, line in enumerate(lines):
        s = line.strip()
        if s.lower().startswith('&ebm'):
            in_ebm = True
            continue
        if in_ebm and s and not s.startswith('!'):
            code = s.split('!')[0]
            if re.search(r"/\s*$", code):
                lines.insert(i, '   ' + switch + ',')
                return '\n'.join(lines)
    raise ValueError('no &ebm terminator found')


def prepare(rundir, modeldir):
    if os.path.isdir(rundir):
        shutil.rmtree(rundir)
    os.makedirs(os.path.join(rundir, 'out'))
    data = os.path.join(rundir, 'data')
    os.makedirs(data)
    for name in os.listdir(os.path.join(modeldir, 'data')):
        if name != 'restart.dat':
            os.symlink(os.path.join(modeldir, 'data', name), os.path.join(data, name))
    os.symlink(os.path.join(modeldir, 'radiation'), os.path.join(rundir, 'radiation'))


def run(case, switch, tag, driver):
    modeldir = os.path.dirname(os.path.abspath(driver))
    name = os.path.splitext(os.path.basename(case))[0]
    rundir = os.path.join(RUNS, 'switch_' + tag, name)
    prepare(rundir, modeldir)
    with open(os.path.join(rundir, 'input.nml'), 'w') as f:
        f.write(inject(open(case).read(), switch))
    cmd = ['bash', '-c', 'unset SETVARS_COMPLETED; source "$1" > /dev/null 2>&1; exec "$2"',
           '_', ENV, os.path.abspath(driver)]
    proc = subprocess.run(cmd, cwd=rundir, text=True, stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT, timeout=3600)
    if proc.returncode != 0:
        return name, tag, None
    return name, tag, diagnostics(rundir)


def diagnostics(rundir):
    out = os.path.join(rundir, 'out')
    rows = [l.split() for l in open(os.path.join(out, 'tempseries.out')) if l.strip()]
    tglob, years = float(rows[-1][1]), int(rows[-1][0])
    lines = open(os.path.join(out, 'model.out')).read().split('\n')
    i = next(k for k, l in enumerate(lines) if l.startswith('ZONAL STATISTICS'))
    belts = []
    for l in lines[i + 2:i + 2 + NBELTS]:
        p = l.split()
        belts.append((float(p[0]), float(p[1]), float(p[7]), float(p[8])))  # lat, T, OLR, ASR
    area = [abs(math.sin(math.radians(lat + 5)) - math.sin(math.radians(lat - 5))) / 2 for lat, _, _, _ in belts]
    imbal = sum(a * (asr - olr) for a, (_, _, olr, asr) in zip(area, belts))
    conv = open(os.path.join(out, 'convergence.out')).read().split('\n')[1].split()
    return dict(tglob=tglob, years=years, t1=belts[0][1], t18=belts[-1][1], imbal=imbal,
                dolr=float(conv[2]), asym=float(conv[4]), orbits=int(conv[0]))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--driver', default=os.path.join(ROOT, 'model', 'driver'))
    ap.add_argument('--switch', action='append', default=[],
                    help="&ebm assignment(s) to add, e.g. 'diffcons = .true.' (repeatable)")
    ap.add_argument('--cases', nargs='*')
    ap.add_argument('--jobs', type=int, default=4)
    args = ap.parse_args()

    cases = sorted(glob.glob(os.path.join(CASES, '*.nml')))
    if args.cases:
        cases = [c for c in cases if os.path.splitext(os.path.basename(c))[0] in args.cases]
    variants = [('published', '')] + [(re.sub(r'[^A-Za-z0-9]+', '_', s).strip('_'), s) for s in args.switch]
    jobs = [(c, s, t) for c in cases for t, s in variants]
    results = {}
    with concurrent.futures.ThreadPoolExecutor(args.jobs) as pool:
        for name, tag, d in pool.map(lambda j: run(j[0], j[1], j[2], args.driver), jobs):
            results[(name, tag)] = d

    print('%-20s %-34s %9s %8s %8s %8s %8s %9s %7s %6s' % (
        'case', 'variant', 'Tglob', 'dTglob', 'T(1)', 'T(18)', 'ASR-OLR', 'dOLR', 'asymNS', 'yrs'))
    for c in cases:
        name = os.path.splitext(os.path.basename(c))[0]
        base = results.get((name, 'published'))
        for tag, s in variants:
            d = results.get((name, tag))
            if d is None:
                print('%-20s %-34s FAILED' % (name, s or 'published'))
                continue
            dT = d['tglob'] - base['tglob'] if base else float('nan')
            print('%-20s %-34s %9.3f %+8.3f %8.2f %8.2f %+8.3f %9.2e %7.3f %6d' % (
                name, s or 'published', d['tglob'], dT, d['t1'], d['t18'], d['imbal'],
                d['dolr'], d['asym'], d['orbits']))


if __name__ == '__main__':
    main()
