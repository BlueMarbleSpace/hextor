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

Method: the two knobs are nearly decoupled, and the search has to respect
that.  cloudir shifts the global mean and barely touches the day-night
contrast; D controls the contrast and barely touches the mean.  So:

  * for each D, interpolate the cloudir that reproduces the THAI 4-GCM
    ensemble global mean;
  * among those solutions, choose the D whose day-night contrast best matches
    the ensemble contrast.

Doing it the other way round -- solving for D at fixed cloudir -- is
ill-conditioned, because the global mean hardly depends on D: the returned D
is then set by interpolation noise, and the contrast it happens to produce is
not meaningful.

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

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import moist_stability as ms

HEXTOR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TEMPLATE = os.path.join(HEXTOR, 'namelists', 'input.nml.thai.hab1')
DRIVER = os.path.join(HEXTOR, 'model', 'driver')

# THAI 4-GCM ensemble means for Hab1 [K]
TARGET = 240.9          # global mean
TARGET_MIN = 194.7      # antistellar minimum
TARGET_MAX = 291.6      # substellar maximum
TARGET_CONTRAST = TARGET_MAX - TARGET_MIN

D0_SWEEP = [2.6, 2.8, 3.0, 3.1, 3.2, 3.3, 3.4, 3.6, 3.8, 4.0]
CLOUDIR_SWEEP = [-60., -55., -50., -45., -40., -35., -30., -25., -20.]

# Diffusion stability.  These mirror the Hab1 template's &ebm values; with
# moistdiff the diffused quantity is h/cp, which responds to T by a factor
# beta > 1, so the explicit limit tightens by beta (see tools/moist_stability).
HEATCAP = 4.0e6
DT_PUBLISHED = 1350.0
T_HOT_GUESS = 330.0      # first guess for the substellar belt; refined per run


def template_composition():
    """(pg0, fco2) as the Hab1 template sets them, for the moist beta."""
    text = open(TEMPLATE).read()
    pg0 = float(re.search(r'[ \t]*pg0[ \t]*=[ \t]*([\d.eE+-]+)', text).group(1))
    fco2 = float(re.search(r'[ \t]*fco2[ \t]*=[ \t]*([\d.eE+-]+)', text).group(1))
    return pg0, fco2


# The three transport treatments, spelled exactly as run_samosa.py spells
# them, so a calibration and the run it feeds cannot disagree about what
# 'perbar' means.  The trap this closes: the Hab1 template used to carry no
# diffadj_rot, so --diffadj picked up the driver default (.true.) and silently
# calibrated the FULL scaling, rotation term included, when perbar was meant.
TRANSPORT = {'constant': dict(diffadj=False, rot=False),
             'perbar':   dict(diffadj=True,  rot=False),
             'diffadj':  dict(diffadj=True,  rot=True)}


def effective_scale(diffadj, rot):
    """Factor the driver applies to d0 at THAI Hab1, so dt can be set from the
    diffusion coefficient the model will actually use rather than from d0."""
    if not diffadj:
        return 1.0
    text = open(TEMPLATE).read()
    def g(key, default=None):
        m = re.search(r'[ \t]*%s[ \t]*=[ \t]*([\d.eE+-]+)' % key, text)
        return float(m.group(1)) if m else default
    pg0, fco2 = g('pg0'), g('fco2')
    rot_p, rot0 = g('rot'), 7.27e-5
    pn2, pco2 = pg0 * (1 - fco2), pg0 * fco2
    avemol = (28.0 * pn2 + 44.0 * pco2) / pg0
    hcp = (0.2484 * pn2 + 0.2105 * pco2) / pg0
    f = pg0 * (28.89 / avemol) ** 2 * (hcp / 0.2401)
    if rot:
        f *= (rot0 / rot_p) ** 2
    return f


def run_dt(d0, moistdiff, rhmoist, t_hot):
    """The timestep this (d0, thermal state) needs, and the beta behind it."""
    if not moistdiff:
        return ms.stable_dt(d0, HEATCAP, DT_PUBLISHED), 1.0
    pg0, fco2 = template_composition()
    dt, b = ms.moist_dt(d0, pg0, pg0 * (1 - fco2), pg0 * fco2, t_hot,
                        rhmoist, HEATCAP, DT_PUBLISHED)
    return dt, b


def prepare_rundir(rundir):
    os.makedirs(os.path.join(rundir, 'out'), exist_ok=True)
    for name, target in (('data', os.path.join(HEXTOR, 'model', 'data')),
                         ('radiation', os.path.join(HEXTOR, 'model', 'radiation')),
                         ('driver', DRIVER)):
        link = os.path.join(rundir, name)
        if not os.path.exists(link):
            os.symlink(target, link)


def patch_nml(text, d0, cloudir, table, diffadj=False, moistdiff=False,
              rhmoist=0.8, dt=DT_PUBLISHED, rot=True):
    text = re.sub(r'([ \t]+d0[ \t]*=[ \t]*)[\d.eE+-]+', r'\g<1>%s' % d0, text)
    text = re.sub(r'([ \t]*cloudir[ \t]*=[ \t]*)[-\d.eE+]+', r'\g<1>%s' % cloudir, text)
    text = re.sub(r"([ \t]*radfile[ \t]*=[ \t]*)'[^']*'", r"\g<1>'%s'" % table, text)
    text = re.sub(r'([ \t]*diffadj[ \t]*=[ \t]*)\.\w+\.',
                  r'\g<1>%s' % ('.true.' if diffadj else '.false.'), text)
    text = re.sub(r'([ \t]*diffadj_rot[ \t]*=[ \t]*)\.\w+\.',
                  r'\g<1>%s' % ('.true.' if rot else '.false.'), text)
    text = re.sub(r'([ \t]*moistdiff[ \t]*=[ \t]*)\.\w+\.',
                  r'\g<1>%s' % ('.true.' if moistdiff else '.false.'), text)
    text = re.sub(r'([ \t]*rhmoist[ \t]*=[ \t]*)[\d.eE+-]+',
                  r'\g<1>%s' % rhmoist, text)
    text = re.sub(r'([ \t]*dt[ \t]*=[ \t]*)[\d.eE+-]+', r'\g<1>%.1f' % dt, text)
    return text


def run_one(task):
    """One (d0, cloudir) point, with the timestep set from its own warmth.

    With moistdiff the stable timestep depends on the substellar temperature
    through beta, which is not known until the run is done.  So the run is
    made at a first-guess t_hot and then repeated if it turned out warmer than
    that guess -- a two-step fixed point, rather than a single conservative
    timestep that would slow every cold point in the sweep.  Without moistdiff
    beta is 1 and nothing is repeated.
    """
    d0, cloudir, table, work, diffadj, moistdiff, rhmoist, t_hot, rot = task
    for _attempt in range(2):
        out = _run_at(d0, cloudir, table, work, diffadj, moistdiff, rhmoist,
                      t_hot, rot)
        tmax = out[4]
        if not moistdiff or tmax != tmax or tmax <= t_hot + 0.5:
            return out
        t_hot = tmax + 1.0        # refine and redo at the tighter limit
    return out


def _run_at(d0, cloudir, table, work, diffadj, moistdiff, rhmoist, t_hot,
            rot=True):
    dt, bet = run_dt(d0 * effective_scale(diffadj, rot), moistdiff, rhmoist,
                     t_hot)
    rundir = os.path.join(work, 'd%.5f_c%+.1f' % (d0, cloudir))
    prepare_rundir(rundir)
    with open(TEMPLATE) as f:
        text = patch_nml(f.read(), d0, cloudir, table, diffadj, moistdiff,
                         rhmoist, dt, rot)
    with open(os.path.join(rundir, 'input.nml'), 'w') as f:
        f.write(text)

    nan = (d0, cloudir, float('nan'), float('nan'), float('nan'), dt, bet)
    try:
        res = subprocess.run(['./driver'], cwd=rundir, capture_output=True,
                             text=True, timeout=7200)
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
    return (d0, cloudir, T, float(zon[0].split()[1]),
            float(zon[-1].split()[1]), dt, bet)


def interpolate_to_target(xs, Ts):
    """The x reproducing TARGET on the first monotone crossing."""
    pairs = [(d, t) for d, t in zip(xs, Ts) if not np.isnan(t)]
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
    ap.add_argument('--transport', default=None,
                    choices=sorted(TRANSPORT),
                    help="which diffusion treatment to calibrate, using the "
                         "same names as run_samosa.py: 'constant' (D = d0), "
                         "'perbar' (D = d0 * p * composition) or 'diffadj' "
                         "(adds the (rot0/rot)^2 rotation factor). Calibrate "
                         "with the SAME option the SAMOSA run will use -- d0 "
                         "means a different D under each.")
    ap.add_argument('--diffadj', action='store_true',
                    help="deprecated alias for --transport diffadj")
    ap.add_argument('--d0-sweep', default=None,
                    help='comma-separated d0 values overriding the default grid')
    ap.add_argument('--cloudir-sweep', default=None,
                    help='comma-separated cloudir values overriding the grid')
    ap.add_argument('--moistdiff', action='store_true',
                    help='diffuse moist static energy (&ebm::moistdiff), the '
                         'configuration the RH < 1 radiation tables are built '
                         'for.  The timestep is then set per run from the '
                         'substellar beta rather than left at 1350 s.')
    ap.add_argument('--rhmoist', type=float, default=0.8,
                    help='relative humidity of the diffused moisture; match '
                         "the table's own RH (default 0.8)")
    args = ap.parse_args()

    mode = TRANSPORT[args.transport or ('diffadj' if args.diffadj
                                       else 'constant')]
    args.diffadj, args.rot = mode['diffadj'], mode['rot']

    global D0_SWEEP, CLOUDIR_SWEEP
    if args.d0_sweep:
        D0_SWEEP = [float(x) for x in args.d0_sweep.split(',')]
    if args.cloudir_sweep:
        CLOUDIR_SWEEP = [float(x) for x in args.cloudir_sweep.split(',')]

    work = args.work or os.path.join('/tmp/claude-1000/-hugespace-models-hextor/'
                                     '380c2fc7-6424-4334-8480-5cc50da5b2ff/'
                                     'scratchpad', 'thaical',
                                     os.path.basename(args.table))
    os.makedirs(work, exist_ok=True)

    print('THAI Hab1 calibration against %s' % args.table)
    print('  targets: global %.1f K   antistellar %.1f K   substellar %.1f K'
          '   contrast %.1f K' % (TARGET, TARGET_MIN, TARGET_MAX, TARGET_CONTRAST))
    scale = effective_scale(args.diffadj, args.rot)
    print('  transport: %s   (D = d0 x %.4f at Hab1)'
          % (args.transport or ('diffadj' if args.diffadj else 'constant'),
             scale))
    print('  moist  : %s' % ('.true., rhmoist = %.2f (dt set per run from '
                              'the substellar beta)' % args.rhmoist
                              if args.moistdiff else '.false. (dry diffusion)'))
    print('  sweep  : %d d0 x %d cloudir = %d runs'
          % (len(D0_SWEEP), len(CLOUDIR_SWEEP), len(D0_SWEEP) * len(CLOUDIR_SWEEP)))
    print()

    tasks = [(d, c, args.table, work, args.diffadj, args.moistdiff,
              args.rhmoist, T_HOT_GUESS, args.rot)
             for c in CLOUDIR_SWEEP for d in D0_SWEEP]
    results = []
    with Pool(args.workers) as pool:
        for r in pool.imap_unordered(run_one, tasks):
            results.append(r)
    results.sort(key=lambda r: (r[1], r[0]))

    if args.out:
        with open(args.out, 'w', newline='') as f:
            w = csv.writer(f)
            w.writerow(['d0', 'cloudir', 'T_global', 'T_antistellar',
                        'T_substellar', 'dt', 'beta'])
            w.writerows(results)

    # For each d0: the cloudir that hits the global-mean target, then verify.
    print('%9s %9s %9s %9s %9s %10s %8s %7s' %
          ('d0', 'cloudir*', 'T_global', 'T_anti', 'T_sub', 'contrast',
           'dt', 'beta'))
    print('-' * 77)
    curve = []
    verify = []
    for d in D0_SWEEP:
        sub = [r for r in results if r[0] == d]
        cs = [r[1] for r in sub]
        Ts = [r[2] for r in sub]
        cstar = interpolate_to_target(cs, Ts)
        if cstar is None:
            print('%9.3f %9s   target not bracketed (T range %.1f-%.1f K)'
                  % (d, '-', np.nanmin(Ts), np.nanmax(Ts)))
            continue
        verify.append((d, cstar, args.table, work, args.diffadj,
                       args.moistdiff, args.rhmoist, T_HOT_GUESS, args.rot))
        curve.append((d, cstar))

    if verify:
        with Pool(args.workers) as pool:
            vres = list(pool.imap_unordered(run_one, verify))
        vres.sort(key=lambda r: r[0])
        best = None
        for d0, c, T, tmin, tmax, dt, bet in vres:
            contrast = tmax - tmin
            print('%9.3f %9.2f %9.2f %9.2f %9.2f %10.2f %8.1f %7.2f'
                  % (d0, c, T, tmin, tmax, contrast, dt, bet))
            if best is None or abs(contrast - TARGET_CONTRAST) < abs(best[5] - TARGET_CONTRAST):
                best = (d0, c, T, tmin, tmax, contrast)
        print('-' * 77)
        print('%9s %9s %9.1f %9.1f %9.1f %10.1f  <- THAI ensemble'
              % ('target', '', TARGET, TARGET_MIN, TARGET_MAX, TARGET_CONTRAST))
        if best:
            print()
            print('best match on day-night contrast:')
            print('   d0 = %.4f   cloudir = %.2f W/m2   (effective D = %.4f)'
                  % (best[0], best[1], best[0] * scale))
            print('   T_global = %.2f K (target %.1f)   contrast = %.2f K (target %.1f)'
                  % (best[2], TARGET, best[5], TARGET_CONTRAST))
    return 0


if __name__ == '__main__':
    sys.exit(main())
