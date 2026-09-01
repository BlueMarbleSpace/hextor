#!/usr/bin/env python3
"""
samosa_submit.py — write HEXTOR's SAMOSA results as the protocol's global
output table.

Produces global_output_HEXTOR.dat from a completed run directory: one row per
case, in the column order of the distributed template, with the model
description and the omitted cases recorded in the header.

Diagnostics that require a vertical dimension, water vapor, clouds or ice
thickness are left as NaN rather than being invented; the protocol asks
participating models to supply what applies to them.

    python tools/samosa_submit.py [--rundir samosa] [--outdir samosa]
"""

import argparse
import csv
import math
import os
import sys

import numpy as np

HEXTOR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MODEL = 'HEXTOR'
CONTACT = 'Jacob Haqq-Misra (jacob@bmsis.org)'

ICETEMP = 263.15          # driver.f:251


def belt_weights(coords_deg, nbelts=18):
    """Fractional area of each HEXTOR belt (equally spaced in the coordinate)."""
    half = math.pi / (2 * nbelts)
    w = [abs(math.sin(math.radians(c) + half) - math.sin(math.radians(c) - half))
         for c in coords_deg]
    t = sum(w)
    return np.array([x / t for x in w])


def ice_fraction(T):
    """Sea-ice fraction from belt temperature, following driver.f:806-810."""
    T = np.asarray(T, dtype=float)
    f = np.where(T >= 273.15, 0.0,
                 np.where(T < ICETEMP, 1.0, 1.0 - np.exp((T - 273.15) / 10.0)))
    return f


def read_case(case_dir):
    """The ZONAL STATISTICS block, ordered from sub-stellar to anti-stellar."""
    path = os.path.join(case_dir, 'out', 'model.out')
    rows, seen = [], False
    for line in open(path):
        if line.startswith('ZONAL STATISTICS'):
            seen = True
            continue
        if seen:
            p = line.split()
            if len(p) != 9:
                if rows:
                    break
                continue
            try:
                rows.append([float(x) for x in p])
            except ValueError:
                if rows:
                    break
    if not rows:
        return None
    z = np.array(rows)
    order = np.argsort(-z[:, 0])          # sub-stellar first
    z = z[order]
    theta = 90.0 - z[:, 0]                # 0 = sub-stellar, 180 = anti-stellar
    return dict(theta=theta, coord=z[:, 0], T=z[:, 1], albedo=z[:, 6],
                olr=z[:, 7], asr=z[:, 8], w=belt_weights(z[:, 0]))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--rundir', default=os.path.join(HEXTOR, 'samosa'))
    ap.add_argument('--outdir', default=os.path.join(HEXTOR, 'samosa'))
    ap.add_argument('--notes', default='')
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    summary = os.path.join(args.rundir, 'samosa_summary.csv')
    if not os.path.exists(summary):
        print('no samosa_summary.csv in %s' % args.rundir)
        return 1
    rows = [r for r in csv.DictReader(open(summary))]
    # One row per case, preferring the warm start where both were run.
    by_case = {}
    for r in rows:
        c = int(r['case'])
        if c not in by_case or r['init'] == 'warm':
            by_case[c] = r

    reported, skipped, cfg = [], [], None
    for c in sorted(by_case):
        r = by_case[c]
        state = r.get('state', '')
        if state != 'equilibrium':
            skipped.append((c, state))
            continue
        case_dir = os.path.join(args.rundir, 'case_%02d_%s' % (c, r['init']))
        d = read_case(case_dir)
        if d is None:
            skipped.append((c, 'no zonal output'))
            continue
        cfg = r
        inst, pres = float(r['instellation']), float(r['ps_bar'])
        fice = ice_fraction(d['T'])
        insol = inst * np.clip(np.cos(np.radians(d['theta'])), 0.0, None)
        w = d['w']
        reported.append(dict(
            sample=c, inst=inst, pres=pres,
            tglob=float(np.sum(w * d['T'])),
            tmax=float(np.max(d['T'])), tmin=float(np.min(d['T'])),
            olr=float(np.sum(w * d['olr'])), asr=float(np.sum(w * d['asr'])),
            fsdn=float(np.sum(w * insol)), fnet=float(np.sum(w * d['olr'])),
            ocnfrac=float(np.sum(w * (1.0 - fice)))))

    # ---- global output table, in the template's column order -----------------
    path = os.path.join(args.outdir, 'global_output_%s.dat' % MODEL)
    with open(path, 'w') as f:
        f.write('# SAMOSA Global Output File\n#\n')
        f.write('# Model name: %s (energy balance model)\n' % MODEL)
        f.write('# Name of primary contact: %s\n' % CONTACT)
        f.write('# Notes:\n')
        f.write('#   One-dimensional EBM in the tidally locked coordinate, 18 belts.\n')
        f.write('#   Radiative transfer by interpolation in lookup tables computed\n')
        f.write('#   with ExoColumn/ExoRT n68equiv for a 3000 K blackbody, indexed by\n')
        f.write('#   dry surface pressure, CO2 mixing ratio, surface temperature,\n')
        f.write('#   zenith angle and surface albedo. Clear sky: HEXTOR has no\n')
        f.write('#   clouds, and no vertical dimension, water vapor or ice thickness,\n')
        f.write('#   so Qstrat, Qmass, Icethick, Cldliq, Cldice and Cldfrac are NaN.\n')
        if cfg is not None:
            f.write('#   Diffusion D = %s W/m^2/K (constant); cloud infrared\n'
                    % cfg.get('d0', '?'))
            f.write('#   correction to OLR = %s W/m^2, both calibrated against the\n'
                    % cfg.get('cloudir', '?'))
            f.write('#   THAI Hab 1 GCM ensemble for TRAPPIST-1 e.\n')
        f.write('#   Fsdn is the downward stellar flux at the top of the atmosphere.\n')
        f.write('#   Fnet is the net longwave flux at the top of the atmosphere,\n')
        f.write('#   which for this model equals the outgoing longwave.\n')
        if args.notes:
            f.write('#   %s\n' % args.notes)
        if skipped:
            f.write('#\n# Cases with no steady climate state, omitted per the\n')
            f.write('# protocol\'s allowance for incipient runaway:\n')
            for c, st in skipped:
                f.write('#   sample %-2d : %s\n' % (c, st))
        f.write('#\n')
        f.write('# Sample Inst Pres Tglob Tmax Tmin OLR ASR Fsdn Fnet Qstrat '
                'Qmass Ocnfrac Icethick Cldliq Cldice Cldfrac\n')
        for r in reported:
            f.write('%d %.1f %.2f %.3f %.3f %.3f %.3f %.3f %.3f %.3f '
                    'NaN NaN %.4f NaN NaN NaN NaN\n'
                    % (r['sample'], r['inst'], r['pres'], r['tglob'], r['tmax'],
                       r['tmin'], r['olr'], r['asr'], r['fsdn'], r['fnet'],
                       r['ocnfrac']))

    print('wrote %s' % path)
    print('  %d cases reported, %d omitted' % (len(reported), len(skipped)))
    for c, st in skipped:
        print('     omitted sample %-2d : %s' % (c, st))
    return 0


if __name__ == '__main__':
    sys.exit(main())
