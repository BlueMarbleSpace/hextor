#!/usr/bin/env python3
"""
samosa_compare.py — score HEXTOR's SAMOSA runs against the other submissions.

Reads one or more `samosa_summary.csv` files written by tools/run_samosa.py and
compares them, case by case, with the models that submitted to the same
protocol.  Two references are built in, and they answer different questions:

  ExoCAM     a 3-D GCM with clouds and dynamics.  The closest thing the
             intercomparison has to a reference answer, and the only one that
             constrains the day-night contrast.  Differences against it mix
             radiative transfer, transport, clouds and sea ice.

  ExoColumn  a 1-D radiative-convective column using the SAME radiation core
             HEXTOR's lookup tables are built from (ExoRT n68equiv).  It has
             no horizontal dimension, so Tmax = Tmin = Tglob and it says
             nothing about contrast -- but for the global mean it isolates
             everything EXCEPT the radiative transfer.  A HEXTOR-ExoColumn gap
             is transport, clouds or surface albedo, not spectroscopy.
             Caveat: it is cloud-free with a uniform surface albedo of 0.2736
             standing in for cloud shortwave, where HEXTOR uses ocnalb 0.06
             plus a uniform `cloudir` offset to OLR.  The two cloud stand-ins
             are not the same approximation, and that is most of what the cold
             cases show.

None of this is a calibration target: HEXTOR's (d0, cloudir) are fitted at THAI
Hab1, around a 2600 K star, so every number here is out of sample.

    python tools/samosa_compare.py samosa [other_dir ...]
    python tools/samosa_compare.py "D constant=samosa" "D ~ p=/tmp/perbar"
"""

import argparse
import csv
import glob
import math
import os
import re
import sys

SAMOSA_DATA = '/models/data/samosa'


def read_exocam(root):
    """Global means from ExoCAM's analysis.py digest, contrast from the files.

    The digest's TS is area-weighted; the extrema come from the maps, where
    weighting does not matter.
    """
    out = {}
    path = os.path.join(root, 'exocam', 'output.txt')
    if os.path.exists(path):
        for ln in open(path):
            m = re.match(r'\s*\d+\s+samosa(\d+)\.cam\S*\s+([\d.]+)\s+([\d.]+)'
                         r'\s+([\d.]+)\s+([\d.]+)', ln)
            if m:
                out[int(m.group(1))] = dict(T=float(m.group(2)),
                                            ice=float(m.group(3)),
                                            alb=float(m.group(4)),
                                            olr=float(m.group(5)))
    try:
        from netCDF4 import Dataset
        import numpy as np
        for f in glob.glob(os.path.join(root, 'exocam', 'samosa*.cam.h0.avg.nc')):
            c = int(re.search(r'samosa(\d+)', os.path.basename(f)).group(1))
            if c not in out:
                continue
            ts = np.squeeze(Dataset(f).variables['TS'][:])
            out[c]['Tmin'], out[c]['Tmax'] = float(ts.min()), float(ts.max())
    except ImportError:
        pass
    return out


def read_protocol_dat(path):
    """A submission in the protocol's own global-output format."""
    out = {}
    if not os.path.exists(path):
        return out
    for ln in open(path):
        if ln.startswith('#') or not ln.strip():
            continue
        p = ln.split()
        if len(p) < 10:
            continue
        try:
            c = int(p[0])
        except ValueError:
            continue
        out[c] = dict(T=float(p[3]), Tmax=float(p[4]), Tmin=float(p[5]),
                      olr=float(p[6]), asr=float(p[7]), fsdn=float(p[8]),
                      alb=1.0 - float(p[7]) / float(p[8]) if float(p[8]) else
                      float('nan'),
                      ice=1.0 - float(p[12]) if p[12] not in ('-999',) else
                      float('nan'))
    return out


def read_hextor(path):
    """Equilibrium cases from a run_samosa summary, warm start where both ran."""
    rows = {}
    with open(os.path.join(path, 'samosa_summary.csv')) as f:
        for r in csv.DictReader(f):
            c = int(r['case'])
            if c in rows and r['init'] != 'warm':
                continue
            rows[c] = r
    out = {}
    for c, r in rows.items():
        if r.get('state') != 'equilibrium':
            out[c] = None
            continue
        out[c] = dict(T=float(r['T_global']), Tmin=float(r['T_min']),
                      Tmax=float(r['T_max']),
                      olr=float(r['OLR_global']), asr=float(r['ASR_global']))
    return out


def score(hex_runs, ref, key):
    """RMSE and bias of each HEXTOR run against a reference, on `key`."""
    out = []
    for label, runs in hex_runs:
        e = []
        for c, rv in sorted(ref.items()):
            hv = runs.get(c)
            if hv is None or (key == 'T' and 'T' not in rv) or \
               (key != 'T' and 'Tmax' not in rv):
                continue
            a = hv['T'] if key == 'T' else hv['Tmax'] - hv['Tmin']
            b = rv['T'] if key == 'T' else rv['Tmax'] - rv['Tmin']
            if b != b:
                continue
            e.append(a - b)
        if e:
            out.append((label, len(e), sum(e) / len(e),
                        math.sqrt(sum(x * x for x in e) / len(e)),
                        max(abs(x) for x in e)))
    return out


def table(title, hex_runs, ref, key, note=''):
    print('=' * 78)
    print(title + ('   [%s]' % note if note else ''))
    print('=' * 78)
    hdr = '%5s %6s %9s' % ('case', 'S', 'ref')
    for label, _ in hex_runs:
        hdr += ' | %-20s' % label[:20]
    print(hdr)
    for c, rv in sorted(ref.items()):
        b = rv['T'] if key == 'T' else rv.get('Tmax', float('nan')) - rv.get('Tmin', float('nan'))
        if b != b:
            continue
        row = '%5d %6s %9.1f' % (c, rv.get('S', ''), b)
        for label, runs in hex_runs:
            hv = runs.get(c)
            if hv is None:
                row += ' | %-20s' % '        runaway'
            else:
                a = hv['T'] if key == 'T' else hv['Tmax'] - hv['Tmin']
                row += ' | %8.1f (%+8.1f)' % (a, a - b)
        print(row)
    print('-' * 78)
    for label, n, bias, rmse, mx in score(hex_runs, ref, key):
        print('  %-24s n=%2d   bias %+7.1f K   RMSE %6.1f K   |max| %6.1f K'
              % (label, n, bias, rmse, mx))
    print()


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('rundirs', nargs='+',
                    help='directories holding a samosa_summary.csv, each '
                         'optionally prefixed "label=" to name it in the table')
    ap.add_argument('--data', default=SAMOSA_DATA,
                    help='root of the other submissions')
    args = ap.parse_args()

    dirs, labels = [], []
    for spec in args.rundirs:
        label, _, d = spec.rpartition('=')
        dirs.append(d)
        labels.append(label or os.path.basename(d.rstrip('/')))
    hex_runs = [(l, read_hextor(d)) for l, d in zip(labels, dirs)]

    inst = {}
    with open(os.path.join(dirs[0], 'samosa_summary.csv')) as f:
        for r in csv.DictReader(f):
            inst[int(r['case'])] = r['instellation']

    cam = read_exocam(args.data)
    col = read_protocol_dat(os.path.join(args.data, 'exocolumn',
                                         'global_output_ExoColumn_a2736.dat'))
    for d in (cam, col):
        for c in d:
            d[c]['S'] = inst.get(c, '')

    if cam:
        table('Global mean surface temperature vs ExoCAM (3-D GCM, clouds)',
              hex_runs, cam, 'T')
        if any('Tmax' in v for v in cam.values()):
            table('Day-night contrast vs ExoCAM', hex_runs, cam, 'contrast',
                  'Tmax - Tmin')
    if col:
        table('Global mean surface temperature vs ExoColumn (1-D, same ExoRT)',
              hex_runs, col, 'T',
              'isolates transport + clouds, not radiative transfer')
    return 0


if __name__ == '__main__':
    sys.exit(main())
