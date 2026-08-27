#!/usr/bin/env python3
"""
check_table_interpolation.py — measure the interpolation error of the
radiation table's angular and surface-albedo axes.

The table is only as good as its grid spacing.  This script computes a fine
reference block of planetary albedo over (zenith x surface albedo) with
ExoColumn, then asks how well bilinear interpolation on the coarse table nodes
(MU_NODES x ALBEDOS from make_radiation_table) reproduces it — including at
the off-node points a running EBM will actually query.

Reported: max and RMS error in planetary albedo, and the equivalent error in
absorbed shortwave for a representative instellation.

    python tools/check_table_interpolation.py
"""

import os
import subprocess
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from make_radiation_table import (build_namelist, EXOCOL_ROOT, MU_NODES,
                                  ALBEDOS, t_strato_for)

WORK = os.environ.get('INTERP_WORK',
                      '/tmp/claude-1000/-hugespace-models-hextor/'
                      '380c2fc7-6424-4334-8480-5cc50da5b2ff/scratchpad/interp')
EXE = os.path.join(EXOCOL_ROOT, 'run', 'exocol_sweepts200.exe')
STAR = 'blackbody_3000K_n68.nc'

# Fine reference axes: the coarse nodes PLUS their midpoints, so the coarse
# table is exact at its own nodes and the error is measured exactly where a
# running EBM interpolates.  (Must stay within the Fortran MAX_SWEEP = 64.)
def _refine(a):
    mid = 0.5 * (a[:-1] + a[1:])
    return np.sort(np.concatenate([a, mid]))[::-1] if a[0] > a[-1] else \
           np.sort(np.concatenate([a, mid]))

FINE_MU = _refine(MU_NODES)
FINE_ALB = _refine(ALBEDOS)

CASES = [(0.1, [240.0, 300.0]), (1.0, [240.0, 300.0, 360.0]), (10.0, [300.0, 360.0])]
FCO2 = 4.0e-4
INSTELLATION = 1000.0   # W/m2, for converting an albedo error into a flux error


def run_fine(pdry, temps):
    os.makedirs(os.path.join(WORK, 'iofiles'), exist_ok=True)
    nml = build_namelist(pdry, FCO2, temps[0], STAR, 'iofiles/sweep.txt',
                         mus=FINE_MU, albs=FINE_ALB)
    ts_line = '  sweep_ts        = ' + ', '.join('%.4f' % t for t in temps) + '\n'
    tst_line = ('  sweep_t_strato  = '
                + ', '.join('%.4f' % t_strato_for(t) for t in temps) + '\n')
    nml = nml.replace('  sweep_outfile', ts_line + tst_line + '  sweep_outfile')
    with open(os.path.join(WORK, 'exocol_config.nml'), 'w') as f:
        f.write(nml)
    res = subprocess.run([EXE], cwd=WORK, capture_output=True, text=True, timeout=3600)
    if res.returncode != 0:
        print(res.stdout[-500:])
        raise RuntimeError('ExoColumn failed for p=%g' % pdry)

    rows = []
    for line in open(os.path.join(WORK, 'iofiles/sweep.txt')):
        if line.startswith('#'):
            continue
        p = line.split()
        if len(p) == 10:
            rows.append([float(x) for x in p])
    return np.array(rows)


def bilinear(xg, yg, zg, x, y):
    """Bilinear interpolation on an ascending grid, clamped at the edges."""
    i = np.clip(np.searchsorted(xg, x) - 1, 0, len(xg) - 2)
    j = np.clip(np.searchsorted(yg, y) - 1, 0, len(yg) - 2)
    tx = np.clip((x - xg[i]) / (xg[i + 1] - xg[i]), 0.0, 1.0)
    ty = np.clip((y - yg[j]) / (yg[j + 1] - yg[j]), 0.0, 1.0)
    return ((1 - tx) * (1 - ty) * zg[i, j] + tx * (1 - ty) * zg[i + 1, j]
            + (1 - tx) * ty * zg[i, j + 1] + tx * ty * zg[i + 1, j + 1])


def main():
    coarse_zen = np.degrees(np.arccos(MU_NODES))
    order = np.argsort(coarse_zen)
    coarse_zen = coarse_zen[order]

    print('coarse grid : %d zenith x %d albedo' % (len(MU_NODES), len(ALBEDOS)))
    print('fine grid   : %d zenith x %d albedo' % (len(FINE_MU), len(FINE_ALB)))
    print()
    print('%-8s %-6s  %-10s %-10s  %-12s' %
          ('p_dry', 'Ts', 'max|dpalb|', 'rms|dpalb|', 'max dASR'))
    print('-' * 56)

    worst = 0.0
    for pdry, temps in CASES:
        rows = run_fine(pdry, temps)
        for ts in temps:
            sub = rows[np.abs(rows[:, 0] - ts) < 1e-6]
            # Fine reference block, indexed [zenith ascending, albedo]
            fine = np.full((len(FINE_MU), len(FINE_ALB)), np.nan)
            for r in sub:
                iz = int(np.argmin(np.abs(FINE_MU - r[4])))
                ia = int(np.argmin(np.abs(FINE_ALB - r[5])))
                fine[iz, ia] = r[9]
            fine_zen = np.degrees(np.arccos(FINE_MU))
            fo = np.argsort(fine_zen)
            fine_zen = fine_zen[fo]
            fine = fine[fo, :]

            # Coarse table built from the fine block at the coarse nodes.
            coarse = np.full((len(coarse_zen), len(ALBEDOS)), np.nan)
            for iz, z in enumerate(coarse_zen):
                for ia, a in enumerate(ALBEDOS):
                    coarse[iz, ia] = bilinear(fine_zen, FINE_ALB, fine, z, a)

            # Error at every fine point that is NOT a coarse node.
            errs = []
            for iz, z in enumerate(fine_zen):
                for ia, a in enumerate(FINE_ALB):
                    on_node = (np.min(np.abs(coarse_zen - z)) < 1e-6 and
                               np.min(np.abs(ALBEDOS - a)) < 1e-9)
                    if on_node:
                        continue
                    errs.append(bilinear(coarse_zen, ALBEDOS, coarse, z, a) - fine[iz, ia])
            errs = np.array(errs)
            mx = np.max(np.abs(errs))
            rms = np.sqrt(np.mean(errs ** 2))
            # An albedo error maps into absorbed shortwave as S*mu*dalb; use the
            # global-mean geometry (mu -> 1/4 of the substellar flux).
            dasr = mx * INSTELLATION * 0.25
            worst = max(worst, mx)
            print('%-8.2f %-6.0f  %-10.5f %-10.5f  %-6.2f W/m2'
                  % (pdry, ts, mx, rms, dasr), flush=True)

    print('-' * 56)
    print('worst-case planetary albedo interpolation error: %.5f' % worst)
    print('  -> %.2f W/m2 of absorbed shortwave at S = %.0f W/m2 (global mean)'
          % (worst * INSTELLATION * 0.25, INSTELLATION))


if __name__ == '__main__':
    main()
