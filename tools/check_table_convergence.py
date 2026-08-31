#!/usr/bin/env python3
"""
check_table_convergence.py — verify that the radiation table's numerical
settings (vertical resolution and model top) are converged before committing
to a full table generation.

Runs a small corner set of (p_dry, Ts) at fixed CO2 through ExoColumn with
different layer counts and p_top ratios, and reports the spread in OLR and
planetary albedo.  Anything much below ~1 W/m^2 in OLR and ~0.002 in albedo is
comfortably inside the physical uncertainty of the table itself.

    python tools/check_table_convergence.py
"""

import os
import subprocess
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from make_radiation_table import build_namelist, parse_sweep, EXOCOL_ROOT

WORK = os.environ.get('CONVERGENCE_WORK',
                      '/tmp/claude-1000/-hugespace-models-hextor/380c2fc7-6424-4334-8480-5cc50da5b2ff/scratchpad/conv')
STAR = 'blackbody_3000K_n68.nc'
FCO2 = 4.0e-4
MU = np.array([0.5])
ALB = np.array([0.3])

PRESSURES = [0.1, 1.0, 10.0]
TEMPS = [200.0, 260.0, 300.0, 360.0, 400.0]

VARIANTS = [
    ('pver70  ptop/1e-5', 'run/exocol_sweep70.exe',  1.0e-5),
    ('pver200 ptop/1e-5', 'run/exocol_sweep200.exe', 1.0e-5),
    ('pver70  ptop/1e-6', 'run/exocol_sweep70.exe',  1.0e-6),
    ('pver70  ptop/1e-4', 'run/exocol_sweep70.exe',  1.0e-4),
]


def run(exe, pdry, ts, ptop_ratio):
    os.makedirs(os.path.join(WORK, 'iofiles'), exist_ok=True)
    nml = build_namelist(pdry, FCO2, ts, STAR, 'iofiles/sweep.txt',
                         mus=MU, albs=ALB, ptop_ratio=ptop_ratio)
    with open(os.path.join(WORK, 'exocol_config.nml'), 'w') as f:
        f.write(nml)
    res = subprocess.run([os.path.join(EXOCOL_ROOT, exe)], cwd=WORK,
                         capture_output=True, text=True, timeout=600)
    if res.returncode != 0:
        return None, None
    head, rows = parse_sweep(os.path.join(WORK, 'iofiles/sweep.txt'))
    if rows.shape[0] != 1:
        return None, None
    return rows[0, 2], rows[0, 5]


def main():
    print('%-8s %-6s | ' % ('p_dry', 'Ts')
          + ' | '.join('%-20s' % v[0] for v in VARIANTS))
    print('%-8s %-6s | ' % ('[bar]', '[K]')
          + ' | '.join('%-9s %-10s' % ('OLR', 'palb') for _ in VARIANTS))
    print('-' * (17 + 23 * len(VARIANTS)))

    max_dolr = 0.0
    max_dalb = 0.0
    for pdry in PRESSURES:
        for ts in TEMPS:
            cells = []
            olrs, albs = [], []
            for _, exe, ratio in VARIANTS:
                olr, palb = run(exe, pdry, ts, ratio)
                if olr is None:
                    cells.append('%-9s %-10s' % ('FAIL', ''))
                else:
                    olrs.append(olr)
                    albs.append(palb)
                    cells.append('%-9.3f %-10.5f' % (olr, palb))
            if len(olrs) > 1:
                dolr = max(olrs) - min(olrs)
                dalb = max(albs) - min(albs)
                max_dolr = max(max_dolr, dolr)
                max_dalb = max(max_dalb, dalb)
            print('%-8.2f %-6.0f | ' % (pdry, ts) + ' | '.join(cells), flush=True)

    print('-' * (17 + 23 * len(VARIANTS)))
    print('max spread across variants:  dOLR = %.3f W/m2   dpalb = %.5f'
          % (max_dolr, max_dalb))


if __name__ == '__main__':
    main()
