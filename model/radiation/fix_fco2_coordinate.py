#!/usr/bin/env python3
"""
fix_fco2_coordinate.py -- one-time repair of the fCO2 coordinate column in
the HEXTOR radiation lookup tables.

The stored fCO2 column was found to be pCO2/(X + pCO2) with X = pi for the
2600 K table and X = pi/2 for the Sun table, instead of the true mixing
ratio pCO2/(1 bar N2 + pCO2). (Verified against the total-pressure column,
which cleanly satisfies P = 1 + pCO2 for pCO2 = 1e-4..10 bar; the factor is
believed to be related to the pi factor removed from the instellation
function, see Haqq-Misra 2026 Appendix B.)

This script recovers pCO2 = X*f/(1-f) from the stored coordinate (full
relative precision at all magnitudes, unlike the rounded P column) and
rewrites the column as f_true = pCO2/(1 + pCO2), equivalently
f_true = X*f / (1 + (X-1)*f), in both the 'olr' and 'palb' datasets.

Originals are backed up as <name>.h5.bak-pi-coordinate. Idempotent: skips
files whose coordinate maximum already exceeds 0.9.
"""
import shutil
import numpy as np
import h5py

TABLES = {
    'radiation_N2_CO2_2600K.h5': np.pi,
    'radiation_N2_CO2_Sun.h5':   np.pi / 2.0,
}

for name, X in TABLES.items():
    with h5py.File(name, 'r') as f:
        fmax = f['olr'][:, 1].max()
    if fmax > 0.9:
        print(f'{name}: coordinate max {fmax:.4f} already > 0.9, skipping')
        continue

    backup = name + '.bak-pi-coordinate'
    shutil.copy2(name, backup)
    print(f'{name}: backed up to {backup}')

    with h5py.File(name, 'r+') as f:
        for ds in ('olr', 'palb'):
            col = f[ds][:, 1]
            f_true = X * col / (1.0 + (X - 1.0) * col)
            data = f[ds][:]
            data[:, 1] = f_true
            f[ds][...] = data
            print(f'  {ds}: f range {col.min():.4e}..{col.max():.4f} -> '
                  f'{f_true.min():.4e}..{f_true.max():.4f}')
