#!/usr/bin/env python3
"""
check_table_sanity.py — physical sanity checks and an inspection figure for a
v2 (pressure-resolved) HEXTOR radiation table.

A finished table holds thousands of independently computed grid points, and a
handful failing quietly would show up in an EBM run only as odd climate.  This
checks the properties the physics guarantees, and flags — rather than fails —
the ones that have legitimate exceptions:

  hard   finite everywhere; OLR > 0; planetary albedo in [0, 1]
  hard   OLR rises with surface temperature at fixed (p, fCO2)
  hard   planetary albedo rises with surface albedo
  soft   OLR falls with CO2 (broken where CO2 condenses out of the profile)
  soft   OLR falls with surface pressure (pressure broadening; weakest, and
         reversible, where the column is already water-dominated)
  stat   how the planetary albedo behaves toward the limb — reported, not
         checked.  For a red host there is almost no Rayleigh scattering to
         brighten the limb, while the slant path absorbs more of the surface
         return, so the albedo commonly FALLS toward the limb; under a solar
         spectrum it rises.  Which way it goes is a property of the star, not
         an error.

Also writes a four-panel inspection figure (PNG + PDF) next to the table.

    python tools/check_table_sanity.py model/radiation/radiation_N2_CO2_3000K_p.h5
"""

import os
import sys

import numpy as np
import h5py


def report(name, bad, total, hard, examples=()):
    tag = 'FAIL' if (hard and bad) else ('note' if bad else 'ok  ')
    print('  [%s] %-52s %d/%d' % (tag, name, bad, total))
    for e in examples[:3]:
        print('         %s' % e)
    return (hard and bad > 0)


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return 2
    path = sys.argv[1]

    with h5py.File(path, 'r') as f:
        pre = f['pressure'][:]
        fc = f['fco2'][:]
        tm = f['temperature'][:]
        zn = f['zenith'][:]
        sa = f['surfalb'][:]
        olr = f['olr'][:] / 1000.0     # -> W/m^2
        palb = f['palb'][:]
        attrs = dict(f.attrs)

    print('table : %s' % path)
    print('grid  : %d p x %d fCO2 x %d T x %d zenith x %d albedo = %d points'
          % (len(pre), len(fc), len(tm), len(zn), len(sa), olr.size))
    print('star  : %s' % attrs.get('star', '?'))
    print()

    failed = False

    # ---- hard checks --------------------------------------------------------
    failed |= report('finite OLR', int(np.sum(~np.isfinite(olr))), olr.size, True)
    failed |= report('finite planetary albedo',
                     int(np.sum(~np.isfinite(palb))), palb.size, True)
    failed |= report('OLR > 0', int(np.sum(olr <= 0)), olr.size, True)
    failed |= report('planetary albedo in [0, 1]',
                     int(np.sum((palb < 0) | (palb > 1))), palb.size, True)

    d = np.diff(olr, axis=2)
    bad = np.argwhere(d <= 0)
    ex = ['p=%.3g fco2=%.3g between T=%.0f and %.0f: %.2f -> %.2f W/m2'
          % (pre[i], fc[j], tm[k], tm[k + 1], olr[i, j, k], olr[i, j, k + 1])
          for i, j, k in bad[:3]]
    failed |= report('OLR increases with surface temperature',
                     len(bad), d.size, True, ex)

    d = np.diff(palb, axis=4)
    failed |= report('planetary albedo increases with surface albedo',
                     int(np.sum(d < -1e-12)), d.size, True)

    # ---- soft checks --------------------------------------------------------
    d = np.diff(olr, axis=1)
    bad = np.argwhere(d > 0)
    ex = ['p=%.3g T=%.0f between fco2=%.2e and %.2e: %.2f -> %.2f W/m2'
          % (pre[i], tm[k], fc[j], fc[j + 1], olr[i, j, k], olr[i, j + 1, k])
          for i, j, k in bad[:3]]
    report('OLR decreases with CO2 (CO2 condensation excepted)',
           len(bad), d.size, False, ex)

    d = np.diff(olr, axis=0)
    bad = np.argwhere(d > 0)
    ex = ['fco2=%.2e T=%.0f between p=%.3g and %.3g bar: %.2f -> %.2f W/m2'
          % (fc[j], tm[k], pre[i], pre[i + 1], olr[i, j, k], olr[i + 1, j, k])
          for i, j, k in bad[:3]]
    report('OLR decreases with surface pressure', len(bad), d.size, False, ex)

    d = np.diff(palb, axis=3)
    frac_up = float(np.mean(d > 0))
    print('  [stat] %-52s %.0f%% rising' % ('planetary albedo trend toward the limb',
                                            100 * frac_up))

    # ---- scale checks against blackbody emission ---------------------------
    sigma = 5.670374419e-8
    ratio = olr / (sigma * tm[None, None, :] ** 4)
    print()
    print('  OLR / sigma T^4 : %.3f .. %.3f  (greenhouse; > 1 would be '
          'unphysical for a grey-ish column)' % (ratio.min(), ratio.max()))
    print('  planetary albedo: %.4f .. %.4f' % (palb.min(), palb.max()))
    print('  OLR             : %.2f .. %.2f W/m2' % (olr.min(), olr.max()))

    make_figure(path, pre, fc, tm, zn, sa, olr, palb, attrs)

    print()
    print('SANITY FAILED' if failed else 'SANITY OK')
    return 1 if failed else 0


def make_figure(path, pre, fc, tm, zn, sa, olr, palb, attrs):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib import cm, colors

    ico2 = int(np.argmin(np.abs(fc - 4.0e-4)))   # SAMOSA composition
    ialb = int(np.argmin(np.abs(sa - 0.06)))     # open ocean
    izen = int(np.argmin(np.abs(zn - 60.0)))
    it = int(np.argmin(np.abs(tm - 280.0)))

    norm = colors.LogNorm(vmin=pre[0], vmax=pre[-1])
    smap = cm.ScalarMappable(norm=norm, cmap='viridis')

    fig, ax = plt.subplots(2, 2, figsize=(11, 8.5))

    for i, p in enumerate(pre):
        c = smap.to_rgba(p)
        ax[0, 0].plot(tm, olr[i, ico2, :], color=c, lw=1.4)
        ax[0, 1].plot(zn, palb[i, ico2, it, :, ialb], color=c, lw=1.4)
        ax[1, 0].plot(sa, palb[i, ico2, it, izen, :], color=c, lw=1.4)

    ax[0, 0].set_xlabel('surface temperature (K)')
    ax[0, 0].set_ylabel('OLR (W m$^{-2}$)')
    ax[0, 0].set_title('Outgoing longwave, fCO$_2$ = %.0e' % fc[ico2])

    ax[0, 1].set_xlabel('solar zenith angle (deg)')
    ax[0, 1].set_ylabel('planetary albedo')
    ax[0, 1].set_title('Planetary albedo vs zenith, T = %.0f K, '
                       r'$\alpha_s$ = %.2f' % (tm[it], sa[ialb]))

    ax[1, 0].set_xlabel('surface albedo')
    ax[1, 0].set_ylabel('planetary albedo')
    ax[1, 0].set_title('Planetary albedo vs surface, T = %.0f K, z = %.0f deg'
                       % (tm[it], zn[izen]))

    # OLR over the pressure x temperature plane at the SAMOSA composition
    im = ax[1, 1].pcolormesh(np.arange(len(tm)), np.arange(len(pre)),
                             olr[:, ico2, :], cmap='magma', shading='nearest')
    ax[1, 1].set_xlabel('temperature index')
    ax[1, 1].set_ylabel('pressure index')
    ax[1, 1].set_title('OLR (W m$^{-2}$) over pressure and temperature')
    fig.colorbar(im, ax=ax[1, 1])

    for a in (ax[0, 0], ax[0, 1], ax[1, 0]):
        a.grid(False)
    cb = fig.colorbar(smap, ax=[ax[0, 0], ax[0, 1], ax[1, 0]], shrink=0.6,
                      pad=0.02)
    cb.set_label('dry surface pressure (bar)')

    fig.suptitle('HEXTOR radiation table: %s' % os.path.basename(path),
                 fontweight='normal')
    base = os.path.splitext(path)[0]
    for ext in ('png', 'pdf'):
        fig.savefig('%s_inspect.%s' % (base, ext), dpi=150,
                    bbox_inches='tight')
    plt.close(fig)
    print()
    print('  wrote %s_inspect.png / .pdf' % base)


if __name__ == '__main__':
    sys.exit(main())
