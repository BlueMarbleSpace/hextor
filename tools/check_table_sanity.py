#!/usr/bin/env python3
"""
check_table_sanity.py — physical sanity checks and an inspection figure for a
v2 (pressure-resolved) or v3 (CH4-resolved) HEXTOR radiation table.

A v2 table is handled by inserting a single degenerate CH4 level, so there is
one code path for both formats.

A finished table holds thousands of independently computed grid points, and a
handful failing quietly would show up in an EBM run only as odd climate.  This
checks the properties the physics guarantees, and flags — rather than fails —
the ones that have legitimate exceptions:

  hard   finite everywhere; OLR > 0; planetary albedo in [0, 1]
  hard   OLR rises with surface temperature at fixed (p, fCO2) BELOW the
         runaway plateau.  Above it, OLR asymptotes to the Simpson-Nakajima
         limit and wobbles by ~1 W/m^2 either way; that is the runaway
         greenhouse, not an error, so only decreases well below a column's own
         maximum OLR are treated as defects.
  hard   planetary albedo rises with surface albedo
  soft   OLR falls with CO2 (broken where CO2 condenses out of the profile)
  soft   OLR falls with CH4, where the table has a CH4 axis (same isothermal
         exception: a cold column at t_strato = Ts has no greenhouse effect
         and so no composition dependence at all)
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
        has_ch4 = 'ch4' in f
        ch4 = f['ch4'][:] if has_ch4 else np.array([0.0])
        # Give a v2 table a degenerate CH4 axis so every check below can
        # assume the rank-4 / rank-6 shapes.
        if not has_ch4:
            olr = olr[:, :, None, :]
            palb = palb[:, :, None, :, :, :]

    print('table : %s' % path)
    print('grid  : %d p x %d fCO2 x %d fCH4 x %d T x %d zenith x %d albedo '
          '= %d points'
          % (len(pre), len(fc), len(ch4), len(tm), len(zn), len(sa), olr.size))
    print('star  : %s' % attrs.get('star', '?'))
    print('CH4   : %s' % ('%.1e .. %.1e (%d levels)'
                          % (ch4[0], ch4[-1], len(ch4)) if has_ch4
                          else 'no CH4 axis'))
    if attrs.get('ch4_caveat'):
        print('        caveat: %s' % attrs['ch4_caveat'])
    print()

    failed = False

    # ---- hard checks --------------------------------------------------------
    failed |= report('finite OLR', int(np.sum(~np.isfinite(olr))), olr.size, True)
    failed |= report('finite planetary albedo',
                     int(np.sum(~np.isfinite(palb))), palb.size, True)
    failed |= report('OLR > 0', int(np.sum(olr <= 0)), olr.size, True)
    failed |= report('planetary albedo in [0, 1]',
                     int(np.sum((palb < 0) | (palb > 1))), palb.size, True)

    # A column is "on the plateau" where OLR has stopped responding to
    # temperature: the local slope is a small fraction of that column's
    # steepest slope.  Defining it by proximity to the column maximum instead
    # fails at low CO2, where OLR peaks near 300 K and then settles slightly
    # BELOW that peak, leaving genuinely saturated points looking like errors.
    # This form also copes with the uneven temperature spacing of an extended
    # table, since it works in dOLR/dT rather than per step.
    PLATEAU_FRAC = 0.05
    dT = np.diff(tm)[None, None, None, :]
    slope = np.diff(olr, axis=3) / dT
    maxslope = np.abs(slope).max(axis=3)[:, :, :, None]
    d = slope
    on_plateau = np.abs(slope) < PLATEAU_FRAC * maxslope
    bad = np.argwhere((d <= 0) & ~on_plateau)
    ex = ['p=%.3g fco2=%.3g fch4=%.3g between T=%.0f and %.0f: '
          '%.2f -> %.2f W/m2'
          % (pre[i], fc[j], ch4[m], tm[k], tm[k + 1],
             olr[i, j, m, k], olr[i, j, m, k + 1])
          for i, j, m, k in bad[:3]]
    failed |= report('OLR increases with T below the runaway plateau',
                     len(bad), int(np.sum(~on_plateau)), True, ex)
    npl = int(np.sum((d <= 0) & on_plateau))
    print('  [stat] %-52s %d, max %.2f W/m2 per step'
          % ('OLR dips on the runaway plateau (expected)', npl,
             -d[(d <= 0) & on_plateau].min() if npl else 0.0))

    # Where the plateau sets in, as a function of pressure: the physics SAMOSA
    # is probing, and a quick way to see the table is resolving it.
    print()
    ic = int(np.argmin(np.abs(fc - 4.0e-4)))
    im = 0     # the CH4 floor: the slice that reproduces CH4-free physics
    print('  runaway plateau at fCO2 = %.1e, fCH4 = %.1e (dOLR/dT below %.0f%% '
          'of the column maximum slope):' % (fc[ic], ch4[im], 100 * PLATEAU_FRAC))
    for i in range(0, len(pre), max(1, len(pre) // 6)):
        col = olr[i, ic, im, :]
        hit = int(np.argmax(on_plateau[i, ic, im, :]))
        print('    p = %6.2f bar : from Ts = %3.0f K, OLR -> %.1f W/m2'
              % (pre[i], tm[hit], col.max()))

    d = np.diff(palb, axis=5)
    failed |= report('planetary albedo increases with surface albedo',
                     int(np.sum(d < -1e-12)), d.size, True)

    # ---- soft checks --------------------------------------------------------
    # Tolerance: in the cold isothermal columns (t_strato = min(200 K, Ts), so
    # the atmosphere is at the surface temperature) there is no greenhouse
    # effect at all and OLR is independent of composition to rounding.
    TOL = 0.01   # W/m^2
    d = np.diff(olr, axis=1)
    bad = np.argwhere(d > TOL)
    ex = ['p=%.3g fch4=%.3g T=%.0f between fco2=%.2e and %.2e: '
          '%.2f -> %.2f W/m2'
          % (pre[i], ch4[m], tm[k], fc[j], fc[j + 1],
             olr[i, j, m, k], olr[i, j + 1, m, k])
          for i, j, m, k in bad[:3]]
    report('OLR decreases with CO2 (CO2 condensation excepted)',
           len(bad), d.size, False, ex)

    if len(ch4) > 1:
        d = np.diff(olr, axis=2)
        bad = np.argwhere(d > TOL)
        ex = ['p=%.3g fco2=%.2e T=%.0f between fch4=%.2e and %.2e: '
              '%.2f -> %.2f W/m2'
              % (pre[i], fc[j], tm[k], ch4[m], ch4[m + 1],
                 olr[i, j, m, k], olr[i, j, m + 1, k])
              for i, j, m, k in bad[:3]]
        report('OLR decreases with CH4', len(bad), d.size, False, ex)

        d = np.diff(palb, axis=2)
        frac_down = float(np.mean(d < 0))
        print('  [stat] %-52s %.0f%% falling'
              % ('planetary albedo trend with CH4', 100 * frac_down))

    d = np.diff(olr, axis=0)
    bad = np.argwhere(d > TOL)
    ex = ['fco2=%.2e fch4=%.3g T=%.0f between p=%.3g and %.3g bar: '
          '%.2f -> %.2f W/m2'
          % (fc[j], ch4[m], tm[k], pre[i], pre[i + 1],
             olr[i, j, m, k], olr[i + 1, j, m, k])
          for i, j, m, k in bad[:3]]
    report('OLR decreases with surface pressure', len(bad), d.size, False, ex)

    d = np.diff(palb, axis=4)
    frac_up = float(np.mean(d > 0))
    print('  [stat] %-52s %.0f%% rising' % ('planetary albedo trend toward the limb',
                                            100 * frac_up))

    # ---- scale checks against blackbody emission ---------------------------
    sigma = 5.670374419e-8
    ratio = olr / (sigma * tm[None, None, None, :] ** 4)
    print()
    print('  OLR / sigma T^4 : %.3f .. %.3f  (greenhouse; > 1 would be '
          'unphysical for a grey-ish column)' % (ratio.min(), ratio.max()))
    print('  planetary albedo: %.4f .. %.4f' % (palb.min(), palb.max()))
    print('  OLR             : %.2f .. %.2f W/m2' % (olr.min(), olr.max()))

    make_figure(path, pre, fc, ch4, tm, zn, sa, olr, palb, attrs, has_ch4)

    print()
    print('SANITY FAILED' if failed else 'SANITY OK')
    return 1 if failed else 0


def make_figure(path, pre, fc, ch4, tm, zn, sa, olr, palb, attrs, has_ch4):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib import cm, colors

    ico2 = int(np.argmin(np.abs(fc - 4.0e-4)))   # SAMOSA composition
    ialb = int(np.argmin(np.abs(sa - 0.06)))     # open ocean
    izen = int(np.argmin(np.abs(zn - 60.0)))
    it = int(np.argmin(np.abs(tm - 280.0)))
    ich4 = 0                                     # the CH4 floor
    ip = int(np.argmin(np.abs(pre - 1.0)))       # 1 bar, for the CH4 panels

    norm = colors.LogNorm(vmin=pre[0], vmax=pre[-1])
    smap = cm.ScalarMappable(norm=norm, cmap='viridis')

    # A CH4-resolved table gets two extra panels for the new axis; the first
    # four are unchanged and drawn at the CH4 floor, so they stay comparable
    # with the CH4-free tables' figures.
    nrow = 3 if has_ch4 else 2
    fig, ax = plt.subplots(nrow, 2, figsize=(11, 4.25 * nrow))

    for i, p in enumerate(pre):
        c = smap.to_rgba(p)
        ax[0, 0].plot(tm, olr[i, ico2, ich4, :], color=c, lw=1.4)
        ax[0, 1].plot(zn, palb[i, ico2, ich4, it, :, ialb], color=c, lw=1.4)
        ax[1, 0].plot(sa, palb[i, ico2, ich4, it, izen, :], color=c, lw=1.4)

    ax[0, 0].set_xlabel('surface temperature (K)')
    ax[0, 0].set_ylabel('OLR (W m$^{-2}$)')
    ax[0, 0].set_title('Outgoing longwave, fCO$_2$ = %.0e, fCH$_4$ = %.0e'
                       % (fc[ico2], ch4[ich4]))

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
                             olr[:, ico2, ich4, :], cmap='magma',
                             shading='nearest')
    ax[1, 1].set_xlabel('temperature index')
    ax[1, 1].set_ylabel('pressure index')
    ax[1, 1].set_title('OLR (W m$^{-2}$) over pressure and temperature')
    fig.colorbar(im, ax=ax[1, 1])

    if has_ch4:
        # The CH4 axis itself: OLR and albedo against fCH4, as a family over
        # surface temperature, at 1 bar.
        tnorm = colors.Normalize(vmin=tm[0], vmax=tm[-1])
        tmap = cm.ScalarMappable(norm=tnorm, cmap='plasma')
        for k, t in enumerate(tm):
            c = tmap.to_rgba(t)
            ax[2, 0].semilogx(ch4, olr[ip, ico2, :, k], color=c, lw=1.2)
            ax[2, 1].semilogx(ch4, palb[ip, ico2, :, k, izen, ialb],
                              color=c, lw=1.2)
        ax[2, 0].set_xlabel('CH$_4$ mixing ratio')
        ax[2, 0].set_ylabel('OLR (W m$^{-2}$)')
        ax[2, 0].set_title('Outgoing longwave vs CH$_4$, p = %.2f bar' % pre[ip])
        ax[2, 1].set_xlabel('CH$_4$ mixing ratio')
        ax[2, 1].set_ylabel('planetary albedo')
        ax[2, 1].set_title('Planetary albedo vs CH$_4$, z = %.0f deg, '
                           r'$\alpha_s$ = %.2f' % (zn[izen], sa[ialb]))
        for a in (ax[2, 0], ax[2, 1]):
            a.grid(False)
        cbt = fig.colorbar(tmap, ax=[ax[2, 0], ax[2, 1]], shrink=0.8, pad=0.02)
        cbt.set_label('surface temperature (K)')

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
