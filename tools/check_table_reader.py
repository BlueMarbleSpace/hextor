#!/usr/bin/env python3
"""
check_table_reader.py — verify that radiation.f90 interpolates a v2 table the
way the table's own axes say it should.

Builds an independent NumPy interpolator straight from the HDF5 file
(trilinear in log10 p / log10 fCO2 / T for OLR, pentalinear with the two extra
axes for the planetary albedo, edge-clamped), runs the same queries through the
Fortran `tabletest` probe, and reports the largest disagreement.

Queries deliberately include on-node points, off-node points, and points
outside every axis so the clamping and the OLR power-law extrapolation are
exercised.

    python tools/check_table_reader.py <table.h5>
"""

import os
import subprocess
import sys

import numpy as np
import h5py

HEXTOR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROBE = os.path.join(HEXTOR, 'model', 'radiation', 'tabletest')


def interp_axis(levels, val, logaxis=False):
    """Bracket + fraction, matching radiation.f90's bracket()/axfrac()."""
    n = len(levels)
    if val <= levels[0]:
        return 0, 0, 0.0
    if val >= levels[-1]:
        return n - 1, n - 1, 0.0
    lo = int(np.searchsorted(levels, val, side='right') - 1)
    hi = lo + 1
    if logaxis:
        frac = (np.log10(val) - np.log10(levels[lo])) / \
               (np.log10(levels[hi]) - np.log10(levels[lo]))
    else:
        frac = (val - levels[lo]) / (levels[hi] - levels[lo])
    return lo, hi, frac


def reference(f, pdry, fco2, tg0, zy, alb):
    """Independent reimplementation of getOLR/getPALB from the HDF5 axes."""
    pre = f['pressure'][:]
    fc = f['fco2'][:]
    tm = f['temperature'][:]
    zn = f['zenith'][:]
    sa = f['surfalb'][:]
    olrt = f['olr'][:]     # (npre, nco2, ntmp), mW/m^2
    palbt = f['palb'][:]   # (npre, nco2, ntmp, nzen, nsab)

    ip0, ip1, fp = interp_axis(pre, pdry, logaxis=True)
    ic0, ic1, fc_ = interp_axis(fc, fco2, logaxis=True)

    # ---- OLR: clamp T for the interpolation, then extrapolate by power law
    tg_eval = min(max(tg0, tm[0]), tm[-1])
    it0, it1, ft = interp_axis(tm, tg_eval, logaxis=False)
    olr = 0.0
    for ip, wp in ((ip0, 1 - fp), (ip1, fp)):
        for ic, wc in ((ic0, 1 - fc_), (ic1, fc_)):
            for it, wt in ((it0, 1 - ft), (it1, ft)):
                w = wp * wc * wt
                if w:
                    olr += w * olrt[ip, ic, it]

    if tg0 < tm[0] or tg0 > tm[-1]:
        if tg0 > tm[-1]:
            ea, eb = -1, -2
            ta, tb = tm[-1], tm[-2]
        else:
            ea, eb = 1, 0
            ta, tb = tm[1], tm[0]
        n_eff = 0.0
        for ip, wp in ((ip0, 1 - fp), (ip1, fp)):
            for ic, wc in ((ic0, 1 - fc_), (ic1, fc_)):
                n = np.log(olrt[ip, ic, ea] / olrt[ip, ic, eb]) / np.log(ta / tb)
                # radiation.f90 floors the upper exponent at zero: on the
                # runaway plateau the fit can come out slightly negative, and
                # extrapolating that would make OLR fall as the model heats.
                # The reference has to mirror the module, not the algebra.
                if tg0 > tm[-1]:
                    n = max(n, 0.0)
                n_eff += wp * wc * n
        olr *= (tg0 / tg_eval) ** n_eff

    # ---- planetary albedo (T is NOT clamped-then-extrapolated here)
    jt0, jt1, gt = interp_axis(tm, tg0, logaxis=False)
    iz0, iz1, fz = interp_axis(zn, zy, logaxis=False)
    ia0, ia1, fa = interp_axis(sa, alb, logaxis=False)
    palb = 0.0
    for ip, wp in ((ip0, 1 - fp), (ip1, fp)):
        for ic, wc in ((ic0, 1 - fc_), (ic1, fc_)):
            for it, wt in ((jt0, 1 - gt), (jt1, gt)):
                for iz, wz in ((iz0, 1 - fz), (iz1, fz)):
                    for ia, wa in ((ia0, 1 - fa), (ia1, fa)):
                        w = wp * wc * wt * wz * wa
                        if w:
                            palb += w * palbt[ip, ic, it, iz, ia]

    return olr / 1000.0, palb


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return 2
    table = sys.argv[1]

    with h5py.File(table, 'r') as f:
        pre, fc = f['pressure'][:], f['fco2'][:]
        tm, zn, sa = f['temperature'][:], f['zenith'][:], f['surfalb'][:]

        rng = np.random.default_rng(20260827)
        queries = []
        # on-node
        for _ in range(60):
            queries.append((rng.choice(pre), rng.choice(fc), rng.choice(tm),
                            rng.choice(zn), rng.choice(sa)))
        # off-node, inside the grid
        for _ in range(160):
            queries.append((10 ** rng.uniform(np.log10(pre[0]), np.log10(pre[-1])),
                            10 ** rng.uniform(np.log10(fc[0]), np.log10(fc[-1])),
                            rng.uniform(tm[0], tm[-1]),
                            rng.uniform(zn[0], zn[-1]),
                            rng.uniform(sa[0], sa[-1])))
        # outside every axis: clamping and OLR extrapolation
        for _ in range(60):
            queries.append((rng.choice([pre[0] * 0.4, pre[-1] * 2.5]),
                            rng.choice([fc[0] * 0.2, min(fc[-1] * 1.5, 0.99)]),
                            rng.choice([tm[0] - 25.0, tm[-1] + 40.0]),
                            rng.choice([0.0, 90.0]),
                            rng.choice([0.0, 1.0])))

        inp = '\n'.join('%.10e %.10e %.10e %.10e %.10e' % q for q in queries) + '\n'
        res = subprocess.run([PROBE, table], input=inp, capture_output=True,
                             text=True, timeout=300)
        if res.returncode != 0:
            print(res.stdout[-2000:])
            print(res.stderr[-2000:])
            return 1

        got = np.array([[float(x) for x in ln.split()]
                        for ln in res.stdout.splitlines() if len(ln.split()) == 7])
        if len(got) != len(queries):
            print('probe returned %d of %d rows' % (len(got), len(queries)))
            return 1

        max_dolr = max_dpalb = 0.0
        worst = None
        for row in got:
            p, c, t, z, a, olr_f, palb_f = row
            olr_r, palb_r = reference(f, p, c, t, z, a)
            dolr = abs(olr_f - olr_r) / max(abs(olr_r), 1e-30)
            dpalb = abs(palb_f - palb_r)
            if dolr > max_dolr:
                max_dolr, worst = dolr, (p, c, t, z, a, olr_f, olr_r)
            max_dpalb = max(max_dpalb, dpalb)

    print('queries            : %d (on-node, off-node, and out-of-range)' % len(got))
    print('max relative dOLR  : %.3e' % max_dolr)
    print('max absolute dpalb : %.3e' % max_dpalb)
    if worst and max_dolr > 1e-6:
        print('worst OLR point    : p=%.4g fco2=%.4g T=%.1f -> fortran %.6f vs ref %.6f'
              % (worst[0], worst[1], worst[2], worst[5], worst[6]))
    ok = max_dolr < 1e-6 and max_dpalb < 1e-9
    print('READER OK' if ok else 'READER MISMATCH')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
