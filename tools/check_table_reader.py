#!/usr/bin/env python3
"""
check_table_reader.py — verify that radiation.f90 interpolates a v2 or v3
table the way the table's own axes say it should.

Builds an independent NumPy interpolator straight from the HDF5 file
(multilinear in log10 p / log10 fCO2 / log10 fCH4 / T for OLR, plus the two
extra axes for the planetary albedo, edge-clamped), runs the same queries
through the Fortran `tabletest` probe, and reports the largest disagreement.

A v2 table (no /ch4 axis) is handled by inserting a single degenerate CH4
level, mirroring what the Fortran reader does, so one code path checks both
formats.

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


def numeric_row(line, nfield):
    """Parse a probe output row, or None if the line is not one.

    radiation_init writes its banner to stdout too, so rows are identified by
    being nfield parseable numbers -- not by field count alone, since the
    CH4-free banner line splits into exactly 8 fields.
    """
    parts = line.split()
    if len(parts) != nfield:
        return None
    try:
        return [float(x) for x in parts]
    except ValueError:
        return None


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


def axes(f):
    """Table axes, with a degenerate CH4 axis inserted for a v2 table."""
    has_ch4 = 'ch4' in f
    ch4 = f['ch4'][:] if has_ch4 else np.array([0.0])
    olrt = f['olr'][:]     # (npre, nco2, [nch4,] ntmp), mW/m^2
    palbt = f['palb'][:]   # (npre, nco2, [nch4,] ntmp, nzen, nsab)
    if not has_ch4:
        olrt = olrt[:, :, None, :]
        palbt = palbt[:, :, None, :, :, :]
    return (f['pressure'][:], f['fco2'][:], ch4, f['temperature'][:],
            f['zenith'][:], f['surfalb'][:], olrt, palbt, has_ch4)


def reference(f, pdry, fco2, fch4, tg0, zy, alb):
    """Independent reimplementation of getOLR/getPALB from the HDF5 axes."""
    pre, fc, ch4, tm, zn, sa, olrt, palbt, has_ch4 = axes(f)

    ip0, ip1, fp = interp_axis(pre, pdry, logaxis=True)
    ic0, ic1, fc_ = interp_axis(fc, fco2, logaxis=True)
    # A degenerate axis never reaches the log10, exactly as axfrac() short
    # circuits on lo == hi; guard it here so a v2 table's 0.0 level is safe.
    if len(ch4) == 1:
        im0, im1, fm = 0, 0, 0.0
    else:
        im0, im1, fm = interp_axis(ch4, fch4, logaxis=True)

    # ---- OLR: clamp T for the interpolation, then extrapolate by power law
    tg_eval = min(max(tg0, tm[0]), tm[-1])
    it0, it1, ft = interp_axis(tm, tg_eval, logaxis=False)
    olr = 0.0
    for ip, wp in ((ip0, 1 - fp), (ip1, fp)):
        for ic, wc in ((ic0, 1 - fc_), (ic1, fc_)):
            for im, wm in ((im0, 1 - fm), (im1, fm)):
                for it, wt in ((it0, 1 - ft), (it1, ft)):
                    w = wp * wc * wm * wt
                    if w:
                        olr += w * olrt[ip, ic, im, it]

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
                for im, wm in ((im0, 1 - fm), (im1, fm)):
                    w = wp * wc * wm
                    if not w:
                        continue
                    nn = (np.log(olrt[ip, ic, im, ea] / olrt[ip, ic, im, eb])
                          / np.log(ta / tb))
                    # radiation.f90 floors the upper exponent at zero: on the
                    # runaway plateau the fit can come out slightly negative,
                    # and extrapolating that would make OLR fall as the model
                    # heats.  The reference has to mirror the module, not the
                    # algebra.
                    if tg0 > tm[-1]:
                        nn = max(nn, 0.0)
                    n_eff += w * nn
        olr *= (tg0 / tg_eval) ** n_eff

    # ---- planetary albedo (T is NOT clamped-then-extrapolated here)
    jt0, jt1, gt = interp_axis(tm, tg0, logaxis=False)
    iz0, iz1, fz = interp_axis(zn, zy, logaxis=False)
    ia0, ia1, fa = interp_axis(sa, alb, logaxis=False)
    palb = 0.0
    for ip, wp in ((ip0, 1 - fp), (ip1, fp)):
        for ic, wc in ((ic0, 1 - fc_), (ic1, fc_)):
            for im, wm in ((im0, 1 - fm), (im1, fm)):
                for it, wt in ((jt0, 1 - gt), (jt1, gt)):
                    for iz, wz in ((iz0, 1 - fz), (iz1, fz)):
                        for ia, wa in ((ia0, 1 - fa), (ia1, fa)):
                            w = wp * wc * wm * wt * wz * wa
                            if w:
                                palb += w * palbt[ip, ic, im, it, iz, ia]

    return olr / 1000.0, palb


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return 2
    table = sys.argv[1]

    with h5py.File(table, 'r') as f:
        pre, fc = f['pressure'][:], f['fco2'][:]
        tm, zn, sa = f['temperature'][:], f['zenith'][:], f['surfalb'][:]
        has_ch4 = 'ch4' in f
        ch4 = f['ch4'][:] if has_ch4 else None
        print('table              : %s' % ('v3 (%d CH4 levels)' % len(ch4)
                                           if has_ch4 else 'v2 (no CH4 axis)'))

        def rand_ch4(mode):
            """A CH4 query value; anything is valid when there is no axis."""
            if not has_ch4:
                return 0.0
            if mode == 'node':
                return rng.choice(ch4)
            if mode == 'inside':
                return 10 ** rng.uniform(np.log10(ch4[0]), np.log10(ch4[-1]))
            return rng.choice([0.0, ch4[0] * 0.01, min(ch4[-1] * 3.0, 0.99)])

        rng = np.random.default_rng(20260827)
        queries = []
        # on-node
        for _ in range(60):
            queries.append((rng.choice(pre), rng.choice(fc), rand_ch4('node'),
                            rng.choice(tm), rng.choice(zn), rng.choice(sa)))
        # off-node, inside the grid
        for _ in range(160):
            queries.append((10 ** rng.uniform(np.log10(pre[0]), np.log10(pre[-1])),
                            10 ** rng.uniform(np.log10(fc[0]), np.log10(fc[-1])),
                            rand_ch4('inside'),
                            rng.uniform(tm[0], tm[-1]),
                            rng.uniform(zn[0], zn[-1]),
                            rng.uniform(sa[0], sa[-1])))
        # outside every axis: clamping and OLR extrapolation
        for _ in range(60):
            queries.append((rng.choice([pre[0] * 0.4, pre[-1] * 2.5]),
                            rng.choice([fc[0] * 0.2, min(fc[-1] * 1.5, 0.99)]),
                            rand_ch4('outside'),
                            rng.choice([tm[0] - 25.0, tm[-1] + 40.0]),
                            rng.choice([0.0, 90.0]),
                            rng.choice([0.0, 1.0])))

        inp = '\n'.join('%.10e %.10e %.10e %.10e %.10e %.10e' % q
                        for q in queries) + '\n'
        res = subprocess.run([PROBE, table], input=inp, capture_output=True,
                             text=True, timeout=300)
        if res.returncode != 0:
            print(res.stdout[-2000:])
            print(res.stderr[-2000:])
            return 1

        got = np.array([r for r in (numeric_row(ln, 8)
                                    for ln in res.stdout.splitlines())
                        if r is not None])
        if len(got) != len(queries):
            print('probe returned %d of %d rows' % (len(got), len(queries)))
            return 1

        max_dolr = max_dpalb = 0.0
        worst = None
        for row in got:
            p, c, m, t, z, a, olr_f, palb_f = row
            olr_r, palb_r = reference(f, p, c, m, t, z, a)
            dolr = abs(olr_f - olr_r) / max(abs(olr_r), 1e-30)
            dpalb = abs(palb_f - palb_r)
            if dolr > max_dolr:
                max_dolr, worst = dolr, (p, c, m, t, z, a, olr_f, olr_r)
            max_dpalb = max(max_dpalb, dpalb)

    print('queries            : %d (on-node, off-node, and out-of-range)' % len(got))
    print('max relative dOLR  : %.3e' % max_dolr)
    print('max absolute dpalb : %.3e' % max_dpalb)
    if worst and max_dolr > 1e-6:
        print('worst OLR point    : p=%.4g fco2=%.4g fch4=%.4g T=%.1f '
              '-> fortran %.6f vs ref %.6f'
              % (worst[0], worst[1], worst[2], worst[3], worst[6], worst[7]))
    ok = max_dolr < 1e-6 and max_dpalb < 1e-9
    print('READER OK' if ok else 'READER MISMATCH')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
