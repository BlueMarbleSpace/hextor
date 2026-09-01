#!/usr/bin/env python3
"""
samosa_submit.py — package HEXTOR results in the format the SAMOSA
intercomparison expects.

Produces, for a completed run directory from tools/run_samosa.py:

  global_output_HEXTOR.dat   the protocol's global summary table, one row per
                             case, in the column order of the distributed
                             template
  samosa<N>_HEXTOR.nc        per case, the Table 5 two-dimensional diagnostics
                             on the same 46 x 72 grid the ExoCAM submission
                             uses, with CESM variable names so the files drop
                             into the existing analysis
  README_HEXTOR.txt          what an energy balance model can and cannot supply

HEXTOR solves one dimension: the angle from the sub-stellar point.  Its fields
are therefore axisymmetric about that point, and the maps here are that profile
rotated onto the comparison grid, with cos(theta) = cos(lat) cos(lon) for a
sub-stellar point at (0, 0).  The maps carry no more information than the
18-belt profile, which is written into the same files as a native
one-dimensional variable so nothing is hidden by the interpolation.

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
from datetime import datetime, timezone

import numpy as np

HEXTOR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MODEL = 'HEXTOR'
CONTACT = 'Jacob Haqq-Misra (jacob@bmsis.org)'

# Comparison grid, matching the ExoCAM submission (4 deg x 5 deg).
NLAT, NLON = 46, 72
ICETEMP = 263.15          # driver.f:251
FILL = float('nan')


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


def to_map(theta, values, lat, lon):
    """Rotate a profile in angle-from-sub-stellar onto a lat/lon grid."""
    LON, LAT = np.meshgrid(np.radians(lon), np.radians(lat))
    cos_theta = np.clip(np.cos(LAT) * np.cos(LON), -1.0, 1.0)
    grid_theta = np.degrees(np.arccos(cos_theta))
    idx = np.argsort(theta)
    return np.interp(grid_theta, theta[idx], np.asarray(values)[idx])


def write_netcdf(path, sample, inst, pres, d, lat, lon):
    import netCDF4 as nc

    fice = ice_fraction(d['T'])
    insol = inst * np.clip(np.cos(np.radians(d['theta'])), 0.0, None)

    with nc.Dataset(path, 'w', format='NETCDF4_CLASSIC') as f:
        f.createDimension('latitude', len(lat))
        f.createDimension('longitude', len(lon))
        f.createDimension('theta', len(d['theta']))

        def var(name, dims, data, units, long_name):
            v = f.createVariable(name, 'f8', dims, fill_value=FILL)
            v.units = units
            v.long_name = long_name
            v[:] = data
            return v

        var('latitude', ('latitude',), lat, 'degrees_north', 'latitude')
        var('longitude', ('longitude',), lon, 'degrees_east', 'longitude')

        # Two-dimensional diagnostics (Table 5), axisymmetric about (0, 0).
        var('TS', ('latitude', 'longitude'), to_map(d['theta'], d['T'], lat, lon),
            'K', 'Surface temperature')
        var('FLUT', ('latitude', 'longitude'),
            to_map(d['theta'], d['olr'], lat, lon),
            'W/m2', 'Upwelling longwave flux at top of model')
        var('FLNT', ('latitude', 'longitude'),
            to_map(d['theta'], d['olr'], lat, lon),
            'W/m2', 'Net longwave flux at top of model')
        var('FSNTOA', ('latitude', 'longitude'),
            to_map(d['theta'], d['asr'], lat, lon),
            'W/m2', 'Net solar flux at top of atmosphere')
        var('SOLIN', ('latitude', 'longitude'), to_map(d['theta'], insol, lat, lon),
            'W/m2', 'Downward solar flux at top of atmosphere')
        var('ALBEDO', ('latitude', 'longitude'),
            to_map(d['theta'], d['albedo'], lat, lon),
            'fraction', 'Planetary albedo')
        var('ICEFRAC', ('latitude', 'longitude'), to_map(d['theta'], fice, lat, lon),
            'fraction', 'Fraction of sfc area covered by sea-ice')
        var('OCNFRAC', ('latitude', 'longitude'),
            to_map(d['theta'], 1.0 - fice, lat, lon),
            'fraction', 'Fraction of sfc area covered by open ocean')

        # The native solution, so the maps can be checked against it.
        var('theta', ('theta',), d['theta'], 'degrees',
            'angle from sub-stellar point')
        var('TS_theta', ('theta',), d['T'], 'K',
            'Surface temperature on the native HEXTOR belts')
        var('FLUT_theta', ('theta',), d['olr'], 'W/m2',
            'Outgoing longwave on the native HEXTOR belts')
        var('FSNTOA_theta', ('theta',), d['asr'], 'W/m2',
            'Absorbed stellar radiation on the native HEXTOR belts')
        var('area_theta', ('theta',), d['w'], 'fraction',
            'Fractional area of each HEXTOR belt')

        f.model = MODEL
        f.contact = CONTACT
        f.sample = sample
        f.instellation_W_m2 = inst
        f.surface_pressure_bar = pres
        f.description = ('HEXTOR is a one-dimensional energy balance model in '
                         'the tidally locked coordinate. The lat/lon fields are '
                         'the 18-belt profile rotated about the sub-stellar '
                         'point at (0, 0) and carry no additional information; '
                         'the native profile is included as the *_theta '
                         'variables.')
        f.omitted = ('No vertical dimension, water vapor, clouds or ice '
                     'thickness: the corresponding Table 5 diagnostics are not '
                     'provided.')
        f.created = datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')


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

    lat = np.linspace(-90.0, 90.0, NLAT)
    lon = np.linspace(0.0, 360.0 - 360.0 / NLON, NLON)

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
        write_netcdf(os.path.join(args.outdir, 'samosa%d_%s.nc' % (c, MODEL)),
                     c, inst, pres, d, lat, lon)

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

    # ---- README --------------------------------------------------------------
    with open(os.path.join(args.outdir, 'README_%s.txt' % MODEL), 'w') as f:
        f.write('SAMOSA submission: %s\n%s\n\n' % (MODEL, '=' * 30))
        f.write('Contact: %s\n\n' % CONTACT)
        f.write('Files\n-----\n')
        f.write('  global_output_%s.dat   global summary, template column order\n' % MODEL)
        f.write('  samosa<N>_%s.nc        per-case fields for the cases that\n' % MODEL)
        f.write('                          reach a steady climate\n\n')
        f.write('What this model provides\n------------------------\n')
        f.write('HEXTOR is a one-dimensional energy balance model. It solves for\n')
        f.write('surface temperature against the angle from the sub-stellar point\n')
        f.write('on 18 belts, so its fields are axisymmetric about that point. The\n')
        f.write('latitude/longitude maps in the netCDF files are that profile\n')
        f.write('rotated onto the 46 x 72 grid used by the ExoCAM submission, with\n')
        f.write('cos(theta) = cos(lat) cos(lon) and the sub-stellar point at (0, 0).\n')
        f.write('They contain no information beyond the profile, which is included\n')
        f.write('in the same files as the *_theta variables.\n\n')
        f.write('Supplied : TS, FLUT, FLNT, FSNTOA, SOLIN, ALBEDO, ICEFRAC, OCNFRAC\n')
        f.write('Not supplied: every diagnostic requiring a vertical dimension,\n')
        f.write('water vapor, clouds or ice thickness. The model has no vertical\n')
        f.write('structure and no cloud physics; cloud radiative effects enter only\n')
        f.write('as a uniform offset to the outgoing longwave, fitted to the THAI\n')
        f.write('ensemble, which is not a cloud field and is not reported as one.\n\n')
        if skipped:
            f.write('Cases without a steady state\n----------------------------\n')
            for c, st in skipped:
                f.write('  sample %-2d : %s\n' % (c, st))
            f.write('\nThe protocol allows incipient runaway cases to be omitted or\n')
            f.write('reported at the last stable state; they are omitted here.\n')

    print('wrote %d case files to %s' % (len(reported), args.outdir))
    print('  global_output_%s.dat  (%d cases reported, %d omitted)'
          % (MODEL, len(reported), len(skipped)))
    for c, st in skipped:
        print('     omitted sample %-2d : %s' % (c, st))
    return 0


if __name__ == '__main__':
    sys.exit(main())
