#!/usr/bin/env python3
"""
merge_radiation_tables.py — join two v2 radiation tables along the temperature
axis.

Extending a table's temperature range does not require regenerating it.  The
pressure, CO2 and CH4 axes are the outer loop of make_radiation_table.py, so
the extra temperatures can be generated on their own and stitched on here,
which costs only the fraction of the work the new levels represent.

Both v2 (no CH4 axis) and v3 (CH4-resolved) tables are accepted, but the two
inputs must be the same format: temperature is a different array axis in each,
so merging across formats would silently misalign the data.

Every other axis, and the physical assumptions recorded in the file
attributes, must agree between the two inputs; a mismatch is an error rather
than something to paper over, since the result would silently mix conventions.

    python tools/merge_radiation_tables.py base.h5 extension.h5 -o merged.h5
"""

import argparse
import sys

import numpy as np
import h5py

# Axes that must agree between the two inputs.  'ch4' is optional, so it is
# checked only when present in both -- see the format comparison below.
AXES = ['pressure', 'fco2', 'zenith', 'surfalb']

# Attributes that describe the physics rather than the bookkeeping.  These must
# match, or the two halves of the table were not computed the same way.
PHYSICAL_ATTRS = ['star', 'water', 't_strato', 'p_top_ratio', 'h2o_eos',
                  'h2o_continuum', 'co2_condense', 'clouds', 'olr_units',
                  'pressure_units']


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('base')
    ap.add_argument('extension')
    ap.add_argument('-o', '--out', required=True)
    args = ap.parse_args()

    with h5py.File(args.base, 'r') as a, h5py.File(args.extension, 'r') as b:
        # Temperature sits at a different array axis in v2 and v3, so the two
        # inputs have to carry the same format before anything is joined.
        ach4, bch4 = 'ch4' in a, 'ch4' in b
        if ach4 != bch4:
            print('format mismatch: %s has %s CH4 axis, %s has %s'
                  % (args.base, 'a' if ach4 else 'no',
                     args.extension, 'one' if bch4 else 'none'))
            return 1
        axes = AXES + (['ch4'] if ach4 else [])
        # A CH4 axis shifts temperature from array axis 2 to axis 3.
        tax = 3 if ach4 else 2

        for ax in axes:
            if a[ax].shape != b[ax].shape or not np.allclose(a[ax][:], b[ax][:]):
                print('axis mismatch: %s differs between the two tables' % ax)
                return 1
        for key in PHYSICAL_ATTRS:
            va, vb = a.attrs.get(key), b.attrs.get(key)
            if va != vb:
                print('assumption mismatch: %s is %r in the base and %r in the '
                      'extension' % (key, va, vb))
                return 1

        ta, tb = a['temperature'][:], b['temperature'][:]
        overlap = np.intersect1d(ta, tb)
        if overlap.size:
            print('temperature axes overlap at %s; refusing to guess which to '
                  'keep' % overlap)
            return 1

        temperature = np.concatenate([ta, tb])
        order = np.argsort(temperature)
        temperature = temperature[order]

        olr = np.concatenate([a['olr'][:], b['olr'][:]], axis=tax)
        palb = np.concatenate([a['palb'][:], b['palb'][:]], axis=tax)
        olr = np.take(olr, order, axis=tax)
        palb = np.take(palb, order, axis=tax)

        with h5py.File(args.out, 'w') as o:
            o.create_dataset('olr', data=olr)
            o.create_dataset('palb', data=palb)
            o.create_dataset('temperature', data=temperature)
            for ax in axes:
                o.create_dataset(ax, data=a[ax][:])
            for k, v in a.attrs.items():
                o.attrs[k] = v
            o.attrs['merged_from'] = '%s + %s' % (args.base, args.extension)
            o.attrs['temperature_note'] = (
                'base %.0f-%.0f K extended to %.0f K'
                % (ta.min(), ta.max(), tb.max()))

    print('wrote %s' % args.out)
    print('  temperature axis: %d levels, %.0f to %.0f K'
          % (len(temperature), temperature[0], temperature[-1]))
    print('  %s' % '  '.join('%.0f' % t for t in temperature))
    return 0


if __name__ == '__main__':
    sys.exit(main())
