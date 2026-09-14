#!/usr/bin/env python3
"""
migrate_cache_names.py — rename an index-named radiation-table cache to the
value-named scheme.

make_radiation_table.py used to name each cached column by its position on the
grid (`col_<ip>_<ic>.txt`).  That is only meaningful alongside the axes it was
built with: regrid the CO2 axis and the same names refer to different mixing
ratios, so a resume would load the wrong columns and report a full cache hit.
Cache files are now named by their physical values instead.

This script performs the one-time rename, taking the axes from the HDF5 table
the cache was built for -- the table is the only record of what the indices
meant, which is why the axes are not guessed from the file names.

    python tools/migrate_cache_names.py model/radiation/radiation_N2_CO2_Sun_p.h5
    python tools/migrate_cache_names.py TABLE.h5 --cache DIR --apply

Without --apply it reports what it would do and changes nothing.
"""

import argparse
import os
import sys

import h5py
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from make_radiation_table import cache_name


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('table', help='HDF5 table whose axes the cache was built on')
    ap.add_argument('--cache', default=None,
                    help='cache directory (default: <table>.cache)')
    ap.add_argument('--apply', action='store_true',
                    help='actually rename; without it, only report')
    args = ap.parse_args()

    cache = args.cache or (args.table + '.cache')
    if not os.path.isdir(cache):
        print('no such cache directory: %s' % cache)
        return 1

    with h5py.File(args.table, 'r') as f:
        pressures = f['pressure'][:]
        fco2 = f['fco2'][:]
        ch4 = f['ch4'][:] if 'ch4' in f else None

    print('%s' % args.table)
    print('  axes: %d pressures, %d fco2%s'
          % (len(pressures), len(fco2),
             ', %d ch4' % len(ch4) if ch4 is not None else ''))

    plan, missing, already = [], [], 0
    for ip, p in enumerate(pressures):
        for ic, c in enumerate(fco2):
            for im, m in enumerate(ch4 if ch4 is not None else [None]):
                if ch4 is None:
                    old = 'col_%03d_%03d.txt' % (ip, ic)
                else:
                    old = 'col_%03d_%03d_%03d.txt' % (ip, ic, im)
                new = cache_name(p, c, m)
                if os.path.exists(os.path.join(cache, new)):
                    already += 1
                elif os.path.exists(os.path.join(cache, old)):
                    plan.append((old, new))
                else:
                    missing.append(old)

    print('  %d to rename, %d already value-named, %d absent'
          % (len(plan), already, len(missing)))

    # Anything left over is a file the table's own axes do not account for --
    # a leftover from a different grid.  Renaming those would be guesswork, so
    # they are listed and left alone.
    known = set(o for o, _ in plan)
    stray = [f for f in sorted(os.listdir(cache))
             if f.startswith('col_') and f.endswith('.txt')
             and f not in known
             and not f.startswith('col_p')]
    if stray:
        print('  %d cache files do not correspond to this table\'s axes and are'
              ' left untouched:' % len(stray))
        for f in stray[:10]:
            print('     %s' % f)
        if len(stray) > 10:
            print('     ... and %d more' % (len(stray) - 10))

    if not args.apply:
        print('  (dry run -- pass --apply to rename)')
        return 0

    for old, new in plan:
        os.rename(os.path.join(cache, old), os.path.join(cache, new))
    print('  renamed %d files' % len(plan))
    return 0


if __name__ == '__main__':
    sys.exit(main())
