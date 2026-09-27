#!/usr/bin/env python3
"""Golden-output regression test for the HEXTOR driver.

HEXTOR has no unit tests; its correctness rests on published runs being
reproducible.  This script protects that.  It runs a fixed set of namelists
through a driver binary and compares every output file byte for byte against
a recorded baseline.  A physics or numerics change that is hidden behind a
namelist switch must leave every case here identical with the switch off, and
a change to a diagnostic must change only the files that carry it.

    # record the baseline from the binary that produced the published results
    python tests/regression/run_regression.py --driver /models/hextor/model/driver --record

    # check a rebuilt driver (default: model/driver in this checkout)
    python tests/regression/run_regression.py

Cases are the namelists in tests/regression/cases/*.nml, copied verbatim from
the published configurations (FILLET benchmarks, pre-industrial Earth, THAI
Hab1, SAMOSA).  Baselines live in tests/regression/baseline/<case>/, scratch
runs in tests/regression/runs/ (gitignored).  The driver needs the Intel
runtime, so each run is started through config/machine.sh with
SETVARS_COMPLETED unset first, or setvars.sh silently does nothing.

Exit status is 0 when every case matches its baseline, 1 otherwise.
"""
import argparse
import concurrent.futures
import glob
import os
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
CASES = os.path.join(HERE, 'cases')
BASELINE = os.path.join(HERE, 'baseline')
RUNS = os.path.join(HERE, 'runs')
ENV = os.path.join(ROOT, 'config', 'machine.sh')

# The driver writes data/restart.dat, so data/ is a directory of links to the
# real inputs rather than a link to the directory: parallel cases must not
# write into a shared file.
def prepare_rundir(rundir, modeldir):
    if os.path.isdir(os.path.join(rundir, 'out')):
        shutil.rmtree(os.path.join(rundir, 'out'))
    os.makedirs(os.path.join(rundir, 'out'))
    data = os.path.join(rundir, 'data')
    os.makedirs(data, exist_ok=True)
    for name in os.listdir(os.path.join(modeldir, 'data')):
        if name == 'restart.dat':
            continue
        link = os.path.join(data, name)
        if not os.path.lexists(link):
            os.symlink(os.path.join(modeldir, 'data', name), link)
    rad = os.path.join(rundir, 'radiation')
    if not os.path.lexists(rad):
        os.symlink(os.path.join(modeldir, 'radiation'), rad)


def run_case(case, driver, timeout):
    tag = os.path.splitext(os.path.basename(case))[0]
    modeldir = os.path.dirname(os.path.abspath(driver))
    rundir = os.path.join(RUNS, tag)
    prepare_rundir(rundir, modeldir)
    shutil.copy(case, os.path.join(rundir, 'input.nml'))
    cmd = ['bash', '-c',
           'unset SETVARS_COMPLETED; source "$1" > /dev/null 2>&1; exec "$2"',
           '_', ENV, os.path.abspath(driver)]
    t0 = time.time()
    try:
        proc = subprocess.run(cmd, cwd=rundir, text=True, timeout=timeout,
                              stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        rc, stdout = proc.returncode, proc.stdout
    except subprocess.TimeoutExpired as exc:
        rc, stdout = -1, (exc.stdout or '') + '\n[timeout]\n'
    with open(os.path.join(rundir, 'out', 'stdout.txt'), 'w') as f:
        f.write(stdout)
    return tag, rc, time.time() - t0, rundir


def output_files(outdir):
    return sorted(f for f in os.listdir(outdir)
                  if os.path.getsize(os.path.join(outdir, f)) > 0)


def first_difference(a, b):
    with open(a, errors='replace') as fa, open(b, errors='replace') as fb:
        for n, (la, lb) in enumerate(zip(fa, fb), 1):
            if la != lb:
                return n, la.rstrip(), lb.rstrip()
    return None


def compare(tag, rundir):
    base = os.path.join(BASELINE, tag)
    if not os.path.isdir(base):
        return False, ['no baseline recorded']
    out = os.path.join(rundir, 'out')
    got, want = output_files(out), output_files(base)
    notes, warnings = [], []
    for f in sorted(set(got) | set(want)):
        if f not in got:
            notes.append('%s: missing from run' % f)
        elif f not in want:
            # A new diagnostic file is not a regression; it is listed so the
            # baseline can be re-recorded once the new file is reviewed.
            warnings.append('%s: new file, not in baseline' % f)
        else:
            d = first_difference(os.path.join(out, f), os.path.join(base, f))
            if d is not None:
                notes.append('%s: first difference at line %d\n      run:  %s\n      base: %s'
                             % (f, d[0], d[1][:110], d[2][:110]))
            elif os.path.getsize(os.path.join(out, f)) != os.path.getsize(os.path.join(base, f)):
                notes.append('%s: one file is longer' % f)
    return not notes, notes + warnings


def record(tag, rundir):
    base = os.path.join(BASELINE, tag)
    if os.path.isdir(base):
        shutil.rmtree(base)
    os.makedirs(base)
    out = os.path.join(rundir, 'out')
    for f in output_files(out):
        shutil.copy(os.path.join(out, f), os.path.join(base, f))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--driver', default=os.path.join(ROOT, 'model', 'driver'))
    ap.add_argument('--record', action='store_true',
                    help='write the baseline from this driver instead of checking it')
    ap.add_argument('--cases', nargs='*', help='case names to run (default: all)')
    ap.add_argument('--jobs', type=int, default=4)
    ap.add_argument('--timeout', type=float, default=3600.0)
    args = ap.parse_args()

    cases = sorted(glob.glob(os.path.join(CASES, '*.nml')))
    if args.cases:
        cases = [c for c in cases
                 if os.path.splitext(os.path.basename(c))[0] in args.cases]
    if not cases:
        sys.exit('no cases selected')
    if not os.path.exists(args.driver):
        sys.exit('driver not found: %s' % args.driver)
    os.makedirs(RUNS, exist_ok=True)

    print('driver: %s' % os.path.abspath(args.driver))
    ok_all = True
    with concurrent.futures.ThreadPoolExecutor(args.jobs) as pool:
        futures = [pool.submit(run_case, c, args.driver, args.timeout) for c in cases]
        for fut in futures:
            tag, rc, secs, rundir = fut.result()
            if rc != 0:
                ok_all = False
                print('%-24s FAILED (rc=%s, %.0f s)' % (tag, rc, secs))
                continue
            if args.record:
                record(tag, rundir)
                print('%-24s recorded (%.0f s, %d files)'
                      % (tag, secs, len(output_files(os.path.join(rundir, 'out')))))
                continue
            ok, notes = compare(tag, rundir)
            ok_all &= ok
            print('%-24s %s (%.0f s)' % (tag, 'identical' if ok else 'DIFFERS', secs))
            for n in notes:
                print('    ' + n)
    sys.exit(0 if ok_all else 1)


if __name__ == '__main__':
    main()
