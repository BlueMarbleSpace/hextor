#!/usr/bin/env python3
"""Run the FILLET benchmarks and experiments with HEXTOR and assemble the archive.

Replaces the csh scripts under fillet/<exp>/ (June 2024), which had no error
checking -- a failed case silently repeated the previous case's line under a
new number -- and ran Experiment 2a on the Experiment 1a grid.  Every case is
run in its own scratch directory, in parallel, from a namelist written here,
so the whole configuration is in one place and is printed into the headers.

    python tools/run_fillet.py --exp all                    # the 4.3.0 submission
    python tools/run_fillet.py --exp ben1 ben2 ben3 --jobs 2
    python tools/run_fillet.py --exp exp4 --no-diffcons --nstepyr 0 \
        --label 'published numerics' --outdir fillet_published

The defaults are the configuration of the September 2026 re-file (HEXTOR
4.3.0): Benchmark 1 tuned through cloudir alone (9.3695 W/m2), untuned
Benchmarks 2/3 and experiments at D = 0.5, flux-conservative diffusion, 730
steps per orbit, fluxcnvg 1e-3 with cnvgcycle 4.  Spelled out in full:

    python tools/run_fillet.py --exp all --label 'FILLET re-file, September 2026' \
        --cloudir-ben1 9.3695 --diffcons --nstepyr 730

About seven minutes at four jobs.  Repack the tarball afterwards:

    tar czf fillet/fillet_hextor.tar.gz -C fillet/Results hextor

Outputs, under --outdir (default fillet/):
    Results/hextor/<exp>/global_output.dat        the projectcuisines/fillet layout
    Results/hextor/ben?/case_0/lat_output.dat     (exp3/exp4 split into _cold/_warm)
    latfiles/<exp>/case_<n>.dat                   per-case latitude files (kept here,
                                                  the archive README asks for none)
    logs/<exp>_convergence.log                    per-case convergence record
    runs/<exp>/case_<n>/                          scratch (gitignored), resume-safe

Grids follow FILLET Protocol v1.1: Experiments 1/2 on 19 instellations x 10
obliquities, 1a/2a on 19 and 15 semi-major axes (1a 0.875-1.1 au, 2a 0.8-0.975
au, year scaling as a^1.5), Experiment 3 on 57 instellations per branch,
Experiment 4 on 50 CO2 values per branch.  Warm starts are a uniform 288 K,
cold starts 233 K, every case from that state (no continuation).

The driver needs the Intel runtime: each run is started through
config/machine.sh with SETVARS_COMPLETED unset first.
"""
import argparse
import concurrent.futures
import json
import math
import os
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
ENV = os.path.join(ROOT, 'config', 'machine.sh')
RADFILE = './radiation/radiation_N2_CO2_Sun.h5'
AU_CM = 1.495978707e13
SOLAR = 1361.0           # W/m2, the protocol's S_earth
SECONDS_PER_DAY = 86400.0

OBLIQUITIES = [0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0]

LONG_NAME = {
    'ben1': 'Benchmark 1', 'ben2': 'Benchmark 2', 'ben3': 'Benchmark 3',
    'exp1': 'Experiment 1', 'exp1a': 'Experiment 1a',
    'exp2': 'Experiment 2', 'exp2a': 'Experiment 2a',
    'exp3_cold': 'Experiment 3 (cold-start branch)',
    'exp3_warm': 'Experiment 3 (warm-start branch)',
    'exp4_cold': 'Experiment 4 (cold-start branch)',
    'exp4_warm': 'Experiment 4 (warm-start branch)',
}
ALL_EXPS = list(LONG_NAME)


def frange(start, stop, step):
    """Inclusive float grid, rounded so 0.8125 prints as 0.8125."""
    n = int(round((stop - start) / step)) + 1
    return [round(start + i * step, 6) for i in range(n)]


def cases_for(exp, twarm, tcold):
    """List of dicts: inst (S_earth), obl (deg), xco2 (ppm), a (au), tinit (K)."""
    cs = []

    def add(inst, obl, xco2, a, tinit):
        cs.append(dict(inst=inst, obl=obl, xco2=xco2, a=a, tinit=tinit))

    if exp == 'ben1':
        add(1.0, 23.5, 280.0, 1.0, twarm)
    elif exp == 'ben2':
        add(1.0, 23.5, 280.0, 1.0, twarm)
    elif exp == 'ben3':
        add(1.0, 60.0, 280.0, 1.0, twarm)
    elif exp in ('exp1', 'exp2'):
        lo, hi = (0.8, 1.25) if exp == 'exp1' else (1.05, 1.5)
        for s in frange(lo, hi, 0.025):
            for obl in OBLIQUITIES:
                add(s, obl, 280.0, 1.0, twarm if exp == 'exp1' else tcold)
    elif exp in ('exp1a', 'exp2a'):
        # semi-major axis descending, so instellation ascends as in exp1/exp2
        lo, hi = (0.875, 1.1) if exp == 'exp1a' else (0.8, 0.975)
        for a in reversed(frange(lo, hi, 0.0125)):
            for obl in OBLIQUITIES:
                add(round(1.0 / a**2, 6), obl, 280.0, a,
                    twarm if exp == 'exp1a' else tcold)
    elif exp == 'exp3_cold':
        for s in frange(0.8, 1.5, 0.0125):
            add(s, 23.5, 280.0, 1.0, tcold)
    elif exp == 'exp3_warm':
        for s in reversed(frange(0.8, 1.5, 0.0125)):
            add(s, 23.5, 280.0, 1.0, twarm)
    elif exp in ('exp4_cold', 'exp4_warm'):
        xs = [10.0 ** (-6.0 + 5.0 * i / 49.0) * 1.0e6 for i in range(50)]  # ppm
        if exp == 'exp4_warm':
            xs = list(reversed(xs))
        for x in xs:
            add(1.0, 23.5, x, 1.0, tcold if exp == 'exp4_cold' else twarm)
    else:
        raise ValueError(exp)
    for i, c in enumerate(cs):
        c['case'] = i
    return cs


def orbital_period_days(a_au, msun=1.9891e33, grav=6.6732e-8):
    a = a_au * AU_CM
    return 2.0 * math.pi / (math.sqrt(grav * msun) * a ** -1.5) / SECONDS_PER_DAY


def namelist(exp, case, args):
    """The namelist text for one case.  Benchmark 1 keeps its published Earth
    configuration (geography, Fresnel ocean, HEXTOR heat capacities, the
    published transport); everything else is the protocol's Table 4."""
    ben1 = exp == 'ben1'
    a_cm = case['a'] * AU_CM
    solarcon = SOLAR / case['a'] ** 2 if exp in ('exp1a', 'exp2a') else SOLAR
    relsolcon = case['inst'] if exp not in ('exp1a', 'exp2a') else 1.0
    lines = ['&ebm',
             '   seasons      = .true.,',
             '   do_manualseasons = .false.,',
             '   tempinit     = %.1f,' % case['tinit'],
             '   resfile      = 0,',
             '   tend         = 7.e11,',
             '   dt           = %.1f,' % args.dt,
             '   rot          = 7.27e-5,',
             '   a            = %.9E,' % a_cm,
             '   ecc          = 0.0,',
             '   peri         = 0.0,',
             '   obl          = %.1f,' % case['obl'],
             '   d0           = %.4f,' % (args.d0_ben1 if ben1 else args.d0),
             '   diffadj      = %s,' % ('.true.' if (ben1 or args.diffadj) else '.false.'),
             '   pg0          = 1.0,',
             '   fco2         = %.8e,' % (case['xco2'] * 1.0e-6),
             '   fh2          = 0.,',
             '   ocean        = %s,' % ('0.7' if ben1 else '0.75'),
             '   igeog        = %d,' % (1 if ben1 else 6),
             '   constheatcap = .false.,',
             '   heatcap      = 2.00e6,',
             '   yrstep       = 1,',
             '   msun         = 1.9891e33,',
             '   iterhalt     = .true.,',
             '   fluxcnvg     = %.3e,' % args.fluxcnvg,
             '   icelinetemp  = %.2f,' % args.icelinetemp,
             '   cnvgcycle    = %d,' % args.cnvgcycle,
             '   do_longitudinal = .false.,']
    if not ben1:
        lines += ['   cl           = 1.e7,',
                  '   cw           = 4.e8,',
                  '   ci           = 1.e7,']
    for extra in switch_lines(args):
        lines.append('   %s,' % extra)
    lines += ['   fillet       = .true. /',
              '',
              '&radiation',
              '   relsolcon    = %.6f,' % relsolcon,
              '   solarcon     = %.6f,' % solarcon,
              '   radparam     = 3,',
              '   groundalb    = 0.30,',
              '   snowalb      = %s,' % ('0.663' if ben1 else '0.6'),
              '   oceanalbconst = %s,' % ('.false.' if ben1 else '.true.'),
              '   ocnalb       = %s,' % ('0.06' if ben1 else '0.2'),
              '   landsnowfrac = 1.00,',
              '   cloudir      = %.4f,' % (args.cloudir_ben1 if ben1 else args.cloudir),
              '   fcloud       = 0.0,',
              '   cloudalb     = .false.,',
              '   soladj       = .false.,',
              '   linalb       = .false.,',
              '   linrad       = .false.,',
              "   radfile      = '%s' /" % RADFILE,
              '',
              '&co2cycle',
              '   do_cs_cycle  = .false.,',
              '   outgassing   = 7.0,',
              '   weathering   = 7.0,',
              '   betaexp      = 0.50,',
              '   kact         = 0.09,',
              '   krun         = 0.045 /',
              '',
              '&h2cycle',
              '   do_h2_cycle  = .false.,',
              '   h2outgas     = 1.0e12 /',
              '',
              '&stochastic',
              '   do_stochastic = .false.,',
              '   noisevar      = 0.334 /',
              '']
    return '\n'.join(lines)


def switch_lines(args):
    """The numerics-switch lines of &ebm: --diffcons and --nstepyr, then any
    --extra-ebm lines.  An extra line that sets the same variable replaces the
    flag, so the spelling the archive was first generated with,
    --extra-ebm 'diffcons = .true.' --extra-ebm 'nstepyr = 730', still works."""
    extras = [e.strip().rstrip(',') for e in args.extra_ebm]

    def given(name):
        return any(e.replace(' ', '').startswith(name + '=') for e in extras)

    lines = []
    if args.diffcons and not given('diffcons'):
        lines.append('diffcons = .true.')
    if args.nstepyr and not given('nstepyr'):
        lines.append('nstepyr = %d' % args.nstepyr)
    return lines + extras


def prepare_rundir(rundir, modeldir):
    os.makedirs(os.path.join(rundir, 'out'), exist_ok=True)
    data = os.path.join(rundir, 'data')
    os.makedirs(data, exist_ok=True)
    for name in os.listdir(os.path.join(modeldir, 'data')):
        if name == 'restart.dat':
            continue          # written by the driver; keep it per case
        link = os.path.join(data, name)
        if not os.path.lexists(link):
            os.symlink(os.path.join(modeldir, 'data', name), link)
    rad = os.path.join(rundir, 'radiation')
    if not os.path.lexists(rad):
        os.symlink(os.path.join(modeldir, 'radiation'), rad)


def series_has_nan(path):
    try:
        with open(path) as f:
            return any('nan' in line.lower() for line in f)
    except OSError:
        return False


def run_case(exp, case, args, modeldir):
    tag = 'case_%d' % case['case']
    rundir = os.path.join(args.outdir, 'runs', exp, tag)
    result_file = os.path.join(rundir, 'result.json')
    if os.path.exists(result_file) and not args.force:
        with open(result_file) as f:
            return json.load(f)
    prepare_rundir(rundir, modeldir)
    for f in os.listdir(os.path.join(rundir, 'out')):
        os.remove(os.path.join(rundir, 'out', f))
    with open(os.path.join(rundir, 'input.nml'), 'w') as f:
        f.write(namelist(exp, case, args))
    cmd = ['bash', '-c',
           'unset SETVARS_COMPLETED; source "$1" > /dev/null 2>&1; exec "$2"',
           '_', ENV, os.path.abspath(args.driver)]
    t0 = time.time()
    try:
        proc = subprocess.run(cmd, cwd=rundir, text=True, timeout=args.timeout,
                              stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        rc, stdout = proc.returncode, proc.stdout
    except subprocess.TimeoutExpired as exc:
        rc, stdout = -1, (exc.stdout or '') + '\n[timeout]\n'
    with open(os.path.join(rundir, 'stdout.txt'), 'w') as f:
        f.write(stdout)
    res = dict(exp=exp, case=case, rc=rc, seconds=round(time.time() - t0, 1))
    out = os.path.join(rundir, 'out')
    if rc == 0 and series_has_nan(os.path.join(out, 'tempseries.out')):
        res['rc'] = 'non-finite'
    if res['rc'] == 0:
        rows = [l for l in open(os.path.join(out, 'fillet_global.out'))
                if l.strip() and not l.startswith('#')]
        lat = [l.rstrip('\n') for l in open(os.path.join(out, 'fillet.out'))
               if l.strip() and not l.startswith('#')]
        conv = [l for l in open(os.path.join(out, 'convergence.out'))
                if l.strip() and not l.startswith('#')]
        if len(rows) != 1 or len(lat) != 18 or len(conv) != 1:
            res['rc'] = 'incomplete output'
        else:
            res['global'] = rows[0].split()[1:]
            res['lat'] = lat
            orbits, flag, dolr, dtg, asym = conv[0].split()
            res['conv'] = dict(orbits=int(orbits), converged=flag != '0', cycle=int(flag),
                               dolr=float(dolr), dtglob=float(dtg), asym=float(asym))
    if res['rc'] == 0:
        with open(result_file, 'w') as f:
            json.dump(res, f)
    return res


def config_lines(exp, args, version):
    ben1 = exp == 'ben1'
    if ben1:
        d_rule = ('d0 = %.2f rescaled at run time by pressure, composition and rotation '
                  '(diffadj; 1.101 for 1 bar N2 at 280 ppm), the published Earth transport; '
                  'the Diff column is the value used' % args.d0_ben1)
        cloud = 'cloudir = %.2f W m^-2 (tuned so that Benchmark 1 reaches 288 K)' % args.cloudir_ben1
        surf = ('Fresnel ocean albedo, land 0.30, snow/ice 0.663; heat capacities land 5.25e6, '
                'water 2.1e8, ice 1.05e7 J m^-2 K^-1; present-day geography, ocean fraction 0.7')
    else:
        d_rule = ('D = %.2f W m^-2 K^-1 constant%s' %
                  (args.d0, ' rescaled by pressure, composition and rotation (diffadj)'
                   if args.diffadj else ' (diffadj off); the Diff column is the value used'))
        cloud = ('cloudir = %.2f W m^-2 (%s)' %
                 (args.cloudir, 'untuned' if args.cloudir == 0.0 else 'declared offset'))
        surf = ('surface albedo ocean 0.20 / land 0.30 / snow-ice 0.60; heat capacities land 1e7, '
                'water 4e8, ice 1e7 J m^-2 K^-1; ocean fraction 0.75 in every belt')
    nstep = args.nstepyr
    for e in args.extra_ebm:
        if e.replace(' ', '').startswith('nstepyr='):
            nstep = int(e.split('=')[1].strip(' ,'))
    period = orbital_period_days(1.0)
    if nstep:
        steps = '%d steps per orbit (dt = %.1f s, set by nstepyr)' % (nstep, period * SECONDS_PER_DAY / nstep)
    else:
        steps = 'dt = %.0f s (%.1f steps per orbit; the model year is the next whole step)' % (
            args.dt, period * SECONDS_PER_DAY / args.dt)
    orbit = 'a = 1 au, e = 0, year %.2f d' % period
    if exp in ('exp1a', 'exp2a'):
        orbit = 'a varied with S = S_earth / a^2, e = 0, year = %.2f d x a^1.5' % period
    switches = [s for s in switch_lines(args) if not s.replace(' ', '').startswith('nstepyr=')]
    numerics = ', '.join(switches) if switches else 'published defaults (no numerics switches set)'
    return [
        '# Model: HEXTOR %s (%s); 18 latitude belts of 10 deg, seasonal cycle, explicit %s.' % (version, args.label, steps),
        '# Radiation: 1 bar lookup table %s (clear sky; CO2 axis 100 ppm - 0.91, T 190-370 K with power-law'
        ' extrapolation outside; CO2 below 100 ppm is clamped to the axis floor); no cloud albedo (fcloud = 0); %s.' % (RADFILE, cloud),
        '# Transport: %s.' % d_rule,
        '# Surface: %s; sea-ice fraction 1 - exp((T - 273.15 K)/10 K) between 273.15 and 263.15 K and 1 below; land snow below 273.15 K.' % surf,
        '# Orbit: %s.' % orbit,
        '# Initial state: uniform %.0f K (warm start) or %.0f K (cold start) as noted per experiment; every case starts from that state (no continuation).' % (args.twarm, args.tcold),
        '# Convergence: halt when the year-to-year change in global-mean OLR falls below %.0e W m^-2%s (cap 5000 orbits); one further orbit provides the annual means. Per-case records in logs/.' % (
            args.fluxcnvg, ', or when it repeats to that tolerance over a cycle of up to %d years (a periodic seasonal oscillation; the case then reports one phase of the cycle)' % args.cnvgcycle if args.cnvgcycle >= 2 else ''),
        '# Numerics switches: %s.' % numerics,
    ]


def global_header(exp, args, version, start):
    lines = ['# FILLET %s with HEXTOR %s' % (LONG_NAME[exp], version)]
    lines += config_lines(exp, args, version)
    lines += [
        '#',
        '# Name of benchmark/experiment: %s' % LONG_NAME[exp],
        '# Initial state for this experiment: uniform %.0f K (%s start)' % (start, 'warm' if start == args.twarm else 'cold'),
        '# Describe how ice line latitude is determined: latitude where the annual-mean belt temperature crosses %.2f K,'
        ' linearly interpolated between belt centres. Ice-free: NMin = 90, SMax = -90; snowball: NMin = 0, SMax = 0.'
        ' A cap coexisting with a belt in one hemisphere is reported as the cap (noted in logs/).' % args.icelinetemp,
        '#',
        '# Columns of data',
        '# Case = case number (0 -> total number of cases in experiment)',
        "# Inst = instellation (S_earth; i.e., relative to Earth's 1361 W m^-2)",
        '# Obl = obliquity (degrees)',
        '# XCO2 = mixing ratio (volume) of CO2 (ppm)',
        '# Tglob = global, annual mean surface temperature (K)',
        '# IceLineNMax = latitude of northern edge of ice in northern hemisphere (deg)',
        '# IceLineNMin = latitude of southern edge of ice in northern hemisphere (deg)',
        '# IceLineSMax = latitude of northern edge of ice in southern hemisphere (deg)',
        '# IceLineSMin = latitude of southern edge of ice in southern hemisphere (deg)',
        '# Diff = diffusion coefficient actually used by the run (W m^-2 K^-1)',
        '# OLR = global, annual mean outgoing longwave radiation (W m^-2)',
        '#',
        '# Examples for IceLine max and min:',
        '# Polar caps: IceLineNMax = 90, IceLineNMin = 60, IceLineSMax = -60, IceLineSMin = -90',
        '# Snowball state: IceLineNMax = 90, IceLineNMin = 0, IceLineSMax = 0, IceLineSMin = -90',
        '# Ice free state: IceLineNMax = 90, IceLineNMin = 90, IceLineSMax = -90, IceLineSMin = -90',
        '# Ice belt: IceLineNMax = 30, IceLineNMin = 0, IceLineSMax = 0, IceLineSMin = -30',
        "# Note ice belt doesn't really have \"edges\" at 0 deg, but record as 0 anyway to",
        '# keep the columns the same length',
        '#',
        '# Case Inst Obl XCO2 Tglob IceLineNMax IceLineNMin IceLineSMax IceLineSMin Diff OLRglob',
    ]
    return '\n'.join(lines) + '\n'


def convergence_text(cv):
    if not cv['converged']:
        return 'ITERATION CAP REACHED'
    if cv.get('cycle', 1) > 1:
        return ('year-to-year OLR test not met, OLR repeats over a %d-year cycle (the means are one '
                'phase of it; the year-to-year change is its amplitude)' % cv['cycle'])
    return 'year-to-year OLR test met'


def lat_file(exp, res, args, version):
    c, cv = res['case'], res['conv']
    lines = ['# FILLET %s with HEXTOR %s' % (LONG_NAME[exp], version)]
    lines += config_lines(exp, args, version)
    lines += [
        '#',
        '# Name of benchmark/experiment: %s' % LONG_NAME[exp],
        '# Case number: %d' % c['case'],
        '# Instellation (S_earth): %.4f' % c['inst'],
        '# XCO2 (ppm): %.4f' % c['xco2'],
        '# Obliquity (degrees): %.1f' % c['obl'],
        '# Initial state: uniform %.0f K' % c['tinit'],
        '# Convergence: %d orbits, %s; final year-to-year |dOLR| = %.3e W m^-2, |dTglob| = %.3e K;'
        ' north-south asymmetry of the annual-mean belt temperatures = %.4f K' % (
            cv['orbits'], convergence_text(cv), cv['dolr'], cv['dtglob'], cv['asym']),
        '#',
        '# Columns of data (annually averaged for last orbit)',
        '# Lat = latitude of belt centre (degrees)',
        '# Tsurf = annual mean surface temperature (K)',
        '# Asurf = insolation-weighted annual mean surface albedo',
        '# ATOA = insolation-weighted annual mean top-of-atmosphere/planetary albedo (1 - ASR/S)',
        '# OLR = annual mean outgoing longwave radiation (W m^-2)',
        '# Lat Tsurf Asurf ATOA OLR',
    ]
    return '\n'.join(lines + res['lat']) + '\n'


def assemble(exp, results, args, version):
    archive = os.path.join(args.outdir, 'Results', 'hextor', exp)
    latdir = os.path.join(args.outdir, 'latfiles', exp)
    logdir = os.path.join(args.outdir, 'logs')
    for d in (archive, latdir, logdir):
        os.makedirs(d, exist_ok=True)
    start = results[0]['case']['tinit']
    with open(os.path.join(archive, 'global_output.dat'), 'w') as g:
        g.write(global_header(exp, args, version, start))
        for r in results:
            g.write('%d %s\n' % (r['case']['case'], ' '.join(r['global'])))
    for r in results:
        text = lat_file(exp, r, args, version)
        with open(os.path.join(latdir, 'case_%d.dat' % r['case']['case']), 'w') as f:
            f.write(text)
        if exp.startswith('ben'):
            os.makedirs(os.path.join(archive, 'case_0'), exist_ok=True)
            with open(os.path.join(archive, 'case_0', 'lat_output.dat'), 'w') as f:
                f.write(text)
    with open(os.path.join(logdir, '%s_convergence.log' % exp), 'w') as f:
        f.write('# case inst obl xco2 tglob orbits flag dOLR dTglob asymNS seconds'
                '   [flag: 1 year-to-year test, p = 2..4 p-year cycle, 0 cap]\n')
        for r in results:
            c, cv = r['case'], r['conv']
            f.write('%d %.4f %.1f %.4f %s %d %d %.3e %.3e %.4f %.1f\n' % (
                c['case'], c['inst'], c['obl'], c['xco2'], r['global'][3], cv['orbits'],
                cv.get('cycle', int(cv['converged'])), cv['dolr'], cv['dtglob'], cv['asym'], r['seconds']))
    n_cap = sum(1 for r in results if not r['conv']['converged'])
    n_cyc = sum(1 for r in results if r['conv'].get('cycle', 1) > 1)
    n_asym = sum(1 for r in results if r['conv']['asym'] > args.asym_warn)
    print('%-10s %3d cases -> %s   (%d hit the iteration cap, %d on a multi-year cycle, '
          '%d with N-S asymmetry > %.2f K)'
          % (exp, len(results), os.path.relpath(archive, ROOT), n_cap, n_cyc, n_asym, args.asym_warn))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--exp', nargs='+', default=['all'],
                    help='experiments to run: ben1 ben2 ben3 exp1 exp1a exp2 exp2a exp3 exp4, or all')
    ap.add_argument('--outdir', default=os.path.join(ROOT, 'fillet'))
    ap.add_argument('--driver', default=os.path.join(ROOT, 'model', 'driver'))
    ap.add_argument('--jobs', type=int, default=4)
    ap.add_argument('--timeout', type=float, default=7200.0)
    ap.add_argument('--force', action='store_true', help='rerun cases that have results')
    ap.add_argument('--dry-run', action='store_true', help='list the cases and stop')
    ap.add_argument('--label', default='FILLET re-file, September 2026',
                    help='short description written into every header')
    # configuration
    ap.add_argument('--cloudir', type=float, default=0.0,
                    help='longwave cloud offset for Ben2/3 and all experiments (W/m2; 0 = untuned)')
    ap.add_argument('--cloudir-ben1', type=float, default=9.3695,
                    help='the tuned Benchmark 1 offset (W/m2; tools/tune_fillet_ben1.py gives '
                         '9.3695 for 288.0 K with the other defaults)')
    ap.add_argument('--d0', type=float, default=0.5, help='D for Ben2/3 and experiments')
    ap.add_argument('--d0-ben1', type=float, default=0.38, help='Benchmark 1 d0 (diffadj on)')
    ap.add_argument('--diffadj', action='store_true',
                    help='rescale D by pressure, composition and rotation in Ben2/3 and experiments '
                         '(the pre-4.3.0 files did this while filing d0)')
    ap.add_argument('--dt', type=float, default=43200.0, help='time step (s) when --nstepyr is 0')
    ap.add_argument('--nstepyr', type=int, default=730,
                    help='steps per orbit (even; 0 = the published clock, dt from --dt and a model '
                         'year of the next whole step)')
    ap.add_argument('--diffcons', action=argparse.BooleanOptionalAction, default=True,
                    help='flux-conservative diffusion on the true belt edges '
                         '(--no-diffcons: the published operator)')
    ap.add_argument('--fluxcnvg', type=float, default=1.0e-3,
                    help='halt when the year-to-year change in global OLR is below this (W/m2)')
    ap.add_argument('--icelinetemp', type=float, default=263.15,
                    help='threshold of the reported ice line (K)')
    ap.add_argument('--cnvgcycle', type=int, default=4,
                    help='also halt when global OLR repeats over a cycle of up to this many years '
                         '(0 = year-to-year test only)')
    ap.add_argument('--twarm', type=float, default=288.0)
    ap.add_argument('--tcold', type=float, default=233.0)
    ap.add_argument('--extra-ebm', action='append', default=[],
                    help="extra &ebm line, e.g. 'icecont = .true.' (repeatable; a line that sets "
                         "diffcons or nstepyr replaces the flag)")
    ap.add_argument('--asym-warn', type=float, default=0.05,
                    help='flag cases whose N-S asymmetry exceeds this (K)')
    args = ap.parse_args()

    exps = []
    for e in args.exp:
        if e == 'all':
            exps = list(ALL_EXPS)
            break
        if e in ('exp3', 'exp4'):
            exps += [e + '_cold', e + '_warm']
        elif e in LONG_NAME:
            exps.append(e)
        else:
            sys.exit('unknown experiment %s' % e)
    if not os.path.exists(args.driver):
        sys.exit('driver not found: %s' % args.driver)
    modeldir = os.path.dirname(os.path.abspath(args.driver))
    try:
        version = subprocess.run(['git', 'describe', '--tags', '--always', '--dirty'],
                                 cwd=ROOT, text=True, stdout=subprocess.PIPE).stdout.strip()
    except OSError:
        version = 'unknown'

    for exp in exps:
        cases = cases_for(exp, args.twarm, args.tcold)
        if args.dry_run:
            print('%-10s %3d cases; first: %s; last: %s' % (exp, len(cases), cases[0], cases[-1]))
            continue
        t0 = time.time()
        with concurrent.futures.ThreadPoolExecutor(args.jobs) as pool:
            results = list(pool.map(lambda c: run_case(exp, c, args, modeldir), cases))
        failed = [r for r in results if r['rc'] != 0]
        for r in failed:
            print('  %s case %d FAILED: %s' % (exp, r['case']['case'], r['rc']))
        if failed:
            print('%s: %d of %d cases failed; not assembling' % (exp, len(failed), len(results)))
            continue
        assemble(exp, results, args, version)
        print('   %.0f s' % (time.time() - t0))


if __name__ == '__main__':
    main()
