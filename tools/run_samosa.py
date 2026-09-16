#!/usr/bin/env python3
"""
run_samosa.py — run HEXTOR over the SAMOSA intercomparison cases.

SAMOSA (Haqq-Misra et al. 2022, PSJ 3, 260) samples a (surface pressure,
instellation) parameter space for a synchronously rotating aquaplanet around a
3000 K blackbody star.  This script fills namelists/input.nml.samosa for each
case, runs each in its own scratch directory (so cases can run in parallel and
nothing touches model/out), and collects the diagnostics.

Participation options from the protocol, selectable with --sequence:

    warm      1 case  — the warm ExoCAM case only (case 4)
    seq2      8 cases — Sequence 2 (cases 9-16), the GCM-stable subspace
    stable   10 cases — every case ExoCAM ran stably (1, 4, 8-12, 14-16)
    all16    16 cases — Sequences 1 and 2 (the primary set)
    all32    32 cases — adds Sequences 1b and 2b
    all64    64 cases — adds Sequence 3

Outputs (in --outdir):
    samosa_summary.csv          one row per case: global means and ice line
    case_NN_<init>/zonal.txt    per-case zonal profile (ASCII, which the
                                protocol accepts from EBMs)
    case_NN_<init>/model.out    the raw HEXTOR report

Usage:
    python tools/run_samosa.py --sequence all16 --d0-ref 3.10
    python tools/run_samosa.py --sequence stable --init both
"""

import argparse
import csv
import os
import re
import shutil
import subprocess
import sys
import time
from multiprocessing import Pool

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import moist_stability as ms

HEXTOR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TEMPLATE = os.path.join(HEXTOR, 'namelists', 'input.nml.samosa')
DRIVER = os.path.join(HEXTOR, 'model', 'driver')
DEFAULT_TABLE = './radiation/radiation_N2_CO2_3000K_p.h5'

# CO2 is a fixed PARTIAL PRESSURE of 400 ubar (Haqq-Misra et al. 2024
# erratum), not a fixed 400 ppm mixing ratio: the CO2 column is the same in
# every case and the mixing ratio falls as the N2 pressure rises.  HEXTOR's
# pg0 is the dry surface pressure (pN2 + pCO2), so each case gets
# pg0 = pN2 + P_CO2 and fco2 = P_CO2 / pg0.
P_CO2 = 4.0e-4           # bar

# (sample, instellation [W/m2], N2 pressure [bar]) from the protocol.
# Table 1 = Sequence 1 (1-8) and Sequence 2 (9-16); Table 4 = Sequence 1b
# (17-24), Sequence 2b (25-32) and Sequence 3 (33-64).
CASES = [
    (1,  500,  0.70), (2,  1900, 7.85), (3,  2400, 0.21), (4,  1200, 2.34),
    (5,  1500, 0.16), (6,  2100, 1.83), (7,  1600, 0.55), (8,  800,  6.16),
    (9,  1100, 0.70), (10, 400,  4.83), (11, 900,  0.10), (12, 1500, 2.98),
    (13, 1600, 0.16), (14, 900,  1.44), (15, 600,  0.43), (16, 1400, 10.0),
    (17, 700,  0.26), (18, 1700, 2.98), (19, 2300, 0.89), (20, 1300, 10.0),
    (21, 900,  0.43), (22, 2600, 4.83), (23, 2000, 0.10), (24, 500,  1.13),
    (25, 1300, 0.16), (26, 500,  2.34), (27, 1000, 0.89), (28, 1700, 3.79),
    (29, 1400, 0.55), (30, 800,  7.85), (31, 500,  0.21), (32, 1200, 1.13),
]

SEQUENCES = {
    'warm':   [4],
    'seq2':   list(range(9, 17)),
    'stable': [1, 4, 8, 9, 10, 11, 12, 14, 15, 16],
    'all16':  list(range(1, 17)),
    'all32':  list(range(1, 33)),
    'all64':  list(range(1, 65)),
}

# Initial surface temperatures.  HEXTOR's ice-albedo feedback is bistable, so
# a warm and a cold start can land in different states at the same forcing;
# running both is how that is detected rather than accidentally sampled.
INITS = {'warm': 300.0, 'cold': 233.0}


# Diffusion stability, matching how driver.f discretises: 18 belts in
# x = sin(lat), explicit in time, heat capacity C.  With moistdiff the
# diffused quantity is h/cp rather than T, which responds to a temperature
# perturbation beta = 1 + (L/cp) dq/dT times faster, so the limit tightens by
# beta -- 3.5 at 300 K, 13 at 343 K.  See tools/moist_stability.py.
NBELTS = 18
HEATCAP = 4.0e6          # &ebm::heatcap in the SAMOSA template
DT_PUBLISHED = 1350.0    # the published THAI timestep
CFL_SAFETY = 0.4         # explicit diffusion needs D*dt/(C*dx^2) below ~1/2

# Ceiling on the temperature used to evaluate beta.  420 K is the top of the
# non-extrapolated part of the radiation table, so a case hotter than this is
# already labelled a runaway and there is nothing to be gained by shrinking dt
# further to resolve it; without a cap, beta at a 3000 K surface would drive dt
# to zero and the run would never finish.
T_HOT_CAP = 420.0
T_HOT_START = 300.0      # first guess for the substellar belt


# How heat transport is allowed to depend on the case.  All four are
# calibrated to the SAME effective D at THAI Hab1 (1 bar), so they differ only
# in how that is carried across the SAMOSA parameter space.
TRANSPORT = {
    # D = d0, identical in every case.  diffadj = .false. means the driver
    # uses d0 literally (driver.f:504), so this is what that flag alone gives.
    'constant': dict(diffadj=False, p_exp=None, rot=False),
    # D = d0 * sqrt(p), applied by this script through d0 with the driver's own
    # scaling off.  Measured against ExoCAM, neither of the two published
    # options is right at both ends of the pressure axis: constant D leaves
    # Case 16 (10 bar) with a 45 K day-night contrast where the GCM has 7, and
    # D ~ p collapses Case 11 (0.1 bar) to a 92 K night side and destroys its
    # climate solution.  A sub-linear exponent is the obvious interpolation
    # between them, and 1/2 needs no recalibration: SAMOSA's reference point is
    # 1 bar, where sqrt(p) = p = 1 and all three give the same D.
    'psqrt': dict(diffadj=False, p_exp=0.5, rot=False),
    # D = d0 * p * composition.  Pressure scaling applied by the driver.
    'perbar': dict(diffadj=True, p_exp=None, rot=False),
    # D = d0 * p * composition * (rot0/rot)^2, HEXTOR's full scaling.
    'diffadj': dict(diffadj=True, p_exp=None, rot=True),
}


# How the cloud correction is allowed to depend on the case.  cloudir stands in
# for the shortwave cooling of the clouds the clear-sky table lacks, and a
# cloud shortwave effect is a fraction of the incident stellar flux, so a
# fixed W/m2 offset is too strong at low instellation (45 W/m2 is 45% of the
# absorbed flux at 500 W/m2) and too weak at high.  'instellation' scales it
# with S, anchored at the instellation where it was calibrated, so both modes
# give the SAME cloudir at THAI Hab1 and the calibration carries over unchanged
# -- the same construction as the TRANSPORT options.  This is equivalent to a
# fixed increment in planetary albedo of -cloudir_ref / (S_REF/4), 0.20 for
# -45 W/m2, still applied uniformly to every belt.
S_REF = 900.0            # THAI Hab1 instellation [W/m2], namelists/input.nml.thai.hab1
CLOUD_SCALING = ('constant', 'instellation')


def effective_cloudir(cloudir_ref, inst, scaling):
    """The cloudir the driver will use for a case at instellation inst."""
    if scaling == 'instellation':
        return cloudir_ref * inst / S_REF
    return cloudir_ref


def effective_D(d0_ref, pn2, diffadj, rot_scaling, pco2=P_CO2,
                rot=4.84813681e-6, rot0=7.27e-5):
    """The diffusion coefficient the driver will end up using."""
    if not diffadj:
        return d0_ref
    pg0 = pn2 + pco2
    avemol = (28.0 * pn2 + 44.0 * pco2) / pg0
    hcp = (0.2484 * pn2 + 0.2105 * pco2) / pg0
    d = d0_ref * pg0 * (28.89 / avemol) ** 2 * (hcp / 0.2401)
    if rot_scaling:
        d *= (rot0 / rot) ** 2
    return d


def stable_dt(D, pg0=1.0, fco2=0.0, t_hot=None, rhmoist=0.0):
    """Timestep satisfying the explicit-diffusion limit, capped at published.

    With rhmoist > 0 (moistdiff) the limit is divided by beta evaluated at
    t_hot, the warmest belt.  Returns (dt, beta).
    """
    if rhmoist <= 0.0 or t_hot is None:
        return ms.stable_dt(D, HEATCAP, DT_PUBLISHED, nbelts=NBELTS,
                            cfl=CFL_SAFETY), 1.0
    return ms.moist_dt(D, pg0, pg0 * (1 - fco2), pg0 * fco2,
                       min(t_hot, T_HOT_CAP), rhmoist, HEATCAP, DT_PUBLISHED,
                       nbelts=NBELTS, cfl=CFL_SAFETY)


def prepare_rundir(rundir):
    """A HEXTOR run directory: the driver plus the relative paths it opens."""
    os.makedirs(os.path.join(rundir, 'out'), exist_ok=True)
    for name, target in (('data', os.path.join(HEXTOR, 'model', 'data')),
                         ('radiation', os.path.join(HEXTOR, 'model', 'radiation')),
                         ('driver', DRIVER)):
        link = os.path.join(rundir, name)
        if not os.path.exists(link):
            os.symlink(target, link)


SECONDS_PER_YEAR = 365.0 * 86400.0


def ffloat(text):
    """float() that tolerates what a runaway does to Fortran's output.

    HEXTOR writes its diagnostics with fixed-width descriptors (f6.2, f8.3),
    so a case that runs away to 10^5 K overflows them and Fortran emits
    '********'.  That is information -- the run left the physical range -- not
    a reason for the harness to die, which is what it used to do, taking the
    whole worker pool with it.
    """
    try:
        return float(text)
    except ValueError:
        return float('nan')


def parse_zonal(model_out):
    """Pull the ZONAL STATISTICS block out of model.out.

    Columns: coordinate(deg), Tave, Tmin, dec@Tmin, Tmax, dec@Tmax, albedo,
    OLR, ASR.  In do_longitudinal mode the first column is the angle from the
    substellar point expressed on HEXTOR's latitude grid.
    """
    rows = []
    if not os.path.exists(model_out):
        return rows
    with open(model_out) as f:
        lines = f.readlines()
    for i, line in enumerate(lines):
        if line.startswith('ZONAL STATISTICS'):
            for ln in lines[i + 2:]:
                parts = ln.split()
                if len(parts) != 9:
                    break
                if not any(ch.isdigit() for ch in parts[0]):
                    break
                rows.append([ffloat(x) for x in parts])
            break
    return rows


def belt_area_weights(coords_deg, nbelts=18):
    """Fractional area of each HEXTOR belt.

    The belts are equally spaced in the coordinate, NOT in area
    (driver.f:376, 490): a belt centred at c spans sin(c + pi/2n) - sin(c -
    pi/2n).  Unweighted means over the 18 belts would over-weight the poles —
    or, in the tidally-locked coordinate, the substellar and antistellar
    points.  The same expression serves both readings, since x = sin(lat) is
    reinterpreted as cos(angle from substellar).
    """
    import math
    half = math.pi / (2 * nbelts)
    w = [abs(math.sin(math.radians(c) + half) - math.sin(math.radians(c) - half))
         for c in coords_deg]
    tot = sum(w)
    return [x / tot for x in w]


def parse_scalars(rundir, stdout=''):
    """Global diagnostics from tempseries.out, icelines.out and model.out."""
    out = {}
    ts = os.path.join(rundir, 'out', 'tempseries.out')
    if os.path.exists(ts):
        lines = [l for l in open(ts) if l.strip()]
        if lines:
            p = lines[-1].split()
            if len(p) >= 8:
                out['T_global'] = ffloat(p[1])
                out['pg0'] = ffloat(p[2])
                out['pco2'] = ffloat(p[3])
                out['q'] = ffloat(p[6])
                out['d'] = ffloat(p[7])
            # HEXTOR halts on the year-to-year change in global OLR, which is
            # a false positive on the runaway-greenhouse plateau: there OLR is
            # genuinely insensitive to surface temperature, so it can go flat
            # while the temperature is still climbing.  Record the temperature
            # trend over the final years so that case can be recognised rather
            # than reported as an equilibrium.
            tail = [ffloat(l.split()[1]) for l in lines[-11:] if len(l.split()) >= 2]
            if len(tail) >= 2:
                out['dT_last'] = tail[-1] - tail[0]
                out['n_years'] = len(lines)
            # The residual drift, per decade.  Taking it as the change over the
            # last ten years (dT_last) is wrong for a run that halts soon after
            # it converges: these cases approach equilibrium exponentially and
            # reach it in ~20 years, so a ten-year window still spans most of
            # the transient.  Case 11 settles to 0.001 K/yr by year 19 and yet
            # scores -0.51 K over its last decade, which labelled a converged
            # run 'drifting' and dropped it from the submission.  The median of
            # the last three annual increments has neither problem: it equals
            # the ten-year change for anything drifting linearly (a runaway, or
            # the slow creep on the OLR plateau) and goes to zero once the
            # exponential has died, and taking the median rather than the last
            # increment alone keeps the +-0.001 K flicker of a moist run out of
            # it.
            steps = [b - a for a, b in zip(tail, tail[1:])]
            if steps:
                steps = sorted(steps[-3:])
                out['dT_rate10'] = 10.0 * steps[len(steps) // 2]
    il = os.path.join(rundir, 'out', 'icelines.out')
    if os.path.exists(il):
        lines = [l for l in open(il) if l.strip()]
        if lines:
            p = lines[-1].split()
            if len(p) >= 3:
                out['icelineS'] = ffloat(p[1])
                out['icelineN'] = ffloat(p[2])
    # Fraction of belts with CO2 condensing at the surface.  SAMOSA's cold
    # corner can reach this: with only 400 ubar of CO2 a cold enough surface
    # drops below the frost point, and HEXTOR then pins it there.
    # That is a distinct climate regime, not an equilibrium and not a runaway,
    # and it needs its own label.
    cc = os.path.join(rundir, 'out', 'co2clouds.out')
    if os.path.exists(cc):
        lines = [l for l in open(cc) if l.strip()]
        if lines:
            p_ = lines[-1].split()
            if len(p_) >= 19:
                flags = [ffloat(x) for x in p_[1:19]]
                good = [x for x in flags if x == x]
                if good:
                    out['co2_condensing'] = sum(1 for x in good if x > 0.5) / len(good)

    # The flux-convergence message goes to stdout, not model.out.
    m = re.search(r'Flux converged at year\s+([0-9.Ee+-]+)', stdout)
    out['converged'] = bool(m)
    m = re.search(r'dOLR =\s*\n?\s*([0-9.Ee+-]+)', stdout)
    if m:
        out['dOLR'] = float(m.group(1))
    return out


def substellar_longitude(iceline_n, iceline_s):
    """Ice-line position as degrees from the substellar point.

    In do_longitudinal mode the belt coordinate is the angle from substellar,
    so 0 deg = substellar, 90 deg = terminator, 180 deg = antistellar.
    Sentinels: (90, -90) means ice free, (0, 0) means fully glaciated.
    """
    if iceline_n is None or iceline_s is None:
        return None
    if iceline_n == 90.0 and iceline_s == -90.0:
        return 180.0
    if iceline_n == 0.0 and iceline_s == 0.0:
        return 0.0
    if iceline_n < 90.0:
        return 90.0 - iceline_n
    return 90.0 - iceline_s


RC_NONFINITE = -2        # the run went non-finite and was stopped
DEFAULT_TIMEOUT = 1800   # a healthy case here finishes in seconds


def _series_has_nan(path):
    """Has the annual temperature series gone non-finite?"""
    try:
        with open(path, 'rb') as f:
            f.seek(0, os.SEEK_END)
            size = f.tell()
            f.seek(max(0, size - 4096))
            tail = f.read()
        return b'NaN' in tail or b'nan' in tail
    except OSError:
        return False


def _integrate(rundir, template_fields, timeout=DEFAULT_TIMEOUT, poll=2.0):
    """Write the namelist and run the driver once, watching for a blow-up.

    An unstable timestep does not make the driver fail: it produces NaN and
    then keeps integrating NaN until the hard-coded niter = 5000 is reached,
    which at these timesteps is two hours of wall time per case (tend cannot
    stop it -- see the --years note).  Worse, the run then returns a timeout,
    which used to abort the caller's dt refinement before it could retry at a
    smaller step, so the case was reported as a runaway.  Polling the annual
    series turns that into a few seconds and a distinguishable return code.
    """
    nml = open(TEMPLATE).read().format(**template_fields)
    with open(os.path.join(rundir, 'input.nml'), 'w') as f:
        f.write(nml)
    series = os.path.join(rundir, 'out', 'tempseries.out')
    if os.path.exists(series):
        os.remove(series)          # never inspect the previous pass's output

    proc = subprocess.Popen(['./driver'], cwd=rundir, text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    t0, stdout, err, rc = time.time(), '', '', None
    while True:
        try:
            stdout, stderr = proc.communicate(timeout=poll)
            rc = proc.returncode
            err = '' if rc == 0 else (stderr or stdout)[-300:].replace('\n', ' ')
            break
        except subprocess.TimeoutExpired:
            pass
        if _series_has_nan(series):
            proc.kill()
            proc.communicate()
            rc, err = RC_NONFINITE, 'non-finite temperature'
            break
        if time.time() - t0 > timeout:
            proc.kill()
            proc.communicate()
            rc, err = -1, 'timeout'
            break
    return rc, stdout, err


def run_case(task):
    (sample, inst, ps, init_name, init_t, d0_ref, cloudir_ref, cloud_scaling,
     table, outdir, diffadj, rot_scaling, years, moistdiff, rhmoist,
     fluxcnvg, p_exp) = task

    tag = 'case_%02d_%s' % (sample, init_name)
    rundir = os.path.join(outdir, tag)
    prepare_rundir(rundir)

    # With diffadj the driver applies the pressure scaling itself (along with
    # the composition and rotation factors), so d0 is passed through as the
    # reference value; without it the runner does the pressure scaling here.
    pg0 = ps + P_CO2
    # With diffadj the driver applies the pressure scaling itself; p_exp is the
    # other route, where this script folds it into d0 and the driver uses that
    # literally.  Only one of the two is ever active.
    d0 = d0_ref if p_exp is None else d0_ref * pg0 ** p_exp
    D = effective_D(d0, ps, diffadj, rot_scaling)
    fco2 = P_CO2 / pg0
    cloudir = effective_cloudir(cloudir_ref, inst, cloud_scaling)
    rh = rhmoist if moistdiff else 0.0

    fields = dict(
        tempinit='%.1f' % init_t, d0='%.6f' % d0, pg0='%.5f' % pg0,
        fco2='%.6e' % fco2,
        solarcon='%.1f' % inst, cloudir='%.2f' % cloudir, radfile=table,
        diffadj='.true.' if diffadj else '.false.',
        diffadj_rot='.true.' if rot_scaling else '.false.',
        moistdiff='.true.' if moistdiff else '.false.',
        rhmoist='%.4f' % rhmoist, fluxcnvg='%.3e' % fluxcnvg,
        tend='%.4e' % (years * SECONDS_PER_YEAR))

    # The moist timestep limit depends on the substellar temperature, which is
    # not known until the case has run.  Iterate: run at a guess, then redo at
    # the beta the run actually produced.  Without moistdiff beta is 1 and this
    # is a single pass.  A first pass that goes non-finite is retried once at
    # T_HOT_CAP, the most conservative in-range beta, so that an unstable
    # timestep is not mistaken for a runaway -- the two look identical in the
    # output, and only the retry tells them apart.
    t_hot = min(max(init_t, T_HOT_START), T_HOT_CAP)
    t0 = time.time()
    passes, dt_capped = 0, False
    while True:
        dt, bet = stable_dt(D, pg0, fco2, t_hot, rh)
        fields['dt'] = '%.1f' % dt
        rc, stdout, err = _integrate(rundir, fields)
        passes += 1
        if not moistdiff or passes >= 4:
            break
        zon = parse_zonal(os.path.join(rundir, 'out', 'model.out'))
        t_obs = max((r[1] for r in zon if r[1] == r[1]), default=float('nan'))
        if rc != 0:
            # A blow-up or a timeout is exactly the case the retry exists for,
            # so it must not break out of the loop before the retry happens.
            t_obs = float('nan')
        if t_obs != t_obs:                     # NaN: unstable, or a real runaway
            if t_hot >= T_HOT_CAP:
                break                          # already at the tightest dt
            t_hot, dt_capped = T_HOT_CAP, True
            continue
        if t_obs <= t_hot + 0.5:
            break                              # dt was already tight enough
        t_next = min(t_obs + 1.0, T_HOT_CAP)
        dt_capped = t_obs + 1.0 > T_HOT_CAP
        if t_next <= t_hot:
            break        # already at the cap: a case this hot is a runaway,
        t_hot = t_next   # and shrinking dt further cannot make it a climate

    rec = {'case': sample, 'init': init_name, 'instellation': inst,
           'ps_bar': ps, 'fco2': fco2, 'd0': d0, 'cloudir': cloudir,
           'cloudir_ref': cloudir_ref, 'cloud_scaling': cloud_scaling,
           'diffadj': diffadj,
           'rot_scaling': rot_scaling, 'D': D, 'dt': dt,
           'moistdiff': moistdiff, 'rhmoist': rhmoist if moistdiff else 0.0,
           'table': table, 'rh_table': table_rh(table),
           'fluxcnvg': fluxcnvg, 'beta': bet, 't_hot': t_hot,
           'dt_passes': passes, 'dt_capped': dt_capped,
           'rc': rc, 'error': err, 'wall_s': round(time.time() - t0, 1)}
    rec.update(parse_scalars(rundir, stdout))

    zonal = parse_zonal(os.path.join(rundir, 'out', 'model.out'))
    if zonal:
        # Report from substellar outward.
        zonal.sort(key=lambda r: -r[0])
        with open(os.path.join(rundir, 'zonal.txt'), 'w') as f:
            f.write('# HEXTOR / SAMOSA case %d (%s start)\n' % (sample, init_name))
            f.write('# instellation = %.1f W/m2   surface pressure = %.3f bar\n'
                    % (inst, ps))
            f.write('# theta = angle from substellar point [deg]; '
                    '0 = substellar, 180 = antistellar\n')
            f.write('# theta  T_ave_K  T_min_K  T_max_K  planetary_albedo  '
                    'OLR_Wm2  ASR_Wm2\n')
            for r in zonal:
                f.write('%7.1f %9.3f %9.3f %9.3f %9.4f %10.3f %10.3f\n'
                        % (90.0 - r[0], r[1], r[2], r[4], r[6], r[7], r[8]))
        w = belt_area_weights([r[0] for r in zonal])
        rec['T_min'] = min(r[1] for r in zonal)
        rec['T_max'] = max(r[1] for r in zonal)
        if any(r[1] != r[1] for r in zonal):
            rec['T_min'] = rec['T_max'] = float('nan')
        rec['OLR_global'] = sum(wi * r[7] for wi, r in zip(w, zonal))
        rec['ASR_global'] = sum(wi * r[8] for wi, r in zip(w, zonal))
        # What the protocol asks models to drive toward +-1 W/m2.
        rec['TOA_imbalance'] = rec['ASR_global'] - rec['OLR_global']

    rec['iceline_lon'] = substellar_longitude(rec.get('icelineN'),
                                              rec.get('icelineS'))
    rec['state'] = classify(rec, table)
    return rec


def classify(rec, table):
    """Label the outcome, since HEXTOR's own 'converged' flag is not enough.

    HEXTOR halts on the year-to-year change in global OLR.  On the runaway
    plateau OLR is insensitive to surface temperature, so that criterion trips
    while the planet is still heating: every case in this parameter space
    reports converged.  What separates a climate from a runaway is whether the
    top of the atmosphere actually balances, and whether the model stayed
    inside the range its radiation table covers.

    The protocol anticipates this — incipient-runaway cases may be reported at
    the last stable state or omitted — so runaways are labelled, not hidden.
    """
    if rec.get('rc') == RC_NONFINITE:
        # Non-finite even after the timestep was tightened to the smallest
        # beta the table covers.  That is a numerical failure, and calling it
        # a runaway would hide it in the omitted-cases list.
        return 'unstable (non-finite, dt at floor)' if rec.get('dt_capped') \
            else 'unstable (non-finite)'
    if rec.get('rc') == -1:
        return 'no result (timeout)'
    T = rec.get('T_global')
    if T is None or T != T:               # NaN, including a format overflow
        return 'runaway (out of range)'
    tmax = table_tmax(table)
    if tmax is not None and T > tmax:
        return 'runaway (beyond table, T > %.0f K)' % tmax
    cond = rec.get('co2_condensing', 0.0)
    if cond and cond > 0.5:
        return 'CO2 condensing'
    imb = rec.get('TOA_imbalance')
    dT = rec.get('dT_rate10', rec.get('dT_last'))
    if imb is None or imb != imb:
        return 'unknown'
    # Classify on the temperature trend, not the imbalance.  HEXTOR closes its
    # global budget only to ~1 W/m2 on the 18-belt grid (pre-industrial Earth
    # sits at -0.91 W/m2, and the discretised dayside insolation carries
    # +0.86 W/m2 by itself), so a 1-2 W/m2 residual is the model's floor rather
    # than a sign of non-convergence.  The protocol allows exactly this: a
    # stable trend suffices where the balance cannot be driven to +-1 W/m2.
    # The threshold on the drift has to be loose enough not to bisect a set of
    # equally converged runs.  Across the cloudir sweep the converged cases sit
    # at |dT| <= 0.13 K per decade while the runaways are at 0.6 K and above,
    # so 0.5 K separates them cleanly; a 0.1 K cut instead flipped cases in and
    # out of 'equilibrium' while their temperatures varied perfectly smoothly.
    # The protocol accepts "a stable trend over at least 10 orbits" where exact
    # balance cannot be reached, and 0.12 K per decade on a 180 K planet is
    # 0.07%.
    if dT is not None and abs(dT) < 0.5 and abs(imb) <= 3.0:
        return 'equilibrium'
    if dT is not None and abs(dT) < 2.0 and abs(imb) <= 10.0:
        return 'drifting'
    return 'runaway'


_TMAX = {}
_TATTR = {}


def table_path(table):
    """The table as this script sees it, given the path the driver sees."""
    if os.path.isabs(table):
        return table
    return os.path.join(HEXTOR, 'model', table.lstrip('./'))


def table_rh(table):
    """The relative humidity the table's columns were built at, or None.

    Recorded per run because the calibration is table-specific: an RH 0.8
    table and an RH 1 table are different radiative transfer, and (d0, cloudir)
    fitted against one does not carry to the other.
    """
    if table in _TATTR:
        return _TATTR[table]
    val = None
    try:
        import h5py
        with h5py.File(table_path(table), 'r') as f:
            val = f.attrs.get('rh', '1.0')
            if isinstance(val, bytes):
                val = val.decode()
            val = str(val)
    except Exception:
        val = None
    _TATTR[table] = val
    return val


def table_tmax(table):
    """Top of a table's temperature axis; None for the legacy 1 bar format."""
    if table in _TMAX:
        return _TMAX[table]
    path = table_path(table)
    val = None
    try:
        import h5py
        with h5py.File(path, 'r') as f:
            if 'temperature' in f:
                val = float(f['temperature'][:].max())
    except Exception:
        val = None
    _TMAX[table] = val
    return val


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--sequence', default='all16', choices=sorted(SEQUENCES))
    ap.add_argument('--init', default='warm', choices=['warm', 'cold', 'both'],
                    help='initial surface temperature (both = test bistability)')
    ap.add_argument('--d0-ref', type=float, default=3.10,
                    help='d0 written to the namelist. With --transport constant '
                         'it is D itself; with perbar/diffadj the driver '
                         'multiplies it by pressure (and composition, and '
                         'rotation). Default 3.10 is the published THAI Hab1 '
                         'value; the SAMOSA submission uses 3.33 constant.')
    ap.add_argument('--cloudir', type=float, default=-35.0,
                    help='cloud correction to OLR [W/m2]; driver.f applies '
                         'ir = ir - cloudir, so the published TRAPPIST-1 value '
                         'of -35 raises OLR and cools. Pass 0 for clear sky. '
                         'With --cloud-scaling instellation this is the value '
                         'at S = %.0f W/m2 (THAI Hab1).' % S_REF)
    ap.add_argument('--cloud-scaling', default='constant', choices=CLOUD_SCALING,
                    help="how cloudir varies across the parameter space: "
                         "'constant' (the same W/m2 in every case) or "
                         "'instellation' (cloudir * S / %.0f, a fixed cloud "
                         "albedo increment). Both give the calibrated value "
                         "at THAI Hab1." % S_REF)
    ap.add_argument('--table', default=DEFAULT_TABLE,
                    help='radiation table, as the driver sees it (relative to '
                         'the run directory)')
    ap.add_argument('--transport', default='perbar', choices=sorted(TRANSPORT),
                    help="how D varies across the parameter space: 'constant' "
                         "(D = d0 everywhere), 'psqrt' (D = d0 * sqrt(p)), "
                         "'perbar' (D = d0 * p * "
                         "composition, no rotation term), or 'diffadj' "
                         "(HEXTOR's full scaling, including (rot0/rot)^2). All "
                         "are calibrated to the same D at THAI Hab1.")
    ap.add_argument('--moistdiff', action='store_true',
                    help='diffuse moist static energy instead of temperature '
                         '(&ebm::moistdiff).  Match --rhmoist to the RH the '
                         'radiation table was built at.')
    ap.add_argument('--rhmoist', type=float, default=0.8,
                    help='relative humidity of the diffused moisture '
                         '(default 0.8, the RH of the _rh0.8 tables)')
    ap.add_argument('--fluxcnvg', type=float, default=5.0e-4,
                    help='OLR convergence threshold [W/m2] written to the '
                         'namelist (default 5e-4)')
    ap.add_argument('--years', type=float, default=200.0,
                    help='years written to &ebm::tend (default 200). NOTE this '
                         'is not a wall on the run: with seasons = .true. the '
                         'halt at tend is commented out in driver.f:593-600, '
                         'so tend only freezes the orbital clock and the run '
                         'continues to flux convergence or the hard-coded '
                         'niter = 5000. Harmless for these cases -- they are '
                         'synchronous at zero obliquity and eccentricity, so '
                         'the insolation does not vary within a year anyway.')
    ap.add_argument('--outdir', default=os.path.join(HEXTOR, 'samosa'))
    ap.add_argument('--workers', type=int, default=8)
    args = ap.parse_args()

    mode = TRANSPORT[args.transport]
    samples = SEQUENCES[args.sequence]
    known = {c[0]: c for c in CASES}
    missing = [s for s in samples if s not in known]
    if missing:
        print('cases not tabulated in this script (Sequence 3 is not '
              'transcribed): %s' % missing)
        samples = [s for s in samples if s in known]

    inits = ['warm', 'cold'] if args.init == 'both' else [args.init]

    os.makedirs(args.outdir, exist_ok=True)
    tasks = []
    for s in samples:
        _, inst, ps = known[s]
        for init_name in inits:
            tasks.append((s, inst, ps, init_name, INITS[init_name],
                          args.d0_ref, args.cloudir, args.cloud_scaling,
                          args.table, args.outdir,
                          mode['diffadj'], mode['rot'], args.years,
                          args.moistdiff, args.rhmoist, args.fluxcnvg,
                          mode['p_exp']))

    print('SAMOSA sequence %s: %d cases x %d start(s) = %d runs'
          % (args.sequence, len(samples), len(inits), len(tasks)))
    desc = {'constant': 'D = %.4f (constant)' % args.d0_ref,
            'psqrt': 'D = %.4f * sqrt(p)' % args.d0_ref,
            'perbar': 'D = %.4f * p * comp' % args.d0_ref,
            'diffadj': 'D = %.5f * p * comp * (rot0/rot)^2' % args.d0_ref}
    cdesc = {'constant': 'cloudir = %.2f W/m2' % args.cloudir,
             'instellation': 'cloudir = %.2f W/m2 * S/%.0f'
                             % (args.cloudir, S_REF)}
    print('  transport: %-42s %s'
          % (desc[args.transport], cdesc[args.cloud_scaling]))
    print('  table    : %s' % args.table)
    print('  diffusion: %s'
          % ('moist static energy, rhmoist = %.2f (dt set per case from the '
             'substellar beta)' % args.rhmoist if args.moistdiff
             else 'dry (temperature)'))
    print()

    results = []
    with Pool(min(args.workers, len(tasks))) as pool:
        for rec in pool.imap_unordered(run_case, tasks):
            results.append(rec)
            print('  case %2d %-4s  S=%4.0f  p=%5.2f bar  ->  '
                  % (rec['case'], rec['init'], rec['instellation'], rec['ps_bar'])
                  + ('T=%7.2f K  Tmin=%7.2f  Tmax=%7.2f  ice@%s  (%.0f s)'
                     % (rec.get('T_global', float('nan')),
                        rec.get('T_min', float('nan')),
                        rec.get('T_max', float('nan')),
                        ('%5.1f deg' % rec['iceline_lon'])
                        if rec.get('iceline_lon') is not None else '   n/a',
                        rec['wall_s'])
                     + '  [%s]' % rec.get('state', '?')
                     if rec['rc'] == 0 else 'FAILED: %s' % rec['error']),
                  flush=True)

    results.sort(key=lambda r: (r['case'], r['init']))
    cols = ['case', 'init', 'instellation', 'ps_bar', 'fco2', 'd0', 'cloudir',
            'cloudir_ref', 'cloud_scaling', 'T_global', 'T_min', 'T_max', 'OLR_global', 'ASR_global',
            'TOA_imbalance', 'dT_last', 'dT_rate10', 'n_years', 'state', 'diffadj',
            'rot_scaling', 'D', 'dt', 'table', 'rh_table', 'moistdiff',
            'rhmoist', 'beta', 't_hot',
            'dt_passes', 'dt_capped', 'fluxcnvg', 'co2_condensing',
            'icelineN', 'icelineS', 'iceline_lon', 'converged', 'dOLR',
            'wall_s', 'rc', 'error']
    path = os.path.join(args.outdir, 'samosa_summary.csv')
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction='ignore')
        w.writeheader()
        for r in results:
            w.writerow(r)

    nfail = sum(1 for r in results if r['rc'] != 0)
    print('\nwrote %s  (%d runs, %d failed)' % (path, len(results), nfail))
    tally = {}
    for r in results:
        tally[r.get('state', '?')] = tally.get(r.get('state', '?'), 0) + 1
    for k in sorted(tally):
        print('  %-38s %d' % (k, tally[k]))
    return 1 if nfail else 0


if __name__ == '__main__':
    sys.exit(main())
