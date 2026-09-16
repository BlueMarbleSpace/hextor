#!/usr/bin/env python3
"""
make_radiation_table.py — build a pressure-resolved HEXTOR radiation lookup
table (v2 format) with ExoColumn / ExoRT.

WHAT THIS PRODUCES
------------------
An HDF5 file readable by model/radiation/radiation.f90 (radparam = 3), holding

    /olr       (npre, nco2, [nch4,] ntmp)             outgoing longwave
    /palb      (npre, nco2, [nch4,] ntmp, nzen, nsab) planetary albedo
    /pressure  (npre)    dry surface pressure pN2 + pCO2 + pCH4  [bar]
    /fco2      (nco2)    CO2 mixing ratio pCO2 / p_dry     [-]
    /ch4       (nch4)    CH4 mixing ratio pCH4 / p_dry     [-]  (--ch4-axis only)
    /temperature (ntmp)  surface temperature               [K]
    /zenith    (nzen)    solar zenith angle                [degrees]
    /surfalb   (nsab)    broadband surface albedo          [-]

Without --ch4-axis (or an explicit --ch4 list) the CH4 axis and its dataset are
omitted entirely, giving a table byte-for-byte in the previous v2 layout.

OLR is stored in mW/m^2 (W/m^2 x 1000) to match the historical v1 convention:
driver.f divides the value returned by getOLR by 1000.

HOW THE TABLE IS COMPUTED
-------------------------
Each grid point is one ExoColumn run in `flux_only` + `sweep_mode`:

  * the column is a prescribed profile, not an RCE solution — a moist adiabat
    integrated upward from the tabulated surface temperature and capped at the
    stratospheric temperature, which is the inverse-mode calculation an EBM
    lookup table needs (an RCE run would solve for Ts, which is the EBM's job);
  * surface pressure is the DRY background (N2 + CO2); `variable_ps` adds the
    surface water vapour pressure on top, so the total surface pressure is
    p_dry + esat(Ts).  This matches the SAMOSA protocol's definition of surface
    pressure and the Kopparapu et al. (2013) convention;
  * neither the solar zenith angle nor the surface albedo feeds back on the
    profile, so ExoColumn builds the column once and loops the radiation call
    over the whole (zenith x albedo) block — ~0.014 s per extra point against
    ~1.6 s of fixed ExoRT initialisation.

Runs are distributed over NWORKERS processes, each in its own scratch directory
(ExoRT resolves its data files from an absolute root, so only the local
`iofiles/` and `exocol_config.nml` need to be per-worker).  Every completed grid
point is cached as a text file, so the script is resume-safe: re-running skips
work already done.

Cache files are named by the PHYSICAL VALUES of the column, not by axis index.
Index names would be reused across a regrid -- change the CO2 axis from 14
nodes to 11 and `col_000_005` still exists but now refers to a different CO2
mixing ratio, so the resume would load the wrong physics into the new grid and
report a full cache hit.  With value names a regridded axis simply misses and
recomputes.  (Caches written before this change use the old index names and
are ignored; tools/migrate_cache_names.py renames them given the table whose
axes they were built on.)

USAGE
-----
    python tools/make_radiation_table.py [--out FILE] [--workers N]
                                         [--exe PATH] [--star FILE]
                                         [--cache DIR] [--dry-run]
                                         [--ch4-axis | --ch4 LIST]

The grid and the physical assumptions are the CONFIG block below.
"""

import argparse
import math
import os
import shutil
import subprocess
import sys
import time
from multiprocessing import Pool

import numpy as np
import h5py

# ---------------------------------------------------------------------------
# CONFIG
# ---------------------------------------------------------------------------

EXOCOL_ROOT = '/models/ExoColumn'
# PVER=200 layers.  Vertical convergence at the demanding corner of the grid
# (10 bar, Ts = 360 K) is 255.6 / 253.1 / 252.5 / 252.1 W/m^2 for 70 / 140 /
# 200 / 400 layers, so 200 sits within ~1 W/m^2 of the converged OLR at ~6 s
# per grid point; see tools/check_table_convergence.py.
#
# exocol_sweepts200_psrh.exe adds &exocol_nml::variable_ps_rh (written below
# only when RH < 1).  With the switch off it reproduces exocol_sweepts200.exe
# bit for bit (all sweep records, RH 1 and 0.8, verified 2026-09-15), so it is a
# drop-in replacement for every build.
DEFAULT_EXE = os.path.join(EXOCOL_ROOT, 'run', 'exocol_sweepts200_psrh.exe')

# Host-star SED, as an ExoRT n68 solar file name (data/solar/ in ExoRT).
# 'blackbody_3000K_n68.nc' is the SAMOSA host star; 'G2V_SUN_n68.nc' is the Sun.
DEFAULT_STAR = 'blackbody_3000K_n68.nc'

# ---- Grid ------------------------------------------------------------------
# Dry surface pressure [bar].  SAMOSA spans 0.1-10 bar; the grid carries a
# margin either side so that the protocol's endpoints are interior points.
PRESSURES = np.array([0.05, 0.1, 0.2, 0.35, 0.5, 0.75, 1.0, 1.5,
                      2.0, 3.0, 5.0, 7.5, 10.0, 15.0, 20.0])

# CO2 mixing ratio of dry air.  SAMOSA fixes CO2 at 400 ppm, but HEXTOR's
# carbonate-silicate cycle needs the full range.  This axis is the most
# forgiving of the five, which is why it carries the coarsest spacing: measured
# with CH4 present, 11 nodes (0.57 dex) leaves 0.57 W/m^2 in OLR and 0.0023 in
# albedo, against 0.375 / 0.0015 for the 14 nodes (0.44 dex) used before the
# CH4 axis existed -- both inside the ~1 W/m^2 and 0.003 budgets.  Since the
# axis multiplies the whole build, the 11 nodes buy back 11 h of the CH4
# axis's cost.  12 nodes (0.52 dex, 0.497 / 0.0020) is the middle option.
FCO2 = np.logspace(-6.0, -0.3, 11)

# CH4 mixing ratio of dry air, used only with --ch4-axis.  Spacing was
# measured the same way as the other axes (see the note in notes/
# ch4_lookup_table.md): over 1e-8 to 1e-1 the PLANETARY ALBEDO, not the OLR, is
# what sets the node count, because CH4 absorbs in the near-infrared.  Against
# a 0.25 dex reference under the solar SED, uniform 0.75 dex leaves 0.0061 in
# albedo -- twice the zenith axis's 0.003 -- while the spacing below (1 dex up
# to 1e-6, where the response is nearly flat, then 0.5 dex through the steep
# part) holds 0.0023 in albedo and 0.50 W/m^2 in OLR on 13 nodes instead of 15.
#
# The floor is 1e-8, not zero, so that the log10 interpolation is well posed.
# That is radiatively indistinguishable from CH4-free: at 1 bar, 400 ppm CO2
# and 280 K, 1e-8 of CH4 moves OLR by 0.08 W/m^2 and the albedo by 1e-4.
CH4_NODES = np.concatenate([
    np.array([1.0e-8, 1.0e-7, 1.0e-6]),
    np.logspace(-5.5, -1.0, 10),
])

# Surface temperature [K].
TEMPERATURES = np.arange(180.0, 421.0, 10.0)

# Solar zenith angle, as cos(z) nodes; stored in the table as degrees
# (ascending) to match the getPALB interface.  Spacing is roughly GEOMETRIC in
# cos(z) rather than uniform: tools/check_table_interpolation.py shows the
# planetary albedo's curvature in this axis is concentrated at the limb, where
# a uniform grid left a ~0.009 albedo interpolation error, against ~0.0005 on
# the surface-albedo axis.  Geometric spacing spreads that error evenly.
MU_NODES = np.array([1.0, 0.9, 0.78, 0.65, 0.52, 0.42, 0.33, 0.25, 0.19,
                     0.13, 0.09, 0.05])

# Broadband surface albedo.  Extends down to open-ocean values (the v1 table
# floored at 0.2, well above HEXTOR's ocnalb = 0.06).  This axis is nearly
# linear, so it needs far fewer nodes than the zenith axis.
ALBEDOS = np.array([0.0, 0.05, 0.1, 0.2, 0.3, 0.45, 0.6, 0.8, 1.0])

# ---- Physical assumptions --------------------------------------------------
T_STRATO_REF = 200.0   # stratospheric temperature cap [K]; see t_strato_for()
RH = 1.0               # tropospheric relative humidity of the prescribed profile
PTOP_RATIO = 1.0e-5    # model top as a fraction of the dry surface pressure
N2_FILL = True         # background gas is N2, filling whatever CO2 leaves
CO2_CONDENSE = True    # cap the adiabat at the CO2 saturation curve (Kasting 1991)
CP_CO2_TDEP = True     # temperature-dependent cp for the CO2 share of dry air
H2O_EOS = 'nonideal'   # Kasting (1988) non-ideal steam adiabat (needed hot-end)
H2O_CONTINUUM = 'mtckd'  # MT_CKD 3.3, as used by ExoCAM in the SAMOSA protocol

# This machine is a 6-core Xeon with hyperthreading (12 logical CPUs).  The
# workload is memory-bandwidth bound — each process carries a ~240 MB ExoRT
# working set — so throughput saturates at ~38 records/s around 3-4 workers
# and DEGRADES beyond that (12 workers measured 9 records/s, worse than a
# single process).  Raise this only on a machine with real cores to spare.
NWORKERS = 4
TIMEOUT = 600          # seconds per ExoColumn run


def t_strato_for(ts):
    """Stratospheric temperature cap for a surface temperature ts.

    A fixed 200 K cap would put the whole atmosphere ABOVE the surface
    temperature for the cold end of the table (ExoColumn then builds an
    isothermal column at t_strato, i.e. a strong inversion over a 180 K
    surface).  Capping at ts instead gives the physically sensible cold-column
    limit: isothermal at the surface temperature.
    """
    return min(T_STRATO_REF, ts)


# ---------------------------------------------------------------------------
# NAMELIST
# ---------------------------------------------------------------------------

NML = """\
&exocol_nml
  flux_only       = .true.
  o3_profile      = 'none'
  msdist          = 1.0
  solar_file      = '{star}'
  h2o_eos         = '{h2o_eos}'
  h2o_continuum   = '{continuum}'
  cold_trap_phase = 'ice'
  variable_ps     = .true.
  co2_condense    = {co2_condense}
  cp_co2_tdep     = {cp_co2_tdep}
/
&exocol_sweep
  sweep_mode      = .true.
  sweep_coszrs    = {mus}
  sweep_albedo    = {albs}
  sweep_ts        = {temps}
  sweep_t_strato  = {tstrats}
  sweep_outfile   = '{sweepfile}'
/
&exocol_init
  input_file = ''
  ts         = {ts:.4f}
  t_strato   = {t_strato:.4f}
  p_top      = {p_top:.6e}
  rh_init    = {rh:.4f}
  coszrs     = 0.5
  asdir      = 0.3
  asdif      = 0.3
  aldir      = 0.3
  aldif      = 0.3
/
&exocol_composition
  ps       = {ps_pa:.8e}
  co2_vmr  = {co2_vmr:.8e}
  n2_vmr   = {n2_vmr:.8e}
  o2_vmr   = 0.0
  ar_vmr   = 0.0
  ch4_vmr  = {ch4_vmr:.8e}
  o3_vmr   = 0.0
/
"""


def build_namelist(pdry_bar, fco2, temps, star, sweepfile,
                   mus=MU_NODES, albs=ALBEDOS, rh=RH, ptop_ratio=PTOP_RATIO,
                   fch4=0.0):
    """Render the ExoColumn namelist for one (p_dry, fCO2, fCH4) column.

    `temps` is the surface-temperature axis swept inside the single ExoColumn
    process (a scalar is accepted for the single-column diagnostic scripts).
    &exocol_init::ts is set to the first entry so that a build without the
    sweep_ts support still produces that column rather than failing silently.

    N2 fills what CO2 and CH4 leave.  ExoColumn only WARNS when the dry VMRs
    fail to sum to one and then uses them as given, so the fill has to be right
    here: leaving CH4 out of it would quietly add its partial pressure on top
    of p_dry and break the table's own pressure coordinate.
    """
    temps = np.atleast_1d(np.asarray(temps, dtype=float))
    ps_pa = pdry_bar * 1.0e5
    n2_vmr = (1.0 - fco2 - fch4) if N2_FILL else 1.0
    nml = NML.format(
        star=star,
        h2o_eos=H2O_EOS,
        continuum=H2O_CONTINUUM,
        co2_condense='.true.' if CO2_CONDENSE else '.false.',
        cp_co2_tdep='.true.' if CP_CO2_TDEP else '.false.',
        mus=', '.join('%.6f' % m for m in mus),
        albs=', '.join('%.6f' % a for a in albs),
        temps=', '.join('%.4f' % t for t in temps),
        tstrats=', '.join('%.4f' % t_strato_for(t) for t in temps),
        sweepfile=sweepfile,
        ts=temps[0],
        t_strato=t_strato_for(temps[0]),
        p_top=ps_pa * ptop_ratio,
        rh=rh,
        ps_pa=ps_pa,
        co2_vmr=fco2,
        ch4_vmr=fch4,
        n2_vmr=n2_vmr,
    )
    if rh < 1.0:
        # Sub-saturated column: carry rh*esat, not esat, on top of p_dry, so
        # the column's dry mass is p_dry (the table's pressure coordinate) and
        # the surface vapour matches HEXTOR's moist static energy diffusion.
        # Without it the pressure coordinate holds (1-rh)*esat(Ts) of extra dry
        # gas: up to 7x the real dry gas in a 0.05 bar column at 390 K.
        nml = nml.replace('  variable_ps     = .true.\n',
                          '  variable_ps     = .true.\n'
                          '  variable_ps_rh  = .true.\n')
        assert 'variable_ps_rh' in nml
    return nml


# ---------------------------------------------------------------------------
# RUNNING
# ---------------------------------------------------------------------------

_WORKER = {}


def _init_worker(workroot, exe):
    """Give each pool process its own scratch directory and ExoColumn copy."""
    wid = os.getpid()
    wdir = os.path.join(workroot, 'w%d' % wid)
    os.makedirs(os.path.join(wdir, 'iofiles'), exist_ok=True)
    # ExoRT finds its data through an absolute root, so a link to the binary
    # is all the worker directory needs beyond its own iofiles/.
    rundir = os.path.join(wdir, 'run')
    if not os.path.exists(rundir):
        os.makedirs(rundir, exist_ok=True)
    exelink = os.path.join(rundir, os.path.basename(exe))
    if not os.path.exists(exelink):
        os.symlink(exe, exelink)
    _WORKER['dir'] = wdir
    _WORKER['exe'] = exelink


def cache_name(pdry, fco2, fch4=None, rh=1.0):
    """Cache file name for one column, keyed by its physical values.

    Formatted to a fixed number of significant digits so the same grid point
    always maps to the same name; see the note on index names in the module
    docstring.  A sub-saturated column carries its RH in the name, so a build
    pointed at an RH = 1 cache cannot silently reuse saturated columns (RH = 1
    names are unchanged, keeping existing caches valid).
    """
    name = 'col_p%.6e_c%.6e' % (pdry, fco2)
    if fch4 is not None:
        name += '_m%.6e' % fch4
    if rh != 1.0:
        name += '_rh%.4f' % rh
    return name + '.txt'


# Record layout written by exocol_sweep::run_flux_sweep.
C_TS, C_PS, C_MWDRY, C_QSRF, C_MU, C_ALB, C_OLR, C_SWDN, C_SWUP, C_PALB = range(10)


def parse_sweep(path):
    """Read an ExoColumn sweep file into (header dict, records array)."""
    head = {}
    rows = []
    with open(path) as f:
        for line in f:
            if line.startswith('#'):
                if '=' in line:
                    key, val = line[1:].split('=', 1)
                    try:
                        head[key.strip()] = float(val)
                    except ValueError:
                        pass
                continue
            parts = line.split()
            if len(parts) == 10:
                rows.append([float(x) for x in parts])
    return head, np.array(rows)


def run_column(task):
    """Run one (p_dry, fCO2, fCH4) column over the whole temperature axis.

    Returns (key, olr[ntmp], palb[ntmp, nzen, nsab], err); on failure the
    arrays are None and err carries the reason.
    """
    key, pdry, fco2, fch4, temps, star, cachefile, rh = task
    wdir = _WORKER['dir']
    exe = _WORKER['exe']
    sweepfile = 'iofiles/sweep.txt'

    nml = build_namelist(pdry, fco2, temps, star, sweepfile, fch4=fch4, rh=rh)
    with open(os.path.join(wdir, 'exocol_config.nml'), 'w') as f:
        f.write(nml)

    try:
        res = subprocess.run([exe], cwd=wdir, capture_output=True, text=True,
                             timeout=TIMEOUT)
    except subprocess.TimeoutExpired:
        return key, None, None, 'timeout'

    if res.returncode != 0:
        tail = (res.stderr or res.stdout or '')[-300:]
        return key, None, None, 'rc=%d %s' % (res.returncode, tail.replace('\n', ' '))

    spath = os.path.join(wdir, sweepfile)
    if not os.path.exists(spath):
        return key, None, None, 'no sweep output'

    # Keep the raw records so a re-run can skip this column.
    shutil.copyfile(spath, cachefile)

    head, rows = parse_sweep(spath)
    olr, palb, err = assemble_column(rows, temps)
    return key, olr, palb, err


def assemble_column(rows, temps):
    """Turn one sweep file's records into (olr[ntmp], palb[ntmp, nzen, nsab]).

    The temperature axis is passed in rather than read from the module global:
    worker processes re-import this module (multiprocessing does not
    necessarily fork), so a --temps override on the command line would not
    reach them otherwise.
    """
    temps = np.asarray(temps, dtype=float)
    ntmp, nzen, nsab = len(temps), len(MU_NODES), len(ALBEDOS)
    if rows.ndim != 2 or rows.shape[0] != ntmp * nzen * nsab:
        got = 0 if rows.ndim != 2 else rows.shape[0]
        return None, None, 'expected %d records, got %d' % (ntmp * nzen * nsab, got)

    olr = np.full(ntmp, np.nan)
    palb = np.full((ntmp, nzen, nsab), np.nan)

    for r in rows:
        it = int(np.argmin(np.abs(temps - r[C_TS])))
        iz = int(np.argmin(np.abs(MU_NODES - r[C_MU])))
        ia = int(np.argmin(np.abs(ALBEDOS - r[C_ALB])))
        palb[it, iz, ia] = r[C_PALB]
        # The longwave cannot depend on zenith angle or surface albedo, so
        # every record of a given Ts must agree; disagreement means the sweep
        # went wrong, not that the physics is subtle.
        if np.isnan(olr[it]):
            olr[it] = r[C_OLR]
        elif abs(olr[it] - r[C_OLR]) > 1e-9 * max(abs(olr[it]), 1.0):
            return None, None, 'OLR varies within a Ts block at %.0f K' % r[C_TS]

    if np.any(~np.isfinite(olr)) or np.any(~np.isfinite(palb)):
        return None, None, 'missing records in sweep block'
    if np.any(palb < 0.0):
        return None, None, 'invalid planetary albedo in sweep block'
    return olr, palb, None


# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--out', default='model/radiation/radiation_N2_CO2_3000K_p.h5',
                    help='output HDF5 table')
    ap.add_argument('--star', default=DEFAULT_STAR, help='ExoRT n68 solar file name')
    ap.add_argument('--exe', default=DEFAULT_EXE, help='ExoColumn executable')
    ap.add_argument('--workers', type=int, default=NWORKERS)
    ap.add_argument('--cache', default=None,
                    help='cache directory (default: <out>.cache)')
    ap.add_argument('--work', default=None,
                    help='scratch directory for worker runs (default: <cache>/work)')
    ap.add_argument('--dry-run', action='store_true',
                    help='report grid size and cost estimate, then exit')
    ap.add_argument('--pressures', default=None,
                    help='comma-separated dry surface pressures [bar], '
                         'overriding the built-in grid (for reduced-grid '
                         'validation tables)')
    ap.add_argument('--fco2', default=None,
                    help='comma-separated CO2 mixing ratios, overriding the grid')
    ap.add_argument('--ch4-axis', action='store_true',
                    help='resolve CH4 as a table axis on the built-in CH4_NODES '
                         'grid, writing a /ch4 dataset.  Without this (and '
                         'without --ch4) the table is CH4-free and carries no '
                         'CH4 axis, exactly as before.')
    ap.add_argument('--ch4', default=None,
                    help='comma-separated CH4 mixing ratios to use as the CH4 '
                         'axis, overriding CH4_NODES (implies --ch4-axis)')
    ap.add_argument('--temps', default=None,
                    help='comma-separated surface temperatures [K], '
                         'overriding the grid')
    ap.add_argument('--rh', type=float, default=RH,
                    help='tropospheric relative humidity of the prescribed '
                         'column (default %.2f).  Match HEXTOR\'s rhmoist when '
                         'the EBM diffuses moist static energy.' % RH)
    args = ap.parse_args()
    if not 0.0 < args.rh <= 1.0:
        print('--rh must be in (0, 1]')
        return 2

    global PRESSURES, FCO2, TEMPERATURES, CH4_NODES
    if args.pressures:
        PRESSURES = np.array([float(x) for x in args.pressures.split(',')])
    if args.fco2:
        FCO2 = np.array([float(x) for x in args.fco2.split(',')])
    if args.temps:
        TEMPERATURES = np.array([float(x) for x in args.temps.split(',')])
    if args.ch4:
        CH4_NODES = np.array([float(x) for x in args.ch4.split(',')])

    # A CH4-resolved table carries the axis; otherwise the run is CH4-free and
    # the axis collapses to a single implicit level that is never written, so
    # the output keeps the previous v2 layout byte for byte.
    ch4_axis = bool(args.ch4_axis or args.ch4)
    if ch4_axis:
        if np.any(CH4_NODES <= 0.0):
            print('CH4 axis levels must be positive: the axis is interpolated '
                  'in log10, so use a small floor (1e-8) instead of zero')
            return 2
        CH4 = np.sort(CH4_NODES)
    else:
        CH4 = np.array([0.0])

    npre, nco2, ntmp = len(PRESSURES), len(FCO2), len(TEMPERATURES)
    nch4 = len(CH4)
    nzen, nsab = len(MU_NODES), len(ALBEDOS)
    ncase = npre * nco2 * nch4 * ntmp

    zenith_deg = np.degrees(np.arccos(np.clip(MU_NODES, -1.0, 1.0)))
    order = np.argsort(zenith_deg)          # ascending zenith angle for getPALB
    zenith_deg = zenith_deg[order]

    print('grid : npre=%d  nco2=%d  nch4=%d  ntmp=%d  nzen=%d  nsab=%d'
          % (npre, nco2, nch4, ntmp, nzen, nsab))
    print('       %d grid points x %d (zenith x albedo) = %d radiation calls'
          % (ncase, nzen * nsab, ncase * nzen * nsab))
    print('       p_dry %.3g..%.3g bar   fco2 %.3g..%.3g   T %.0f..%.0f K'
          % (PRESSURES[0], PRESSURES[-1], FCO2[0], FCO2[-1],
             TEMPERATURES[0], TEMPERATURES[-1]))
    if ch4_axis:
        print('       fch4  %.3g..%.3g (%d levels)' % (CH4[0], CH4[-1], nch4))
    else:
        print('       fch4  none (CH4-free table, no /ch4 axis)')
    print('       RH    %.2f' % args.rh)
    # One ExoColumn process per (p_dry, fCO2, fCH4) column, sweeping Ts
    # internally.  0.101 s per radiation call and 1.6 s of fixed ExoRT
    # initialisation are what the 210-column Sun build actually measured
    # (4.0 h of wall time on 4 workers); the older 0.049 s figure was optimistic.
    ncol = npre * nco2 * nch4
    est = ncol * (1.6 + ntmp * nzen * nsab * 0.101) / max(args.workers, 1)
    print('       %d ExoColumn processes (one per p_dry x fCO2 x fCH4 column)'
          % ncol)
    print('       estimated wall time on %d workers: %.1f h (before contention)'
          % (args.workers, est / 3600.0))
    print('       estimated table size: %.0f MB'
          % (npre * nco2 * nch4 * ntmp * (1 + nzen * nsab) * 8 / 1024.0**2))
    if args.dry_run:
        return 0

    cache = os.path.abspath(args.cache or (args.out + '.cache'))
    work = os.path.abspath(args.work or os.path.join(cache, 'work'))
    os.makedirs(cache, exist_ok=True)
    os.makedirs(work, exist_ok=True)

    olr = np.full((npre, nco2, nch4, ntmp), np.nan)
    palb = np.full((npre, nco2, nch4, ntmp, nzen, nsab), np.nan)

    # Collect the work, reusing any cached columns.
    tasks = []
    ncached = 0
    for ip in range(npre):
        for ic in range(nco2):
            for im in range(nch4):
                key = (ip, ic, im)
                cachefile = os.path.join(cache, cache_name(
                    PRESSURES[ip], FCO2[ic], CH4[im] if ch4_axis else None,
                    rh=args.rh))
                if os.path.exists(cachefile):
                    head, rows = parse_sweep(cachefile)
                    o, p, err = assemble_column(rows, TEMPERATURES)
                    if err is None:
                        olr[ip, ic, im] = o
                        palb[ip, ic, im] = p
                        ncached += 1
                        continue
                tasks.append((key, float(PRESSURES[ip]), float(FCO2[ic]),
                              float(CH4[im]), TEMPERATURES, args.star,
                              cachefile, args.rh))

    print('       %d columns cached, %d to run' % (ncached, len(tasks)))

    failures = []
    if tasks:
        t0 = time.time()
        done = 0
        with Pool(args.workers, initializer=_init_worker,
                  initargs=(work, os.path.abspath(args.exe))) as pool:
            for key, o, p, err in pool.imap_unordered(run_column, tasks, chunksize=1):
                done += 1
                if err is not None:
                    failures.append((key, err))
                else:
                    olr[key] = o
                    palb[key] = p
                if done % 5 == 0 or done == len(tasks):
                    el = time.time() - t0
                    rate = done / el
                    print('  %4d/%4d columns  %.2f col/s  eta %.1f min  (%d failed)'
                          % (done, len(tasks), rate,
                             (len(tasks) - done) / max(rate, 1e-9) / 60.0,
                             len(failures)), flush=True)

    nbad = int(np.sum(~np.isfinite(olr)))
    if nbad:
        print('WARNING: %d of %d grid points have no OLR' % (nbad, olr.size))
        for key, err in failures[:20]:
            ip, ic, im = key
            print('   p=%.3g fco2=%.3g fch4=%.3g : %s'
                  % (PRESSURES[ip], FCO2[ic], CH4[im], err))
        if len(failures) > 20:
            print('   ... and %d more' % (len(failures) - 20))

    # Reorder the zenith axis to ascending degrees and convert to table units.
    palb = palb[:, :, :, :, order, :]
    olr_mw = olr * 1000.0     # W/m^2 -> mW/m^2 (v1 convention; see radiation.f90)

    # A CH4-free build drops the degenerate CH4 dimension so the file is the
    # v2 layout the reader already knows, with no /ch4 dataset to detect.
    if not ch4_axis:
        olr_mw = olr_mw[:, :, 0, :]
        palb = palb[:, :, 0, :, :, :]

    with h5py.File(args.out, 'w') as f:
        f.create_dataset('olr', data=olr_mw)
        f.create_dataset('palb', data=palb)
        f.create_dataset('pressure', data=PRESSURES)
        f.create_dataset('fco2', data=FCO2)
        if ch4_axis:
            f.create_dataset('ch4', data=CH4)
        f.create_dataset('temperature', data=TEMPERATURES)
        f.create_dataset('zenith', data=zenith_deg)
        f.create_dataset('surfalb', data=ALBEDOS)
        f.attrs['format'] = ('hextor-radiation-v3' if ch4_axis
                             else 'hextor-radiation-v2')
        f.attrs['source'] = 'ExoColumn (flux_only + sweep_mode) / ExoRT n68equiv'
        f.attrs['star'] = args.star
        f.attrs['olr_units'] = 'mW m-2 (divide by 1000 for W m-2)'
        f.attrs['pressure_units'] = \
            'bar, DRY surface pressure (pN2 + pCO2 + pCH4)'
        if ch4_axis:
            f.attrs['ch4'] = ('CH4 mixing ratio of dry air; log10-interpolated. '
                              'The %.0e floor stands in for CH4-free (worth '
                              '0.08 W/m2 in OLR).' % CH4[0])
            # Organic haze is the known gap at the top of this axis: ExoRT's
            # calc_opd_mod.F90 carries no haze optics (only a note to hook up
            # CARMA aerosols), so a real Archean atmosphere above CH4/CO2 ~ 0.1
            # would have an anti-greenhouse the table cannot represent.
            f.attrs['ch4_caveat'] = ('clear sky, no organic haze: biased warm '
                                     'where CH4/CO2 exceeds ~0.1')
        if args.rh < 1.0:
            f.attrs['water'] = ('variable_ps + variable_ps_rh: total ps = p_dry '
                                '+ RH*esat(Ts), dry mass = p_dry; saturated moist '
                                'adiabat in the non-condensable pressure, '
                                'vapour RH*esat(T), RH = %.2f' % args.rh)
        else:
            f.attrs['water'] = ('variable_ps: total ps = p_dry + esat(Ts); '
                                'prescribed moist adiabat, RH = %.2f' % args.rh)
        f.attrs['rh'] = args.rh
        f.attrs['exocolumn_exe'] = os.path.basename(args.exe)
        f.attrs['t_strato'] = 'min(%.1f K, Ts)' % T_STRATO_REF
        f.attrs['p_top_ratio'] = PTOP_RATIO
        f.attrs['h2o_eos'] = H2O_EOS
        f.attrs['h2o_continuum'] = H2O_CONTINUUM
        f.attrs['co2_condense'] = str(CO2_CONDENSE)
        f.attrs['clouds'] = 'none (clear sky)'
        f.attrs['generated'] = time.strftime('%Y-%m-%d %H:%M:%S')

    print('wrote %s' % args.out)
    return 1 if nbad else 0


if __name__ == '__main__':
    sys.exit(main())
