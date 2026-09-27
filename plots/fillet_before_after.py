#!/usr/bin/env python3
"""Before/after figure for HEXTOR's FILLET benchmarks and experiments.

Three states of the submission are compared:
  archived  the files in projectcuisines/fillet Results/hextor, byte-identical
            to commit 852fd03 (June 2024), read from git history;
  4.2.0     the July 2026 regeneration under fillet/<exp>/ (zenith, CO2-lookup
            and CO2-coordinate fixes, Benchmark 1 tuning carried into the
            untuned runs), kept in the tree for comparison;
  now       fillet/Results/hextor/ written by tools/run_fillet.py.

    python plots/fillet_before_after.py            # run from the repository root

Writes plots/fillet_before_after.png and .pdf.  Panels whose files are
missing are skipped, so the script also runs before the new set exists.
"""
import os
import subprocess
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

sys.dont_write_bytecode = True

ARCHIVED_COMMIT = "852fd03"
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "plots", "fillet_before_after")

# Palette: the documented steps used for the AVALON figure.  check_ramp()
# below verifies the ordinal ramp (monotone OKLab lightness, steps >= 0.06,
# light end >= 2:1 on white); text wears ink, never a series colour.
INK, INK2, MUTED, AXIS = "#0b0b0b", "#52514e", "#898781", "#c3c2b7"
BLUE, ORANGE = "#2a78d6", "#eb6834"
RAMP = ["#86b6ef", "#3987e5", "#1c5cab", "#0d366b"]
STATES = ["ice-free", "polar caps", "ice belt", "snowball"]
BEFORE_ALPHA = 0.45


def _oklab_L(hexcol):
    r, g, b = (int(hexcol[i:i + 2], 16) / 255 for i in (1, 3, 5))
    lin = [c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4 for c in (r, g, b)]
    l = 0.4122214708 * lin[0] + 0.5363325363 * lin[1] + 0.0514459929 * lin[2]
    m = 0.2119034982 * lin[0] + 0.6806995451 * lin[1] + 0.1073969566 * lin[2]
    s = 0.0883024619 * lin[0] + 0.2817188376 * lin[1] + 0.6299787005 * lin[2]
    l, m, s = np.cbrt(l), np.cbrt(m), np.cbrt(s)
    return 0.2104542553 * l + 0.7936177850 * m - 0.0040720468 * s


def _rel_lum(hexcol):
    r, g, b = (int(hexcol[i:i + 2], 16) / 255 for i in (1, 3, 5))
    lin = [c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4 for c in (r, g, b)]
    return 0.2126 * lin[0] + 0.7152 * lin[1] + 0.0722 * lin[2]


def check_ramp():
    L = [_oklab_L(c) for c in RAMP]
    dL = [L[i] - L[i + 1] for i in range(len(L) - 1)]
    contrast = 1.05 / (_rel_lum(RAMP[0]) + 0.05)
    ok = all(d >= 0.06 for d in dL) and contrast >= 2.0
    print(f"ramp L = {np.round(L, 3).tolist()}  dL = {np.round(dL, 3).tolist()}  "
          f"light-end contrast = {contrast:.2f}:1  -> {'PASS' if ok else 'FAIL'}")
    if not ok:
        sys.exit("ramp check failed")


# ---- data ------------------------------------------------------------------
def parse(text):
    """FILLET .dat -> dict of column arrays, from the last '# Case ...' header."""
    headers, rows = None, []
    for line in text.splitlines():
        s = line.strip()
        if not s:
            continue
        if s.startswith("#"):
            tokens = s.lstrip("# ").split()
            if tokens and tokens[0] == "Case":
                headers = tokens
        else:
            rows.append(s.split())
    if headers is None or not rows:
        return None
    cols = {h: [] for h in headers}
    for r in rows:
        for h, v in zip(headers, r):
            cols[h].append(float(v))
    return {k: np.array(v) for k, v in cols.items()}


def parse_lat(text):
    lat, T = [], []
    for line in text.splitlines():
        s = line.strip()
        if s and not s.startswith("#"):
            p = s.split()
            lat.append(float(p[0])); T.append(float(p[1]))
    return np.array(lat), np.array(T)


def git_show(path):
    r = subprocess.run(["git", "-C", ROOT, "show", f"{ARCHIVED_COMMIT}:{path}"],
                       capture_output=True, text=True)
    return r.stdout if r.returncode == 0 else None


def read_file(path):
    full = os.path.join(ROOT, path)
    if not os.path.exists(full):
        return None
    with open(full) as f:
        return f.read()


def old_layout(exp):
    """archived and 4.2.0 files: one file per experiment, exp3/exp4 with the
    cold branch first (57 or 50 rows) then the warm branch."""
    base = exp.split("_")[0]
    return f"fillet/{base}/global_output_HEXTOR_{base}.dat"


def load_global(exp, version):
    if version == "archived":
        text = git_show(old_layout(exp))
    elif version == "4.2.0":
        text = read_file(old_layout(exp))
    else:
        text = read_file(f"fillet/Results/hextor/{exp}/global_output.dat")
    g = parse(text) if text else None
    if g is None:
        return None
    if version != "now" and "_" in exp:
        n = len(g["Case"]) // 2
        sel = slice(0, n) if exp.endswith("cold") else slice(n, 2 * n)
        g = {k: v[sel] for k, v in g.items()}
    return g


def load_lat(exp, version):
    if version == "archived":
        text = git_show(f"fillet/{exp}/lat_output_HEXTOR_{exp}.dat")
    elif version == "4.2.0":
        text = read_file(f"fillet/{exp}/lat_output_HEXTOR_{exp}.dat")
    else:
        text = read_file(f"fillet/Results/hextor/{exp}/case_0/lat_output.dat")
    return parse_lat(text) if text else None


def hemi_state(pmax, pmin):
    if pmax == 90.0 and pmin == 90.0:
        return 0      # ice-free
    if pmax == 90.0 and pmin == 0.0:
        return 3      # ice to the equator
    if pmax == 90.0:
        return 1      # cap
    return 2          # belt


def climate_state(g, i):
    nh = hemi_state(g["IceLineNMax"][i], g["IceLineNMin"][i])
    sh = hemi_state(-g["IceLineSMin"][i], -g["IceLineSMax"][i])
    if nh == 2 or sh == 2:
        return 2
    if nh == 3 and sh == 3:
        return 3
    if 1 in (nh, sh) or 3 in (nh, sh):
        return 1
    return 0


def state_grid(g):
    S = np.unique(np.round(g["Inst"], 4))
    O = np.unique(g["Obl"])
    grid = np.full((len(O), len(S)), -1)
    for i in range(len(g["Inst"])):
        grid[np.searchsorted(O, g["Obl"][i]), np.searchsorted(S, round(g["Inst"][i], 4))] = climate_state(g, i)
    return S, O, grid


def edges(v):
    m = 0.5 * (v[1:] + v[:-1])
    return np.concatenate([[v[0] - (m[0] - v[0])], m, [v[-1] + (v[-1] - m[-1])]])


def describe_ice(g, i=0):
    s = climate_state(g, i)
    if s == 0:
        return "ice-free"
    if s == 3:
        return "snowball"
    if s == 1:
        return f"caps to {g['IceLineNMin'][i]:.1f}N / {abs(g['IceLineSMax'][i]):.1f}S"
    return f"belt {g['IceLineNMin'][i]:.0f}-{g['IceLineNMax'][i]:.0f}"


# ---- style -----------------------------------------------------------------
def style_axis(ax):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(AXIS)
        ax.spines[side].set_linewidth(0.8)
    ax.tick_params(colors=INK2, labelsize=8, width=0.8, length=3)
    ax.xaxis.label.set_color(INK2); ax.yaxis.label.set_color(INK2)
    ax.grid(False)


def note(ax, text, loc="lower left"):
    x, y, ha, va = {"lower left": (0.03, 0.04, "left", "bottom"),
                    "upper left": (0.03, 0.96, "left", "top"),
                    "lower right": (0.97, 0.04, "right", "bottom"),
                    "upper right": (0.97, 0.96, "right", "top")}[loc]
    ax.text(x, y, text, transform=ax.transAxes, ha=ha, va=va, fontsize=7.4, color=INK2,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="none", alpha=0.85))


# ---- panels ----------------------------------------------------------------
def panel_benchmark(ax, tag, title):
    series = [("archived", BEFORE_ALPHA, "--", 1.5), ("4.2.0", 0.75, ":", 1.5), ("now", 1.0, "-", 2.0)]
    ax.axhline(273.15, color=AXIS, linewidth=0.8)
    lines = []
    for version, alpha, ls, lw in series:
        lat = load_lat(tag, version)
        g = load_global(tag, version)
        if lat is None or g is None:
            continue
        ax.plot(lat[0], lat[1], color=BLUE, alpha=alpha, linestyle=ls, linewidth=lw, label=version)
        lines.append(f"{version}: {g['Tglob'][0]:.1f} K, {describe_ice(g)}")
    ax.set_xlim(-90, 90); ax.set_xticks([-90, -60, -30, 0, 30, 60, 90])
    ax.set_ylim(230, 305)
    ax.set_xlabel("latitude"); ax.set_ylabel("annual-mean T (K)")
    ax.set_title(title, fontsize=10, color=INK, loc="left")
    ax.text(-88, 274.2, "273.15 K", fontsize=7, color=MUTED, ha="left", va="bottom")
    note(ax, "\n".join(lines))
    style_axis(ax)


def panel_sweep(ax_b, ax_n, tag, title):
    g_b, g_n = load_global(tag, "archived"), load_global(tag, "now")
    cmap = ListedColormap(RAMP); norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5], cmap.N)
    grids = {}
    for ax, g, sub in ((ax_b, g_b, f"{title}, archived"), (ax_n, g_n, f"{title}, now")):
        if g is None:
            ax.text(0.5, 0.5, "not available", transform=ax.transAxes, ha="center", color=MUTED)
            style_axis(ax); continue
        S, O, grid = state_grid(g)
        grids[ax] = (S, O, grid)
        ax.pcolormesh(edges(S), edges(O), grid, cmap=cmap, norm=norm, edgecolors="white", linewidth=0.6)
        if ax is ax_n and g_b is not None:
            S_b0 = np.unique(np.round(g_b["Inst"], 4))
            if len(S_b0) == len(S):
                S = S_b0 if np.allclose(S_b0, S, atol=0.006) else S
        ax.set_xlim(edges(S)[0], edges(S)[-1]); ax.set_ylim(-5, 95)
        ax.set_yticks([0, 30, 60, 90])
        ax.set_xlabel("instellation (S_earth)")
        ax.set_title(sub, fontsize=9.5, color=INK, loc="left")
        style_axis(ax)
    ax_b.set_ylabel("obliquity")
    if ax_b in grids and ax_n in grids:
        S_b, O_b, grid_b = grids[ax_b]; S_n, O_n, grid_n = grids[ax_n]
        if grid_b.shape == grid_n.shape and np.allclose(S_b, S_n, atol=0.006):
            changed = grid_b != grid_n
            yy, xx = np.nonzero(changed)
            ax_n.scatter(S_n[xx], O_n[yy], s=9, color=INK, edgecolors="white", linewidths=0.8, zorder=3)
            note(ax_n, f"{changed.sum()} of {changed.size} states changed", loc="upper right")
            return int(changed.sum()), int(changed.size)
        note(ax_n, "grids differ: the archived run\nused the Experiment 1a axis range", loc="upper right")
    return None


def branch(g, xkey):
    x, T = g[xkey], g["Tglob"]
    snow = np.array([climate_state(g, i) == 3 for i in range(len(x))])
    o = np.argsort(x)
    return x[o], T[o], snow[o]


def thresholds(warm, cold):
    """Glaciation: largest x at which the warm-start branch is a snowball;
    deglaciation: smallest x at which the cold-start branch is not (the
    code comparison's Table 10 convention)."""
    x, _, snow = warm
    glac = x[np.max(np.nonzero(snow)[0])] if snow.any() and not snow.all() else None
    x, _, snow = cold
    deglac = x[np.min(np.nonzero(~snow)[0])] if snow.any() and not snow.all() else None
    return glac, deglac


def panel_hysteresis(ax, base, title, xkey, xlabel, log=False):
    lines = []
    for version, alpha, ls, lw in (("archived", BEFORE_ALPHA, "--", 1.5), ("now", 1.0, "-", 2.0)):
        gw, gc = load_global(base + "_warm", version), load_global(base + "_cold", version)
        if gw is None or gc is None:
            continue
        warm, cold = branch(gw, xkey), branch(gc, xkey)
        ax.plot(warm[0], warm[1], color=ORANGE, alpha=alpha, linestyle=ls, linewidth=lw,
                label=f"warm-start branch, {version}")
        ax.plot(cold[0], cold[1], color=BLUE, alpha=alpha, linestyle=ls, linewidth=lw,
                label=f"cold-start branch, {version}")
        gl, dg = thresholds(warm, cold)
        fmt = (lambda v: "none" if v is None else (f"{v:.3g} ppm" if log else f"{v:.4g}"))
        s = f"{version}: glaciation {fmt(gl)}, deglaciation {fmt(dg)}"
        if not log and gl is not None and dg is not None:
            s += f", width {dg - gl:.3f}"
        lines.append(s)
    if log:
        ax.set_xscale("log")
    ax.set_xlabel(xlabel); ax.set_ylabel("Tglob (K)")
    ax.set_title(title, fontsize=10, color=INK, loc="left")
    style_axis(ax)
    ax.legend(fontsize=7.5, frameon=False, loc="upper left", labelcolor=INK2)
    if lines:
        note(ax, "\n".join(lines), loc="lower right")


# ---- figure ----------------------------------------------------------------
def main():
    check_ramp()
    plt.rcParams.update({"font.family": "sans-serif", "font.size": 8.5, "axes.titleweight": "normal",
                         "text.color": INK, "axes.labelsize": 8.5})
    fig = plt.figure(figsize=(13, 14.2), facecolor="white")
    gs = fig.add_gridspec(4, 12, height_ratios=[2.9, 2.5, 2.5, 3.0], hspace=0.55, wspace=0.75,
                          left=0.055, right=0.985, top=0.893, bottom=0.05)
    fig.suptitle("HEXTOR: FILLET benchmarks and experiments, archived submission against the re-file",
                 fontsize=13, color=INK, x=0.055, ha="left", y=0.975)
    fig.text(0.055, 0.952,
             f"Archived = FILLET Results/hextor (June 2024, commit {ARCHIVED_COMMIT}); 4.2.0 = the July 2026 regeneration "
             "kept under fillet/; now = this release.  Dashed and light: archived; dotted: 4.2.0; solid: now.\n"
             "Changes in between: zenith argument, CO2 lookup and coordinate, ghost cells, effective D filed, untuned "
             "Benchmarks 2/3 (cloudir 0, D = 0.5), Experiment 2a grid,\nconvergence 0.1 to 0.001 W/m2, and the numerics "
             "switches named in the file headers.",
             fontsize=8.5, color=INK2, va="top", linespacing=1.4)

    for k, (tag, title) in enumerate((("ben1", "Benchmark 1, tuned Earth (obliquity 23.5)"),
                                      ("ben2", "Benchmark 2, Table 4 defaults (obliquity 23.5)"),
                                      ("ben3", "Benchmark 3, Table 4 defaults (obliquity 60)"))):
        ax = fig.add_subplot(gs[0, 4 * k:4 * k + 4])
        panel_benchmark(ax, tag, title)
        if k == 0:
            ax.legend(fontsize=7.5, frameon=False, loc="upper right", labelcolor=INK2)

    changes, sweep_axes = {}, {1: [], 2: []}
    for row, pair in ((1, (("exp1", "Exp 1 (warm start)"), ("exp1a", "Exp 1a (warm start, a varies)"))),
                      (2, (("exp2", "Exp 2 (cold start)"), ("exp2a", "Exp 2a (cold start, a varies)")))):
        for k, (tag, title) in enumerate(pair):
            ax_b = fig.add_subplot(gs[row, 6 * k:6 * k + 3])
            ax_n = fig.add_subplot(gs[row, 6 * k + 3:6 * k + 6], sharey=ax_b)
            ax_n.tick_params(labelleft=False, left=False)
            changes[tag] = panel_sweep(ax_b, ax_n, tag, title)
            sweep_axes[row] += [ax_b, ax_n]

    handles = [Patch(facecolor=c, edgecolor="none", label=s) for c, s in zip(RAMP, STATES)]
    handles.append(Line2D([], [], marker="o", color="none", markerfacecolor=INK, markeredgecolor="white",
                          markersize=5, label="state changed since the archived submission"))
    hyst_axes = [fig.add_subplot(gs[3, 0:6]), fig.add_subplot(gs[3, 6:12])]
    panel_hysteresis(hyst_axes[0], "exp3", "Exp 3, instellation hysteresis (Benchmark 2 configuration)",
                     "Inst", "instellation (S_earth)")
    panel_hysteresis(hyst_axes[1], "exp4", "Exp 4, CO2 hysteresis (Benchmark 2 configuration)",
                     "XCO2", "CO2 (ppm)", log=True)

    fig.canvas.draw()
    r, inv = fig.canvas.get_renderer(), fig.transFigure.inverted()
    row3_bottom = min(inv.transform(ax.get_tightbbox(r))[0, 1] for ax in sweep_axes[2])
    row4_top = max(inv.transform(ax.get_tightbbox(r))[1, 1] for ax in hyst_axes)
    fig.legend(handles=handles, loc="center", bbox_to_anchor=(0.52, 0.5 * (row3_bottom + row4_top)), ncol=5,
               fontsize=8, frameon=False, labelcolor=INK2, handlelength=1.4, columnspacing=1.6)

    fig.savefig(OUT + ".png", dpi=200, facecolor="white")
    fig.savefig(OUT + ".pdf", facecolor="white")
    print("state changes:", {k: (f"{v[0]}/{v[1]}" if v else "n/a") for k, v in changes.items()})
    print("wrote", OUT + ".png", "and .pdf")


if __name__ == "__main__":
    main()
