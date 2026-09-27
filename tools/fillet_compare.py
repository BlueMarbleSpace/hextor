#!/usr/bin/env python3
"""Compare two states of HEXTOR's FILLET submission, experiment by experiment.

    python tools/fillet_compare.py archived now
    python tools/fillet_compare.py 4.2.0 now --dir fillet
    python tools/fillet_compare.py now experiments/fillet_eval/clock_730

A state is one of
    archived   the files in projectcuisines/fillet, byte-identical to commit
               852fd03, read from git history;
    4.2.0      the July 2026 regeneration kept under fillet/<exp>/;
    now        fillet/Results/hextor/ (or another --dir written by run_fillet.py);
    <path>     any directory holding Results/hextor/ in the archive layout.

For every experiment present in both: mean, RMS and largest change in Tglob,
the number of cases whose climate state changed, the state census
(ice-free / caps / belt / snowball), the largest north-south ice-edge
asymmetry, Experiment 3/4 thresholds and bistable widths, and the
warm-start-never-colder-than-cold-start check on matched grids.  Cases are
matched by (Inst, Obl, XCO2), so grids that differ (the archived
Experiment 2a) compare only where they overlap.
"""
import argparse
import os
import subprocess
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ARCHIVED_COMMIT = "852fd03"
EXPS = ["ben1", "ben2", "ben3", "exp1", "exp1a", "exp2", "exp2a",
        "exp3_cold", "exp3_warm", "exp4_cold", "exp4_warm"]
STATES = ["ice-free", "caps", "belt", "snowball"]


def parse(text):
    headers, rows = None, []
    for line in text.splitlines():
        s = line.strip()
        if not s:
            continue
        if s.startswith("#"):
            tok = s.lstrip("# ").split()
            if tok and tok[0] == "Case":
                headers = tok
        else:
            rows.append(s.split())
    if headers is None or not rows:
        return None
    cols = {h: np.array([float(r[i]) for r in rows]) for i, h in enumerate(headers)}
    return cols


def load(state, exp):
    base = exp.split("_")[0]
    if state == "archived":
        r = subprocess.run(["git", "-C", ROOT, "show", f"{ARCHIVED_COMMIT}:fillet/{base}/global_output_HEXTOR_{base}.dat"],
                           capture_output=True, text=True)
        text = r.stdout if r.returncode == 0 else None
    elif state == "4.2.0":
        p = os.path.join(ROOT, "fillet", base, f"global_output_HEXTOR_{base}.dat")
        text = open(p).read() if os.path.exists(p) else None
    else:
        d = os.path.join(ROOT, "fillet") if state == "now" else state
        p = os.path.join(d, "Results", "hextor", exp, "global_output.dat")
        text = open(p).read() if os.path.exists(p) else None
    g = parse(text) if text else None
    if g is None:
        return None
    if state in ("archived", "4.2.0") and "_" in exp:
        n = len(g["Case"]) // 2
        sel = slice(0, n) if exp.endswith("cold") else slice(n, 2 * n)
        g = {k: v[sel] for k, v in g.items()}
    return g


def hemi(pmax, pmin):
    if pmax == 90.0 and pmin == 90.0:
        return 0
    if pmax == 90.0 and pmin == 0.0:
        return 3
    if pmax == 90.0:
        return 1
    return 2


def state(g, i):
    nh = hemi(g["IceLineNMax"][i], g["IceLineNMin"][i])
    sh = hemi(-g["IceLineSMin"][i], -g["IceLineSMax"][i])
    if nh == 2 or sh == 2:
        return 2
    if nh == 3 and sh == 3:
        return 3
    if 1 in (nh, sh) or 3 in (nh, sh):
        return 1
    return 0


def key(g, i):
    return (round(g["Inst"][i], 3), round(g["Obl"][i], 1), round(g["XCO2"][i], 2))


def census(g):
    c = [0, 0, 0, 0]
    for i in range(len(g["Case"])):
        c[state(g, i)] += 1
    return c


def asym(g):
    return float(np.max(np.abs(g["IceLineNMin"] + g["IceLineSMax"])))


def thresholds(gw, gc, xkey):
    """Glaciation: largest x at which the warm-start branch is a snowball;
    deglaciation: smallest x at which the cold-start branch is not."""
    def snow(g):
        return np.array([state(g, i) == 3 for i in range(len(g["Case"]))])
    xw, sw = gw[xkey], snow(gw)
    xc, sc = gc[xkey], snow(gc)
    glac = xw[sw].max() if sw.any() and not sw.all() else None
    deglac = xc[~sc].min() if sc.any() and not sc.all() else None
    return glac, deglac


def pairs(a, b):
    """Matched case indices.  The old files print Inst with two decimals and
    XCO2 with one, so cases are paired by position when the grids have the
    same length and order (every experiment but the archived 2a), checked
    against the printed values, and by nearest forcing otherwise."""
    na, nb = len(a["Case"]), len(b["Case"])

    def same(i, j, tol):
        # two-decimal printing can put 0.825 at 0.82 in one file and 0.83 in the other
        return (abs(a["Inst"][i] - b["Inst"][j]) < tol and a["Obl"][i] == b["Obl"][j]
                and abs(a["XCO2"][i] - b["XCO2"][j]) <= 0.06 * max(a["XCO2"][i], 1.0))
    if na == nb and all(same(i, i, 0.011) for i in range(na)):
        return [(i, i) for i in range(na)]
    out = []
    for i in range(na):
        js = [j for j in range(nb) if same(i, j, 0.006)]
        if js:
            out.append((i, min(js, key=lambda j: abs(a["Inst"][i] - b["Inst"][j]))))
    return out


def compare(a, b, exp):
    pp = pairs(a, b)
    if not pp:
        return None
    dT = np.array([b["Tglob"][j] - a["Tglob"][i] for i, j in pp])
    changed = sum(1 for i, j in pp if state(a, i) != state(b, j))
    return dict(n=len(pp), na=len(a["Case"]), nb=len(b["Case"]), mean=dT.mean(),
                rms=np.sqrt((dT ** 2).mean()), imax=int(np.argmax(np.abs(dT))), dTmax=dT[np.argmax(np.abs(dT))],
                changed=changed, cens_a=census(a), cens_b=census(b), asym_a=asym(a), asym_b=asym(b))


def branch_order(gw, gc):
    """Warm start never colder than cold start at the same forcing."""
    viol, worst, n = 0, 0.0, 0
    for i, j in pairs(gc, gw):
        n += 1
        d = gc["Tglob"][i] - gw["Tglob"][j]
        if d > 1e-6:
            viol += 1
            worst = max(worst, d)
    return n, viol, worst


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("a")
    ap.add_argument("b")
    args = ap.parse_args()

    print("%-10s %5s %8s %8s %8s %8s %-22s %-22s %7s %7s" % (
        "exp", "n", "meandT", "rmsdT", "maxdT", "changed", "census A (free/cap/belt/snow)",
        "census B", "asymA", "asymB"))
    for exp in EXPS:
        a, b = load(args.a, exp), load(args.b, exp)
        if a is None or b is None:
            print("%-10s %s" % (exp, "missing in " + ("A" if a is None else "B")))
            continue
        c = compare(a, b, exp)
        if c is None:
            print("%-10s no matching cases (A %d, B %d)" % (exp, len(a["Case"]), len(b["Case"])))
            continue
        print("%-10s %5d %+8.2f %8.2f %+8.2f %8d %-22s %-22s %7.1f %7.1f" % (
            exp, c["n"], c["mean"], c["rms"], c["dTmax"], c["changed"],
            "/".join(map(str, c["cens_a"])), "/".join(map(str, c["cens_b"])), c["asym_a"], c["asym_b"]))

    for base, xkey, unit in (("exp3", "Inst", "S_earth"), ("exp4", "XCO2", "ppm")):
        for label, st in (("A", args.a), ("B", args.b)):
            gw, gc = load(st, base + "_warm"), load(st, base + "_cold")
            if gw is None or gc is None:
                continue
            gl, dg = thresholds(gw, gc, xkey)
            n, viol, worst = branch_order(gw, gc)
            fmt = lambda v: "none" if v is None else ("%.4g" % v)
            width = ("%.4f" % (dg - gl)) if (gl is not None and dg is not None and xkey == "Inst") else "-"
            print("%s %s: glaciation %s, deglaciation %s %s, width %s; branch order: %d of %d matched cases "
                  "have warm start colder than cold start (worst %.3f K)"
                  % (base, label, fmt(gl), fmt(dg), unit, width, viol, n, worst))
    for label, st in (("A", args.a), ("B", args.b)):
        g1, g2 = load(st, "exp1"), load(st, "exp2")
        if g1 is None or g2 is None:
            continue
        n, viol, worst = branch_order(g1, g2)
        print("exp1/exp2 %s: %d matched cases, %d with warm start colder than cold start (worst %.3f K)"
              % (label, n, viol, worst))


if __name__ == "__main__":
    main()
