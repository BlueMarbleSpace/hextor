#!/usr/bin/env python3
"""
Compare OLR and PALB output between two radiation module versions.
Prints statistics and saves a multi-panel comparison figure.

Usage:
    python compare_radiation.py out_old.txt out_new.txt [figure.png]
"""

import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

# ── Load ──────────────────────────────────────────────────────────────────────

def load(path):
    return np.loadtxt(path, comments='#', usecols=(0, 1, 2, 3, 4, 5))

old = load(sys.argv[1])
new = load(sys.argv[2])
fig_path = sys.argv[3] if len(sys.argv) > 3 else "radiation_comparison.png"

assert old.shape == new.shape, "Output files have different row counts"

fco2    = old[:,0]
tg0     = old[:,1]
zy      = old[:,2]
surfalb = old[:,3]
olr_old,  palb_old  = old[:,4], old[:,5]
olr_new,  palb_new  = new[:,4], new[:,5]

dolr  = olr_new  - olr_old
dpalb = palb_new - palb_old

valid_olr  = (olr_old  != 1.0)  & (olr_new  != 1.0)  & np.isfinite(dolr)
valid_palb = (palb_old != -1.0) & (palb_new != -1.0) & np.isfinite(dpalb)

# ── Text statistics ───────────────────────────────────────────────────────────

def stats(label, diff, mask):
    d = diff[mask]
    print(f"  {label}:")
    print(f"    max |diff|  = {np.max(np.abs(d)):.4f}")
    print(f"    mean |diff| = {np.mean(np.abs(d)):.4f}")
    print(f"    RMS diff    = {np.sqrt(np.mean(d**2)):.4f}")
    n_exact = np.sum(np.abs(d) < 1e-10)
    print(f"    exact match = {n_exact}/{len(d)} ({100*n_exact/len(d):.1f}%)")

def top_rows(label, diff, n=10):
    print(f"\nTop 10 largest {label} differences:")
    idx = np.argsort(np.abs(diff))[::-1][:n]
    print(f"  {'fco2':>12}  {'T(K)':>7}  {'zen':>6}  {'salb':>6}  {'diff':>14}")
    for i in idx:
        print(f"  {fco2[i]:12.3e}  {tg0[i]:7.1f}  {zy[i]:6.1f}  {surfalb[i]:6.2f}  {diff[i]:14.4f}")

print("=" * 60)
print("Overall statistics")
print("=" * 60)
stats("OLR  (raw units)", dolr,  valid_olr)
stats("PALB (0–1)",       dpalb, valid_palb)

top_rows("OLR",  dolr)
top_rows("PALB", dpalb)

print("\nMean |PALB diff| by zenith angle:")
for z in np.unique(zy):
    mask = valid_palb & (zy == z)
    print(f"  zy = {z:4.0f}°:  {np.mean(np.abs(dpalb[mask])):.4f}")

print("\nMean |OLR diff| by temperature:")
for t in np.unique(tg0):
    mask = valid_olr & (tg0 == t)
    print(f"  T = {t:5.0f} K:  {np.mean(np.abs(dolr[mask])):.4f}")

n_miss_old = np.sum(~valid_olr)
n_miss_new = np.sum(olr_new == 1.0) + np.sum(palb_new == -1.0)
if n_miss_old or n_miss_new:
    print(f"\nSentinel / hash-miss counts:  old = {n_miss_old}   new = {n_miss_new}")

# ── Reshape into 4D ───────────────────────────────────────────────────────────
# Loop order in compare.f90: surfalb (outer) → zenith → T → CO2 (inner)

u_co2 = np.sort(np.unique(fco2))     # 28
u_tmp = np.sort(np.unique(tg0))      # 8
u_zen = np.sort(np.unique(zy))       # 7
u_alb = np.sort(np.unique(surfalb))  # 7
nco2, ntmp, nzen, nalb = len(u_co2), len(u_tmp), len(u_zen), len(u_alb)

# Shape: (nalb, nzen, ntmp, nco2)
adolr  = np.where(valid_olr,  np.abs(dolr),  np.nan).reshape(nalb, nzen, ntmp, nco2)
adpalb = np.where(valid_palb, np.abs(dpalb), np.nan).reshape(nalb, nzen, ntmp, nco2)

# Reduce axes for each panel
olr_T_co2    = np.nanmean(adolr,  axis=(0, 1))  # (ntmp, nco2)
palb_T_co2   = np.nanmean(adpalb, axis=(0, 1))  # (ntmp, nco2)
palb_alb_zen = np.nanmean(adpalb, axis=(2, 3))  # (nalb, nzen)

# ── Style ─────────────────────────────────────────────────────────────────────
plt.rcParams.update({
    'font.family'      : 'sans-serif',
    'font.size'        : 10,
    'axes.labelsize'   : 11,
    'axes.titlesize'   : 11,
    'axes.titleweight' : 'bold',
    'axes.spines.top'  : False,
    'axes.spines.right': False,
    'axes.grid'        : False,
    'xtick.labelsize'  : 9,
    'ytick.labelsize'  : 9,
    'figure.facecolor' : 'white',
})

CMAP     = 'viridis'
log10co2 = np.log10(u_co2)

# ── Helpers ───────────────────────────────────────────────────────────────────

def add_colorbar(fig, ax, im, label=''):
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.06)
    cb = fig.colorbar(im, cax=cax)
    cb.ax.tick_params(labelsize=8)
    if label:
        cb.set_label(label, fontsize=9)
    return cb

def co2_xticks(ax):
    ticks = np.arange(np.floor(log10co2.min()), np.ceil(log10co2.max()) + 1)
    ax.set_xticks(ticks)
    ax.set_xticklabels([f'$10^{{{int(v)}}}$' for v in ticks],
                       rotation=35, ha='right', rotation_mode='anchor')
    ax.set_xlabel('CO$_2$ fraction')

def panel_label(ax, letter):
    ax.text(-0.13, 1.04, f'({letter})', transform=ax.transAxes,
            fontsize=11, fontweight='bold', va='bottom', ha='left')

# ── Figure ────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))

# Integer coordinates give uniform cell heights; relabel with actual values
tmp_idx = np.arange(ntmp)
zen_idx = np.arange(nzen)
alb_idx = np.arange(nalb)

# ── (a)  |ΔOLR| over CO₂ × T ─────────────────────────────────────────────────
ax = axes[0]
im = ax.pcolormesh(log10co2, tmp_idx, olr_T_co2, cmap=CMAP, shading='auto')
co2_xticks(ax)
ax.set_yticks(tmp_idx)
ax.set_yticklabels([f'{int(v)}' for v in u_tmp])
ax.set_ylabel('Temperature (K)')
ax.set_title('Mean |ΔOLR| (W m$^{-2}$)')
add_colorbar(fig, ax, im)
panel_label(ax, 'a')

# ── (b)  |ΔPALB| over CO₂ × T ────────────────────────────────────────────────
ax = axes[1]
im = ax.pcolormesh(log10co2, tmp_idx, palb_T_co2, cmap=CMAP, shading='auto')
co2_xticks(ax)
ax.set_yticks(tmp_idx)
ax.set_yticklabels([f'{int(v)}' for v in u_tmp])
ax.set_ylabel('Temperature (K)')
ax.set_title('Mean |ΔPALB|')
add_colorbar(fig, ax, im)
panel_label(ax, 'b')

# ── (c)  |ΔPALB| over zenith × surface albedo ────────────────────────────────
ax = axes[2]
im = ax.pcolormesh(zen_idx, alb_idx, palb_alb_zen, cmap=CMAP, shading='auto')
ax.set_xticks(zen_idx)
ax.set_xticklabels([f'{int(v)}°' for v in u_zen])
ax.set_yticks(alb_idx)
ax.set_yticklabels([f'{v:.2f}' for v in u_alb])
ax.set_xlabel('Zenith angle')
ax.set_ylabel('Surface albedo')
ax.set_title('Mean |ΔPALB|')
add_colorbar(fig, ax, im)
panel_label(ax, 'c')

# ── Save ──────────────────────────────────────────────────────────────────────
fig.tight_layout(rect=[0, 0, 1, 0.93])
eps_path = fig_path.rsplit('.', 1)[0] + '.eps'
fig.savefig(fig_path, dpi=150, bbox_inches='tight')
fig.savefig(eps_path, bbox_inches='tight')
print(f"\nFigure saved: {fig_path}")
print(f"Figure saved: {eps_path}")
