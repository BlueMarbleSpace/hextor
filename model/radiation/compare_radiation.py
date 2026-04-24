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
olr_T_co2    = np.nanmean(adolr,  axis=(0, 1))   # (ntmp, nco2) — mean over zen & alb
palb_T_co2   = np.nanmean(adpalb, axis=(0, 1))   # (ntmp, nco2)
palb_alb_zen = np.nanmean(adpalb, axis=(2, 3))   # (nalb, nzen) — mean over T & CO2
palb_by_zen  = np.nanmean(adpalb, axis=(0, 2, 3))  # (nzen,)
palb_by_alb  = np.nanmean(adpalb, axis=(1, 2, 3))  # (nalb,)

# ── Figure ────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 3, figsize=(16, 9))
fig.suptitle("Radiation module: improved vs original interpolation",
             fontsize=13, fontweight='bold')

CMAP = 'plasma'
log10co2 = np.log10(u_co2)

def co2_xticks(ax):
    ticks = log10co2[::4]
    ax.set_xticks(ticks)
    ax.set_xticklabels([f'$10^{{{v:.1f}}}$' for v in ticks], fontsize=7)
    ax.set_xlabel('CO₂ fraction')

def tmp_yticks(ax):
    ax.set_yticks(u_tmp)
    ax.set_yticklabels([f'{int(v)}' for v in u_tmp], fontsize=8)
    ax.set_ylabel('Temperature (K)')

# ── (0,0)  |ΔOLR| over CO2 × T ───────────────────────────────────────────────
ax = axes[0, 0]
im = ax.pcolormesh(log10co2, u_tmp, olr_T_co2, cmap=CMAP, shading='auto')
co2_xticks(ax); tmp_yticks(ax)
ax.set_title('Mean |ΔOLR| (raw units)\naveraged over zenith & surface albedo')
fig.colorbar(im, ax=ax).ax.tick_params(labelsize=8)

# ── (0,1)  |ΔPALB| over CO2 × T ──────────────────────────────────────────────
ax = axes[0, 1]
im = ax.pcolormesh(log10co2, u_tmp, palb_T_co2, cmap=CMAP, shading='auto')
co2_xticks(ax); tmp_yticks(ax)
ax.set_title('Mean |ΔPALB|\naveraged over zenith & surface albedo')
fig.colorbar(im, ax=ax).ax.tick_params(labelsize=8)

# ── (0,2)  |ΔPALB| over zenith × surface albedo ──────────────────────────────
ax = axes[0, 2]
# palb_alb_zen shape: (nalb, nzen) → rows = surfalb, cols = zenith
im = ax.pcolormesh(u_zen, u_alb, palb_alb_zen, cmap=CMAP, shading='auto')
ax.set_xticks(u_zen)
ax.set_xticklabels([f'{int(v)}°' for v in u_zen], fontsize=9)
ax.set_yticks(u_alb)
ax.set_yticklabels([f'{v:.2f}' for v in u_alb], fontsize=9)
ax.set_xlabel('Zenith angle')
ax.set_ylabel('Surface albedo')
ax.set_title('Mean |ΔPALB|\naveraged over CO₂ & temperature')
fig.colorbar(im, ax=ax).ax.tick_params(labelsize=8)

# ── (1,0)  ΔOLR histogram ─────────────────────────────────────────────────────
ax = axes[1, 0]
d_olr = dolr[valid_olr]
ax.hist(d_olr, bins=60, color='steelblue', edgecolor='none', alpha=0.85)
ax.axvline(0, color='k', lw=0.8, ls='--')
ax.set_xlabel('ΔOLR (new − old, raw units)')
ax.set_ylabel('Count')
ax.set_title('Distribution of OLR differences')
ax.text(0.97, 0.95,
        f'max = {np.max(np.abs(d_olr)):.1f}\nRMS = {np.sqrt(np.mean(d_olr**2)):.1f}',
        transform=ax.transAxes, ha='right', va='top', fontsize=9,
        bbox=dict(boxstyle='round', fc='white', alpha=0.8))

# ── (1,1)  mean |ΔPALB| vs zenith, with per-T curves ─────────────────────────
ax = axes[1, 1]
for it in range(ntmp):
    curve = np.nanmean(adpalb[:, :, it, :], axis=(0, 2))  # (nzen,)
    ax.plot(u_zen, curve, '-', color='grey', lw=0.7, alpha=0.5)
ax.plot(u_zen, palb_by_zen, 'o-', color='darkorange', lw=2, ms=7,
        zorder=3, label='overall mean')
ax.set_xticks(u_zen)
ax.set_xticklabels([f'{int(v)}°' for v in u_zen])
ax.set_xlabel('Zenith angle (°)')
ax.set_ylabel('Mean |ΔPALB|')
ax.set_title('PALB difference vs zenith\n(bold = overall mean; grey = per-T)')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

# ── (1,2)  mean |ΔPALB| vs surface albedo, with per-T curves ─────────────────
ax = axes[1, 2]
for it in range(ntmp):
    curve = np.nanmean(adpalb[:, :, it, :], axis=(1, 2))  # (nalb,)
    ax.plot(u_alb, curve, '-', color='grey', lw=0.7, alpha=0.5)
ax.plot(u_alb, palb_by_alb, 's-', color='forestgreen', lw=2, ms=7,
        zorder=3, label='overall mean')
ax.set_xticks(u_alb)
ax.set_xlabel('Surface albedo')
ax.set_ylabel('Mean |ΔPALB|')
ax.set_title('PALB difference vs surface albedo\n(bold = overall mean; grey = per-T)')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(fig_path, dpi=150, bbox_inches='tight')
print(f"\nFigure saved: {fig_path}")
