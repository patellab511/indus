#!/usr/bin/env python3
"""Plot Ntilde(t) against the moving restraint target for the HFBII dewetting run.

usage: analyze.py [plumed.out] [out.png]
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
matplotlib.rc_file_defaults()
import matplotlib.pyplot as plt

src = sys.argv[1] if len(sys.argv) > 1 else "plumed.out"
out = sys.argv[2] if len(sys.argv) > 2 else "ntilde_vs_time.png"

d = np.loadtxt(src, comments="#")
t, n, nt, bias, work = d[:, 0], d[:, 1], d[:, 2], d[:, 3], d[:, 5]
dt_ps, ramp_ps, hold_ps, at0 = 0.002, 2000.0, 1000.0, 542.0
target = np.where(t <= ramp_ps, at0 * (1 - t / ramp_ps), 0.0)

blue, orange, ink, muted, grid = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e", "#e6e5e1"
plt.rcParams.update({"font.family": "sans-serif", "text.usetex": False, "font.size": 10, "axes.titlesize": 10, "axes.labelsize": 10, "xtick.labelsize": 9, "ytick.labelsize": 9, "legend.fontsize": 9, "axes.linewidth": 0.8, "xtick.major.width": 0.8, "ytick.major.width": 0.8, "axes.titleweight": "semibold",
                     "axes.edgecolor": muted, "xtick.color": muted, "ytick.color": muted})
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 5.6), sharex=True, facecolor="#fcfcfb",
                               gridspec_kw={"height_ratios": [3, 1.3], "hspace": 0.08})
ax1.plot(t / 1000, target, "--", color=orange, lw=1.5, label="restraint center N*")
ax1.plot(t / 1000, nt, "-", color=blue, lw=1.2, label=r"$\tilde N_v$ (INDUS)")
ax1.set_ylabel("waters in the 0.6 nm shell")
ax1.set_title("HFBII: driving water out of the hydration shell of 245 exposed heavy atoms", loc="left")
ax1.legend(frameon=False, loc="upper right")
ax2.plot(t / 1000, nt - target, "-", color=muted, lw=1)
ax2.axhline(0, color=grid, lw=1)
ax2.set_ylabel(r"$\tilde N_v - N^*$"); ax2.set_xlabel("t (ns)")
for ax in (ax1, ax2):
    ax.grid(True, color=grid, lw=0.8); ax.set_axisbelow(True); ax.spines[["top", "right"]].set_visible(False)
fig.savefig(out, dpi=160, bbox_inches="tight", facecolor=fig.get_facecolor())

hold = t > ramp_ps
print(f"frames: {len(t)}, t_end = {t[-1]:.0f} ps")
print(f"start:          N = {n[0]:.0f}, Ntilde = {nt[0]:.1f}")
if hold.any():
    print(f"hold (N* = 0):  Ntilde mean {nt[hold].mean():.1f}, min {nt[hold].min():.1f}, final {nt[-1]:.1f};  N final {n[-1]:.0f}")
    print(f"lag during ramp: mean (Ntilde - N*) = {(nt - target)[~hold].mean():.1f} waters")
print(f"accumulated restraint work: {work[-1]:.0f} kJ/mol = {work[-1] / 2.494:.0f} kT at 300 K")
print("wrote", out)
