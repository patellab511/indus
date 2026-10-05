#!/usr/bin/env python3
"""Compose the four-state HFBII hydration figure with the hydrophobic patch outlined.

Inputs (same viewpoint, rendered by render_count4.pml and render_patchmask.pml):
  start_count.png mid1000_count.png mid1500_count.png final_count.png   surface colored by water count
  patchmask_view1.png   protein grey, patch residues red -> outline
usage: compose_hydration_figure.py <out.png>
"""
import sys
import numpy as np
import matplotlib; matplotlib.use("Agg"); matplotlib.rc_file_defaults()
import matplotlib.pyplot as plt, matplotlib.image as mpimg
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable

out = sys.argv[1] if len(sys.argv) > 1 else "hfbii_union_of_spheres_hydration.png"
plt.rcParams.update({"font.family": "sans-serif", "text.usetex": False, "font.size": 11, "xtick.labelsize": 10, "axes.titlesize": 14})
panels = [("start_count.png", 542), ("mid1000_count.png", 282), ("mid1500_count.png", 151), ("final_count.png", 38)]
crop = lambda im: im[int(0.08 * im.shape[0]):int(0.95 * im.shape[0]), int(0.06 * im.shape[1]):int(0.94 * im.shape[1])]

mask_im = crop(mpimg.imread("patchmask_view1.png")[:, :, :3])
patch = (mask_im[:, :, 0] > 0.6) & (mask_im[:, :, 1] < 0.4) & (mask_im[:, :, 2] < 0.4)
ys, xs = np.where(patch); cx, cy = xs.mean(), ys.mean()
orange, ink = "#eb6834", "#0b0b0b"

fig, axes = plt.subplots(1, 4, figsize=(14, 4.3), facecolor="white")
for i, (ax, (f, nt)) in enumerate(zip(axes, panels)):
    ax.imshow(crop(mpimg.imread(f))); ax.set_axis_off(); ax.set_title(rf"$\tilde N_v$ = {nt}", pad=4)
    ax.contour(patch.astype(float), levels=[0.5], colors=[orange], linewidths=1.3, linestyles="--")
    if i == 0:
        ax.annotate("hydrophobic patch", xy=(cx - 0.18 * patch.shape[1], cy - 0.12 * patch.shape[0]),
                    xytext=(0.03, 0.04), textcoords="axes fraction", color=orange, fontsize=11, fontweight="semibold",
                    arrowprops=dict(arrowstyle="-", color=orange, lw=1.2))
cmap = LinearSegmentedColormap.from_list("w2b", ["#ffffff", "#2a78d6"])
sm = ScalarMappable(norm=Normalize(0, 12), cmap=cmap); sm.set_array([])
cax = fig.add_axes([0.35, 0.12, 0.30, 0.037])
cb = fig.colorbar(sm, cax=cax, orientation="horizontal"); cb.outline.set_edgecolor("#aaaaaa")
cb.set_ticks([0, 3, 6, 9, 12]); cb.ax.tick_params(labelsize=10, length=3)
cb.set_label("water oxygens within 0.6 nm of each protein atom", fontsize=11, labelpad=5)
fig.suptitle("Hydrophobin HFBII: local hydration as INDUS $\\tilde N_v$ in a union of 245 spheres (r = 0.6 nm) is driven from 542 to 0", fontsize=12, y=0.99)
plt.subplots_adjust(left=0.005, right=0.995, top=0.85, bottom=0.21, wspace=0.01)
fig.savefig(out, dpi=150, facecolor="white"); print("wrote", out)
