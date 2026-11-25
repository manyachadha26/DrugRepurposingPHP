"""
Figure 3 - PDE isoform expression across the nephron.
Data: Human Protein Atlas single cell type RNA (nCPM), rna_single_cell_type.tsv.
Input: pde_proximal_tubule.csv, written by extract_pde.py.

Run:  conda activate phppde4 && python figure3.py
Outputs: figure3_pde_expression.pdf (vector, for submission)
         figure3_pde_expression.png (300 dpi, for drafts)
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LogNorm
from matplotlib.ticker import FixedLocator

# ---- EDIT THIS before submitting: the HPA release your numbers came from ----
HPA_VERSION = "HPA version 25.1"
ACCESSED = "accessed 7 October 2026"

PDE4_C = "#2a78d6"   # categorical slot 1
COMP_C = "#eb6834"   # categorical slot 2
INK = "#0b0b0b"
MUTED = "#52514e"
GRID = "#d8d7d2"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
    "font.size": 8,
    "axes.edgecolor": MUTED,
    "axes.linewidth": 0.6,
    "text.color": INK,
    "axes.labelcolor": INK,
    "xtick.color": MUTED,
    "ytick.color": MUTED,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
})

df = pd.read_csv("pde_proximal_tubule.csv").set_index("Gene name")

# anatomical order: glomerulus first, then down the tubule
SEGMENTS = [
    ("podocytes", "Podocytes\n(glomerulus)"),
    ("proximal tubule cells", "Proximal\ntubule"),
    ("loop of henle epithelial cells", "Loop of\nHenle"),
    ("distal convoluted tubule cells", "Distal\nconvoluted"),
    ("renal connecting tubule cells", "Connecting\ntubule"),
    ("renal collecting duct principal cells", "Collecting duct\nprincipal"),
    ("renal collecting duct intercalated cells", "Collecting duct\nintercalated"),
]
KEY = "proximal tubule cells"

fig = plt.figure(figsize=(7.2, 7.4))
gs = fig.add_gridspec(2, 1, height_ratios=[1.05, 1.0], hspace=0.42,
                      left=0.145, right=0.965, top=0.915, bottom=0.115)

# ---------------- Panel A: proximal tubule, the decisive comparison ----------
axA = fig.add_subplot(gs[0])
a = df.sort_values(KEY, ascending=True)
colors = [PDE4_C if g == "PDE4" else COMP_C for g in a["Group"]]
y = np.arange(len(a))
axA.barh(y, a[KEY], color=colors, height=0.68, zorder=3)
axA.set_yticks(y)
axA.set_yticklabels(a.index, fontsize=8)
for t, g in zip(axA.get_yticklabels(), a["Group"]):
    if g == "PDE4":
        t.set_fontweight("bold")
        t.set_color(INK)

for i, v in enumerate(a[KEY]):
    axA.text(v + max(a[KEY]) * 0.012, i, f"{v:,.1f}" if v > 0 else "0.0",
             va="center", ha="left", fontsize=7.2, color=MUTED, zorder=4)

axA.set_xlim(0, max(a[KEY]) * 1.17)
axA.set_xlabel("Expression in proximal tubule cells (nCPM)", fontsize=8.5)
axA.xaxis.grid(True, color=GRID, linewidth=0.6, zorder=0)
axA.set_axisbelow(True)
for s in ("top", "right", "left"):
    axA.spines[s].set_visible(False)
axA.tick_params(axis="y", length=0)

handles = [plt.Rectangle((0, 0), 1, 1, color=PDE4_C),
           plt.Rectangle((0, 0), 1, 1, color=COMP_C)]
axA.legend(handles, ["PDE4 family", "Other cAMP-hydrolysing PDEs"],
           loc="lower right", frameon=False, fontsize=7.6, handlelength=1.1,
           handleheight=1.1, borderpad=0.2)
axA.set_title("A   PDE4D dominates the proximal tubule, the site of PTH action",
              loc="left", fontsize=9.5, fontweight="bold", pad=8)

ratio = a[KEY].max() / a[a["Group"] != "PDE4"][KEY].max()
axA.annotate(f"{ratio:.1f}x the nearest competitor",
             xy=(a[KEY].max() * 0.80, len(a) - 1.32),
             xytext=(a[KEY].max() * 0.46, len(a) - 2.9),
             fontsize=7.6, color=INK,
             arrowprops=dict(arrowstyle="->", color=MUTED, linewidth=0.7))

# ---------------- Panel B: the whole nephron --------------------------------
axB = fig.add_subplot(gs[1])
cols = [c for c, _ in SEGMENTS]
m = df.sort_values(KEY, ascending=False)[cols]
vals = m.to_numpy(dtype=float)
plot = np.where(vals <= 0, 0.05, vals)      # zeros sit at the floor of the log scale

im = axB.imshow(plot, aspect="auto", cmap="Blues",
                norm=LogNorm(vmin=0.1, vmax=np.nanmax(plot)))
axB.set_xticks(np.arange(len(cols)))
axB.set_xticklabels([lab for _, lab in SEGMENTS], fontsize=7.2)
axB.set_yticks(np.arange(len(m)))
axB.set_yticklabels(m.index, fontsize=7.6)
for t, g in zip(axB.get_yticklabels(), df.loc[m.index, "Group"]):
    if g == "PDE4":
        t.set_fontweight("bold")
axB.tick_params(length=0)
for s in axB.spines.values():
    s.set_visible(False)

# value in every cell - the heatmap is read, not just looked at
for i in range(vals.shape[0]):
    for j in range(vals.shape[1]):
        v = vals[i, j]
        shade = im.cmap(im.norm(plot[i, j]))
        lum = 0.299 * shade[0] + 0.587 * shade[1] + 0.114 * shade[2]
        axB.text(j, i, f"{v:,.0f}" if v >= 1 else ("0" if v == 0 else f"{v:.1f}"),
                 ha="center", va="center", fontsize=6.1,
                 color="white" if lum < 0.55 else MUTED)

axB.set_xticks(np.arange(-0.5, len(cols), 1), minor=True)
axB.set_yticks(np.arange(-0.5, len(m), 1), minor=True)
axB.grid(which="minor", color="white", linewidth=1.4)
axB.tick_params(which="minor", length=0)

cb = fig.colorbar(im, ax=axB, pad=0.015, fraction=0.028)
cb.set_label("nCPM (log scale)", fontsize=7.6)
cb.ax.tick_params(labelsize=6.8, length=2)
cb.outline.set_visible(False)
cb.ax.yaxis.set_major_locator(FixedLocator([0.1, 1, 10, 100, 1000, 10000]))
cb.ax.set_yticklabels(["0.1", "1", "10", "100", "1,000", "10,000"])

axB.set_title("B   PDE4D is expressed throughout the nephron, not only the proximal tubule",
              loc="left", fontsize=9.5, fontweight="bold", pad=8)

fig.text(0.145, 0.038,
         f"Human Protein Atlas single cell type RNA data ({HPA_VERSION}, {ACCESSED}). "
         "Values are nCPM per cell type.\nCells showing 0 were reported as zero, "
         "not as missing data.",
         fontsize=6.8, color=MUTED, va="top")

fig.savefig("figure3_pde_expression.pdf", bbox_inches="tight")
fig.savefig("figure3_pde_expression.png", dpi=300, bbox_inches="tight")
print("wrote figure3_pde_expression.pdf and .png")
print(f"PDE4D / nearest competitor in proximal tubule = {ratio:.2f}x")
