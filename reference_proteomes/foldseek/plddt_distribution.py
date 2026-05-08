import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import seaborn as sns
from scipy.stats import gaussian_kde

pdb_dir   = '/home/angelis/thesis/reference_proteomes/foldseek/structures'
save_path = '/home/angelis/thesis/reference_proteomes/foldseek/pLDDT_distribution_thesis.png'

all_plddts = []
files = [f for f in os.listdir(pdb_dir) if f.endswith('.pdb')]
total = len(files)
print(f"Found {total} PDB files. Starting extraction")

for i, pdb_file in enumerate(files):
    if i % 500 == 0:
        print(f"  {i}/{total}")
    plddts = []
    with open(os.path.join(pdb_dir, pdb_file), 'r') as fh:
        for line in fh:
            if line.startswith('ATOM'):
                try:
                    plddts.append(float(line[60:66].strip()))
                except ValueError:
                    continue
    if plddts:
        all_plddts.append(sum(plddts) / len(plddts))

passed = sum(1 for p in all_plddts if p >= 80)
print(f"Total: {len(all_plddts)}  |  Passed (≥80): {passed}")

BAR_COLOR   = "#5B8DB8"
BAR_EDGE    = "#FFFFFF"
KDE_COLOR   = "#1A3A5C"
THRESH_COLOR= "#C0392B"
GRID_COLOR  = "#E8ECF0"
TEXT_COLOR  = "#1C2833"

plt.rcParams.update({
    "font.family"       : "serif",
    "font.serif"        : ["Georgia", "Times New Roman", "DejaVu Serif"],
    "axes.spines.top"   : False,
    "axes.spines.right" : False,
    "axes.spines.left"  : True,
    "axes.spines.bottom": True,
    "axes.linewidth"    : 0.8,
    "xtick.direction"   : "out",
    "ytick.direction"   : "out",
    "xtick.major.size"  : 4,
    "ytick.major.size"  : 4,
    "xtick.minor.size"  : 2,
    "ytick.minor.size"  : 2,
    "xtick.minor.visible": True,
    "ytick.minor.visible": True,
})

fig, ax = plt.subplots(figsize=(8.5, 5.5), dpi=300)
fig.patch.set_facecolor("white")
ax.set_facecolor("white")

n_bins = 120
counts, bin_edges = np.histogram(all_plddts, bins=n_bins)
ax.bar(
    bin_edges[:-1], counts,
    width=np.diff(bin_edges),
    align='edge',
    color=BAR_COLOR,
    edgecolor=BAR_EDGE,
    linewidth=0.4,
    alpha=0.75,
    zorder=2,
)

x_range = np.linspace(min(all_plddts) - 1, max(all_plddts) + 1, 500)
kde      = gaussian_kde(all_plddts, bw_method=0.15)
kde_vals = kde(x_range)

# scale KDE to histogram counts
bin_width = bin_edges[1] - bin_edges[0]
scale     = len(all_plddts) * bin_width
ax.plot(x_range, kde_vals * scale,
        color=KDE_COLOR, linewidth=0.9, zorder=3)

ax.axvline(x=80, color=THRESH_COLOR, linestyle='--',
           linewidth=1.2, zorder=4,
           label='Quality threshold (pLDDT ≥ 80)')

#grid
ax.yaxis.grid(True, color=GRID_COLOR, linewidth=0.6, zorder=0)
ax.set_axisbelow(True)

ax.set_xlabel('Average pLDDT Score', fontsize=11, color=TEXT_COLOR, labelpad=6)
ax.set_ylabel('Number of Structures', fontsize=11, color=TEXT_COLOR, labelpad=6)
ax.set_title(
    'Distribution of AlphaFold2 Confidence Scores\nin Plant Methyltransferases',
    fontsize=13, color=TEXT_COLOR, fontweight='bold', pad=12, linespacing=1.4
)

ax.tick_params(colors=TEXT_COLOR, labelsize=9)
for spine in ax.spines.values():
    spine.set_color("#AABBCC")

ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(5))
ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(5))

leg = ax.legend(
    fontsize=9, frameon=True, framealpha=0.9,
    edgecolor="#CCCCCC", fancybox=False,
    loc='upper left',
)
leg.get_frame().set_linewidth(0.6)

plt.tight_layout(pad=1.5)
plt.savefig(save_path, dpi=300, bbox_inches='tight', facecolor='white')
print(f"Saved {save_path}")