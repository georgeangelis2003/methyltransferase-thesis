import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from collections import Counter

UNIPROT_CSV = '/home/angelis/thesis/reference_proteomes/active_sites.csv'
MAFFT_CSV   = '/home/angelis/thesis/reference_proteomes/active_sites_mafft.csv'
MSTA_CSV    = '/home/angelis/thesis/reference_proteomes/foldseek/fixed_shifts_msta.csv'
SAVE_PATH   = '/home/angelis/thesis/reference_proteomes/active_site_bar_chart.png'

uniprot = pd.read_csv(UNIPROT_CSV)
mafft   = pd.read_csv(MAFFT_CSV)
msta    = pd.read_csv(MSTA_CSV, sep='\t')

uniprot.columns = uniprot.columns.str.strip()
mafft.columns   = mafft.columns.str.strip()
msta.columns    = msta.columns.str.strip()

uniprot_aa = uniprot['Amino_Acid'].astype(str).str.strip()
mafft_aa   = mafft['Residue_at_Site'].astype(str).str.strip()
msta_aa    = msta['Amino_Acid'].astype(str).str.strip()

total_uniprot = len(uniprot_aa)
total_mafft   = len(mafft_aa)
total_msta    = len(msta_aa)

uniprot_counts = Counter(aa for aa in uniprot_aa if aa.lower() not in ('noinfo', 'nan', ''))
mafft_counts   = Counter(aa for aa in mafft_aa   if aa.lower() not in ('nan', ''))
msta_counts    = Counter(aa for aa in msta_aa    if aa.lower() not in ('nan', ''))

all_aas = sorted(
    set(uniprot_counts) | set(mafft_counts) | set(msta_counts),
    key=lambda aa: -(uniprot_counts.get(aa, 0) + mafft_counts.get(aa, 0) + msta_counts.get(aa, 0))
)
if '-' in all_aas:
    all_aas.remove('-')
    all_aas.append('-')

def pct(counts, aa, total):
    return counts.get(aa, 0) / total * 100

uniprot_pct = [pct(uniprot_counts, aa, total_uniprot) for aa in all_aas]
mafft_pct   = [pct(mafft_counts,   aa, total_mafft)   for aa in all_aas]
msta_pct    = [pct(msta_counts,    aa, total_msta)    for aa in all_aas]

COLORS = {'uniprot': '#A8C5DA', 'mafft': '#0a4b7a', 'msta': '#C0392B'}
TEXT_COLOR = '#1C2833'
GRID_COLOR = '#E8ECF0'
EDGE_COLOR = '#AABBCC'

plt.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['Georgia', 'Times New Roman', 'DejaVu Serif'],
    'axes.spines.top': False,
    'axes.spines.right': False,
    'axes.linewidth': 0.8,
})

fig, ax = plt.subplots(figsize=(14, 7), dpi=300)
fig.patch.set_facecolor('white')
ax.set_facecolor('white')

n = len(all_aas)
x = np.arange(n)
w = 0.30
gap = 0.01

bar1 = ax.bar(x - w - gap, uniprot_pct, width=w, label='UniProt annotation', color=COLORS['uniprot'], edgecolor='white', linewidth=0.4, zorder=2)
bar2 = ax.bar(x, mafft_pct, width=w, label='MAFFT (sequence MSA)', color=COLORS['mafft'], edgecolor='white', linewidth=0.4, zorder=2)
bar3 = ax.bar(x + w + gap, msta_pct, width=w, label='FoldMason ±3 (structural MSA)', color=COLORS['msta'], edgecolor='white', linewidth=0.4, zorder=2)

target_aas = {'H', 'D', 'E', 'C'}

def add_labels(rects):
    for i, rect in enumerate(rects):
        if all_aas[i] in target_aas:
            height = rect.get_height()
            if height > 0:
                ax.text(
                    rect.get_x() + rect.get_width()/2., height * 1.05,
                    f'{height:.1f}%',
                    ha='center', va='bottom', rotation=90,
                    fontsize=9, color=TEXT_COLOR, fontweight='bold'
                )

add_labels(bar1)
add_labels(bar2)
add_labels(bar3)

ax.set_yscale('log')
ax.set_ylim(bottom=0.01, top=200)
ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda val, _: f'{val:g}%'))

ax.set_xticks(x)
ax.set_xticklabels([aa if aa != '-' else 'gap (–)' for aa in all_aas], fontsize=11, color=TEXT_COLOR)
ax.set_ylabel('Percentage of proteins (%, log scale)', fontsize=12, color=TEXT_COLOR)
ax.set_xlabel('Amino acid at active site column', fontsize=12, color=TEXT_COLOR)
ax.set_title('Active Site Residue Distribution Across Alignment Methods', fontsize=14, color=TEXT_COLOR, fontweight='bold', pad=15)

ax.yaxis.grid(True, which='both', color=GRID_COLOR, linewidth=0.6, zorder=0)
ax.tick_params(colors=TEXT_COLOR, labelsize=10)

for spine in ax.spines.values():
    spine.set_color(EDGE_COLOR)

n_noinfo = (uniprot_aa.str.lower() == 'noinfo').sum()
noinfo_pct = n_noinfo / total_uniprot * 100
ax.annotate(f'UniProt: {noinfo_pct:.1f}% of proteins\nhave no active site annotation', xy=(0.98, 0.97), xycoords='axes fraction', ha='right', va='top', fontsize=9, color='#666666', style='italic', bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor=EDGE_COLOR, linewidth=0.6))

ax.legend(fontsize=10, frameon=True, framealpha=0.9, edgecolor='#CCCCCC', loc='upper right', bbox_to_anchor=(0.98, 0.85))

plt.tight_layout(pad=2.5)
plt.savefig(SAVE_PATH, dpi=300, bbox_inches='tight', facecolor='white')