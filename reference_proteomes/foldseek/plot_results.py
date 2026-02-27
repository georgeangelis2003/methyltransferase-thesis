import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import gaussian_kde
import numpy as np

cols = ["query", "target", "fident", "alnlen", "mismatch", "gapopen", 
        "qstart", "qend", "tstart", "tend", "evalue", "bits", "lddt", "alntmscore"]

print('loading results')
df = pd.read_csv('/home/angelis/thesis/reference_proteomes/foldseek/results.tsv', sep='\t', names=cols)

df = df[df['query'] != df['target']]
df['fident_pct'] = df['fident'] * 100

print('calculating density...')
sample = df.sample(n=min(len(df), 50000))
x = sample['fident_pct']
y = sample['alntmscore']
xy = np.vstack([x, y])
z = gaussian_kde(xy)(xy)

plt.figure(figsize=(12, 8))
plt.scatter(x, y, c=z, s=5, cmap='viridis', alpha=0.5, edgecolor=None)

plt.axhline(y=0.5, color='red', linestyle='--', label='Structural Similarity > 0.5')
plt.axvline(x=30, color='orange', linestyle='--', label='Sequence ID < 30%')

plt.title('Foldseek Results: alntmscore vs Sequence Identity (Density Weighted)')
plt.xlabel('Sequence Identity (%)')
plt.ylabel('alntmscore')
plt.colorbar(label='Point Density')
plt.legend()
plt.grid(True, linestyle=':', alpha=0.6)

plt.savefig('/home/angelis/thesis/reference_proteomes/foldseek/foldseek_plot_density.png', dpi=300)
print('done, plot saved as foldseek_plot_density.png')