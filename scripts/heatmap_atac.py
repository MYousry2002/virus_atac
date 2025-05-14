import pandas as pd
import pyBigWig
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# === Load BED-like tile file ===
df = pd.read_csv("../MPRA/HSV1_KOS_MPRA_with_signals.bed", sep="\t", header=None)

# Assign column names
df.columns = [
    'chrom', 'start', 'end', 'tile_id', 'score', 'strand',
    'GM12878', 'Jurkat', 'MRC5', 'A549', 'HEK293', 'K562',
    'GM12878_active', 'Jurkat_active', 'MRC5_active', 'A549_active', 'HEK293_active', 'K562_active'
]

# Sort by start coordinate
df = df.sort_values(by='start').reset_index(drop=True)

# === Load bigWig ===
bw = pyBigWig.open("../workdir/SRR10176079/SRR10176079_HSV1_KOS.bw")

# Compute signal for each tile
def get_signal(row):
    signal = bw.values(row['chrom'], int(row['start']), int(row['end']), numpy=True)
    signal = signal[~np.isnan(signal)]
    return np.mean(signal) if len(signal) > 0 else 0

df['ATAC_signal'] = df.apply(get_signal, axis=1)
bw.close()

# === Prepare matrix ===
heatmap_matrix = df[['ATAC_signal']]

# === Generate y-tick positions for every 10kb ===
max_kb = (df['start'].max() // 10000) * 10 + 10  # round up
tick_labels_kb = list(range(0, 151, 10))  # 0, 10, 20, ..., max_kb
tick_positions = []

# For each 10kb label, find the closest tile index
for kb in tick_labels_kb:
    pos = kb * 1000
    # find the index of the tile with closest start
    closest_idx = (df['start'] - pos).abs().idxmin()
    tick_positions.append(closest_idx)

# === Plot heatmap ===
plt.figure(figsize=(6, 12))
sns.heatmap(heatmap_matrix, cmap="viridis", yticklabels=False, cbar_kws={'label': 'ATAC Signal'})

plt.yticks(ticks=tick_positions, labels=tick_labels_kb, fontsize=8)
plt.ylabel("Genomic Position (kb)")
plt.xlabel("ATAC Signal")
plt.title("ATAC Signal per Tile (HSV1 KOS)")

plt.tight_layout()
plt.savefig("../results/atac_signal_heatmap_10kb_ticks.png", dpi=300)
plt.show()