import pandas as pd
import pyBigWig
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import gridspec

# === Load tile BED-like file ===
df = pd.read_csv("../MPRA/HSV1_KOS_MPRA_with_signals.bed", sep="\t", header=None)
df.columns = [
    'chrom', 'start', 'end', 'tile_id', 'score', 'strand',
    'GM12878', 'Jurkat', 'MRC5', 'A549', 'HEK293', 'K562',
    'GM12878_active', 'Jurkat_active', 'MRC5_active', 'A549_active', 'HEK293_active', 'K562_active'
]
df = df.sort_values(by='start').reset_index(drop=True)

# === Load ATAC signal from bigWig ===
bw = pyBigWig.open("../workdir/SRR10176079/SRR10176079_HSV1_KOS.bw")
def get_signal(row):
    values = bw.values(row['chrom'], int(row['start']), int(row['end']), numpy=True)
    values = values[~np.isnan(values)]
    return np.mean(values) if len(values) > 0 else 0
df['ATAC_signal'] = df.apply(get_signal, axis=1)
bw.close()

# === Y-ticks every 10kb ===
tick_labels_kb = list(range(0, 151 + 1, 10))
tick_positions = [(df['start'] - kb * 1000).abs().idxmin() for kb in tick_labels_kb]

# === Loop over all cell types ===
cell_types = ['GM12878', 'Jurkat', 'MRC5', 'A549', 'HEK293', 'K562']

for cell in cell_types:
    matrix = df[['ATAC_signal', cell]].copy()
    matrix.columns = ['ATAC', 'MPRA']

    # Plot
    fig = plt.figure(figsize=(4.5, 10))
    gs = gridspec.GridSpec(ncols=3, nrows=1, width_ratios=[1, 1, 0.1], wspace=0.05)

    # ATAC heatmap
    ax1 = fig.add_subplot(gs[0])
    heat_atac = sns.heatmap(matrix[['ATAC']], cmap='viridis', cbar=False, yticklabels=True, ax=ax1)
    ax1.set_title("ATAC")
    ax1.set_ylabel("Tiles")
    ax1.set_xticklabels(['ATAC Signal'])
    ax1.set_yticks(tick_positions)
    ax1.set_yticklabels(tick_labels_kb, fontsize=8)

    # MPRA heatmap
    ax2 = fig.add_subplot(gs[1], sharey=ax1)
    heat_mpra = sns.heatmap(matrix[['MPRA']], cmap='magma', cbar=False, yticklabels=False, ax=ax2)
    ax2.set_title("MPRA")
    ax2.set_ylabel("")
    ax2.set_xticklabels([f"MPRA Activity ({cell})"])

    # Colorbars
    cax1 = fig.add_axes([0.92, 0.55, 0.015, 0.35])
    cbar1 = plt.colorbar(heat_atac.collections[0], cax=cax1)
    cbar1.set_label("ATAC Signal", rotation=270, labelpad=15)

    cax2 = fig.add_axes([0.92, 0.1, 0.015, 0.35])
    cbar2 = plt.colorbar(heat_mpra.collections[0], cax=cax2)
    cbar2.set_label("MPRA Activity", rotation=270, labelpad=15)

    # Save
    plt.suptitle(f"ATAC Signal vs MPRA Activity ({cell}) per Tile", y=0.98, fontsize=14)
    plt.savefig(f"../results/atac_vs_mpra_{cell}_heatmap.png", dpi=300, bbox_inches='tight')
    plt.close()