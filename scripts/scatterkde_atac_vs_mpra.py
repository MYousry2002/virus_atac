import pandas as pd
import pyBigWig
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# === Load tile BED-like file ===
df = pd.read_csv("../MPRA/HSV1_KOS_MPRA_with_signals.bed", sep="\t", header=None)
df.columns = [
    'chrom', 'start', 'end', 'tile_id', 'score', 'strand',
    'GM12878', 'Jurkat', 'MRC5', 'A549', 'HEK293', 'K562',
    'GM12878_active', 'Jurkat_active', 'MRC5_active', 'A549_active', 'HEK293_active', 'K562_active'
]
df = df.sort_values(by='start').reset_index(drop=True)

# === Load ATAC signal from bigWig ===
bw = pyBigWig.open("../workdir/HSV1_KOS/SRR10176079/SRR10176079_HSV1_KOS.bw")
def get_signal(row):
    values = bw.values(row['chrom'], int(row['start']), int(row['end']), numpy=True)
    values = values[~np.isnan(values)]
    return np.mean(values) if len(values) > 0 else 0
df['ATAC_signal'] = df.apply(get_signal, axis=1)
bw.close()

# === Plot for each cell type ===
cell_types = ['GM12878', 'Jurkat', 'MRC5', 'A549', 'HEK293', 'K562']

for cell in cell_types:
    active_col = f"{cell}_active"
    df[active_col] = df[active_col].astype(bool)

    plt.figure(figsize=(7, 6))
    
    # Plot scatter with transparent dots and no white edges
    sns.scatterplot(
        data=df, x='ATAC_signal', y=cell,
        hue=active_col, palette={True: "orange", False: "steelblue"},
        alpha=0.3, edgecolor="none", s=20
    )
    
    # KDE Contours
    sns.kdeplot(
        data=df, x='ATAC_signal', y=cell,
        levels=5, linewidths=1, color='black', alpha=0.5
    )

    plt.title(f"ATAC Signal vs MPRA Activity – {cell}")
    plt.xlabel("ATAC Signal (mean per tile)")
    plt.ylabel("MPRA Activity")
    plt.legend(title="MPRA Active")
    plt.tight_layout()
    plt.savefig(f"../results/HSV1_KOS/atac_vs_mpra_{cell}_scatter_kde.png", dpi=300)
    plt.close()