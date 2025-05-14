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

# === Load ATAC signal from bigWig ===
bw = pyBigWig.open("../workdir/HSV1_KOS/SRR10176079/SRR10176079_HSV1_KOS.bw")
def get_signal(row):
    values = bw.values(row['chrom'], int(row['start']), int(row['end']), numpy=True)
    values = values[~np.isnan(values)]
    return np.mean(values) if len(values) > 0 else 0
df['ATAC_signal'] = df.apply(get_signal, axis=1)
bw.close()

# === Plot: MPRA distributions grouped by binarized ATAC signal ===
accessibility_threshold = 400  # adjust based on distribution
df['ATAC_binary'] = df['ATAC_signal'] > accessibility_threshold

cell_types = ['GM12878', 'Jurkat', 'MRC5', 'A549', 'HEK293', 'K562']
for cell in cell_types:
    plt.figure(figsize=(6, 5))
    
    sns.violinplot(data=df, x='ATAC_binary', y=cell, inner=None, palette='pastel')
    sns.boxplot(data=df, x='ATAC_binary', y=cell, width=0.2, showcaps=False, boxprops={'facecolor':'none'}, whiskerprops={'linewidth':1})
    
    plt.xticks([0, 1], ['Inaccessible', 'Accessible'])
    plt.title(f"MPRA Activity by ATAC Accessibility – {cell}")
    plt.xlabel("ATAC Accessibility (Binarized)")
    plt.ylabel("MPRA Activity")
    plt.tight_layout()
    plt.savefig(f"../results/HSV1_KOS/mpra_by_accessibility_{cell}.png", dpi=300)
    plt.close()



#### opposite ###
# === Load tile BED-like file ===
df = pd.read_csv("../MPRA/HSV1_KOS_MPRA_with_signals.bed", sep="\t", header=None)
df.columns = [
    'chrom', 'start', 'end', 'tile_id', 'score', 'strand',
    'GM12878', 'Jurkat', 'MRC5', 'A549', 'HEK293', 'K562',
    'GM12878_active', 'Jurkat_active', 'MRC5_active', 'A549_active', 'HEK293_active', 'K562_active'
]

# === Load ATAC signal from bigWig ===
bw = pyBigWig.open("../workdir/HSV1_KOS/SRR10176079/SRR10176079_HSV1_KOS.bw")
def get_signal(row):
    values = bw.values(row['chrom'], int(row['start']), int(row['end']), numpy=True)
    values = values[~np.isnan(values)]
    return np.mean(values) if len(values) > 0 else 0
df['ATAC_signal'] = df.apply(get_signal, axis=1)
bw.close()

# === Plot: ATAC signal grouped by MPRA binary status ===
cell_types = ['GM12878', 'Jurkat', 'MRC5', 'A549', 'HEK293', 'K562']
for cell in cell_types:
    active_col = f"{cell}_active"
    df[active_col] = df[active_col].astype(bool)

    plt.figure(figsize=(6, 5))
    sns.violinplot(data=df, x=active_col, y='ATAC_signal', inner=None, palette='pastel')
    sns.boxplot(data=df, x=active_col, y='ATAC_signal', width=0.2, showcaps=False, boxprops={'facecolor':'none'}, whiskerprops={'linewidth':1})

    plt.xticks([0, 1], ['Inactive', 'Active'])
    plt.title(f"ATAC Signal by MPRA Activity – {cell}")
    plt.xlabel("MPRA Activity (Binarized)")
    plt.ylabel("ATAC Signal (mean per tile)")
    plt.tight_layout()
    plt.savefig(f"../results/HSV1_KOS/atac_by_mpra_{cell}.png", dpi=300)
    plt.close()