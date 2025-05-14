import pandas as pd
import pyBigWig
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os

# === Config ===
mpra_file = "../MPRA/HSV1_KOS_MPRA_with_signals.bed"
bigwig_file = "../workdir/HSV1_KOS/SRR10176079/SRR10176079_HSV1_KOS.bw"
output_table = "../results/mpra_atac_category_table.csv"
output_plot = "../results/mpra_atac_category_heatmap.png"

cell_types = ['GM12878', 'Jurkat', 'MRC5', 'A549', 'HEK293', 'K562']
accessibility_threshold = 650.0  # you can tune this based on signal distribution

# === Load data ===
df = pd.read_csv(mpra_file, sep="\t", header=None)
df.columns = [
    'chrom', 'start', 'end', 'tile_id', 'score', 'strand',
    'GM12878', 'Jurkat', 'MRC5', 'A549', 'HEK293', 'K562',
    'GM12878_active', 'Jurkat_active', 'MRC5_active', 'A549_active', 'HEK293_active', 'K562_active'
]

# === Load bigWig signal ===
bw = pyBigWig.open(bigwig_file)

def average_signal(row):
    values = bw.values(row['chrom'], int(row['start']), int(row['end']), numpy=True)
    values = values[~np.isnan(values)]
    return np.mean(values) if len(values) > 0 else 0

df['ATAC_signal'] = df.apply(average_signal, axis=1)
bw.close()

# === Categorize each tile per cell type ===
category_counts = {}

for cell in cell_types:
    active_col = f"{cell}_active"
    df[active_col] = df[active_col].astype(bool)
    
    categories = []
    for i, row in df.iterrows():
        active = row[active_col]
        accessible = row['ATAC_signal'] > accessibility_threshold
        
        if active and accessible:
            categories.append("Active & Accessible")
        elif active and not accessible:
            categories.append("Active only")
        elif not active and accessible:
            categories.append("Accessible only")
        else:
            categories.append("Neither")
    
    df[f"{cell}_category"] = categories
    counts = pd.Series(categories).value_counts()
    for cat in ['Active & Accessible', 'Active only', 'Accessible only', 'Neither']:
        category_counts.setdefault(cat, {})[cell] = counts.get(cat, 0)

# === Save full table with per-cell-type category ===
df.to_csv(output_table, index=False)
print(f"Saved categorized MPRA+ATAC table: {output_table}")

# === Plot heatmap ===
heatmap_df = pd.DataFrame(category_counts).T[cell_types]
plt.figure(figsize=(10, 6))
sns.heatmap(heatmap_df, annot=True, fmt='d', cmap="YlGnBu")
plt.title("MPRA Activity vs ATAC Accessibility per Cell Type")
plt.xlabel("Cell Type")
plt.ylabel("Tile Category")
plt.tight_layout()
plt.savefig(output_plot)
print(f"Saved heatmap: {output_plot}")