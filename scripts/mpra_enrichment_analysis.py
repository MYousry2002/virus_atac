
import pandas as pd
from scipy.stats import fisher_exact
import matplotlib.pyplot as plt
import seaborn as sns

# === Config ===
input_table = "../results/mpra_atac_category_table.csv"
output_stats = "../results/mpra_enrichment_fisher_results.csv"
output_plot = "../results/mpra_enrichment_odds_ratio.png"

cell_types = ['GM12878', 'Jurkat', 'MRC5', 'A549', 'HEK293', 'K562']


####### Do active tiles fall into accessible regions more than inactive tiles do? ####


# === Load data ===
df = pd.read_csv(input_table)

results = []

for cell in cell_types:
    category_col = f"{cell}_category"
    counts = df[category_col].value_counts()

    active_accessible = counts.get("Active & Accessible", 0)
    active_only = counts.get("Active only", 0)
    accessible_only = counts.get("Accessible only", 0)
    neither = counts.get("Neither", 0)

    # Contingency table:
    #                Accessible     Not Accessible
    # Active             a                b
    # Inactive           c                d
    a = active_accessible
    b = active_only
    c = accessible_only   
    d = neither

    table = [[a, b], [c, d]]
    odds_ratio, p_value = fisher_exact(table, alternative='greater')

    results.append({
        "Cell Type": cell,
        "Active & Accessible": a,
        "Active only": b,
        "Accessible only": c,
        "Neither": d,
        "Odds Ratio": odds_ratio,
        "P-value": p_value
    })

# === Save results ===
res_df = pd.DataFrame(results)
res_df.to_csv(output_stats, index=False)
print(f"Saved Fisher test results: {output_stats}")

# === Plot odds ratios ===
plt.figure(figsize=(8, 5))
sns.barplot(x="Cell Type", y="Odds Ratio", data=res_df)
plt.axhline(1, ls="--", color="gray")
plt.ylabel("Odds Ratio (Active vs Inactive in Accessible regions)")
plt.title("Enrichment of MPRA Activity in Accessible Chromatin")
plt.tight_layout()
plt.savefig(output_plot)
print(f"Saved odds ratio plot: {output_plot}")



###### Are active tiles enriched in accessible regions compared to the background genome (all other tiles)? #####

# === Load data ===
df = pd.read_csv(input_table)

results = []

for cell in cell_types:
    category_col = f"{cell}_category"
    active_col = f"{cell}_active"
    
    df[f"{cell}_accessible"] = df[category_col].isin(["Active & Accessible", "Accessible only"])
    
    active = df[active_col].astype(bool)
    accessible = df[f"{cell}_accessible"]

    # a: Active & Accessible
    a = ((active) & (accessible)).sum()
    # c: Active & Not Accessible
    c = ((active) & (~accessible)).sum()
    
    # b: Inactive & Accessible
    b = ((~active) & (accessible)).sum()
    # d: Inactive & Not Accessible
    d = ((~active) & (~accessible)).sum()
    
    table = [[a, c], [b, d]]
    odds_ratio, p_value = fisher_exact(table, alternative="greater")
    
    results.append({
        "Cell Type": cell,
        "Active & Accessible": a,
        "Active & Not Accessible": c,
        "Background Accessible": b,
        "Background Not Accessible": d,
        "Odds Ratio": odds_ratio,
        "P-value": p_value
    })

# === Save results ===
res_df = pd.DataFrame(results)
res_df.to_csv(output_stats, index=False)
print(f"Saved Fisher enrichment test results: {output_stats}")

# === Plot ===
plt.figure(figsize=(10, 5))
sns.barplot(x="Cell Type", y="Odds Ratio", data=res_df)
plt.axhline(1, ls="--", color="gray")
plt.title("Enrichment of Active Tiles in Accessible Chromatin")
plt.ylabel("Odds Ratio (Active vs Background)")
plt.tight_layout()
plt.savefig(output_plot)
print(f"Saved plot: {output_plot}")