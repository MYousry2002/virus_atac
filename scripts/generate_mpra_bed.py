import pandas as pd

# Load CSV files
activity_df = pd.read_csv("../MPRA/tiles_cell_specific_activity.csv")   # contains signal and activity columns
coords_df = pd.read_csv("../MPRA/genbank_ids_and_tile_positions.csv")       # contains tile_id, genbank_id, start_position, alignment_identity

# Filter for HSV1_KOS tiles
activity_df = activity_df[activity_df['tile_id'].str.contains("Herpesvirus:Herpes_Simplex_1_KOS")]
coords_df = coords_df[coords_df['tile_id'].str.contains("Herpesvirus:Herpes_Simplex_1_KOS")]

# Merge both tables on tile_id
merged = pd.merge(coords_df, activity_df, on="tile_id")

# Create BED fields
merged['chrom'] = merged['genbank_id']
merged['start'] = merged['start_position']
merged['end'] = merged['start'] + merged['alignment_identity']
merged['name'] = merged['tile_id']
merged['score'] = 0
merged['strand'] = merged['tile_id'].str.extract(r':(\+|\-)$')[0]

# Select columns for BED format with signal + activity
signals = ['GM12878', 'Jurkat', 'MRC5', 'A549', 'HEK293', 'K562']
activities = [f"{cell}_is_active" for cell in signals]
bed_columns = ['chrom', 'start', 'end', 'name', 'score', 'strand'] + signals + activities

# Save to file
merged[bed_columns].to_csv("../MPRA/HSV1_KOS_MPRA_with_signals.bed", sep="\t", index=False, header=False)