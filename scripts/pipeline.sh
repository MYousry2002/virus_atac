#!/bin/bash
#$ -N virus_atac_pipeline
#$ -cwd
#$ -o ../logs/pipeline_$TASK_ID.out
#$ -e ../logs/pipeline_$TASK_ID.err
#$ -pe smp 8
#$ -l h_vmem=64G
#$ -l h_rt=8:00:00

# === Load conda environment ===
source /projectnb/bioinfor/myousry/miniconda3/etc/profile.d/conda.sh
conda activate virus_atac

set -euo pipefail

# === Input variables ===
SAMPLES_FILE="/projectnb/vtrs/myousry/virus_atac/data/HCMV/samples_list.txt"
VIRUS_NAME="HCMV"
VIRAL_GENOME_SIZE=152000
OUTDIR="/projectnb/vtrs/myousry/virus_atac/workdir/HCMV"

if [ ! -f "$SAMPLES_FILE" ]; then
  echo "Sample list not found: $SAMPLES_FILE"
  exit 1
fi

cd "$OUTDIR"

while read -r SAMPLE_ID; do
  [ -z "$SAMPLE_ID" ] && continue
  # echo "==== Processing $SAMPLE_ID ===="

  # Alignment
  bash /projectnb/vtrs/myousry/virus_atac/scripts/02_alignment.sh "$SAMPLE_ID" "$VIRUS_NAME" SE

  # Generate bigWig
  bash /projectnb/vtrs/myousry/virus_atac/scripts/03_generate_bigwig.sh "$SAMPLE_ID" "$VIRUS_NAME" "$VIRAL_GENOME_SIZE"

  echo "==== Completed $SAMPLE_ID ===="
  echo

done < "$SAMPLES_FILE"

echo "All samples processed."