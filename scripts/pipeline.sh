#!/bin/bash
#$ -N virus_atac_pipeline
#$ -cwd
#$ -o ../logs/pipeline_$TASK_ID.out
#$ -e ../logs/pipeline_$TASK_ID.err
#$ -pe smp 8
#$ -l h_vmem=64G
#$ -l h_rt=8:00:00
#$ -t 1-12

# === Load conda environment ===
source /projectnb/bioinfor/myousry/miniconda3/etc/profile.d/conda.sh
conda activate virus_atac

set -euo pipefail

# === Input variables ===
SAMPLES_FILE="/projectnb/vtrs/myousry/virus_atac/data/HSV1_17/samples_list.txt"
VIRUS_NAME="HSV1_17"
VIRAL_GENOME_SIZE=152000
OUTDIR="/projectnb/vtrs/myousry/virus_atac/workdir/HSV1_17"

# Get the sample ID for this task
SAMPLE_ID=$(sed -n "${SGE_TASK_ID}p" "$SAMPLES_FILE")
[ -z "$SAMPLE_ID" ] && { echo "Sample ID not found for task $SGE_TASK_ID"; exit 1; }

cd "$OUTDIR"

# Alignment
bash /projectnb/vtrs/myousry/virus_atac/scripts/02_alignment.sh "$SAMPLE_ID" "$VIRUS_NAME" PE

# Generate bigWig
bash /projectnb/vtrs/myousry/virus_atac/scripts/03_generate_bigwig.sh "$SAMPLE_ID" "$VIRUS_NAME" "$VIRAL_GENOME_SIZE"

echo "==== Completed $SAMPLE_ID ===="