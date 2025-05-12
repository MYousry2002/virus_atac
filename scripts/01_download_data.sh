#!/bin/bash
#$ -N download_data
#$ -cwd
#$ -o ../logs/download_data_$TASK_ID.out
#$ -e ../logs/download_data_$TASK_ID.err
#$ -pe smp 4
#$ -l h_vmem=32G
#$ -l h_rt=2:00:00

# Load conda environment
source /projectnb/bioinfor/myousry/miniconda3/etc/profile.d/conda.sh
conda activate virus_atac

set -euo pipefail

SAMPLES_FILE="../data/samples_list.txt"
OUTDIR="../data"

if [ ! -f "$SAMPLES_FILE" ]; then
  echo "Sample list not found: $SAMPLES_FILE"
  exit 1
fi

echo "Starting download and FASTQ conversion for samples listed in $SAMPLES_FILE."

while read -r SRR_ID; do
  [ -z "$SRR_ID" ] && continue

  SAMPLE_DIR="${OUTDIR}/${SRR_ID}_fastq"
  mkdir -p "$SAMPLE_DIR"

  echo "[$SRR_ID] Downloading and converting to FASTQ..."
  fasterq-dump "$SRR_ID" -O "$SAMPLE_DIR" --split-files -e 4

  echo "[$SRR_ID] Compressing FASTQ files..."
  pigz "$SAMPLE_DIR/${SRR_ID}_1.fastq"
  pigz "$SAMPLE_DIR/${SRR_ID}_2.fastq"

  echo "[$SRR_ID] Completed. Files stored in $SAMPLE_DIR"
  echo
done < "$SAMPLES_FILE"

echo "All downloads and conversions completed."