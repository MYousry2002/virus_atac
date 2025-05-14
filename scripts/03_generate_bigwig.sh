#!/bin/bash

# === Inputs ===
SAMPLE_ID="$1"
VIRUS_NAME="$2"
VIRAL_GENOME_SIZE="$3"     # e.g., 152000

if [ $# -ne 3 ]; then
  echo "Usage: $0 <sample_id> <virus_name> <viral_genome_size>"
  exit 1
fi

# === Directories ===
OUTDIR="/projectnb/vtrs/myousry/virus_atac/workdir/${VIRUS_NAME}/${SAMPLE_ID}"
cd "$OUTDIR"

# === Files ===
VIRAL_BAM="${SAMPLE_ID}_${VIRUS_NAME}.bam"
SORTED_BAM="${SAMPLE_ID}_${VIRUS_NAME}_sorted.bam"
VIRAL_BW="${SAMPLE_ID}_${VIRUS_NAME}.bw"


# === Sort and index viral BAM ===
samtools sort -o "$SORTED_BAM" "$VIRAL_BAM"
samtools index "$SORTED_BAM"


# === Generate bigWig ===
bamCoverage \
  -b "$SORTED_BAM" \
  -o "$VIRAL_BW" \
  --binSize 1 \
  --outFileFormat bigwig \
  --normalizeUsing CPM
  
echo "Viral bigWig file created: $VIRAL_BW"