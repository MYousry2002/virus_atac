#!/bin/bash

# === Inputs ===
SAMPLE_ID="$1"
VIRUS_NAME="$2"
VIRAL_GENOME_SIZE="$3"     # e.g., 152000
COPY_NUMBER="$4"           # e.g., 300 viral genomes per nucleus

if [ $# -ne 4 ]; then
  echo "Usage: $0 <sample_id> <virus_name> <viral_genome_size> <genome_copy_number>"
  exit 1
fi

# === Filenames ===
VIRAL_BAM="${SAMPLE_ID}_${VIRUS_NAME}.bam"
SORTED_BAM="${SAMPLE_ID}_${VIRUS_NAME}_sorted.bam"
VIRAL_BW="${SAMPLE_ID}_${VIRUS_NAME}.bw"

# === Sort and index viral BAM ===
samtools sort -o "$SORTED_BAM" "$VIRAL_BAM"
samtools index "$SORTED_BAM"

# === Compute scaling factor ===
# Scale = 1e9 / genome copy number (mapped reads per billion per copy)
SCALE=$(echo "scale=10; 1000000000 / $COPY_NUMBER" | bc)

# === Generate bigWig (deepTools) ===
bamCoverage \
  -b "$SORTED_BAM" \
  -o "$VIRAL_BW" \
  --binSize 1 \
  --scaleFactor "$SCALE" \
  --normalizeUsing None \
  --effectiveGenomeSize "$VIRAL_GENOME_SIZE" \
  --extendReads \
  --centerReads \
  --outFileFormat bigwig

echo "Viral bigWig file created: $VIRAL_BW"