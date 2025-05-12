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
OUTDIR="/projectnb/vtrs/myousry/virus_atac/workdir/${SAMPLE_ID}"
cd "$OUTDIR"

# === Files ===
VIRAL_BAM="${SAMPLE_ID}_${VIRUS_NAME}.bam"
SORTED_BAM="${SAMPLE_ID}_${VIRUS_NAME}_sorted.bam"
VIRAL_BW="${SAMPLE_ID}_${VIRUS_NAME}.bw"
HUMAN_BAM="${SAMPLE_ID}_human.bam"

# === Sort and index viral BAM ===
samtools sort -o "$SORTED_BAM" "$VIRAL_BAM"
samtools index "$SORTED_BAM"

# === Get total number of input read pairs from host BAM ===
TOTAL_READS=$(samtools view -c "$HUMAN_BAM")

if [ "$TOTAL_READS" -eq 0 ]; then
  echo "Error: total read count is zero. Check $HUMAN_BAM"
  exit 1
fi

# === Compute scale factor ===
SCALE=$(echo "scale=10; 1000000000 / $TOTAL_READS" | bc)

echo "Sample: $SAMPLE_ID"
echo "Total reads (from $HUMAN_BAM): $TOTAL_READS"
echo "Scale factor for normalization: $SCALE"

# === Generate bigWig ===
bamCoverage \
  -b "$SORTED_BAM" \
  -o "$VIRAL_BW" \
  --binSize 1 \
  --outFileFormat bigwig
  
echo "Viral bigWig file created: $VIRAL_BW"