#!/bin/bash

# USAGE: ./estimate_copy_number.sh sample_id virus_name
# Example: ./estimate_copy_number.sh SRR1234567 HSV1_KOS

set -e

SAMPLE_ID="$1"
VIRUS_NAME="$2"

if [ -z "$SAMPLE_ID" ] || [ -z "$VIRUS_NAME" ]; then
  echo "Usage: $0 <sample_id> <virus_name>"
  exit 1
fi

# === Input BAM files ===
HOST_BAM="${SAMPLE_ID}_human.bam"
VIRAL_BAM="${SAMPLE_ID}_${VIRUS_NAME}.bam"

# === Genome sizes ===
HUMAN_GENOME_SIZE=3200000000    # Diploid human genome ~3.2 Gb
VIRAL_GENOME_SIZE=152000        # HSV1 KOS genome length

# === Count mapped reads ===
echo "Counting mapped reads..."

HOST_MAPPED=$(samtools view -c -F 4 "$HOST_BAM")
VIRAL_MAPPED=$(samtools view -c -F 4 "$VIRAL_BAM")

# === Calculate estimated viral copy number per nucleus ===
COPY_NUM=$(echo "scale=2; ($VIRAL_MAPPED / $HOST_MAPPED) * ($HUMAN_GENOME_SIZE / $VIRAL_GENOME_SIZE)" | bc)

echo "==============================="
echo "Sample ID: $SAMPLE_ID"
echo "Virus: $VIRUS_NAME"
echo "Host mapped reads: $HOST_MAPPED"
echo "Virus mapped reads: $VIRAL_MAPPED"
echo "Estimated viral genome copies per nucleus: $COPY_NUM"
echo "==============================="