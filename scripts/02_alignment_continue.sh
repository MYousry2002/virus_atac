#!/bin/bash

# USAGE: ./align_to_virus_only.sh sample_id virus_name PE|SE
# Example: ./align_to_virus_only.sh SRR1234567 HSV1_KOS PE

set -euo pipefail

# Inputs
SAMPLE_ID="$1"
VIRUS_NAME="$2"
MODE="${3:-PE}"  # Default to PE

if [[ -z "$SAMPLE_ID" || -z "$VIRUS_NAME" || ! "$MODE" =~ ^(PE|SE)$ ]]; then
  echo "Usage: $0 <sample_id> <virus_name> PE|SE"
  exit 1
fi

# Paths
DATADIR="/projectnb/vtrs/myousry/virus_atac/data/${VIRUS_NAME}/"
OUTDIR="/projectnb/vtrs/myousry/virus_atac/workdir/${VIRUS_NAME}/${SAMPLE_ID}"
mkdir -p "$OUTDIR"
cd "$OUTDIR"

# Viral genome index
VIRAL_INDEX="/projectnb/vtrs/myousry/genomes/herpesviruses/bowtie2_index/${VIRUS_NAME}/${VIRUS_NAME}"

# Output
VIRAL_SAM="${SAMPLE_ID}_${VIRUS_NAME}.sam"
VIRAL_BAM="${SAMPLE_ID}_${VIRUS_NAME}.bam"

if [ "$MODE" == "PE" ]; then
  UNMAPPED_R1="${SAMPLE_ID}_unmapped_R1.fastq"
  UNMAPPED_R2="${SAMPLE_ID}_unmapped_R2.fastq"

  if [[ ! -f "$UNMAPPED_R1" || ! -f "$UNMAPPED_R2" ]]; then
    echo "Missing unmapped FASTQ files: $UNMAPPED_R1 and/or $UNMAPPED_R2"
    exit 1
  fi

  echo "Aligning PE unmapped reads to $VIRUS_NAME..."

  bowtie2 --local --very-sensitive-local --no-discordant --no-mixed --contain --overlap --dovetail --phred33 \
    -x "$VIRAL_INDEX" -1 "$UNMAPPED_R1" -2 "$UNMAPPED_R2" -S "$VIRAL_SAM"

else
  UNMAPPED_SE="${SAMPLE_ID}_unmapped_SE.fastq"

  if [[ ! -f "$UNMAPPED_SE" ]]; then
    echo "Missing unmapped FASTQ file: $UNMAPPED_SE"
    exit 1
  fi

  echo "Aligning SE unmapped reads to $VIRUS_NAME..."

  bowtie2 --local --very-sensitive-local --phred33 \
    -x "$VIRAL_INDEX" -U "$UNMAPPED_SE" -S "$VIRAL_SAM"
fi

# Convert SAM to BAM
samtools view -bS "$VIRAL_SAM" > "$VIRAL_BAM"
rm "$VIRAL_SAM"

echo "Viral alignment complete: $VIRAL_BAM"