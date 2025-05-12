#!/bin/bash

# USAGE: ./02_alignment.sh sample_id virus_name
# Example: ./02_alignment.sh SRR1234567 HSV1_KOS

set -euo pipefail

# ===== Inputs ========
SAMPLE_ID="$1"
VIRUS_NAME="$2"

if [ -z "$SAMPLE_ID" ] || [ -z "$VIRUS_NAME" ]; then
  echo "Usage: $0 <sample_id> <virus_name>"
  exit 1
fi

# Directories
DATADIR="/projectnb/vtrs/myousry/virus_atac/data/"
OUTDIR="/projectnb/vtrs/myousry/virus_atac/workdir/${SAMPLE_ID}"
mkdir -p "$OUTDIR"
cd "$OUTDIR"

# FASTQ files
READ1="${DATADIR}/${SAMPLE_ID}_1.fastq.gz"
READ2="${DATADIR}/${SAMPLE_ID}_2.fastq.gz"

# Genome indices
HG38_INDEX="/projectnb/vtrs/myousry/genomes/human/bowtie2_index/hg38"
VIRAL_INDEX="/projectnb/vtrs/myousry/genomes/herpesviruses/bowtie2_index/${VIRUS_NAME}/${VIRUS_NAME}"

# Output names
HUMAN_SAM="${SAMPLE_ID}_human.sam"
HUMAN_BAM="${SAMPLE_ID}_human.bam"
UNMAPPED_BAM="${SAMPLE_ID}_unmapped.bam"
UNMAPPED_SORTED="${SAMPLE_ID}_unmapped_sorted.bam"
UNMAPPED_R1="${SAMPLE_ID}_unmapped_R1.fastq"
UNMAPPED_R2="${SAMPLE_ID}_unmapped_R2.fastq"
VIRAL_SAM="${SAMPLE_ID}_${VIRUS_NAME}.sam"
VIRAL_BAM="${SAMPLE_ID}_${VIRUS_NAME}.bam"

# ==== Align to human genome =====
bowtie2 \
  --local \
  --very-sensitive-local \
  --no-discordant \
  --no-mixed \
  --contain \
  --overlap \
  --dovetail \
  --phred33 \
  -x "$HG38_INDEX" \
  -1 "$READ1" -2 "$READ2" \
  -S "$HUMAN_SAM"

# Convert SAM to BAM
samtools view -bS "$HUMAN_SAM" > "$HUMAN_BAM"

# Extract unmapped reads
samtools view -b -f 12 -F 256 "$HUMAN_BAM" > "$UNMAPPED_BAM"

# Sort for FASTQ conversion
samtools sort -n "$UNMAPPED_BAM" -o "$UNMAPPED_SORTED"

# Convert to FASTQ
bedtools bamtofastq -i "$UNMAPPED_SORTED" -fq "$UNMAPPED_R1" -fq2 "$UNMAPPED_R2"

# ==== Align to viral genome ====
bowtie2 \
  --local \
  --very-sensitive-local \
  --no-discordant \
  --no-mixed \
  --contain \
  --overlap \
  --dovetail \
  --phred33 \
  -x "$VIRAL_INDEX" \
  -1 "$UNMAPPED_R1" -2 "$UNMAPPED_R2" \
  -S "$VIRAL_SAM"

# Convert viral SAM to BAM
samtools view -bS "$VIRAL_SAM" > "$VIRAL_BAM"

# Cleanup
rm "$HUMAN_SAM" "$VIRAL_SAM" "$UNMAPPED_BAM"