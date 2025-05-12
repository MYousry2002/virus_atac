#!/bin/bash

# USAGE: ./02_alignment_continue.sh sample_id virus_name
# Example: ./02_alignment_continue.sh SRR1234567 HSV1_KOS

set -euo pipefail

# ===== Inputs ========
SAMPLE_ID="$1"
VIRUS_NAME="$2"

if [ -z "$SAMPLE_ID" ] || [ -z "$VIRUS_NAME" ]; then
  echo "Usage: $0 <sample_id> <virus_name>"
  exit 1
fi

# === Directories ===
DATADIR="/projectnb/vtrs/myousry/virus_atac/data"
OUTDIR="/projectnb/vtrs/myousry/virus_atac/workdir/${SAMPLE_ID}"
mkdir -p "$OUTDIR"
cd "$OUTDIR"

# === Files and Indices ===
HUMAN_BAM="${SAMPLE_ID}_human.bam"
UNMAPPED_BAM="${SAMPLE_ID}_unmapped.bam"
UNMAPPED_SORTED="${SAMPLE_ID}_unmapped_sorted.bam"
UNMAPPED_R1="${SAMPLE_ID}_unmapped_R1.fastq"
UNMAPPED_R2="${SAMPLE_ID}_unmapped_R2.fastq"
VIRAL_SAM="${SAMPLE_ID}_${VIRUS_NAME}.sam"
VIRAL_BAM="${SAMPLE_ID}_${VIRUS_NAME}.bam"
VIRAL_INDEX="/projectnb/vtrs/myousry/genomes/herpesviruses/bowtie2_index/${VIRUS_NAME}/${VIRUS_NAME}"

# === Extract unmapped read pairs ===
echo "Extracting unmapped reads from: $HUMAN_BAM"
samtools view -b -f 4 -F 264 "$HUMAN_BAM" > "$UNMAPPED_BAM"

echo "Checking BAM line count:"
samtools view "$UNMAPPED_BAM" | wc -l

# === Sort BAM for FASTQ conversion ===
echo "Sorting BAM by read name..."
samtools sort -n "$UNMAPPED_BAM" -o "$UNMAPPED_SORTED"

# === Convert to FASTQ ===
echo "Converting to FASTQ..."
bedtools bamtofastq \
  -i "$UNMAPPED_SORTED" \
  -fq "$UNMAPPED_R1" \
  -fq2 "$UNMAPPED_R2"

echo "FASTQ file line counts:"
wc -l "$UNMAPPED_R1" "$UNMAPPED_R2"

# === Align to viral genome ===
echo "Aligning to viral genome: $VIRUS_NAME"
bowtie2 \
  --no-unal \
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

# === Convert viral SAM to BAM ===
samtools view -bS "$VIRAL_SAM" > "$VIRAL_BAM"

# === Cleanup ===
rm "$VIRAL_SAM" "$UNMAPPED_BAM"

echo "Viral alignment complete for $SAMPLE_ID"