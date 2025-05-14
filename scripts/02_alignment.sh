#!/bin/bash

# USAGE: ./02_alignment.sh sample_id virus_name PE|SE
# Example: ./02_alignment.sh SRR1234567 HSV1_KOS PE

set -euo pipefail

# ===== Inputs ========
SAMPLE_ID="$1"
VIRUS_NAME="$2"
MODE="${3:-PE}"  # Default is PE

if [[ -z "$SAMPLE_ID" || -z "$VIRUS_NAME" || ! "$MODE" =~ ^(PE|SE)$ ]]; then
  echo "Usage: $0 <sample_id> <virus_name> PE|SE"
  exit 1
fi

# Directories
DATADIR="/projectnb/vtrs/myousry/virus_atac/data/${VIRUS_NAME}/"
OUTDIR="/projectnb/vtrs/myousry/virus_atac/workdir/${VIRUS_NAME}/${SAMPLE_ID}"
mkdir -p "$OUTDIR"
cd "$OUTDIR"

# FASTQ files
READ1="${DATADIR}/${SAMPLE_ID}_1.fastq.gz"
READ2="${DATADIR}/${SAMPLE_ID}_2.fastq.gz"
READ_SE="${DATADIR}/${SAMPLE_ID}.fastq"  # Expected for SE

# Genome indices
HG38_INDEX="/projectnb/vtrs/myousry/genomes/human/bowtie2_index/hg38"
VIRAL_INDEX="/projectnb/vtrs/myousry/genomes/herpesviruses/bowtie2_index/${VIRUS_NAME}/${VIRUS_NAME}"

# Output names
HUMAN_SAM="${SAMPLE_ID}_human.sam"
HUMAN_BAM="${SAMPLE_ID}_human.bam"
UNMAPPED_BAM="${SAMPLE_ID}_unmapped.bam"
UNMAPPED_SORTED="${SAMPLE_ID}_unmapped_sorted.bam"
VIRAL_SAM="${SAMPLE_ID}_${VIRUS_NAME}.sam"
VIRAL_BAM="${SAMPLE_ID}_${VIRUS_NAME}.bam"

if [ "$MODE" == "PE" ]; then
  UNMAPPED_R1="${SAMPLE_ID}_unmapped_R1.fastq"
  UNMAPPED_R2="${SAMPLE_ID}_unmapped_R2.fastq"

  # Align to human genome (PE)
  bowtie2 --local --very-sensitive-local --no-discordant --no-mixed --contain --overlap --dovetail --phred33 \
    -x "$HG38_INDEX" -1 "$READ1" -2 "$READ2" -S "$HUMAN_SAM"

  samtools view -bS "$HUMAN_SAM" > "$HUMAN_BAM"
  samtools view -b -f 12 -F 256 "$HUMAN_BAM" > "$UNMAPPED_BAM"
  samtools sort -n "$UNMAPPED_BAM" -o "$UNMAPPED_SORTED"
  bedtools bamtofastq -i "$UNMAPPED_SORTED" -fq "$UNMAPPED_R1" -fq2 "$UNMAPPED_R2"

  # Align to viral genome (PE)
  bowtie2 --local --very-sensitive-local --no-discordant --no-mixed --contain --overlap --dovetail --phred33 \
    -x "$VIRAL_INDEX" -1 "$UNMAPPED_R1" -2 "$UNMAPPED_R2" -S "$VIRAL_SAM"

else
  """
  UNMAPPED_SE="${SAMPLE_ID}_unmapped_SE.fastq"

  # Align to human genome (SE)
  bowtie2 --local --very-sensitive-local --phred33 \
    -x "$HG38_INDEX" -U "$READ_SE" -S "$HUMAN_SAM"

  samtools view -bS "$HUMAN_SAM" > "$HUMAN_BAM"
  samtools view -b -f 4 "$HUMAN_BAM" > "$UNMAPPED_BAM"
  samtools sort -n "$UNMAPPED_BAM" -o "$UNMAPPED_SORTED"
  bedtools bamtofastq -i "$UNMAPPED_SORTED" -fq "$UNMAPPED_SE"
  """
  
  # Align to viral genome (SE)
  bowtie2 --local --very-sensitive-local --phred33 \
    -x "$VIRAL_INDEX" -U "$UNMAPPED_SE" -S "$VIRAL_SAM"
fi

samtools view -bS "$VIRAL_SAM" > "$VIRAL_BAM"

# Cleanup
rm "$HUMAN_SAM" "$VIRAL_SAM" "$UNMAPPED_BAM"