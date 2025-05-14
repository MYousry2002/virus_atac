#!/bin/bash
#$ -N genome_bowtie2
#$ -cwd
#$ -o ../logs/genome_bowtie2$TASK_ID.out
#$ -e ../logs/genome_bowtie2$TASK_ID.err
#$ -pe smp 8
#$ -l h_vmem=64G
#$ -l h_rt=2:00:00

# conda env
source /projectnb/bioinfor/myousry/miniconda3/etc/profile.d/conda.sh
conda activate virus_atac

# Human genome
# cd /projectnb/vtrs/myousry/genomes/human/
# mkdir -p bowtie2_index
# bowtie2-build Homo_sapiens.GRCh38.dna.primary_assembly.fa bowtie2_index/hg38

# Herpesvirus genome (e.g., HSV1_KOS)
cd /projectnb/vtrs/myousry/genomes/herpesviruses/
mkdir -p bowtie2_index/HCMV
bowtie2-build HCMV.fasta bowtie2_index/HCMV