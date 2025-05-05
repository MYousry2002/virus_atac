#!/bin/bash
#$ -N genome_bowtie2
#$ -cwd
#$ -o ../logs/genome_bowtie2$TASK_ID.out
#$ -e ../logs/genome_bowtie2$TASK_ID.err
#$ -pe smp 4
#$ -l h_vmem=64G
#$ -l h_rt=2:00:00

# Human genome
cd ../../genomes/human/
mkdir -p bowtie2_index
bowtie2-build Homo_sapiens.GRCh38.dna.primary_assembly.fa bowtie2_index/hg38

# Herpesvirus genome (e.g., HSV1_KOS)
cd ../herpesvirus/
mkdir -p bowtie2_index/HSV1_KOS
bowtie2-build HSV1_KOS.fasta bowtie2_index/HSV1_KOS