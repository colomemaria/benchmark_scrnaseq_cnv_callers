#!/bin/bash
#SBATCH --job-name=mouse_pileup
#SBATCH --mem-per-cpu=32G
#SBATCH --output=logs/mouse_sortindex_%j.out
#SBATCH --error=logs/mouse_sortindex_%j.err
#SBATCH -c 1
#SBATCH --time=48:00:00

#Stop script once an error occurs
set -e

module load ngs/samtools/1.9

cd snakemake_pipeline/data/input_mouse

samtools sort T989Aligned.out.bam -o possorted_genome_bam.bam
samtools index possorted_genome_bam.bam

