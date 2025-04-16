#!/bin/bash
#SBATCH --job-name=bamtools
#SBATCH --mem-per-cpu=32G
#SBATCH --output=bamtools_%j.out
#SBATCH --error=bamtools_%j.err
#SBATCH -c 1
#SBATCH --time=48:00:00

#Stop script once an error occurs
set -e

#Load the bamtools modul
module load ngs/bamtools/2.5.1
module load ngs/samtools/1.9

#For the HCT116 dataset
bamtools merge -list snakemake_pipeline/data/input_HCT116/list_bam_files_filtered.txt -out HCT116_DNTRseq_merged_unsorted.bam
samtools sort HCT116_DNTRseq_merged_unsorted.bam -o HCT116_DNTRseq_merged.bam
samtools index HCT116_DNTRseq_merged.bam

#For the A375 dataset
bamtools merge -list snakemake_pipeline/data/input_A375/list_bam_files_filtered.txt -out A375_DNTRseq_merged_unsorted.bam 
samtools sort A375_DNTRseq_merged_unsorted.bam -o A375_DNTRseq_merged.bam
samtools index A375_DNTRseq_merged.bam


