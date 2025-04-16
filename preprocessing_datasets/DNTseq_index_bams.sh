#!/bin/bash
#SBATCH --job-name=index
#SBATCH --mem-per-cpu=32G
#SBATCH --output=index_bam%j.out
#SBATCH --error=index_bam%j.err
#SBATCH -c 1
#SBATCH --time=48:00:00

#Stop script once an error occurs
set -e

#Load the samtools modul
module load ngs/samtools/1.9

#For the HCT116 dataset
#cd scRNA_benchmark/revisions/DNTRseq/mRNA/aligned_dedup

#For the A375 dataset
cd scRNA_benchmark/revisions/DNTRseq/A375/mRNA/aligned_dedup

for fl in *.bam
do 
    samtools index $fl
done


