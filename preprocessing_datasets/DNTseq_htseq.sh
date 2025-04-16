#!/bin/bash
#SBATCH --job-name=htseq
#SBATCH --mem-per-cpu=32G
#SBATCH --output=htseq_%j.out
#SBATCH --error=htseq_%j.err
#SBATCH -c 1
#SBATCH --time=48:00:00

#Stop script once an error occurs
set -e

#Load the htseq modul
#module load ngs/HTSeq/0.10.0

#Activate the right conda environment (this contains a newer htseq version than module)
source miniconda3/etc/profile.d/conda.sh
conda activate emt

#For the HCT116 dataset
#cd scRNA_benchmark/revisions/DNTRseq/mRNA/aligned_dedup

#For the A375 dataset
cd scRNA_benchmark/revisions/DNTRseq/A375/mRNA/aligned_dedup

#Create one count matrix
htseq-count -m intersection-nonempty -s no --additional-attr=gene_name *.bam \
 STAR_reference/Homo_sapiens.GRCh38.113.gtf > ../count_matrix_DNT.txt

#Save the header
echo *.bam > ../count_matrix_DNT_header.txt


