#!/bin/bash
#SBATCH --error=filter_bam_%J.err
#SBATCH --output=filter_bam_%J.out
#SBATCH -J filter_bam
#SBATCH -p slim16
#SBATCH -c 1 #number of CPUs needed 
#SBATCH --mem-per-cpu=32G 

#Get samtools from cluster
module load ngs/samtools/1.9

#Define input files
BAM_FILE="pbmc_unsorted_10k_gex_possorted_bam.bam"
FILTER_FILE="selected_test_cells.txt"
 
# Save the header lines
samtools view -H $BAM_FILE > SAM_header

# Filter alignments using filter_file. Use LC_ALL=C to set C locale instead of UTF-8
samtools view $BAM_FILE | LC_ALL=C grep -F -f $FILTER_FILE > filtered_SAM_body
 
# Combine header and body
cat SAM_header filtered_SAM_body > filtered.sam

# Convert filtered.sam to BAM format
samtools view -b filtered.sam > filtered.bam
 
# Delete intermediate files
rm filtered_SAM_body filtered.sam SAM_header

# Sort and index the result file (not sure if sorting is strictly necessary)
samtools sort -o filtered_sorted.bam filtered.bam
samtools index filtered_sorted.bam

# Rename files to match the snakemake pipeline
mv filtered_sorted.bam possorted_genome_bam.bam
mv filtered_sorted.bam.bai possorted_genome_bam.bam.bai