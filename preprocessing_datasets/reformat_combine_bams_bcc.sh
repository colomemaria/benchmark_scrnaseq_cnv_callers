#!/bin/bash

#SBATCH -o reformat_bams.txt
#SBATCH -e reformat_bams.txt
#SBATCH -J reformat_bam
#SBATCH -c 1
#SBATCH --mem=64G
#SBATCH -t 48:00:00
#SBATCH -p slim18

#Load the required moduls
module load ngs/samtools/1.9

# Save the header lines
#samtools view -H data/input_BCC/BCC05_hg38_mapped/outs/possorted_genome_bam.bam > data/input_BCC/SAM_header_su005
#samtools view -H data/input_BCC/BCC06_hg38_mapped/outs/possorted_genome_bam.bam > data/input_BCC/SAM_header_su006
#samtools view -H data/input_BCC/BCC07_hg38_mapped/outs/possorted_genome_bam.bam > data/input_BCC/SAM_header_su007
#samtools view -H data/input_BCC08/BCC08_hg38_mapped/outs/possorted_genome_bam.bam > data/input_BCC/SAM_header_su008

#Filter each bam files for the selected barcodes
#samtools view data/input_BCC/BCC05_hg38_mapped/outs/possorted_genome_bam.bam | LC_ALL=C grep -F -f data/input_BCC/filtered_barcodes/barcodes_su005.tsv > data/input_BCC/filtered_SAM_body_su005
#samtools view data/input_BCC/BCC06_hg38_mapped/outs/possorted_genome_bam.bam | LC_ALL=C grep -F -f data/input_BCC/filtered_barcodes/barcodes_su006.tsv > data/input_BCC/filtered_SAM_body_su006
#samtools view data/input_BCC/BCC07_hg38_mapped/outs/possorted_genome_bam.bam | LC_ALL=C grep -F -f data/input_BCC/filtered_barcodes/barcodes_su007.tsv > data/input_BCC/filtered_SAM_body_su007
#samtools view data/input_BCC08/BCC08_hg38_mapped/outs/possorted_genome_bam.bam | LC_ALL=C grep -F -f data/input_BCC/filtered_barcodes/barcodes_su008.tsv > data/input_BCC/filtered_SAM_body_su008

#Rename barcode tag in filename to represent also the sample name
#sed -e 's/CB:Z:/CB:Z:bcc.su005.pre.tumor_/' data/input_BCC/filtered_SAM_body_su005 > data/input_BCC/renamed_SAM_body_su005
#sed -e 's/CB:Z:/CB:Z:bcc.su006.pre.tumor_/' data/input_BCC/filtered_SAM_body_su006 > data/input_BCC/renamed_SAM_body_su006
#sed -e 's/CB:Z:/CB:Z:bcc.su007.pre.tumor.cd45_/' data/input_BCC/filtered_SAM_body_su007 > data/input_BCC/renamed_SAM_body_su007
#sed -e 's/CB:Z:/CB:Z:bcc.su008.pre.tumor_/' data/input_BCC/filtered_SAM_body_su008 > data/input_BCC/renamed_SAM_body_su008

# Combine header and body back into a full sam file
#cat data/input_BCC/SAM_header_su005 data/input_BCC/renamed_SAM_body_su005 > data/input_BCC/filtered_su005.sam
#cat data/input_BCC/SAM_header_su006 data/input_BCC/renamed_SAM_body_su006 > data/input_BCC/filtered_su006.sam
#cat data/input_BCC/SAM_header_su007 data/input_BCC/renamed_SAM_body_su007 > data/input_BCC/filtered_su007.sam
#cat data/input_BCC/SAM_header_su008 data/input_BCC/renamed_SAM_body_su008 > data/input_BCC/filtered_su008.sam

# Delete intermediate files
#rm data/input_BCC/SAM_header_su005 data/input_BCC/filtered_SAM_body_su005 data/input_BCC/renamed_SAM_body_su005 
#rm data/input_BCC/SAM_header_su006 data/input_BCC/filtered_SAM_body_su006 data/input_BCC/renamed_SAM_body_su006
#rm data/input_BCC/SAM_header_su007 data/input_BCC/filtered_SAM_body_su007 data/input_BCC/renamed_SAM_body_su007
#rm data/input_BCC/SAM_header_su008 data/input_BCC/filtered_SAM_body_su008 data/input_BCC/renamed_SAM_body_su008

# Convert filtered.sam to BAM format
#samtools view -b data/input_BCC/filtered_su005.sam > data/input_BCC/filtered_su005.bam
#samtools view -b data/input_BCC/filtered_su006.sam > data/input_BCC/filtered_su006.bam
#samtools view -b data/input_BCC/filtered_su007.sam > data/input_BCC/filtered_su007.bam
#samtools view -b data/input_BCC/filtered_su008.sam > data/input_BCC/filtered_su008.bam

# Delete intermediate files
#rm data/input_BCC/filtered_su005.sam data/input_BCC/filtered_su006.sam data/input_BCC/filtered_su007.sam data/input_BCC/filtered_su008.sam

# Index the filtered result files
#samtools index data/input_BCC/filtered_su005.bam
#samtools index data/input_BCC/filtered_su006.bam
#samtools index data/input_BCC/filtered_su007.bam
#samtools index data/input_BCC/filtered_su008.bam

#Combine bam files to one large
samtools merge data/input_BCC/combined_samples.bam data/input_BCC/filtered_su005.bam data/input_BCC/filtered_su006.bam data/input_BCC/filtered_su007.bam data/input_BCC/filtered_su008.bam
samtools index data/input_BCC/combined_samples.bam