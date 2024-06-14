#!/bin/bash
#SBATCH --error=%J.err
#SBATCH --output=%J.out
#SBATCH -J cellranger_BCC
#SBATCH -p fat
#SBATCH -c 6 #number of CPUs needed 
#SBATCH --mem-per-cpu=32G 

#Process sample 05
/work/project/ladcol_005/tools/cellranger-7.0.0/bin/cellranger count \
    --transcriptome=/work/project/ladcol_005/genomes/refdata-gex-GRCh38-2020-A \
    --fastqs=/work/project/ladcol_010/benchmark_scrna_cnv_caller/data/input_BCC \
    --sample=SRR8315751 --localcores 6 --localmem 64 --id=BCC05_hg38_mapped

#Process sample 06
#/work/project/ladcol_005/tools/cellranger-7.0.0/bin/cellranger count \
#    --transcriptome=/work/project/ladcol_005/genomes/refdata-gex-GRCh38-2020-A \
#    --fastqs=/work/project/ladcol_010/benchmark_scrna_cnv_caller/data/input_BCC \
#    --sample=SRR8315756 --localcores 6 --localmem 64 --id=BCC06_hg38_mapped

#Process sample 07
/work/project/ladcol_005/tools/cellranger-7.0.0/bin/cellranger count \
    --transcriptome=/work/project/ladcol_005/genomes/refdata-gex-GRCh38-2020-A \
    --fastqs=/work/project/ladcol_010/benchmark_scrna_cnv_caller/data/input_BCC \
    --sample=SRR8315761 --localcores 6 --localmem 64 --id=BCC07_hg38_mapped
