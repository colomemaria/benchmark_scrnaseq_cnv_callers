#!/bin/bash
#SBATCH --error=%J.err
#SBATCH --output=%J.out
#SBATCH -J cellranger_MCF7
#SBATCH -p fat
#SBATCH -c 6 #number of CPUs needed 
#SBATCH --mem-per-cpu=64G 

/work/project/ladcol_005/tools/cellranger-7.0.0/bin/cellranger count \
    --transcriptome=/work/project/ladcol_005/genomes/refdata-gex-GRCh38-2020-A \
    --fastqs=/work/project/ladcol_010/benchmark_scrna_cnv_caller/data/input_MCF7 \
    --sample=SRR10018149 --localcores 6 --localmem 64 --id=MCF7_hg38_mapped