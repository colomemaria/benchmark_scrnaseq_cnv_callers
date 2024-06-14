#!/bin/bash
#SBATCH --error=%J.err
#SBATCH --output=%J.out
#SBATCH -J cellranger_SNU601_subsampled
#SBATCH -p slim16
#SBATCH -c 4 #number of CPUs needed 
#SBATCH --mem-per-cpu=32G  

/work/project/ladcol_005/tools/cellranger-7.0.0/bin/cellranger count \
    --transcriptome=/work/project/ladcol_005/genomes/refdata-gex-GRCh38-2020-A \
    --fastqs=/work/project/ladcol_010/benchmark_scrna_cnv_caller/data/input_SNU601_sample60 \
    --sample=SRR10805153_sample60 --localcores 4 --localmem 32 --id=SNU601_sample60_hg38_mapped