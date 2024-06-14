#!/bin/bash
#SBATCH --error=%J.err
#SBATCH --output=%J.out
#SBATCH -J cellranger_SNU601
#SBATCH -p fat
#SBATCH -c 6 #number of CPUs needed 
#SBATCH --mem-per-cpu=64G 

/work/project/ladcol_005/tools/cellranger-7.0.0/bin/cellranger count \
    --transcriptome=/work/project/ladcol_005/genomes/refdata-gex-GRCh38-2020-A \
    --fastqs=/work/project/ladcol_005/epiAneufinder/data/SNU601 \
    --sample=SRR10805153 --localcores 6 --localmem 64 --id=SNU601_hg38_mapped