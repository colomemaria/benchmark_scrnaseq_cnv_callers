#!/bin/bash
#SBATCH --error=%J.err
#SBATCH --output=%J.out
#SBATCH -J cellranger_KATOIII
#SBATCH -p slim18
#SBATCH -c 6 #number of CPUs needed 
#SBATCH --mem=128G 
#SBATCH -t 48:00:00

/work/project/ladcol_005/tools/cellranger-7.0.0/bin/cellranger count \
    --transcriptome=/work/project/ladcol_005/genomes/refdata-gex-GRCh38-2020-A \
    --fastqs=/work/project/ladcol_004/gastricCellLines/scRNA \
    --sample=SRR10805147 --localcores 6 --localmem 64 --id=KATOIII_hg38_mapped
