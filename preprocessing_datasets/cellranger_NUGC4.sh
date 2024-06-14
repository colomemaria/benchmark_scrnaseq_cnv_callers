#!/bin/bash
#SBATCH --error=%J.err
#SBATCH --output=%J.out
#SBATCH -J cellranger_NUGC4
#SBATCH -p slim16
#SBATCH -c 4 #number of CPUs needed 
#SBATCH --mem-per-cpu=32G 

/store24/project24/ladcol_005/tools/cellranger-7.0.0/bin/cellranger count \
    --transcriptome=/store24/project24/ladcol_005/genomes/refdata-gex-GRCh38-2020-A \
    --fastqs=/work/project/ladcol_010/benchmark_scrna_cnv_caller/data/input_NUGC4 \
    --sample=SRR10805150 --localcores 4 --localmem 32 --id=NUGC4_hg38_mapped