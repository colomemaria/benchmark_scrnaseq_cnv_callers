#!/bin/bash
#SBATCH --error=%J.err
#SBATCH --output=%J.out
#SBATCH -J cellranger_pbmcext
#SBATCH -p slim16
#SBATCH -c 4 #number of CPUs needed 
#SBATCH --mem=128G 
#SBATCH -t 48:00:00

/work/project/ladcol_005/tools/cellranger-7.0.0/bin/cellranger count \
    --transcriptome=/work/project/ladcol_005/genomes/refdata-gex-GRCh38-2020-A \
    --fastqs=/home/kschmid/benchmark_scrna_cnv_caller/snakemake_pipeline/data/pbmc_ext_reference/A.merged/gemgroup001 \
    --sample=bamtofastq --localcores 4 --localmem 32 --id=pbmc_extref_hg38_mapped
