#!/bin/bash
#SBATCH --error=%J.err
#SBATCH --output=%J.out
#SBATCH -J cellranger_iAMP21
#SBATCH -p slim16
#SBATCH -c 4 #number of CPUs needed 
#SBATCH --mem-per-cpu=32G 
#SBATCH -t 48:00:00

/store24/project24/ladcol_005/tools/cellranger-7.0.0/bin/cellranger count \
    --transcriptome=/store24/project24/ladcol_005/genomes/refdata-gex-GRCh38-2020-A \
    --fastqs=/work/project/ladcol_010/benchmark_scrna_cnv_caller/data/input_iAMP21/raw_data/fastq/SJBALL021901_D1_0_1_H2YWTDMXY \
    --sample=bamtofastq --localcores 4 --localmem 32 --id=iAMP21_hg38_mapped
