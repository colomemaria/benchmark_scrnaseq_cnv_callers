#!/bin/bash
#SBATCH --error=%J.err
#SBATCH --output=%J.out
#SBATCH -J bamtofastq_iAMP21
#SBATCH -p slim16
#SBATCH -c 4 #number of CPUs needed 
#SBATCH --mem-per-cpu=32G 
#SBATCH -t 48:00:00

/store24/project24/ladcol_005/tools/cellranger-7.0.0/lib/bin/bamtofastq \
    /work/project/ladcol_010/benchmark_scrna_cnv_caller/data/input_iAMP21/raw_data/SJBALL021901_D1.bam \
    /work/project/ladcol_010/benchmark_scrna_cnv_caller/data/input_iAMP21/raw_data/fastq

