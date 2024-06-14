#!/bin/bash
#SBATCH --error=bamtofastq_%J.err
#SBATCH --output=bamtofastq_%J.out
#SBATCH -J bamtofastq
#SBATCH -p slim16
#SBATCH -c 1 #number of CPUs needed 
#SBATCH --mem-per-cpu=64G 

# Download fastq files of Kang et al, 2018 sample A

# Problem that fastq file on SRA is corrupt
#module load ngs/sratoolkit/2.10.0
#fastq-dump --split-files --origfmt --gzip SRR5398235

# Download instead from EBI the bam file and convert it back
#wget ftp://ftp.sra.ebi.ac.uk/vol1/SRA550/SRA550660/bam/A.merged.bam

#Convert bam back to fastq files
/work/project/ladcol_005/tools/cellranger-7.0.0/lib/bin/bamtofastq --cr11 data/pbmc_ext_reference/A.merged.bam data/pbmc_ext_reference/A.merged
