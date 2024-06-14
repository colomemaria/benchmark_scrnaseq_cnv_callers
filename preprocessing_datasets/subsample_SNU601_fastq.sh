#!/bin/bash
#SBATCH --error=subsampleSNU601_%J.err
#SBATCH --output=subsampleSNU601_%J.out
#SBATCH -J subsampleSNU601
#SBATCH -p slim16
#SBATCH -c 1 #number of CPUs needed 
#SBATCH --mem-per-cpu=32G 

#Load conda environment
source /home/kschmid/miniconda3/etc/profile.d/conda.sh
conda activate diverse

# #Need to unzip fastq files for that
#gzip -d data/input_SNU601/SRR10805153_S1_L001_R1_001.fastq.gz
#gzip -d data/input_SNU601/SRR10805153_S1_L001_R2_001.fastq.gz

#Subsample SNU601 reads to 80/60/40/20% (40% was already run before)
for cutoff in {80,60,20}
do

echo $cutoff

mkdir data/input_SNU601_sample$cutoff

fastq-sample -p 0.$cutoff -s 1 -o data/input_SNU601_sample$cutoff"/SRR10805153_sample"$cutoff"_S1_L001_R1_001" data/input_SNU601/SRR10805153_S1_L001_R1_001.fastq
fastq-sample -p 0.$cutoff -s 1 -o data/input_SNU601_sample$cutoff"/SRR10805153_sample"$cutoff"_S1_L001_R2_001" data/input_SNU601/SRR10805153_S1_L001_R2_001.fastq

#Zip resulting fastq files
gzip data/input_SNU601_sample$cutoff"/SRR10805153_sample"$cutoff"_S1_L001_R1_001.fastq"
gzip data/input_SNU601_sample$cutoff"/SRR10805153_sample"$cutoff"_S1_L001_R2_001.fastq"

done



