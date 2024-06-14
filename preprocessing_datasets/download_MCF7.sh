#Data downloaded from SRA
module load ngs/sratoolkit/2.10.0
fastq-dump --split-files --origfmt --gzip SRR10018149

#Renaming files to run cellranger
mv SRR10018149_1.fastq.gz SRR10018149_S1_L001_R1_001.fastq.gz
mv SRR10018149_2.fastq.gz SRR10018149_S1_L001_R2_001.fastq.gz