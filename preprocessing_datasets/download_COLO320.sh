#Data downloaded from SRA
module load ngs/sratoolkit/2.10.0
fastq-dump --split-files --origfmt --gzip SRR12900668

#Renaming files to run cellranger
mv SRR12900668_1.fastq.gz SRR12900668_S14_L001_R1_001.fastq.gz
mv SRR12900668_2.fastq.gz SRR12900668_S14_L001_R2_001.fastq.gz
