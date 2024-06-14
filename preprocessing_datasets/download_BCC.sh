#Download BCC from SRA (all samples pre treatment with at least 20 annotated cancer cells in GEO)
module load ngs/sratoolkit/2.10.0

#Sample 05 (tumor pre treatment)
#fastq-dump --split-files --origfmt --gzip SRR8315751

#Sample 06 (tumor pre treatment)
#fastq-dump --split-files --origfmt --gzip SRR8315756

#Sample 07 (tumor pre treatment)
#fastq-dump --split-files --origfmt --gzip SRR8315761

#Sample 08 (tumor pre treatment)
#fastq-dump --split-files --origfmt --gzip SRR8315765

#Renaming files for cell ranger
#mv SRR8315751_1.fastq.gz SRR8315751_S1_L001_R1_001.fastq.gz
#mv SRR8315751_2.fastq.gz SRR8315751_S1_L001_R2_001.fastq.gz

#mv SRR8315761_1.fastq.gz SRR8315761_S1_L001_R1_001.fastq.gz
#mv SRR8315761_2.fastq.gz SRR8315761_S1_L001_R2_001.fastq.gz


#Sample 06 (tumor post treatment)
fastq-dump --split-files --origfmt --gzip SRR8315759

#Sample 08 (tumor post treatment)
fastq-dump --split-files --origfmt --gzip SRR8315766

#Renaming files for cell ranger
mv SRR8315759_1.fastq.gz SRR8315759_S1_L001_R1_001.fastq.gz
mv SRR8315759_2.fastq.gz SRR8315759_S1_L001_R2_001.fastq.gz

mv SRR8315766_1.fastq.gz SRR8315766_S1_L001_R1_001.fastq.gz
mv SRR8315766_2.fastq.gz SRR8315766_S1_L001_R2_001.fastq.gz