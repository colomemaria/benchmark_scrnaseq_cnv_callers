#Data downloaded from SRA
module load ngs/sratoolkit/2.10.0


fastq-dump --split-files --origfmt --gzip SRR10574405

#Renaming files to run cellranger
mv SRR10574405_1.fastq.gz RTBM199T_S14_L001_R1_001.fastq.gz
mv SRR10574405_2.fastq.gz RTBM199T_S14_L001_R2_001.fastq.gz
mv SRR10574405_3.fastq.gz RTBM199T_S14_L001_I1_001.fastq.gz

fastq-dump --split-files --origfmt --gzip SRR10574406

#Renaming files to run cellranger
mv SRR10574406_1.fastq.gz RTBM199T_S14_L002_R1_001.fastq.gz
mv SRR10574406_2.fastq.gz RTBM199T_S14_L002_R2_001.fastq.gz
mv SRR10574406_3.fastq.gz RTBM199T_S14_L002_I1_001.fastq.gz

fastq-dump --split-files --origfmt --gzip SRR10574407

#Renaming files to run cellranger
mv SRR10574407_1.fastq.gz RTBM199T_S14_L003_R1_001.fastq.gz
mv SRR10574407_2.fastq.gz RTBM199T_S14_L003_R2_001.fastq.gz
mv SRR10574407_3.fastq.gz RTBM199T_S14_L003_I1_001.fastq.gz

fastq-dump --split-files --origfmt --gzip SRR10574408

#Renaming files to run cellranger
mv SRR10574408_1.fastq.gz RTBM199T_S14_L004_R1_001.fastq.gz
mv SRR10574408_2.fastq.gz RTBM199T_S14_L004_R2_001.fastq.gz
mv SRR10574408_3.fastq.gz RTBM199T_S14_L004_I1_001.fastq.gz

fastq-dump --split-files --origfmt --gzip SRR10574409

#Renaming files to run cellranger
mv SRR10574409_1.fastq.gz RTBM199T_S14_L005_R1_001.fastq.gz
mv SRR10574409_2.fastq.gz RTBM199T_S14_L005_R2_001.fastq.gz
mv SRR10574409_3.fastq.gz RTBM199T_S14_L005_I1_001.fastq.gz

fastq-dump --split-files --origfmt --gzip SRR10574410

#Renaming files to run cellranger
mv SRR10574410_1.fastq.gz RTBM199T_S14_L006_R1_001.fastq.gz
mv SRR10574410_2.fastq.gz RTBM199T_S14_L006_R2_001.fastq.gz
mv SRR10574410_3.fastq.gz RTBM199T_S14_L006_I1_001.fastq.gz

fastq-dump --split-files --origfmt --gzip SRR10574411

#Renaming files to run cellranger
mv SRR10574411_1.fastq.gz RTBM199T_S14_L007_R1_001.fastq.gz
mv SRR10574411_2.fastq.gz RTBM199T_S14_L007_R2_001.fastq.gz
mv SRR10574411_3.fastq.gz RTBM199T_S14_L007_I1_001.fastq.gz

fastq-dump --split-files --origfmt --gzip SRR10574412

#Renaming files to run cellranger
mv SRR10574412_1.fastq.gz RTBM199T_S14_L008_R1_001.fastq.gz
mv SRR10574412_2.fastq.gz RTBM199T_S14_L008_R2_001.fastq.gz
mv SRR10574412_3.fastq.gz RTBM199T_S14_L008_I1_001.fastq.gz