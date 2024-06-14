#Download more gastric cancer cellines from SRA
module load ngs/sratoolkit/2.10.0

#Cell line  NUGC-4 (SRR10805150)
#fastq-dump --split-files --origfmt --gzip SRR10805150

#Cell line SNU-638 (SRR10805154)
#fastq-dump --split-files --origfmt --gzip SRR10805154

#Cell line  HGC-27 (SRR10805145)
fastq-dump --split-files --origfmt --gzip SRR10805145

#Cell line SNU-16 (SRR10805151)
fastq-dump --split-files --origfmt --gzip SRR10805151

#Cell line SNU-668 (SRR10805155)
fastq-dump --split-files --origfmt --gzip SRR10805155

#Renaming files for cell ranger
#mv SRR10805150_1.fastq.gz SRR10805150_S1_L001_R1_001.fastq.gz
#mv SRR10805150_2.fastq.gz SRR10805150_S1_L001_R2_001.fastq.gz

#mv SRR10805154_1.fastq.gz SRR10805154_S1_L001_R1_001.fastq.gz
#mv SRR10805154_2.fastq.gz SRR10805154_S1_L001_R2_001.fastq.gz

mv SRR10805145_1.fastq.gz SRR10805145_S1_L001_R1_001.fastq.gz
mv SRR10805145_2.fastq.gz SRR10805145_S1_L001_R2_001.fastq.gz

mv SRR10805151_1.fastq.gz SRR10805151_S1_L001_R1_001.fastq.gz
mv SRR10805151_2.fastq.gz SRR10805151_S1_L001_R2_001.fastq.gz

mv SRR10805155_1.fastq.gz SRR10805155_S1_L001_R1_001.fastq.gz
mv SRR10805155_2.fastq.gz SRR10805155_S1_L001_R2_001.fastq.gz