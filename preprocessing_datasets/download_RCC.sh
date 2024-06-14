#Download RCC sample
module load ngs/sratoolkit/2.10.0
module load ngs/samtools/1.9

#Sample 86
#fastq-dump --split-files --origfmt --gzip SRR19987215

#Download bam files directly (only one fastq file available)
sam-dump SRR19987215 | samtools view -bS - > possorted_genome_bam.bam