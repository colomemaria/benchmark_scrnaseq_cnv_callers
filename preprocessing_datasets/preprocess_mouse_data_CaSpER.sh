#!/bin/bash
#SBATCH --job-name=mouse_pileup
#SBATCH --mem-per-cpu=32G
#SBATCH --output=logs/mouse_pileup_%j.out
#SBATCH --error=logs/mouse_pileup_%j.err
#SBATCH -c 1
#SBATCH --time=48:00:00

#Stop script once an error occurs
set -e

#Download mouse genome (first tested mm10, then mm39)
#wget -c http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/chromFa.tar.gz
#wget -c http://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.chromFa.tar.gz

#Generate the pileup
tar -xvzf mm39.chromFa.tar.gz
mkdir ../mm39
FILES=./*fa
for f in $FILES
do
  echo "Processing $f file..."
  tools/BAFExtract/bin/BAFExtract -preprocess_FASTA $f ../mm39
done

#To get chromosome sizes
#fetchChromSizes mm10 > mm10.list
