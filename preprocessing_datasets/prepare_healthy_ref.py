##################################################################################
# Preprocess dataset for CNV pipeline (testing a healthy diploid sample)
# Generating two versions:
# 1) once splitting the CD4+ naive T cells in half (using half as the reference)
# 2) once using instead monocytes (CD14 mono) as reference
#
# Dataset downloaded from 10X Genomics: 
# https://www.10xgenomics.com/resources/datasets/pbmc-from-a-healthy-donor-no-cell-sorting-10-k-1-standard-2-0-0 
##################################################################################

import muon as mu
import pandas as pd
import os
import random
import scanpy as sc

#Set the working directory
os.chdir('/Users/kschmid/Documents/CNV_RNAseq_benchmark')

#File paths
path_pbmc_annotated='data/pbmc_healthy_multiome/pbmc_10k_RNA_annotated.h5mu'
path_pbmc_raw='data/pbmc_healthy_multiome/pbmc_unsorted_10k_filtered_feature_bc_matrix.h5'

#Read input files (annotated cells)
pbmcs = mu.read(path_pbmc_annotated)

#Extract the RNA part
pbmc_rna = pbmcs.mod['rna']

#Check cell type distribution
pbmc_rna.obs.celltype.value_counts()

#Save the CD14 Monocytes for later (alternative reference dataset)
pbmc_mono = pbmc_rna[pbmc_rna.obs.celltype == "CD14 mono",]

#Filter for CD4+ naive T cells
pbmc_rna = pbmc_rna[pbmc_rna.obs.celltype == "CD4+ naïve T",]

#Load raw count matrix (as raw counts are required for CNV methods later)
raw_counts = mu.read_10x_h5(path_pbmc_raw)
raw_counts.var_names_make_unique()

#Filter genes (taking gene filter from intial annotation: expressed in at least 5 cells)
raw_counts_rna = raw_counts.mod['rna']
sc.pp.calculate_qc_metrics(raw_counts_rna, percent_top=None, log1p=False, inplace=True)
mu.pp.filter_var(raw_counts_rna, 'n_cells_by_counts', lambda x: x >= 5)

#Filter for the respective cells
raw_counts_mono = raw_counts_rna[pbmc_mono.obs.index,]
raw_counts_rna = raw_counts_rna[pbmc_rna.obs.index,]

##################### Save the first dataset with CD4 T cell reference ###########################

count_mat = pd.DataFrame.sparse.from_spmatrix(raw_counts_rna.X.transpose(),
            index=raw_counts_rna.var.index,
            columns=raw_counts_rna.obs.index)
count_mat.to_csv('data/pbmc_healthy_multiome/input_pbmc/count_matrix.txt',sep="\t")

#Randomly sample 50% of cells to be "cases" and 50% to be controls
random.seed(10)
cancer_cells = random.sample(list(pbmc_rna.obs.index),int(len(pbmc_rna.obs)/2))

#Generate sample annotation
grouping = ["pbmc" if cell in cancer_cells else "control" for cell in list(pbmc_rna.obs.index)]
df = pd.DataFrame(data={'cells':list(pbmc_rna.obs.index),
                        'group':grouping})
df.to_csv('data/pbmc_healthy_multiome/input_pbmc/sample_annotation.txt',sep="\t",
          header=False,index=False)

#Get ref file
file=open('data/pbmc_healthy_multiome/input_pbmc/ref_groups.txt','w')
file.write('ref_groups\ncontrol')
file.close()

#Get barcode file to filter bam file (with prefix as recommended by 10X Genomics)
file = open('data/pbmc_healthy_multiome/selected_test_cells.txt','w')
for cc in cancer_cells:
    file.write("CB:Z:"+ cc+"\n")
file.close()

#Get barcode file to filter bam file (without prefix for numbat)
file = open('data/pbmc_healthy_multiome/barcodes.tsv','w')
for cc in cancer_cells:
    file.write(cc+"\n")
file.close()

##################### Save the second dataset with CD14 Mono reference ###########################

#Filter for CD4 T cell (assigned cancer cell subset)
count_mat = count_mat[cancer_cells]
count_mono = pd.DataFrame.sparse.from_spmatrix(raw_counts_mono.X.transpose(),
            index=raw_counts_mono.var.index,
            columns=raw_counts_mono.obs.index)
count_combined = pd.concat([count_mat,count_mono],axis=1)
count_combined.to_csv('data/pbmc_healthy_multiome/input_pbmc_monoref/count_matrix.txt',sep="\t")

#Generate sample annotation
grouping = ["pbmc_monoref" if cell in cancer_cells else "control" for cell in list(count_combined.columns)]
df = pd.DataFrame(data={'cells':list(count_combined.columns),
                        'group':grouping})
df.to_csv('data/pbmc_healthy_multiome/input_pbmc_monoref/sample_annotation.txt',sep="\t",
          header=False,index=False)

#Get ref file
file=open('data/pbmc_healthy_multiome/input_pbmc_monoref/ref_groups.txt','w')
file.write('ref_groups\ncontrol')
file.close()
