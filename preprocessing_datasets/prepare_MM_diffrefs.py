##################################################################################
# Combine MM dataset with healthy PBMCs (different cell types:
# T cells / B cells / Monocytes )
##################################################################################

import muon as mu
import pandas as pd
import os
import random
import scanpy as sc

#Set the working directory
os.chdir('/Users/kschmid/Documents/CNV_RNAseq_benchmark')

#Load MM dataset
path_MM='data/MM_sciCNV/input_MM/count_matrix.txt'
path_MM_samples='data/MM_sciCNV/input_MM/sample_annotation.txt'
counts_MM = pd.read_csv(path_MM,sep="\t")
annot_MM = pd.read_csv(path_MM_samples,sep="\t",header=None)
annot_MM.columns=["barcode","type"]

#Filter count matrix for cancer cells
cancer_cells = list(annot_MM.barcode[annot_MM.type=="MM"])
counts_MM = counts_MM[cancer_cells]

#File paths
path_pbmc_annotated='data/pbmc_healthy_multiome/pbmc_10k_RNA_annotated.h5mu'
path_pbmc_raw='data/pbmc_healthy_multiome/pbmc_unsorted_10k_filtered_feature_bc_matrix.h5'

#Read input files (annotated cells)
pbmcs = mu.read(path_pbmc_annotated)

#Extract the RNA part
pbmc_rna = pbmcs.mod['rna']

#Check cell type distribution
pbmc_rna.obs.celltype.value_counts()

#Load raw count matrix (as raw counts are required for CNV methods later)
raw_counts = mu.read_10x_h5(path_pbmc_raw)
raw_counts.var_names_make_unique()

#Filter genes (taking gene filter from intial annotation: expressed in at least 5 cells)
raw_counts_rna = raw_counts.mod['rna']
sc.pp.calculate_qc_metrics(raw_counts_rna, percent_top=None, log1p=False, inplace=True)
mu.pp.filter_var(raw_counts_rna, 'n_cells_by_counts', lambda x: x >= 5)

#Match the genes between MM and PBMC dataset
common_genes = set(raw_counts_rna.var.index).intersection(list(counts_MM.index))
counts_MM=counts_MM.filter(items=common_genes,axis=0)
raw_counts_rna=raw_counts_rna[:,list(common_genes)]

##################################################################################
# Select T cells as ref
##################################################################################

#Combine count matrices and save them
pbmc_cell_subset = list(pbmc_rna.obs.index[pbmc_rna.obs.celltype == "CD4+ naïve T"])
raw_counts_subset=raw_counts_rna[pbmc_cell_subset]
counts_pbmc = pd.DataFrame.sparse.from_spmatrix(raw_counts_subset.X.transpose(),
            index=raw_counts_subset.var.index,
            columns=["pbmc_"+bp for bp in raw_counts_subset.obs.index])
combined_mat = pd.concat([counts_MM,counts_pbmc], axis=1)
combined_mat.to_csv('data/MM_sciCNV/input_MM_tcell/count_matrix.txt',sep="\t")

#Create a sample annotation file
df = pd.DataFrame(data={'cells':list(combined_mat.columns),
                        'group': len(counts_MM.columns)*['MM_tcell']+len(pbmc_cell_subset)*['Tcell']})
df.to_csv('data/MM_sciCNV/input_MM_tcell/sample_annotation.txt',sep="\t",
          header=False,index=False)

#Create a ref file
file=open('data/MM_sciCNV/input_MM_tcell/ref_groups.txt','w')
file.write('ref_groups\nTcell')
file.close()

##################################################################################
# Select B cells as ref
##################################################################################

#Combine count matrices and save them
pbmc_cell_subset = list(pbmc_rna.obs.index[pbmc_rna.obs.celltype.isin(["memory B","naïve B"])])
raw_counts_subset=raw_counts_rna[pbmc_cell_subset]
counts_pbmc = pd.DataFrame.sparse.from_spmatrix(raw_counts_subset.X.transpose(),
            index=raw_counts_subset.var.index,
            columns=["pbmc_"+bp for bp in raw_counts_subset.obs.index])
combined_mat = pd.concat([counts_MM,counts_pbmc], axis=1)
combined_mat.to_csv('data/MM_sciCNV/input_MM_bcell/count_matrix.txt',sep="\t")

#Create a sample annotation file
df = pd.DataFrame(data={'cells':list(combined_mat.columns),
                        'group': len(counts_MM.columns)*['MM_bcell']+len(pbmc_cell_subset)*['Bcell']})
df.to_csv('data/MM_sciCNV/input_MM_bcell/sample_annotation.txt',sep="\t",
          header=False,index=False)

#Create a ref file
file=open('data/MM_sciCNV/input_MM_bcell/ref_groups.txt','w')
file.write('ref_groups\nBcell')
file.close()

##################################################################################
# Select CD14 Monocytes cells as ref
##################################################################################

#Combine count matrices and save them
pbmc_cell_subset = list(pbmc_rna.obs.index[pbmc_rna.obs.celltype == "CD14 mono"])
raw_counts_subset=raw_counts_rna[pbmc_cell_subset]
counts_pbmc = pd.DataFrame.sparse.from_spmatrix(raw_counts_subset.X.transpose(),
            index=raw_counts_subset.var.index,
            columns=["pbmc_"+bp for bp in raw_counts_subset.obs.index])
combined_mat = pd.concat([counts_MM,counts_pbmc], axis=1)
combined_mat.to_csv('data/MM_sciCNV/input_MM_mono/count_matrix.txt',sep="\t")

#Create a sample annotation file
df = pd.DataFrame(data={'cells':list(combined_mat.columns),
                        'group': len(counts_MM.columns)*['MM_mono']+len(pbmc_cell_subset)*['Mono']})
df.to_csv('data/MM_sciCNV/input_MM_mono/sample_annotation.txt',sep="\t",
          header=False,index=False)

#Create a ref file
file=open('data/MM_sciCNV/input_MM_mono/ref_groups.txt','w')
file.write('ref_groups\nMono')
file.close()