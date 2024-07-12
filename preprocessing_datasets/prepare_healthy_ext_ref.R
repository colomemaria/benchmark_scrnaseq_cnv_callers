# ------------------------------------------------------------------------------
# Combine the PBMC dataset with an external reference,
# again CD4T cells and Monocytes (taking the first matrix, sample A from Kang et al. 2018)
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96583
# ------------------------------------------------------------------------------

library(data.table)
library(Matrix)
library(ggplot2)

# ------------------------------------------------------------------------------
# Pre-check: how much does the external count matrix differ dependent on the
# cellranger version (and genome version)
# ------------------------------------------------------------------------------

# #Old version
# matrix_ext<-readMM("data/pbmc_Kang/GSE96583_RAW/GSM2560245_A.mat.gz")
# barcodes<-fread("data/pbmc_Kang/GSE96583_RAW/GSM2560245_barcodes.tsv.gz",
#                 header=FALSE)
# colnames(matrix_ext)<-barcodes$V1
# genes<-fread("data/pbmc_Kang/GSE96583_genes.txt",
#              header=TRUE)
# genes$x<-gsub(".*\"\"","",genes$x)
# rownames(matrix_ext)<-genes$x
# 
# #New version
# matrix_ext_new<-readMM("data/pbmc_Kang/remapped_with_cr700/filtered_feature_bc_matrix/matrix.mtx.gz")
# barcodes<-fread("data/pbmc_Kang/remapped_with_cr700/filtered_feature_bc_matrix/barcodes.tsv.gz",
#                 header=FALSE)
# colnames(matrix_ext_new)<-barcodes$V1
# genes<-fread("data/pbmc_Kang/remapped_with_cr700/filtered_feature_bc_matrix/features.tsv.gz",
#              header=FALSE)
# rownames(matrix_ext_new)<-genes$V1
# 
# #Filter for the same genes and cells
# common_cells<-intersect(colnames(matrix_ext),colnames(matrix_ext_new))
# common_genes<-intersect(rownames(matrix_ext),rownames(matrix_ext_new))
# matrix_ext<-matrix_ext[common_genes,common_cells]
# matrix_ext_new<-matrix_ext_new[common_genes,common_cells]
# 
# #Compare the count differences
# matrix_ext_melted<-melt(as.matrix(matrix_ext))
# matrix_ext_new_melted<-melt(as.matrix(matrix_ext_new))
# colnames(matrix_ext_melted)<-c("gene","cell","count_c1")
# matrix_ext_melted$count_c7<-matrix_ext_new_melted$value
# rm(matrix_ext_new_melted)
# 
# #Remove points that are both 0
# matrix_ext_melted<-matrix_ext_melted[matrix_ext_melted$count_c1>0 &
#                                       matrix_ext_melted$count_c7>0,]
# 
# g<-ggplot(matrix_ext_melted,aes(x=count_c1,y=count_c7))+
#   geom_point()+geom_abline()+scale_x_log10()+scale_y_log10()+
#   theme_bw()
# ggsave(g,filename="~/Desktop/deviation_counts.png",
#        width=10,hight=6)
# 
# # Calculate deviation
# matrix_ext_melted$dev<-matrix_ext_melted$count_c1-matrix_ext_melted$count_c7
# ggplot(matrix_ext_melted,aes(x=dev))+
#   geom_histogram()+xlim(-5,5)

# ------------------------------------------------------------------------------
# Process dataset
# ------------------------------------------------------------------------------

#Load current PBMC dataset
count_matrix<-fread("data/pbmc_healthy_multiome/input_pbmc/count_matrix.txt")
genes<-count_matrix$V1
count_matrix$V1<-NULL
count_matrix<-as.matrix(count_matrix)
rownames(count_matrix)<-genes

#Annotation file
samples<-fread("data/pbmc_healthy_multiome/input_pbmc/sample_annotation.txt",
               header=FALSE)
samples<-samples[samples$V2=="pbmc",]
count_matrix<-count_matrix[,samples$V1]

#Load new reference
matrix_ext<-readMM("data/pbmc_Kang/GSE96583_RAW/GSM2560245_A.mat.gz")
barcodes<-fread("data/pbmc_Kang/GSE96583_RAW/GSM2560245_barcodes.tsv.gz",
                header=FALSE)
colnames(matrix_ext)<-barcodes$V1
genes<-fread("data/pbmc_Kang/GSE96583_genes.txt",
             header=TRUE)
genes$x<-gsub(".*\"\"","",genes$x)
rownames(matrix_ext)<-genes$x

#Get all cell annotations
cell_types<-fread("data/pbmc_Kang/GSE96583_batch1.total.tsne.df.tsv",
                header=FALSE)
cell_types<-cell_types[cell_types$V4=="A",]

#Remove doublets
matrix_ext<-matrix_ext[,cell_types$V1[cell_types$V6 == "singlet" &
                                        cell_types$V7 %in% c("CD4 T cells",
                                                             "CD14+ Monocytes")]]

#Merge genes to symbols
gene_map<-fread("data/pbmc_Kang/GSE96583_batch1.genes.tsv",
                header=FALSE)

#Check that sorting is the same as in the data frame
all(rownames(matrix_ext)==gene_map$V1)

matrix_reduced<-apply(matrix_ext, 2, tapply, as.factor(gene_map$V2),
           mean)

#Filter for common gene names
gene_names<-intersect(rownames(matrix_reduced),rownames(count_matrix))
print(paste("Common gene names:",length(gene_names)))
matrix_reduced<-matrix_reduced[gene_names,]
count_matrix<-count_matrix[gene_names,]

# ------------------------------------------------------------------------------
# Dataset 1 - CD4 T cells with ext CD4 T cells
# ------------------------------------------------------------------------------

#Define output directory for the files
sample_name<-"pbmc_ext"
output_dir<-paste0("data/pbmc_healthy_multiome/input_",sample_name,"/")
dir.create(output_dir, showWarnings = FALSE)

#Filter matrix for CD4 T cells
matrix_tcells<-matrix_reduced[,cell_types$V1[cell_types$V6 == "singlet" &
                                               cell_types$V7 == "CD4 T cells"]]

#Rename barcodes to not get duplicated barcode names
colnames(matrix_tcells)<-paste0("ext_",colnames(matrix_tcells))

#Combine both matrices
counts_combined<-cbind(count_matrix,matrix_tcells)

#Generate a meta_data object for the combined matrix
meta_data_combined<-data.frame(barcode=colnames(counts_combined),
                               sample=c(rep(sample_name,ncol(count_matrix)),
                                        rep("cd4t_ext",ncol(matrix_tcells))))

#Save result data frames
write.table(as.matrix(counts_combined),
            file=paste0(output_dir,"count_matrix.txt"),
            sep="\t",quote=FALSE,row.names=TRUE)

write.table(meta_data_combined, 
            file=paste0(output_dir,"sample_annotation.txt"),
            sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)

#Save the reference groupfs
df_refs<-data.frame(ref_groups=c("cd4t_ext"))
write.table(df_refs, 
            file=paste0(output_dir,"ref_groups.txt"),
            sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)

# ------------------------------------------------------------------------------
# Dataset 2 - CD4 T cells with ext Monocytes
# ------------------------------------------------------------------------------

#Define output directory for the files
sample_name<-"pbmc_ext_mono"
output_dir<-paste0("data/pbmc_healthy_multiome/input_",sample_name,"/")
dir.create(output_dir, showWarnings = FALSE)

#Filter matrix for CD4 T cells
matrix_mono<-matrix_reduced[,cell_types$V1[cell_types$V6 == "singlet" &
                                               cell_types$V7 == "CD14+ Monocytes"]]

#Rename barcodes to not get duplicated barcode names
colnames(matrix_mono)<-paste0("ext_",colnames(matrix_mono))

#Combine both matrices
counts_combined<-cbind(count_matrix,matrix_mono)

#Generate a meta_data object for the combined matrix
meta_data_combined<-data.frame(barcode=colnames(counts_combined),
                               sample=c(rep(sample_name,ncol(count_matrix)),
                                        rep("mono_ext",ncol(matrix_mono))))

#Save result data frames
write.table(as.matrix(counts_combined),
            file=paste0(output_dir,"count_matrix.txt"),
            sep="\t",quote=FALSE,row.names=TRUE)

write.table(meta_data_combined, 
            file=paste0(output_dir,"sample_annotation.txt"),
            sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)

#Save the reference groups
df_refs<-data.frame(ref_groups=c("mono_ext"))
write.table(df_refs, 
            file=paste0(output_dir,"ref_groups.txt"),
            sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)

# ------------------------------------------------------------------------------
# Dataset 3 - CD4T cells vs ext CD4 T cells (mapped with cellranger 7)
# ------------------------------------------------------------------------------

#Load current PBMC dataset
count_matrix<-fread("data/pbmc_healthy_multiome/input_pbmc/count_matrix.txt")
genes<-count_matrix$V1
count_matrix$V1<-NULL
count_matrix<-as.matrix(count_matrix)
rownames(count_matrix)<-genes

#Annotation file
samples<-fread("data/pbmc_healthy_multiome/input_pbmc/sample_annotation.txt",
               header=FALSE)
samples<-samples[samples$V2=="pbmc",]
count_matrix<-count_matrix[,samples$V1]

#New version of the second PBMC dataset (mapped with CR7)
matrix_ext_new<-readMM("data/pbmc_Kang/remapped_with_cr700/filtered_feature_bc_matrix/matrix.mtx.gz")
barcodes<-fread("data/pbmc_Kang/remapped_with_cr700/filtered_feature_bc_matrix/barcodes.tsv.gz",
                header=FALSE)
colnames(matrix_ext_new)<-barcodes$V1
genes<-fread("data/pbmc_Kang/remapped_with_cr700/filtered_feature_bc_matrix/features.tsv.gz",
             header=FALSE)
rownames(matrix_ext_new)<-genes$V1

#Get cell annotation (filter for only CD4T cell singlets)
cell_types<-fread("data/pbmc_Kang/GSE96583_batch1.total.tsne.df.tsv",
                  header=FALSE)
cell_types<-cell_types[cell_types$V4=="A" &
                         cell_types$V6 == "singlet" &
                         cell_types$V7 == "CD4 T cells",]
matrix_ext_new<-matrix_ext_new[,cell_types$V1]

#Merge genes to symbols
gene_map<-fread("data/pbmc_Kang/GSE96583_batch1.genes.tsv",
                header=FALSE)
gene_map<-as.data.frame(gene_map)
rownames(gene_map)<-gene_map$V1

matching_genes<-intersect(gene_map$V1,rownames(matrix_ext_new))
matrix_ext_new<-matrix_ext_new[matching_genes,]

#Check that sorting is the same as in the data frame
all(rownames(matrix_ext_new)==gene_map$V1)

matrix_reduced<-apply(matrix_ext_new, 2, tapply, as.factor(gene_map$V2),
                      mean)

#Filter for common gene names
gene_names<-intersect(rownames(matrix_reduced),rownames(count_matrix))
print(paste("Common gene names:",length(gene_names)))
matrix_reduced<-matrix_reduced[gene_names,]
count_matrix<-count_matrix[gene_names,]

#Define output directory for the files
sample_name<-"pbmc_ext_cr7"
output_dir<-paste0("data/pbmc_healthy_multiome/input_",sample_name,"/")
dir.create(output_dir, showWarnings = FALSE)

#Rename barcodes to not get duplicated barcode names
colnames(matrix_reduced)<-paste0("ext_",colnames(matrix_reduced))

#Combine both matrices
counts_combined<-cbind(count_matrix,matrix_reduced)

#Generate a meta_data object for the combined matrix
meta_data_combined<-data.frame(barcode=colnames(counts_combined),
                               sample=c(rep(sample_name,ncol(count_matrix)),
                                        rep("cd4t_ext",ncol(matrix_reduced))))

#Save result data frames
write.table(as.matrix(counts_combined),
            file=paste0(output_dir,"count_matrix.txt"),
            sep="\t",quote=FALSE,row.names=TRUE)

write.table(meta_data_combined, 
            file=paste0(output_dir,"sample_annotation.txt"),
            sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)

#Save the reference groups
df_refs<-data.frame(ref_groups=c("cd4t_ext"))
write.table(df_refs, 
            file=paste0(output_dir,"ref_groups.txt"),
            sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)
