# ------------------------------------------------------------------------------
# Prepare the BCC data (analyse both samples separately)
# Current mode of action: analyse the samples separately, 
#                         focusing on pre-treatment (start with sample08)
#                         use non-immune, non-malignant cells as reference
#                         (as done by the original publication)
#
# Include both pre and post treatment cells (separate datasets)
#
# Extension: checking how much the performance drops with other reference datasets
# ------------------------------------------------------------------------------

library(data.table)
library(zellkonverter)
library(Matrix)

sample_name<-"su008"

#Load dataset
meta_data <- fread("data/BCC_processed_matrices/GSE123813_bcc_all_metadata.txt.gz")

count_matrix<-fread("data/BCC_processed_matrices/GSE123813_bcc_scRNA_counts.txt.gz")
gene_names<-count_matrix$V1
count_matrix$V1<-NULL
count_matrix<-as.matrix(count_matrix)
row.names(count_matrix)<-gene_names

#Check distributions
#table(meta_data$patient,meta_data$treatment)
#table(meta_data$cluster)

#Filter data for sample, treatment and cell type
meta_data<-meta_data[meta_data$patient==sample_name &
                       meta_data$treatment == "pre",]

#Save complete data for other reference types (downstream)
meta_data_complete<-meta_data
count_matrix_complete<-count_matrix

# ------------------------------------------------------------------------------
# Dataset with recommended reference (as used in the paper)
# ------------------------------------------------------------------------------

meta_data<-meta_data[meta_data$cluster %in% c("Endothelial","Melanocytes",
                                              "Myofibroblasts","Tumor_1","Tumor_2"),]

#Filter count data respectively
count_matrix<-count_matrix[,meta_data$cell.id]

#Reformat cell barcode to cause no problems with numbat later
meta_data$cell.id<-gsub(".*_","",meta_data$cell.id)
meta_data$cell.id<-paste0(meta_data$cell.id,"-1")
colnames(count_matrix)<-meta_data$cell.id

#Rename "Tumor" to match Casper pipeline later
meta_data$cluster[startsWith(meta_data$cluster,"Tumor")]<-"BCC08"

#Save results
write.table(count_matrix,
            file="data/BCC_processed_matrices/input_BCC08/count_matrix.txt",
            sep="\t",quote=FALSE,row.names=TRUE)

write.table(meta_data[,c("cell.id","cluster")], 
            file="data/BCC_processed_matrices/input_BCC08/sample_annotation.txt",
            sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)

df_refs<-data.frame(ref_groups=setdiff(unique(meta_data$cluster),"BCC08"))
write.table(df_refs, 
            file="data/BCC_processed_matrices/input_BCC08/ref_groups.txt",
            sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)

# ------------------------------------------------------------------------------
# Dataset using other reference cell types (as used in the paper)
# ------------------------------------------------------------------------------

output_name<-"BCC08_immune"

meta_data<-meta_data_complete[!meta_data_complete$cluster %in% 
                                c("CAFs","Endothelial","Melanocytes","Myofibroblasts"),]

#Filter count data respectively
count_matrix<-count_matrix_complete[,meta_data$cell.id]

#Reformat cell barcode to cause no problems with numbat later
shortened_cellid<-gsub(".*_","",meta_data$cell.id)
meta_data$cell.id<-ifelse(meta_data$sort == "CD45- CD3-",shortened_cellid,meta_data$cell.id)
meta_data$cell.id<-paste0(meta_data$cell.id,"-1")
colnames(count_matrix)<-meta_data$cell.id

#Rename "Tumor" to match Casper pipeline later
meta_data$cluster[startsWith(meta_data$cluster,"Tumor")]<-output_name

#Save results
write.table(count_matrix,
            file=paste0("data/BCC_processed_matrices/input_",output_name,"/count_matrix.txt"),
            sep="\t",quote=FALSE,row.names=TRUE)

write.table(meta_data[,c("cell.id","cluster")], 
            file=paste0("data/BCC_processed_matrices/input_",output_name,"/sample_annotation.txt"),
            sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)

df_refs<-data.frame(ref_groups=setdiff(unique(meta_data$cluster),output_name))
write.table(df_refs, 
            file=paste0("data/BCC_processed_matrices/input_",output_name,"/ref_groups.txt"),
            sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)

# ------------------------------------------------------------------------------
# Dataset using an external reference dataset (skin dataset from Tabular Sapiens)
# ------------------------------------------------------------------------------

output_name<-"BCC08_TS"

#Keep only the cancer cells
meta_data<-meta_data_complete[meta_data_complete$cluster %in% 
                                c("Tumor_1","Tumor_2"),]
count_matrix<-count_matrix_complete[,meta_data$cell.id]

#Reformat cell barcode to cause no problems with numbat later
meta_data$cell.id<-gsub(".*_","",meta_data$cell.id)
meta_data$cell.id<-paste0(meta_data$cell.id,"-1")
colnames(count_matrix)<-meta_data$cell.id

#Add the healthy skin cells from the Tabular sapiens dataset
counts_ref<-readH5AD("data/Tabular_sapiens/TS_Skin.h5ad")

#Get technology and cell type annotation
ref_meta_data<-counts_ref@colData
table(ref_meta_data$method)
table(ref_meta_data$cell_ontology_class,ref_meta_data$free_annotation)

#Get raw counts
raw_counts_ref<-counts_ref@assays@data$raw_counts

#Add the column and row names (missing in this matrix for some reason)
colnames(raw_counts_ref)<-colnames(counts_ref@assays@data$X)
rownames(raw_counts_ref)<-rownames(counts_ref@assays@data$X)

#Delete large HDF5 object
rm(counts_ref)

#Check which genes are part of both count matrices
intersect_genes<-intersect(rownames(raw_counts_ref),rownames(count_matrix))
print(paste("Intersect genes:",length(intersect_genes)))

raw_counts_ref<-raw_counts_ref[intersect_genes,]
count_matrix<-count_matrix[intersect_genes,]

#Filter for cells measured with 10X and selected cell types
selected_cts<-c("endothelial cells","epithelial cells",
                "melanocytes","muscle cells", "skeletal muscle cells",
                "smooth muscle cells","stromal cells")
raw_counts_ref<-raw_counts_ref[,ref_meta_data$method== "10X" &
                                 ref_meta_data$free_annotation %in% selected_cts]
ref_meta_data<-ref_meta_data[ref_meta_data$method== "10X" &
                               ref_meta_data$free_annotation %in% selected_cts,]

#Combine both matrices
counts_combined<-cbind(raw_counts_ref,count_matrix)

#Remove white space from cell type names
ref_meta_data$free_annotation<-gsub(" ","_",as.character(ref_meta_data$free_annotation))

#Generate a meta_data object for the combined matrix
meta_data_combined<-data.frame(barcode=colnames(counts_combined),
                               sample=c(as.character(ref_meta_data$free_annotation),
                                        rep(output_name,ncol(count_matrix))))

#Save result data frames
write.table(as.matrix(counts_combined),
            file=paste0("data/BCC_processed_matrices/input_",output_name,"/count_matrix.txt"),
            sep="\t",quote=FALSE,row.names=TRUE)

write.table(meta_data_combined, 
            file=paste0("data/BCC_processed_matrices/input_",output_name,"/sample_annotation.txt"),
            sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)

#Save the reference groupfs
df_refs<-data.frame(ref_groups=as.character(unique(ref_meta_data$free_annotation)))
write.table(df_refs, 
            file=paste0("data/BCC_processed_matrices/input_",output_name,"/ref_groups.txt"),
            sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)


# ------------------------------------------------------------------------------
# Dataset using an external reference dataset (SNU601)
# ------------------------------------------------------------------------------

output_name<-"BCC08_SNU601"

#Keep only the cancer cells
meta_data<-meta_data_complete[meta_data_complete$cluster %in% 
                                c("Tumor_1","Tumor_2"),]
count_matrix<-count_matrix_complete[,meta_data$cell.id]

#Reformat cell barcode to cause no problems with numbat later
meta_data$cell.id<-gsub(".*_","",meta_data$cell.id)
meta_data$cell.id<-paste0(meta_data$cell.id,"-1")
colnames(count_matrix)<-meta_data$cell.id

#Load SNU601 data-matrix
counts_SNU601<-readMM("data/SNU601_processed_matrices/GSM4238687_SNU-601_matrix.mtx.gz")
genenames<-fread("data/SNU601_processed_matrices/GSM4238687_SNU-601_genes.tsv.gz",
                 header=FALSE)
cells<-fread("data/SNU601_processed_matrices/GSM4238687_SNU-601_barcodes.tsv.gz",
             header=FALSE)

counts_SNU601<-as.matrix(counts_SNU601)
colnames(counts_SNU601)<-paste0("SNU601_",cells$V1)
rownames(counts_SNU601)<-genenames$V1

#Check which genes are part of both count matrices
intersect_genes<-intersect(rownames(count_matrix),rownames(counts_SNU601))
print(paste("Intersect genes:",length(intersect_genes)))

count_matrix<-count_matrix[intersect_genes,]
counts_SNU601<-counts_SNU601[intersect_genes,]

#Define output directory for the files
output_dir<-paste0("data/BCC_processed_matrices/input_",output_name,"/")
dir.create(output_dir, showWarnings = FALSE)

#Combine both matrices
counts_combined<-cbind(count_matrix,counts_SNU601)

#Generate a meta_data object for the combined matrix
meta_data_combined<-data.frame(barcode=colnames(counts_combined),
                               sample=c(rep(output_name,ncol(count_matrix)),
                                        rep("SNU601",ncol(counts_SNU601))))

#Save result data frames
write.table(as.matrix(counts_combined),
            file=paste0(output_dir,"count_matrix.txt"),
            sep="\t",quote=FALSE,row.names=TRUE)

write.table(meta_data_combined, 
            file=paste0(output_dir,"sample_annotation.txt"),
            sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)

#Save the reference groupfs
df_refs<-data.frame(ref_groups=c("SNU601"))
write.table(df_refs, 
            file=paste0(output_dir,"ref_groups.txt"),
            sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)

# ------------------------------------------------------------------------------
# BCC08 - post treatment
# ------------------------------------------------------------------------------

#Load dataset
meta_data <- fread("data/BCC_processed_matrices/GSE123813_bcc_all_metadata.txt.gz")

count_matrix<-fread("data/BCC_processed_matrices/GSE123813_bcc_scRNA_counts.txt.gz")
gene_names<-count_matrix$V1
count_matrix$V1<-NULL
count_matrix<-as.matrix(count_matrix)
row.names(count_matrix)<-gene_names

#Check distributions
#table(meta_data$patient,meta_data$treatment)
#table(meta_data$cluster)

#Filter data for sample, treatment and cell type
meta_data<-meta_data[meta_data$patient==sample_name &
                       meta_data$treatment == "post",]

meta_data<-meta_data[meta_data$cluster %in% c("Endothelial","Melanocytes",
                                              "Myofibroblasts","Tumor_1","Tumor_2"),]

#Filter count data respectively
count_matrix<-count_matrix[,meta_data$cell.id]

#Reformat cell barcode to cause no problems with numbat later
meta_data$cell.id<-gsub(".*_","",meta_data$cell.id)
meta_data$cell.id<-paste0(meta_data$cell.id,"-1")
colnames(count_matrix)<-meta_data$cell.id

dataset_name<-"BCC08post"

#Rename "Tumor" to match Casper pipeline later
meta_data$cluster[startsWith(meta_data$cluster,"Tumor")]<-dataset_name

#Save results
write.table(count_matrix,
            file=paste0("data/BCC_processed_matrices/input_",
                        dataset_name,"/count_matrix.txt"),
            sep="\t",quote=FALSE,row.names=TRUE)

write.table(meta_data[,c("cell.id","cluster")], 
            file=paste0("data/BCC_processed_matrices/input_",
                        dataset_name,"/sample_annotation.txt"),
            sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)

df_refs<-data.frame(ref_groups=setdiff(unique(meta_data$cluster),dataset_name))
write.table(df_refs, 
            file=paste0("data/BCC_processed_matrices/input_",
                        dataset_name,"/ref_groups.txt"),
            sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)

