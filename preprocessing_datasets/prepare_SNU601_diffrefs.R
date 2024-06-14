# ------------------------------------------------------------------------------
# Test SNU601 with different reference datasets
# - 
# - the MM dataset (only cancer cells)
#
# Other dataset tested and the later discarded due to 
# - different cell types from the gastric cell atlas
# https://www.cell.com/cell-reports/pdf/S2211-1247(23)01248-2.pdf
# ------------------------------------------------------------------------------

library(Matrix)
library(data.table)

#Load SNU601 matrix we get with cellranger
data_path<-"data/SNU601_processed_matrices/matrix_cellranger70"
counts_SNU601<-readMM(paste0(data_path,
                         "/filtered_feature_bc_matrix/matrix.mtx.gz"))
genenames<-fread(paste0(data_path,
                        "/filtered_feature_bc_matrix/features.tsv.gz"),
                 header=FALSE)
cells<-fread(paste0(data_path,
                    "/filtered_feature_bc_matrix/barcodes.tsv.gz"),
             header=FALSE)
colnames(counts_SNU601)<-cells$V1
rownames(counts_SNU601)<-genenames$V2

# ------------------------------------------------------------------------------
# Gastric cell atlas - load and categorize all cells
# ------------------------------------------------------------------------------

data_path<-"data/gastric_cell_atlas_II/"

counts<-fread(paste0(data_path,"GSE159929_RAW/GSM4850590_Stomach_Counts.csv.gz"))
genes<-counts$V1
counts$V1<-NULL
counts<-as.matrix(counts)
rownames(counts)<-genes

meta_data<-fread(paste0(data_path,
                        "scRNA-AHCA/Cell_barcode_and_corresponding_cell_types_of_AHCA/",
                        "Annotation_AHCA_alltissues_meta.data_84363_cell.txt"))
meta_data<-meta_data[meta_data$orig.ident == "Stomach_cDNA",]
meta_data$barcode<-gsub("Stomach_cDNA_","",meta_data$V1)

group_celltypes<-data.frame(ct_detailed=c("Absorptive Cell","B Cell CD79A",
                                          "Endothelial Cell ACKR1","Epithelial Cell TFF3",
                                          "Fibroblast C7","Fibroblast PLA2G2A","Fibroblast PTGDS",
                                          "FibSmo Cell","High Proliferation Erythrocyte",
                                          "Macrophage C1QB","Macrophage FCN3","Monocyte",
                                          "Mucosal Epithelial Cell","NK/T Cell GNLY","Plasma Cell JCHAIN",
                                          "Secretory Cell","Smooth Muscle Cell",
                                          "T Cell CCL5","T Cell GZMA","T Cell GZMK",
                                          "T Cell IL7R","T Cell XCL1","Tuft Cell"),
                            ct_merged=c("epithelial","immune",
                                        "endothelial","epithelial",
                                        "fibroblast","fibroblast","fibroblast",
                                        "fibsmo","erythrocyte",
                                        "immune","immune","immune",
                                        "epithelial","immune","immune",
                                        "epithelial","smuscle",
                                        "immune","immune","immune",
                                        "immune","immune","epithelial"))

meta_data<-merge(meta_data,group_celltypes,by.x="Cell_type_in_merged_data",
                 by.y="ct_detailed",all.x=TRUE)
meta_data<-meta_data[order(meta_data$barcode),]

#Check the metadata and count matrix are ordered the same way
all(meta_data$barcode == colnames(counts))

#Filter count matrices for common genes
common_genes<-intersect(rownames(counts),rownames(counts_SNU601))
print(paste0("Common genes to the gastric dataset:",length(common_genes)))

counts<-counts[common_genes,]
counts_SNU601_filtered<-counts_SNU601[common_genes,]

# ------------------------------------------------------------------------------
# Gastric cell atlas - epi & endothelial cells
# ------------------------------------------------------------------------------

#Define output directory for the files
dataset_name<-"SNU601_epiendo"
print(dataset_name)
output_dir<-paste0("data/SNU601_processed_matrices/input_",dataset_name,"/")
dir.create(output_dir, showWarnings = FALSE)
  
#Filter gastric cancer dataset
meta_data_subset<-meta_data[meta_data$ct_merged %in% c("epithelial","endothelial"),]
counts_subset<-counts[,meta_data_subset$barcode]

#Combine both matrices
counts_combined<-cbind(counts_SNU601_filtered,counts_subset)

#Generate a meta_data object for the combined matrix
meta_data_combined<-data.frame(barcode=colnames(counts_combined),
                               sample=c(rep(dataset_name,ncol(counts_SNU601_filtered)),
                                        meta_data_subset$ct_merged))

#Save result data frames
write.table(as.matrix(counts_combined),
            file=paste0(output_dir,"count_matrix.txt"),
            sep="\t",quote=FALSE,row.names=TRUE)

write.table(meta_data_combined, 
            file=paste0(output_dir,"sample_annotation.txt"),
            sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)

#Save the reference groupfs
df_refs<-data.frame(ref_groups=unique(meta_data_subset$ct_merged))
write.table(df_refs, 
            file=paste0(output_dir,"ref_groups.txt"),
            sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)

# ------------------------------------------------------------------------------
# Gastric cell atlas - fibroblasts & smooth muscle cells
# ------------------------------------------------------------------------------

#Define output directory for the files
dataset_name<-"SNU601_fibsom"
print(dataset_name)
output_dir<-paste0("data/SNU601_processed_matrices/input_",dataset_name,"/")
dir.create(output_dir, showWarnings = FALSE)

#Filter gastric cancer dataset
meta_data_subset<-meta_data[meta_data$ct_merged %in% c("fibroblast","fibsmo","smuscle"),]
counts_subset<-counts[,meta_data_subset$barcode]

#Combine both matrices
counts_combined<-cbind(counts_SNU601_filtered,counts_subset)

#Generate a meta_data object for the combined matrix
meta_data_combined<-data.frame(barcode=colnames(counts_combined),
                               sample=c(rep(dataset_name,ncol(counts_SNU601_filtered)),
                                        meta_data_subset$ct_merged))

#Save result data frames
write.table(as.matrix(counts_combined),
            file=paste0(output_dir,"count_matrix.txt"),
            sep="\t",quote=FALSE,row.names=TRUE)

write.table(meta_data_combined, 
            file=paste0(output_dir,"sample_annotation.txt"),
            sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)

#Save the reference groupfs
df_refs<-data.frame(ref_groups=unique(meta_data_subset$ct_merged))
write.table(df_refs, 
            file=paste0(output_dir,"ref_groups.txt"),
            sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)

# ------------------------------------------------------------------------------
# Gastric cell atlas - immune cells
# ------------------------------------------------------------------------------

#Define output directory for the files
dataset_name<-"SNU601_immune"
print(dataset_name)
output_dir<-paste0("data/SNU601_processed_matrices/input_",dataset_name,"/")
dir.create(output_dir, showWarnings = FALSE)

#Filter gastric cancer dataset
meta_data_subset<-meta_data[meta_data$ct_merged %in% c("immune"),]
counts_subset<-counts[,meta_data_subset$barcode]

#Combine both matrices
counts_combined<-cbind(counts_SNU601_filtered,counts_subset)

#Generate a meta_data object for the combined matrix
meta_data_combined<-data.frame(barcode=colnames(counts_combined),
                               sample=c(rep(dataset_name,ncol(counts_SNU601_filtered)),
                                        meta_data_subset$ct_merged))

#Save result data frames
write.table(as.matrix(counts_combined),
            file=paste0(output_dir,"count_matrix.txt"),
            sep="\t",quote=FALSE,row.names=TRUE)

write.table(meta_data_combined, 
            file=paste0(output_dir,"sample_annotation.txt"),
            sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)

#Save the reference groupfs
df_refs<-data.frame(ref_groups=unique(meta_data_subset$ct_merged))
write.table(df_refs, 
            file=paste0(output_dir,"ref_groups.txt"),
            sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)

# ------------------------------------------------------------------------------
# MM dataset
# ------------------------------------------------------------------------------

#Load the original MM dataset
count_mat<-fread("data/MM_sciCNV/input_MM/count_matrix.txt")
genes<-count_mat$V1
count_mat$V1<-NULL
count_mat<-as.matrix(count_mat)
rownames(count_mat)<-genes

#Filter for the cancer cells
sample_annot<-fread("data/MM_sciCNV/input_MM/sample_annotation.txt",
                    header=FALSE)
count_mat<-count_mat[,sample_annot$V2 == "MM"]

#Filter count matrices for common genes
common_genes<-intersect(rownames(count_mat),rownames(counts_SNU601))
print(paste0("Common genes to the gastric dataset:",length(common_genes)))

count_mat<-count_mat[common_genes,]
counts_SNU601_filtered<-counts_SNU601[common_genes,]

#Define output directory for the files
dataset_name<-"SNU601_MM"
print(dataset_name)
output_dir<-paste0("data/SNU601_processed_matrices/input_",dataset_name,"/")
dir.create(output_dir, showWarnings = FALSE)

#Combine both matrices
counts_combined<-cbind(counts_SNU601_filtered,count_mat)

#Generate a meta_data object for the combined matrix
meta_data_combined<-data.frame(barcode=colnames(counts_combined),
                               sample=c(rep(dataset_name,ncol(counts_SNU601_filtered)),
                                        rep("MM",ncol(count_mat))))

#Save result data frames
write.table(as.matrix(counts_combined),
            file=paste0(output_dir,"count_matrix.txt"),
            sep="\t",quote=FALSE,row.names=TRUE)

write.table(meta_data_combined, 
            file=paste0(output_dir,"sample_annotation.txt"),
            sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)

#Save the reference groupfs
df_refs<-data.frame(ref_groups=unique(meta_data_subset$ct_merged))
write.table(df_refs, 
            file=paste0(output_dir,"ref_groups.txt"),
            sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)


