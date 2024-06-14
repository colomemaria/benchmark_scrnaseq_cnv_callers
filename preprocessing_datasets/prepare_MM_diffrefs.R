# ------------------------------------------------------------------------------
# Test the performance of the MM dataset when using an external reference
# dataset
# 1) The SNU601 cell line (something completely off)
# 2) Gastric cells (reference for SNU601)
# 3) Additional a python script prepare_MM_diffrefs.py for different PBMC refs
# ------------------------------------------------------------------------------

library(data.table)
library(Matrix)

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

# ------------------------------------------------------------------------------
# SNU601 data matrix
# ------------------------------------------------------------------------------

#Load SNU601 data matrix
# counts_SNU601<-readMM("data/SNU601_processed_matrices/GSM4238687_SNU-601_matrix.mtx.gz")
# genenames<-fread("data/SNU601_processed_matrices/GSM4238687_SNU-601_genes.tsv.gz",
#                  header=FALSE)
# cells<-fread("data/SNU601_processed_matrices/GSM4238687_SNU-601_barcodes.tsv.gz",
#              header=FALSE)
# counts_SNU601<-as.matrix(counts_SNU601)
# colnames(counts_SNU601)<-paste0("SNU601_",cells$V1)
# rownames(counts_SNU601)<-genenames$V1

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
counts_SNU601<-as.matrix(counts_SNU601)
colnames(counts_SNU601)<-paste0("SNU601_",cells$V1)
rownames(counts_SNU601)<-genenames$V2

#Check which genes are part of both count matrices
intersect_genes<-intersect(rownames(count_mat),rownames(counts_SNU601))
print(paste("Intersect genes:",length(intersect_genes)))

count_mat_filtered<-count_mat[intersect_genes,]
counts_SNU601<-counts_SNU601[intersect_genes,]

#Define output directory for the files
output_dir<-"data/MM_sciCNV/input_MM_SNU601/"

#Combine both matrices
counts_combined<-cbind(count_mat_filtered,counts_SNU601)

#Generate a meta_data object for the combined matrix
meta_data_combined<-data.frame(barcode=colnames(counts_combined),
                               sample=c(rep("MM_SNU601",ncol(count_mat_filtered)),
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
# Gastric cels (reference for SNU601)
# ------------------------------------------------------------------------------

#Load reference dataset (focus on control samples 25-29)
patients<-25:29
count_ref_combined<-NULL
sample_vector<-NULL
for(pat in patients){
  
  suppressWarnings(counts_ref<-fread(paste0("data/gastic_cells_SNU601_reference/GSE150290_RAW/",
                           "GSM45463",(pat+22),"_Pat",pat,"-A.txt.gz")))
  gene_names_ref<-counts_ref$V1
  counts_ref$V1<-NULL
  counts_ref<-as.matrix(counts_ref)
  rownames(counts_ref)<-gene_names_ref
  
  if(is.null(count_ref_combined)){
    count_ref_combined<-counts_ref
    sample_vector<-rep(paste0("Patient",pat),ncol(counts_ref))
  } else{
    print(all(rownames(counts_ref)==rownames(count_ref_combined)))
    count_ref_combined<-cbind(count_ref_combined,counts_ref)
    sample_vector<-c(sample_vector,
                     rep(paste0("Patient",pat),ncol(counts_ref)))
  }
  
}

#Check which genes are part of both count matrices
intersect_genes<-intersect(rownames(count_mat),rownames(count_ref_combined))
print(paste("Intersect genes:",length(intersect_genes)))

count_mat_filtered<-count_mat[intersect_genes,]
count_ref_combined<-count_ref_combined[intersect_genes,]

#Define output directory for the files
output_dir<-"data/MM_sciCNV/input_MM_gastric/"
dir.create(output_dir, showWarnings = FALSE)

#Combine both matrices
counts_combined<-cbind(count_mat_filtered,count_ref_combined)

#Generate a meta_data object for the combined matrix
meta_data_combined<-data.frame(barcode=colnames(counts_combined),
                               sample=c(rep("MM_gastric",ncol(count_mat_filtered)),
                                        sample_vector))

#Save result data frames
write.table(as.matrix(counts_combined),
            file=paste0(output_dir,"count_matrix.txt"),
            sep="\t",quote=FALSE,row.names=TRUE)

write.table(meta_data_combined, 
            file=paste0(output_dir,"sample_annotation.txt"),
            sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)

#Save the reference groupfs
df_refs<-data.frame(ref_groups=unique(sample_vector))
write.table(df_refs, 
            file=paste0(output_dir,"ref_groups.txt"),
            sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)

