# ------------------------------------------------------------------------------
# Explore MCF7 dataset and combine it with 
# 1) different cell types from an external reference, also breast from tabular sapiens
# 2) the SNU601 cell line (to get an extreme case)
# 3) gastric cells (reference for the SNU601 cell line)
# ------------------------------------------------------------------------------

library(data.table)
library(zellkonverter)
library(Matrix)

# ------------------------------------------------------------------------------
# New approach: get directly count from cell ranger 
# Old approach: Load MCF7 cell line dataset and filter for timepoint 0
# ------------------------------------------------------------------------------

#Load matrix we get with cellranger
counts<-readMM("data/MCF7_breast_cancer_celline/filtered_feature_bc_matrix/matrix.mtx.gz")
genenames<-fread("data/MCF7_breast_cancer_celline/filtered_feature_bc_matrix/features.tsv.gz",
                 header=FALSE)
cells<-fread("data/MCF7_breast_cancer_celline/filtered_feature_bc_matrix/barcodes.tsv.gz",
             header=FALSE)
colnames(counts)<-cells$V1
rownames(counts)<-genenames$V2
print(dim(counts))

# counts<-as.matrix(read.table("data/MCF7_breast_cancer_celline/Bortezomib.csv",
#                              row.names=1,sep=","))
# 
# #Check cancer cell line
# annotation<-data.frame(barcode=colnames(counts))
# annotation$cell_line<-gsub(".*_","",annotation$barcode)
# table(annotation$cell_line)
# 
# #Filter for AA cell line (time point 0)
# counts<-counts[,annotation$cell_line=="MCF7.AA.t0"]
# annotation<-annotation[annotation$cell_line=="MCF7.AA.t0",]
# 
# #Rename barcode names to match later the pileup files from Numbat
# colnames(counts)<-gsub("_MCF7.AA.t0","-1",colnames(counts))
# annotation$barcode<-gsub("_MCF7.AA.t0","-1",annotation$barcode)
# 
# #Save count matrix separately (for processing the second reference dataset)
# counts_MCF7<-counts

# ------------------------------------------------------------------------------
# Load reference dataset with healthy breast samples
# ------------------------------------------------------------------------------

counts_ref<-readH5AD("data/Tabular_sapiens/TS_Mammary.h5ad")

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
intersect_genes<-intersect(rownames(raw_counts_ref),rownames(counts))
print(paste("Intersect genes:",length(intersect_genes)))

raw_counts_ref<-raw_counts_ref[intersect_genes,]
counts<-counts[intersect_genes,]

#Filter for cells measured with 10X
raw_counts_ref<-raw_counts_ref[,ref_meta_data$method== "10X"]
ref_meta_data<-ref_meta_data[ref_meta_data$method=="10X",]

#Combine the MCF7 counts with different reference cells 
#to evaluate the effect of the reference
#for(comp in unique(ref_meta_data$compartment)){
for(comp in c("epithelial")){
    
  #Define output directory for the files
  if(comp=="epithelial"){
    sample_name<-"MCF7"
  } else {
    sample_name<-paste0("MCF7_",comp)
  }
  output_dir<-paste0("data/MCF7_breast_cancer_celline/input_",sample_name,"/")
  dir.create(output_dir, showWarnings = FALSE)
  
  #Filter for epithelial cells (potentially test other references later)
  ref_counts_filtered<-raw_counts_ref[,ref_meta_data$compartment== comp]
  ref_meta_data_filtered<-ref_meta_data[ref_meta_data$compartment== comp,]
  
  #Combine both matrices
  counts_combined<-cbind(ref_counts_filtered,counts)
  
  #Remove white space from cell type names
  ref_meta_data_filtered$celltype_formated<-gsub(" ","_",
                      as.character(ref_meta_data_filtered$cell_ontology_class))
  
  #Generate a meta_data object for the combined matrix
  meta_data_combined<-data.frame(barcode=colnames(counts_combined),
                                 sample=c(as.character(ref_meta_data_filtered$celltype_formated),
                                          rep(sample_name,ncol(counts))))
  
  #Save result data frames
  write.table(as.matrix(counts_combined),
              file=paste0(output_dir,"count_matrix.txt"),
              sep="\t",quote=FALSE,row.names=TRUE)
  
  write.table(meta_data_combined, 
              file=paste0(output_dir,"sample_annotation.txt"),
              sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
  
  #Save the reference groupfs
  df_refs<-data.frame(ref_groups=as.character(unique(ref_meta_data_filtered$celltype_formated)))
  write.table(df_refs, 
              file=paste0(output_dir,"ref_groups.txt"),
              sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)
}

# ------------------------------------------------------------------------------
# Load second reference dataset with SNU601 cells 
# (to have something completely off)
# ------------------------------------------------------------------------------

# #Load SNU601 data-matrix
# counts_SNU601<-readMM("data/SNU601_processed_matrices/GSM4238687_SNU-601_matrix.mtx.gz")
# genenames<-fread("data/SNU601_processed_matrices/GSM4238687_SNU-601_genes.tsv.gz",
#                  header=FALSE)
# cells<-fread("data/SNU601_processed_matrices/GSM4238687_SNU-601_barcodes.tsv.gz",
#              header=FALSE)
# 
# counts_SNU601<-as.matrix(counts_SNU601)
# colnames(counts_SNU601)<-paste0("SNU601_",cells$V1)
# rownames(counts_SNU601)<-genenames$V1
# 
# #Check which genes are part of both count matrices
# intersect_genes<-intersect(rownames(counts_MCF7),rownames(counts_SNU601))
# print(paste("Intersect genes:",length(intersect_genes)))
# 
# counts_MCF7_filtered<-counts_MCF7[intersect_genes,]
# counts_SNU601<-counts_SNU601[intersect_genes,]
# 
# #Define output directory for the files
# sample_name<-"MCF7_SNU601"
# output_dir<-paste0("data/MCF7_breast_cancer_celline/input_",sample_name,"/")
# dir.create(output_dir, showWarnings = FALSE)
# 
# #Combine both matrices
# counts_combined<-cbind(counts_SNU601,counts_MCF7_filtered)
# 
# #Generate a meta_data object for the combined matrix
# meta_data_combined<-data.frame(barcode=colnames(counts_combined),
#                                sample=c(rep("SNU601",ncol(counts_SNU601)),
#                                         rep(sample_name,ncol(counts_MCF7_filtered))))
# 
# #Save result data frames
# write.table(as.matrix(counts_combined),
#             file=paste0(output_dir,"count_matrix.txt"),
#             sep="\t",quote=FALSE,row.names=TRUE)
# 
# write.table(meta_data_combined, 
#             file=paste0(output_dir,"sample_annotation.txt"),
#             sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
# 
# #Save the reference groupfs
# df_refs<-data.frame(ref_groups=c("SNU601"))
# write.table(df_refs, 
#             file=paste0(output_dir,"ref_groups.txt"),
#             sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)


# ------------------------------------------------------------------------------
# Gastric cells (reference for the SNU601 cell line)
# ------------------------------------------------------------------------------

# #Load reference dataset (focus on control samples 25-29)
# patients<-25:29
# count_ref_combined<-NULL
# sample_vector<-NULL
# for(pat in patients){
#   
#   suppressWarnings(counts_ref<-fread(paste0("data/gastic_cells_SNU601_reference/GSE150290_RAW/",
#                                             "GSM45463",(pat+22),"_Pat",pat,"-A.txt.gz")))
#   gene_names_ref<-counts_ref$V1
#   counts_ref$V1<-NULL
#   counts_ref<-as.matrix(counts_ref)
#   rownames(counts_ref)<-gene_names_ref
#   
#   if(is.null(count_ref_combined)){
#     count_ref_combined<-counts_ref
#     sample_vector<-rep(paste0("Patient",pat),ncol(counts_ref))
#   } else{
#     print(all(rownames(counts_ref)==rownames(count_ref_combined)))
#     count_ref_combined<-cbind(count_ref_combined,counts_ref)
#     sample_vector<-c(sample_vector,
#                      rep(paste0("Patient",pat),ncol(counts_ref)))
#   }
# }
# 
# 
# #Check which genes are part of both count matrices
# intersect_genes<-intersect(rownames(counts_MCF7),rownames(count_ref_combined))
# print(paste("Intersect genes:",length(intersect_genes)))
# 
# counts_MCF7_filtered<-counts_MCF7[intersect_genes,]
# count_ref_combined<-count_ref_combined[intersect_genes,]
# 
# #Define output directory for the files
# sample_name<-"MCF7_gastric"
# output_dir<-paste0("data/MCF7_breast_cancer_celline/input_",sample_name,"/")
# dir.create(output_dir, showWarnings = FALSE)
# 
# #Combine both matrices
# counts_combined<-cbind(counts_MCF7_filtered,count_ref_combined)
# 
# #Generate a meta_data object for the combined matrix
# meta_data_combined<-data.frame(barcode=colnames(counts_combined),
#                                sample=c(rep(sample_name, ncol(counts_MCF7_filtered)),
#                                         sample_vector))
# 
# #Save result data frames
# write.table(as.matrix(counts_combined),
#             file=paste0(output_dir,"count_matrix.txt"),
#             sep="\t",quote=FALSE,row.names=TRUE)
# 
# write.table(meta_data_combined, 
#             file=paste0(output_dir,"sample_annotation.txt"),
#             sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
# 
# #Save the reference groupfs
# df_refs<-data.frame(ref_groups=unique(sample_vector))
# write.table(df_refs, 
#             file=paste0(output_dir,"ref_groups.txt"),
#             sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)
# 
# 
