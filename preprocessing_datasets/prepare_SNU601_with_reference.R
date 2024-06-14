# ------------------------------------------------------------------------------
# Preprocess SNU-601 data
# Data from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4238687
# Combine with count matrix of a second dataset to get a "heathy" reference
# Data from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE150290 
#
# Prepare also the subsampled SNU601 datasets in the same way
#
# ------------------------------------------------------------------------------

library(Matrix)
library(data.table)
library(ggplot2)

# ------------------------------------------------------------------------------
# Process the original SNU601 dataset - using the published count matrix
# ------------------------------------------------------------------------------

# data_path<-"data/"
# result_path<-"data/input_SNU601_withref/"
# 
# #Load SNU601 data-matrix
# counts<-readMM(paste0(data_path,
#                       "SNU601_processed_matrices/GSM4238687_SNU-601_matrix.mtx.gz"))
# genenames<-fread(paste0(data_path,
#                         "SNU601_processed_matrices/GSM4238687_SNU-601_genes.tsv.gz"),
#                  header=FALSE)
# cells<-fread(paste0(data_path,
#                     "SNU601_processed_matrices/GSM4238687_SNU-601_barcodes.tsv.gz"),
#              header=FALSE)
# 
# counts<-as.matrix(counts)
# colnames(counts)<-cells$V1
# rownames(counts)<-genenames$V1
# 
# counts[1:4,1:4]
# 
# #Load reference dataset (focus on control samples 25-29)
# patients<-25:29
# count_ref_combined<-NULL
# sample_vector<-NULL
# for(pat in patients){
#   
# 
#   counts_ref<-fread(paste0(data_path,
#                            "gastic_cells_SNU601_reference/GSE150290_RAW/",
#                            "GSM45463",(pat+22),"_Pat",pat,"-A.txt.gz"))
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
#   
# }
# 
# #Filter both for a common set of genes
# common_gene_names<-intersect(rownames(counts),rownames(count_ref_combined))
# counts_filtered<-counts[common_gene_names,]
# count_ref_filtered<-count_ref_combined[common_gene_names,]
# counts_filtered<-cbind(counts_filtered,count_ref_filtered)
# 
# #Save result matrices
# write.table(counts_filtered,sep="\t",quote=FALSE,
#             file=paste0(result_path,"count_matrix.txt"))
# 
# #Generate a sample annotation matrix
# sample_annot<-data.table(barcode=colnames(counts_filtered),
#                          samples=c(rep("SNU601",nrow(cells)),sample_vector))
# 
# write.table(sample_annot, sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE,
#             file=paste0(result_path,"sample_annotation.txt"))

# ------------------------------------------------------------------------------
# Process the original SNU601 dataset - using our cell ranger results
# ------------------------------------------------------------------------------

#Load SNU601 matrix we get with cellranger
data_path<-"data/SNU601_processed_matrices/matrix_cellranger70"
counts_cr<-readMM(paste0(data_path,
                      "/filtered_feature_bc_matrix/matrix.mtx.gz"))
genenames<-fread(paste0(data_path,
                        "/filtered_feature_bc_matrix/features.tsv.gz"),
                 header=FALSE)
cells<-fread(paste0(data_path,
                    "/filtered_feature_bc_matrix/barcodes.tsv.gz"),
             header=FALSE)
colnames(counts_cr)<-cells$V1
rownames(counts_cr)<-genenames$V2


# #Filter for the same genes and cells
# dim(counts)
# dim(counts_cr)
# 
# common_cells<-intersect(colnames(counts),colnames(counts_cr))
# common_genes<-intersect(rownames(counts),rownames(counts_cr))
# counts_filtered<-counts[common_genes,common_cells]
# counts_cr_filtered<-counts_cr[common_genes,common_cells]
# 
# 
# # Compare the count differences
# counts_filtered<-melt(as.matrix(counts_filtered))
# counts_cr_filtered<-melt(as.matrix(counts_cr_filtered))
# colnames(counts_filtered)<-c("gene","cell","count_geo")
# counts_filtered$count_cr<-counts_cr_filtered$value
# rm(counts_cr_filtered)
#  
# #Remove points that are both 0
# counts_filtered<-counts_filtered[counts_filtered$count_geo>0 &
#                                    counts_filtered$count_cr>0,]
# 
# g<-ggplot(counts_filtered,aes(x=count_geo,y=count_cr))+
#   geom_point()+geom_abline()+scale_x_log10()+scale_y_log10()+
#   theme_bw()
# ggsave(g,filename="~/Desktop/deviation_counts.png",
#        width=10,height=6)
# 
# # Calculate deviation
# counts_filtered$dev<-counts_filtered$count_geo-counts_filtered$count_cr
# ggplot(counts_filtered,aes(x=dev))+
#   geom_histogram()+xlim(-5,5)

#Filter both for a common set of genes
common_gene_names<-intersect(rownames(counts_cr),rownames(count_ref_combined))
counts_filtered<-counts_cr[common_gene_names,]
count_ref_filtered<-count_ref_combined[common_gene_names,]
counts_filtered<-cbind(counts_filtered,count_ref_filtered)

#Save result matrices
write.table(as.matrix(counts_filtered),sep="\t",quote=FALSE,
            file="data/SNU601_processed_matrices/input_SNU601/count_matrix.txt")

#Generate a sample annotation matrix
sample_annot<-data.table(barcode=colnames(counts_filtered),
                         samples=c(rep("SNU601",ncol(counts_cr)),sample_vector))

write.table(sample_annot, sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE,
            file="data/SNU601_processed_matrices/input_SNU601/sample_annotation.txt")

#Save the reference groupfs
df_refs<-data.frame(ref_groups=unique(sample_vector))
write.table(df_refs, 
            file="data/SNU601_processed_matrices/input_SNU601/ref_groups.txt",
            sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)

# ------------------------------------------------------------------------------
# Process the subsampled SNU601 datasets
# ------------------------------------------------------------------------------

data_path<-"data/"

#Load reference dataset (focus on control samples 25-29)
patients<-25:29
count_ref_combined<-NULL
sample_vector<-NULL
for(pat in patients){
  
  
  counts_ref<-fread(paste0(data_path,
                           "gastic_cells_SNU601_reference/GSE150290_RAW/",
                           "GSM45463",(pat+22),"_Pat",pat,"-A.txt.gz"))
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

data_path<-"code/benchmark_scrna_cnv_caller/snakemake_pipeline/data/"
for(subs in c(20,40,60,80)){
  
  #Load SNU601 data-matrix
  counts<-readMM(paste0(data_path,"input_SNU601_sample",subs,
                        "/filtered_feature_bc_matrix/matrix.mtx.gz"))
  genenames<-fread(paste0(data_path,"input_SNU601_sample",subs,
                         "/filtered_feature_bc_matrix/features.tsv.gz"),
                   header=FALSE)
  cells<-fread(paste0(data_path,"input_SNU601_sample",subs,
                      "/filtered_feature_bc_matrix/barcodes.tsv.gz"),
               header=FALSE)
  colnames(counts)<-cells$V1
  rownames(counts)<-genenames$V2
  
  #Filter genes to the ones also existing in the reference
  common_genes<-intersect(rownames(count_ref_combined),rownames(counts))
  counts<-counts[common_genes,]
  counts_ref_filtered<-count_ref_combined[common_genes,]
  
  #Combine both matrices
  counts_combined<-cbind(counts,counts_ref_filtered)
  
  #Generate a meta_data object for the combined matrix
  sample_name<-paste0("SNU601_sample",subs)
  meta_data_combined<-data.frame(barcode=colnames(counts_combined),
                                 sample=c(rep(sample_name,ncol(counts)),
                                          sample_vector))
  
  #Save result data frames
  output_dir<-paste0(data_path,"input_SNU601_sample",subs,"/")
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
}
