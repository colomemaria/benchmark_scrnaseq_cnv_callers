# ------------------------------------------------------------------------------
# Prepare the input files for the other gastric cancer cell lines
# 1) KATOIII
# 2) MKN-45
# 3) NCI-N87
# Use always the reference from the SNU601 as reference
# and the cellranger count matrices
# ------------------------------------------------------------------------------

library(Matrix)
library(data.table)

#Load reference dataset (focus on control samples 25-29)
patients<-25:29
count_ref_combined<-NULL
sample_vector<-NULL
for(pat in patients){
  
  
  counts_ref<-fread(paste0("data/gastic_cells_SNU601_reference/GSE150290_RAW/",
                           "GSM45463",(pat+22),"_Pat",pat,"-A.txt.gz"))
  gene_names_ref<-counts_ref$V1
  counts_ref$V1<-NULL
  counts_ref<-as.matrix(counts_ref)
  rownames(counts_ref)<-gene_names_ref
  
  if(is.null(count_ref_combined)){
    count_ref_combined<-counts_ref
    sample_vector<-rep(paste0("Patient",pat),ncol(counts_ref))
  } else{
    #print(all(rownames(counts_ref)==rownames(count_ref_combined)))
    count_ref_combined<-cbind(count_ref_combined,counts_ref)
    sample_vector<-c(sample_vector,
                     rep(paste0("Patient",pat),ncol(counts_ref)))
  }
  
}

# ------------------------------------------------------------------------------
# Process all cell lines
# ------------------------------------------------------------------------------

#all_datasets<-c("KATOIII","MKN45","NCIN87","NUGC4","SNU638")
all_datasets<-c("HGC27","SNU16","SNU668")
for(dataset_name in all_datasets){

  #Load matrix we get with cellranger
  data_path<-paste0("data/gastric_cell_lines/",dataset_name)
  counts<-readMM(paste0(data_path,
                           "/filtered_feature_bc_matrix/matrix.mtx.gz"))
  genenames<-fread(paste0(data_path,
                          "/filtered_feature_bc_matrix/features.tsv.gz"),
                   header=FALSE)
  cells<-fread(paste0(data_path,
                      "/filtered_feature_bc_matrix/barcodes.tsv.gz"),
               header=FALSE)
  colnames(counts)<-cells$V1
  rownames(counts)<-genenames$V2
  
  print(paste0("Cell ranger matrix for ", dataset_name,":"))
  print(dim(counts))
  
  #Filter both for a common set of genes
  common_gene_names<-intersect(rownames(counts),rownames(count_ref_combined))
  print(paste0("Common genes with ref: ",length(common_gene_names)))
  counts<-counts[common_gene_names,]
  count_ref_filtered<-count_ref_combined[common_gene_names,]
  counts<-cbind(counts,count_ref_filtered)
  
  dir_path<-paste0("data/gastric_cell_lines/",dataset_name,"/input_",dataset_name)
  dir.create(dir_path)
  
  #Save result matrices
  write.table(as.matrix(counts),sep="\t",quote=FALSE,
              file=paste0(dir_path,"/count_matrix.txt"))
  
  #Generate a sample annotation matrix
  sample_annot<-data.table(barcode=colnames(counts),
                           samples=c(rep(dataset_name,nrow(cells)),sample_vector))
  
  write.table(sample_annot, sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE,
              file=paste0(dir_path,"/sample_annotation.txt"))
  
  #Save the reference groupfs
  df_refs<-data.frame(ref_groups=unique(sample_vector))
  write.table(df_refs, 
              file=paste0(dir_path,"/ref_groups.txt"),
              sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)
}

