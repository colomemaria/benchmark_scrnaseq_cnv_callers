# ------------------------------------------------------------------------------
# Prepare BCC sample 05 & 07 - pre-treatment
# ------------------------------------------------------------------------------

library(data.table)
library(zellkonverter)
library(Matrix)

# ------------------------------------------------------------------------------
# Run sample 5
# ------------------------------------------------------------------------------

sample_name<-"su005"
new_name<-"BCC05"

#Load dataset
meta_data <- fread("data/BCC_processed_matrices/GSE123813_bcc_all_metadata.txt.gz")

count_matrix<-fread("data/BCC_processed_matrices/GSE123813_bcc_scRNA_counts.txt.gz")
gene_names<-count_matrix$V1
count_matrix$V1<-NULL
count_matrix<-as.matrix(count_matrix)
row.names(count_matrix)<-gene_names

#Filter data for sample, treatment and cell type
meta_data<-meta_data[meta_data$patient==sample_name &
                       meta_data$treatment == "pre",]

meta_data<-meta_data[meta_data$cluster %in% c("Endothelial","Melanocytes",
                                              "Myofibroblasts","Tumor_1","Tumor_2"),]

#Filter count data respectively
count_matrix<-count_matrix[,meta_data$cell.id]

#Reformat cell barcode to cause no problems with numbat later
meta_data$cell.id<-gsub(".*_","",meta_data$cell.id)
meta_data$cell.id<-paste0(meta_data$cell.id,"-1")
colnames(count_matrix)<-meta_data$cell.id

#Rename "Tumor" to match Casper pipeline later
meta_data$cluster[startsWith(meta_data$cluster,"Tumor")]<-new_name

#Create directory
dir.create(paste0("data/BCC_processed_matrices/input_",new_name),
           showWarnings = FALSE)

#Save results
write.table(count_matrix,
            file=paste0("data/BCC_processed_matrices/input_",new_name,
                        "/count_matrix.txt"),
            sep="\t",quote=FALSE,row.names=TRUE)

write.table(meta_data[,c("cell.id","cluster")], 
            file=paste0("data/BCC_processed_matrices/input_",new_name,
                        "/sample_annotation.txt"),
            sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)

df_refs<-data.frame(ref_groups=setdiff(unique(meta_data$cluster),new_name))
write.table(df_refs, 
            file=paste0("data/BCC_processed_matrices/input_",new_name,
                        "/ref_groups.txt"),
            sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)

# ------------------------------------------------------------------------------
# Run sample 7
# ------------------------------------------------------------------------------

sample_name<-"su007"
new_name<-"BCC07"

#Load dataset
meta_data <- fread("data/BCC_processed_matrices/GSE123813_bcc_all_metadata.txt.gz")

count_matrix<-fread("data/BCC_processed_matrices/GSE123813_bcc_scRNA_counts.txt.gz")
gene_names<-count_matrix$V1
count_matrix$V1<-NULL
count_matrix<-as.matrix(count_matrix)
row.names(count_matrix)<-gene_names

#Filter data for sample, treatment and cell type
meta_data<-meta_data[meta_data$patient==sample_name &
                       meta_data$treatment == "pre",]

meta_data<-meta_data[meta_data$cluster %in% c("Endothelial","Melanocytes",
                                              "Myofibroblasts","Tumor_1","Tumor_2"),]

#Filter count data respectively
count_matrix<-count_matrix[,meta_data$cell.id]

#Reformat cell barcode to cause no problems with numbat later
meta_data$cell.id<-gsub(".*_","",meta_data$cell.id)
meta_data$cell.id<-paste0(meta_data$cell.id,"-1")
colnames(count_matrix)<-meta_data$cell.id

#Rename "Tumor" to match Casper pipeline later
meta_data$cluster[startsWith(meta_data$cluster,"Tumor")]<-new_name

#Create directory
dir.create(paste0("data/BCC_processed_matrices/input_",new_name),
           showWarnings = FALSE)

#Save results
write.table(count_matrix,
            file=paste0("data/BCC_processed_matrices/input_",new_name,
                        "/count_matrix.txt"),
            sep="\t",quote=FALSE,row.names=TRUE)

write.table(meta_data[,c("cell.id","cluster")], 
            file=paste0("data/BCC_processed_matrices/input_",new_name,
                        "/sample_annotation.txt"),
            sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)

df_refs<-data.frame(ref_groups=setdiff(unique(meta_data$cluster),new_name))
write.table(df_refs, 
            file=paste0("data/BCC_processed_matrices/input_",new_name,
                        "/ref_groups.txt"),
            sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)


