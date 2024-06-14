# ------------------------------------------------------------------------------
# Prepare the BCC data (combine samples to identify substructures)
# Current mode of action: analyse the samples together 
#                         (not labeling them independently), 
#                         focusing on pre-treatment (start with sample08)
#                         use non-immune, non-malignant cells as reference
#                         (as done by the original publication)
# ------------------------------------------------------------------------------

library(data.table)
library(dplyr)

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

#Filter data for treatment and cell type (merge Tumor1 & Tumor2)
meta_data<-meta_data[meta_data$treatment == "pre",]
meta_data$cluster[meta_data$cluster %in% c("Tumor_1","Tumor_2")]<-"Tumor"
meta_data<-meta_data[meta_data$cluster %in% c("Endothelial","Melanocytes",
                                              "Myofibroblasts","Tumor"),]

#Keep only patients with at least 10 tumor cells
num_tumor_cells<-meta_data%>%
  filter(cluster=="Tumor")%>%
  group_by(patient)%>%
  summarize(n_cells=n())%>%
  filter(n_cells>=20)

meta_data<-meta_data[meta_data$patient %in% num_tumor_cells$patient,]

#Add -1 to barcode (as in original bam files)
meta_data$cell.id<-paste0(meta_data$cell.id,"-1")
colnames(count_matrix)<-paste0(colnames(count_matrix),"-1")

#Check barcode prefix of each sample
meta_data$prefix_barcode<-sapply(strsplit(meta_data$cell.id,"_"), `[[`, 1)
meta_data$raw_barcode<-sapply(strsplit(meta_data$cell.id,"_"), `[[`, 2)

table(meta_data$prefix_barcode[meta_data$cluster=="Tumor"])
table(meta_data$patient[meta_data$cluster=="Tumor"],meta_data$sort[meta_data$cluster=="Tumor"])

#Save selected barcodes for each sample to filter bam files and reduce runtime
mapped_samples<-setNames(c("bcc.su005.pre.tumor","bcc.su006.pre.tumor",
                           "bcc.su007.pre.tumor.cd45","bcc.su008.pre.tumor"),
                         c("su005","su006","su007","su008"))
for(pat in names(mapped_samples)){
  print(pat)
  barcodes<-data.frame(barcodes=paste0("CB:Z:",
                                       meta_data$raw_barcode[meta_data$prefix_barcode == 
                                                               mapped_samples[pat]]))
  write.table(barcodes,quote=FALSE,row.names=FALSE,col.names=FALSE,
              file=paste0("data/BCC_processed_matrices/filtered_barcodes/barcodes_",
                          pat,".tsv"))
}

#Save a general barcode file (for numbat)
write.table(meta_data[,c("cell.id")],quote=FALSE,row.names=FALSE,col.names=FALSE,
            file=paste0("data/BCC_processed_matrices/barcodes.tsv"))

#Filter count data respectively
count_matrix<-count_matrix[,meta_data$cell.id]

#Rename "Tumor" to match Casper pipeline later
meta_data$cluster[startsWith(meta_data$cluster,"Tumor")]<-"BCC"

#Save results
write.table(count_matrix,
            file="data/BCC_processed_matrices/input_BCC/count_matrix.txt",
            sep="\t",quote=FALSE,row.names=TRUE)

write.table(meta_data[,c("cell.id","cluster")], 
            file="data/BCC_processed_matrices/input_BCC/sample_annotation.txt",
            sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)

df_refs<-data.frame(ref_groups=setdiff(unique(meta_data$cluster),"BCC"))
write.table(df_refs, 
            file="data/BCC_processed_matrices/input_BCC/ref_groups.txt",
            sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)

#Save true distribution of donors (focus on tumor cells)
patient_clusters<-meta_data[meta_data$cluster=="BCC",c("cell.id","patient")]
colnames(patient_clusters)<-c("cell","cluster")
write.table(patient_clusters, 
            file="data/BCC_processed_matrices/input_BCC/BCC_metadata.txt",
            sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)

