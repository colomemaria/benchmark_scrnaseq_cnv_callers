# ------------------------------------------------------------------------------
# Prepare dataset for COLO320 CNV analysis
# 1) Filter COLO320 to contain only this one cell line and only replicate 1
# 2) Add an external reference dataset (gut cell atlas)
#
# Remark: file is too large to process all samples together with scRNA callers
# ------------------------------------------------------------------------------

library(Matrix)
library(data.table)
library(Seurat)
library(zellkonverter)

# ------------------------------------------------------------------------------
# 1) New approach: use Cellranger count matrix instead of GEO count matrix
# ------------------------------------------------------------------------------

#Load matrix we get with cellranger
counts_colo320<-readMM("filtered_feature_bc_matrix/matrix.mtx.gz")
genenames<-fread("filtered_feature_bc_matrix/features.tsv.gz",
                 header=FALSE)
cells<-fread("filtered_feature_bc_matrix/barcodes.tsv.gz",
             header=FALSE)
colnames(counts_colo320)<-cells$V1
rownames(counts_colo320)<-genenames$V2
print(dim(counts_colo320))

# Previously done: 1) Filter COLO320 to contain only this one cell line
# counts_colo320<-readMM("GSE160148_scRNA_COLO320.mtx.gz")
# dim(counts_colo320)
# 
# barcodes<-fread("GSE160148_scRNA_COLO320_barcodes.tsv.gz")
# barcodes<-barcodes[barcodes$V2!="x",]
# genes<-fread("GSE160148_scRNA_COLO320_features.tsv.gz")
# genes<-genes[genes$V2!="x",]
# 
# colnames(counts_colo320)<-barcodes$V2
# rownames(counts_colo320)<-genes$V2
# 
# #Reformat barcodes
# barcodes$cellline<-gsub("#.*","",barcodes$V2)
# table(barcodes$cellline)
# 
# #Filter for the right cell line and replicate 1
# counts_colo320<-counts_colo320[,startsWith(barcodes$V2,"COLO320HSR_5K_rep1")]
# barcodes<-barcodes[startsWith(barcodes$V2,"COLO320HSR_5K_rep1"),]
# 
# #Rename barcode names to match cell-ranger output later
# barcodes$barcode_reduced<-gsub(".*#","",barcodes$V2)
# 
# colnames(counts_colo320)<-barcodes$barcode_reduced

# ------------------------------------------------------------------------------
# 2) Add an external reference dataset
# ------------------------------------------------------------------------------

counts_ref<-readH5AD("Adult_healthy_large_intestine.h5ad")

#Save cell type annotation
ref_meta_data<-counts_ref@colData

#Filter both matrices
raw_counts_ref<-counts_ref@assays@data$X
rm(counts_ref)
intersect_genes<-intersect(rownames(raw_counts_ref),rownames(counts_colo320))
print(paste("Intersect genes:",length(intersect_genes)))

raw_counts_ref<-raw_counts_ref[intersect_genes,]
counts_colo320<-counts_colo320[intersect_genes,]

#Filter for epithelial cells (take only 5000 of them to not exceed memory limits)
raw_counts_ref<-raw_counts_ref[,ref_meta_data$category == "Epithelial"]
set.seed(1)
raw_counts_ref<-raw_counts_ref[,sample(1:ncol(raw_counts_ref),5000)]

#Combine both matrices
raw_counts_ref<-cbind(raw_counts_ref,counts_colo320)

#Generate a meta_data object for the combined matrix
meta_data_combined<-data.frame(barcode=colnames(raw_counts_ref),
                               sample=c(rep("Epithelial",5000),
                                        rep("COLO320",ncol(counts_colo320))))

#Save results
dir.create("input_COLO320")
write.table(as.matrix(raw_counts_ref),file="input_COLO320/count_matrix.txt",
            sep="\t",quote=FALSE,row.names=TRUE)

write.table(meta_data_combined, 
            file="input_COLO320/sample_annotation.txt",
            sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)

df_refs<-data.frame(ref_groups=setdiff(unique(meta_data_combined$sample),"COLO320"))
write.table(df_refs, file="input_COLO320/ref_groups.txt",
            sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)


