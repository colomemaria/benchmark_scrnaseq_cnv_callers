# ------------------------------------------------------------------------------
# Check the count matrix and add an external reference
# HCT116: colon cancer cell line
# Processing the cells as described in the paper:
# - filter cells with < 20,000 counts
# - skipped filtering low ACTB counts
# Take reference colon sample from here:
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE95435 
# ------------------------------------------------------------------------------

library(data.table)
library(ggplot2)
library(GenomicFeatures) #To get TPM counts in the end

mat<-fread("count_matrix_DNT.txt")

mat<-mat[!startsWith(mat$V1,"__"),]

ensg_genes<-mat$V1
gene_symbols<-mat$V2

#To match in the end the ENSG symbols to gene names
mapping<-data.frame(ensg_genes,gene_symbols)

mat$V1<-NULL
mat$V2<-NULL

mat_cols<-fread("count_matrix_DNT_header.txt")
tmp_samples<-gsub(".aligned.sortedByCoord.dedup.bam","",colnames(mat_cols))

mat<-as.matrix(mat)
rownames(mat)<-ensg_genes
colnames(mat)<-tmp_samples
dim(mat)

#Filter for cells with less than < 20,000 counts
mat<-mat[,colSums(mat)>20000]

# #Create a Seurat object
# seurat_obj <- CreateSeuratObject(counts=mat,min.cells=3,min.features=200)
# 
# #Calculate the percentage of mitochondrial counts
# VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2,
#         layer="count")
# seurat_obj <- subset(seurat_obj, subset = nCount_RNA > 20000)
# 
# seurat_obj<-NormalizeData(seurat_obj,normalization.method = "LogNormalize",
#                           scale.factor=10000)
# 
# VlnPlot(seurat_obj, ensg_genes[gene_names=="ACTB"])

# ------------------------------------------------------------------------------
# Process the healthy controls
# In the count matrix are 92 single cell, 3 bulk and one empty well
# => got the respective numbers for bulk and empty from the GEO description
# GSM2511236
# GSM2511244, GSM2511252, GSM2511260
# ------------------------------------------------------------------------------

counts_ref<-fread("GSE95435_P150057_full_gene_count_table.txt")

gene_names<-counts_ref$V1
counts_ref$V1<-NULL
counts_ref<-as.matrix(counts_ref)
rownames(counts_ref)<-gene_names

bulk_headers<-c("WTCHG_168020_710508","WTCHG_168020_711508","WTCHG_168020_712508")
empty_site_control<-"WTCHG_168020_709508"

counts_ref<-counts_ref[,! colnames(counts_ref) %in% c(bulk_headers,empty_site_control)]
dim(counts_ref)
summary(colSums(counts_ref))

#Check matching ENSG names
matched_ensg<-intersect(rownames(mat),rownames(counts_ref))
print(paste("Matching ENSG IDs",length(matched_ensg)))

#Filter both
mat<-mat[matched_ensg,]
counts_ref<-counts_ref[matched_ensg,]

#Check combined coverage
tmp<-data.frame(coverage=c(colSums(mat),colSums(counts_ref)),
           type=c(rep("HCT",ncol(mat)),rep("ref",ncol(counts_ref))))

ggplot(tmp,aes(x=coverage,fill=type))+
  geom_histogram(bins=50)#+
  #geom_vline(xintercept=20000)

combined_counts<-cbind(mat,counts_ref)

#Rename the ENSG to gene symbols
mapping<-mapping[mapping$gene_symbols!="",]
combined_counts<-combined_counts[rownames(combined_counts) %in% mapping$ensg_genes,]
print(paste("Matching ENSG IDs with mapping gene symbol",nrow(combined_counts)))

rownames(mapping)<-mapping$ensg_genes
mapping<-mapping[rownames(combined_counts),]

#Combine duplicate gene names to one row
combined_counts<-apply(combined_counts, 2, tapply, as.factor(mapping$gene_symbols),sum, na.rm=T)

tmp<-data.frame(coverage=colSums(combined_counts),
                type=c(rep("HCT",ncol(mat)),rep("ref",ncol(counts_ref))))

ggplot(tmp,aes(x=coverage,fill=type))+
  geom_histogram(bins=50)

# ------------------------------------------------------------------------------
# Create a TPM normalized matrix
# Transcripts per million:  (count / gene length) / total_count * 1 million
# ------------------------------------------------------------------------------

gtf_file<-"../annotations/Homo_sapiens.GRCh38.113.gtf"
txdb <- makeTxDbFromGFF(gtf_file, format="gtf")
exonic <- exonsBy(txdb, by="gene")
red.exonic <- reduce(exonic)
exon.lengths <- vapply(width(red.exonic), sum, numeric(1))

#Merge to gene symbols
#overlapping_ensid<-intersect(names(exon.lengths),mapping$ensg_genes)
#length(overlapping_ensid)==nrow(mapping)

gene_length<-exon.lengths[mapping$ensg_genes]
gene_symbol_length<-tapply(gene_length,as.factor(mapping$gene_symbols),sum)

#all(names(gene_symbol_length)==rownames(combined_counts))
x <- combined_counts / as.vector(gene_symbol_length)
tpm_mat <- t( t(x) * 1e6 / colSums(x) )

#Save the reduced count matrix
write.table(tpm_mat,
            file="input_HCT116/tpm_count_matrix.txt",
            sep="\t",quote=FALSE,row.names=TRUE)

# ------------------------------------------------------------------------------
# Save the result matrix
# ------------------------------------------------------------------------------

#Save sample annotation
sample_annot<-data.frame(barcode=colnames(combined_counts),
                         sample=c(rep("HCT116",ncol(mat)),rep("healthy",ncol(counts_ref))))

write.table(sample_annot, file="input_HCT116/sample_annotation.txt",
            sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)


df_refs<-data.frame(ref_groups=setdiff(unique(sample_annot$sample),"HCT116"))
write.table(df_refs, file="input_HCT116/ref_groups.txt",
            sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)

#Save the reduced count matrix
write.table(combined_counts,
            file="input_HCT116/count_matrix.txt",
            sep="\t",quote=FALSE,row.names=TRUE)

# ------------------------------------------------------------------------------
# Match bam file list (remove cells filtered out)
# find . -name "*.bam" > list_bam_files.txt‚
# ------------------------------------------------------------------------------

bam_files<-fread("list_bam_files.txt",header=FALSE)

#Extract sample name from file name
bam_files$sample_name<-gsub("./","",bam_files$V1)
bam_files$sample_name<-gsub("\\.aligned.*","",bam_files$sample_name)

#Filter for cells that we keep in the matrix
bam_files<-bam_files[bam_files$sample_name %in% colnames(mat),]

#Add the full path
path_prefix<-"/home/kschmid/benchmark_scrna_cnv_caller/snakemake_pipeline/data/input_HCT116/DNTRseq/HTC116/mRNA/aligned_dedup/"
bam_files$full_path<-gsub("./","",bam_files$V1)
bam_files$full_path<-paste0(path_prefix,bam_files$full_path)

#Save both as separate lists
write.table(data.frame(X=bam_files$V1),file="input_HCT116/list_bam_files_filtered.txt",
            quote=FALSE,row.names = FALSE,col.names = FALSE)

write.table(data.frame(X=bam_files$sample_name),file="input_HCT116/matching_sample_names.txt",
            quote=FALSE,row.names = FALSE,col.names = FALSE)

write.table(data.frame(X=bam_files$full_path),file="input_HCT116/list_bam_files_fullpath.txt",
            quote=FALSE,row.names = FALSE,col.names = FALSE)


