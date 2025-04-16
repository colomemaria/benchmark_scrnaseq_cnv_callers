# ------------------------------------------------------------------------------
# Check the count matrix and add an external reference
# A375: malignant melanoma cell line
# Processing the cells as described in the paper:
# - filter cells with < 20,000 counts
# - skipped filtering low ACTB counts
# Take reference sample from here:
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151091
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

#Combine duplicate gene names to one row (as second sample is based on gene names)
mat<-apply(mat, 2, tapply, as.factor(mapping$gene_symbols),sum, na.rm=T)

# ------------------------------------------------------------------------------
# Adding healthy controls as reference
# Subsample to a matching control (A375 froma 54-year-old female patient)
# ------------------------------------------------------------------------------

counts_ref<-fread("GSE151091/raw_matrix.csv")
gene_names<-counts_ref$V1
counts_ref$V1<-NULL
counts_ref<-as.matrix(counts_ref)
rownames(counts_ref)<-gene_names

#Filter for adult samples
metadata_ref<-fread("GSE151091/GSE151091_Metadata.csv")
metadata_ref<-metadata_ref[metadata_ref$age>=20,]

counts_ref<-counts_ref[,metadata_ref$cell_id]

#Check matching ENSG names
matched_symbols<-intersect(rownames(mat),rownames(counts_ref))
print(paste("Matching gene symbols",length(matched_symbols)))

#Filter both
mat<-mat[matched_symbols,]
counts_ref<-counts_ref[matched_symbols,]

# #Check combined coverage
# tmp<-data.frame(coverage=c(colSums(mat),colSums(counts_ref)),
#                 type=c(rep("A375",ncol(mat)),rep("ref",ncol(counts_ref))))
# 
# ggplot(tmp,aes(x=coverage,fill=type))+
#   geom_histogram(bins=50)#+

#Combine both
combined_counts<-cbind(mat,counts_ref)

# ------------------------------------------------------------------------------
# Create a TPM normalized matrix
# Transcripts per million:  (count / gene length) / total_count * 1 million
# ------------------------------------------------------------------------------

gtf_file<-"../annotations/Homo_sapiens.GRCh38.113.gtf"
txdb <- makeTxDbFromGFF(gtf_file, format="gtf")
exonic <- exonsBy(txdb, by="gene")
red.exonic <- reduce(exonic)
exon.lengths <- vapply(width(red.exonic), sum, numeric(1))

#Filter mapping file
mapping<-mapping[mapping$gene_symbols %in% matched_symbols,]

gene_length<-exon.lengths[mapping$ensg_genes]
gene_symbol_length<-tapply(gene_length,as.factor(mapping$gene_symbols),sum)

#all(names(gene_symbol_length)==rownames(combined_counts))
x <- combined_counts / as.vector(gene_symbol_length)
tpm_mat <- t( t(x) * 1e6 / colSums(x) )

#Save the reduced count matrix
write.table(tpm_mat,
            file="input_A375/tpm_count_matrix.txt",
            sep="\t",quote=FALSE,row.names=TRUE)

# ------------------------------------------------------------------------------
# Save the result matrix
# ------------------------------------------------------------------------------

#Save sample annotation
sample_annot<-data.frame(barcode=colnames(combined_counts),
                         sample=c(rep("A375",ncol(mat)),rep("healthy",ncol(counts_ref))))

write.table(sample_annot, file="input_A375/sample_annotation.txt",
            sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)


df_refs<-data.frame(ref_groups=setdiff(unique(sample_annot$sample),"A375"))
write.table(df_refs, file="input_A375/ref_groups.txt",
            sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)

#Save the reduced count matrix
write.table(combined_counts,
            file="input_A375/count_matrix.txt",
            sep="\t",quote=FALSE,row.names=TRUE)

# ------------------------------------------------------------------------------
# Match bam file list (remove cells filtered out)
# find . -name "*.bam" > list_bam_files.txt
# ------------------------------------------------------------------------------

bam_files<-fread("list_bam_files.txt",header=FALSE)

#Extract sample name from file name
bam_files$sample_name<-gsub("./","",bam_files$V1)
bam_files$sample_name<-gsub("\\.aligned.*","",bam_files$sample_name)

#Filter for cells that we keep in the matrix
bam_files<-bam_files[bam_files$sample_name %in% colnames(mat),]

#Add the full path
path_prefix<-"/home/kschmid/benchmark_scrna_cnv_caller/snakemake_pipeline/data/input_A375/A375/mRNA/aligned_dedup/"
bam_files$full_path<-gsub("./","",bam_files$V1)
bam_files$full_path<-paste0(path_prefix,bam_files$full_path)

#Save both as separate lists
write.table(data.frame(X=bam_files$V1),file="input_A375/list_bam_files_filtered.txt",
            quote=FALSE,row.names = FALSE,col.names = FALSE)

write.table(data.frame(X=bam_files$sample_name),file="input_A375/matching_sample_names.txt",
            quote=FALSE,row.names = FALSE,col.names = FALSE)

write.table(data.frame(X=bam_files$full_path),file="input_A375/list_bam_files_fullpath.txt",
            quote=FALSE,row.names = FALSE,col.names = FALSE)



