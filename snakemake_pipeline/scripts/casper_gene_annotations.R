# ------------------------------------------------------------------------------
# Get gene annotations to run Casper on a specific data set
# This part requires internet access (to get gene annotations) and is run
# locally for this reason
# ------------------------------------------------------------------------------

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# ------------------------------------------------------------------------------
print("Load libraries")
# ------------------------------------------------------------------------------

library(CaSpER)
library(data.table)

# ------------------------------------------------------------------------------
print("Get input parameters from snakemake")
# ------------------------------------------------------------------------------

input_file<-snakemake@input$matrix

output_gene_annot<-snakemake@output$gene_annot

# ------------------------------------------------------------------------------
print("Get gene annotations to run casper")
# ------------------------------------------------------------------------------

#Load count matrix
count_matrix<-fread(input_file)

#Get centromere information
data(hg38_cytoband)

# Get annotation (works with ensembl_gene_id and hgnc_symbol)
annotation <- generateAnnotation(id_type="hgnc_symbol",
                                 genes=count_matrix$V1,
                                 ishg19=F, centromere_hg38)
write.table(annotation,file=output_gene_annot,
            sep="\t",quote=FALSE,row.names = FALSE)

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()

