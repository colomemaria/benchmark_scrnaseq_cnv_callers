# ------------------------------------------------------------------------------
# Get gene annotations to run Casper on a specific data set
# This part requires internet access (to get gene annotations) and is run
# locally for this reason
#
# Adaptions made for mouse based on:
# https://github.com/akdess/CaSpER/blob/master/demo/CaSpER_for_Mouse.R
# ------------------------------------------------------------------------------

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# ------------------------------------------------------------------------------
print("Load libraries")
# ------------------------------------------------------------------------------

library(CaSpER)
library(data.table)

source("scripts/casper_functions_for_mouse.R")

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


# Get annotation (works with ensembl_gene_id and hgnc_symbol)
annotation <- generateAnnotationMouse(id_type="mgi_symbol",
                                 genes=count_matrix$V1)

write.table(annotation,file=output_gene_annot,
            sep="\t",quote=FALSE,row.names = FALSE)

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()

