# ------------------------------------------------------------------------------
# Run CONICSmat - Part 1: obtain gene positions
# This part requires internet access (to get gene annotations) and is run
# locally for this reason
# ------------------------------------------------------------------------------

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# ------------------------------------------------------------------------------
print("Load libraries")
# ------------------------------------------------------------------------------

library(CONICSmat)

# ------------------------------------------------------------------------------
print("Get input parameters from snakemake")
# ------------------------------------------------------------------------------

input_expression_df <- snakemake@input$expression_df

output_gene_pos_df <- snakemake@output$gene_pos_df

# ------------------------------------------------------------------------------
print("Execute CONICSmat  import script")
# ------------------------------------------------------------------------------
expr_df <- as.matrix(read.table(input_expression_df,sep="\t",header=T,row.names=1,check.names=F))

#Source gene positions from ensembl and save them for future use
gene_pos=getGenePositions(rownames(expr_df), ensembl_version = "https://grch37.ensembl.org/")
write.table(gene_pos, file = output_gene_pos_df, sep = "\t", row.names = FALSE, eol = "\r")

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()