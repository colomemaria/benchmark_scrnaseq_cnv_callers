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

input_species<-snakemake@params$species
if(is.null(input_species)){
  input_species<-"human"
} 
print(paste("Running CONICSmat for the following organism:",input_species))

output_gene_pos_df <- snakemake@output$gene_pos_df

# ------------------------------------------------------------------------------
print("Execute CONICSmat  import script")
# ------------------------------------------------------------------------------
expr_df <- as.matrix(read.table(input_expression_df,sep="\t",header=T,row.names=1,check.names=F))

#Source gene positions from ensembl and save them for future use
if (input_species=="human"){
  gene_pos<-getGenePositions(rownames(expr_df), species = "human",
                            ensembl_version = "https://www.ensembl.org") #currently: GRCh38.p14
} else if (input_species == "mouse"){
  gene_pos<-getGenePositions(rownames(expr_df), species = "mouse",
                             ensembl_version = "https://nov2020.archive.ensembl.org") #mm10 (GRCm38.p6)
} else {
  stop("Species definition not known. Only human and mouse implemented currently.")
}

write.table(gene_pos, file = output_gene_pos_df, sep = "\t", row.names = FALSE, eol = "\r")

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()
