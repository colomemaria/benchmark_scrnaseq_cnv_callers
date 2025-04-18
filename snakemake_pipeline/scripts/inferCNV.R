# ------------------------------------------------------------------------------
# Run inferCNV using 6 state HMM and the analysis method "subclusters"
# ------------------------------------------------------------------------------

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# ------------------------------------------------------------------------------
print("Load libraries")
# ------------------------------------------------------------------------------

library(infercnv)

# ------------------------------------------------------------------------------
print("Get input parameters from snakemake")
# ------------------------------------------------------------------------------

input_matrix<-snakemake@input$matrix
input_annotations<-snakemake@input$annot
input_ref_groups<-snakemake@input$ref_groups
input_gene_annotation<-snakemake@input$gene_pos

#Set the threshold for filtering expression dependent on 10X (0.1) vs Smartseq2 (1)
#Default the 10X threshold
param_expr_filter<-snakemake@params$expr_filter
if(is.null(param_expr_filter)){
  param_expr_filter<-0.1
}
param_expr_filter<-as.numeric(param_expr_filter)
  
output_file<-snakemake@output$cnv_file

threads<-as.numeric(snakemake@threads)
print(paste("Running inferCNV with",threads,"threads."))

# ------------------------------------------------------------------------------
print("Execute inferCNV (6 state model and subcluster mode).")
# ------------------------------------------------------------------------------

#Create output directory in case it doesn't exist yet
output_dir<-dirname(output_file)
dir.create(output_dir,recursive=TRUE) # gives a warning if the directory exists already

#Read reference groups (saved in one tsv file)
ref_groups<-read.table(input_ref_groups,header=TRUE)

infercnv_obj = CreateInfercnvObject(raw_counts_matrix=input_matrix,
                                    annotations_file=input_annotations,
                                    delim="\t",
                                    gene_order_file=input_gene_annotation,
                                    ref_group_names=ref_groups$ref_groups)

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=param_expr_filter, # dependent on technology
                             out_dir=output_dir, 
                             cluster_by_groups=TRUE, 
                             denoise=TRUE,
                             HMM=TRUE,
			                       HMM_type="i6",
			                       analysis_mode='subclusters',
			                       num_threads = threads)

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()
