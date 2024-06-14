# ------------------------------------------------------------------------------
# Evaluate the results of the different CNV callers compared to a ground truth
# dataset (typically WGS)
# Variant: evaluate the results separately for each method (not only looking 
# at the overlap between methods)
# ------------------------------------------------------------------------------

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# ------------------------------------------------------------------------------
print("Load libraries")
# ------------------------------------------------------------------------------

library(data.table)
library(dplyr)
library(GenomicRanges)
library(reshape2)
library(ROCR)
library(pROC)
library(numbat)
library(crfsuite)

#Source script with help functions
source("scripts/lib.R")

# ------------------------------------------------------------------------------
print("Get input parameters from snakemake")
# ------------------------------------------------------------------------------

#Input files
input_annot<-snakemake@input$annot
input_ref_groups<-snakemake@input$ref_groups
input_wgs<-snakemake@input$genomic_ground_truth

input_infercnv_cnv<-snakemake@input$infercnv_cnv
input_infercnv_expr<-snakemake@input$infercnv_expr
input_infercnv_gene_pos<-snakemake@input$infercnv_gene_pos

input_casper<-snakemake@input$casper

input_copykat<-snakemake@input$copykat

input_scevan_expr<-snakemake@input$scevan
input_scevan_gene_pos<-snakemake@input$scevan_gene_pos
input_scevan_clone1<-snakemake@input$scevan_clone1
input_scevan_annot<-snakemake@input$scevan_annot

input_numbat_obj<-snakemake@input$numbat_obj

# input_CONICSmat_chrom_pos<-snakemake@input$CONICSmat_chr_pos
# input_CONICSmat_cnv<-snakemake@input$CONICSmat_cnv

#Parameter how to interpret the WGS/WES data (default set to WGS)
param_genomic_format<-snakemake@params$genomic_format
if(is.null(param_genomic_format)){
  param_genomic_format<-"WGS"
}

#File for the output
output_res<-snakemake@output$res

# ------------------------------------------------------------------------------
print("Read WGS results")
# ------------------------------------------------------------------------------

print("Load genomic data (scWGS/WGS/WES)")
if(param_genomic_format == "WGS"){
  wgs_results<-read_wgs_results(input_wgs)
  
  print(paste("Reading (sc)WGS results with in total",
              length(wgs_results),"bins!"))
} else if (param_genomic_format == "CNVkit"){
  wgs_results<-read_results_cnvkit_segs(input_wgs)
  
  print(paste("Reading CNVkit results with in total",
              length(wgs_results),"bins!")) 
} else if (param_genomic_format == "GATK"){
  wgs_results<-read_wes_results_gatk(input_wgs)
  
  print(paste("Reading GATK results with in total",
              length(wgs_results),"bins!"))
} else {
  stop(paste("Format of the genomic data is unknown!",
             "Need to be WGS, CNVkit or GATK, not",param_genomic_format))
}

#Initialize all output values
res_per_method<-NULL

# ------------------------------------------------------------------------------
print("Evaluate inferCNV (cnv model)")
# ------------------------------------------------------------------------------

#InferCNV (CNV model)
infercnv_results<-read_infercnv_6state_model(input_infercnv_cnv,
                                             input_infercnv_gene_pos)
  
combined_range<-combine_range_objects(wgs_results,infercnv_results[[1]],
                                      method_colname="infercnv_cnv")

res<-get_all_metrics(combined_range,"infercnv_cnv")

res_per_method<-rbind(res_per_method,res)
                      
                      
# ------------------------------------------------------------------------------
print("Evaluate inferCNV (expr model)")
# ------------------------------------------------------------------------------

#InferCNV (normalized expression)
infercnv_results<-read_infercnv_expr(input_infercnv_expr,
                                     input_infercnv_gene_pos)

combined_range<-combine_range_objects(wgs_results,infercnv_results[[1]],
                                      method_colname="infercnv_expr")

res<-get_all_metrics(combined_range,"infercnv_expr")

res_per_method<-rbind(res_per_method,res)

# ------------------------------------------------------------------------------
print("Evaluate Casper")
# ------------------------------------------------------------------------------

casper_results<-read_casper(input_casper)

combined_range<-combine_range_objects(wgs_results,casper_results,
                                      method_colname="casper")

res<-get_all_metrics(combined_range,"casper")

res_per_method<-rbind(res_per_method,res)

# ------------------------------------------------------------------------------
print("Evaluate copyKat")
# ------------------------------------------------------------------------------

copykat_results<-read_copykat(input_copykat,input_annot,input_ref_groups)

combined_range<-combine_range_objects(wgs_results,copykat_results[[1]],
                                      method_colname="copykat")

res<-get_all_metrics(combined_range,"copykat")

res_per_method<-rbind(res_per_method,res)

# ------------------------------------------------------------------------------
print("Evaluate SCEVAN (expr)")
# ------------------------------------------------------------------------------

scevan_results<-read_scevan(input_scevan_expr,input_scevan_gene_pos,
                            input_annot,input_ref_groups)

combined_range<-combine_range_objects(wgs_results,scevan_results[[1]],
                                      method_colname="scevan")

res<-get_all_metrics(combined_range,"scevan")

res_per_method<-rbind(res_per_method,res)

# ------------------------------------------------------------------------------
print("Evaluate SCEVAN (CNV)")
# ------------------------------------------------------------------------------

binned_genome<-combined_range
elementMetadata(binned_genome)<-NULL
scevan_results<-read_scevan_cn_status(input_scevan_clone1, input_scevan_annot,
                                      binned_genome)

combined_range<-combine_range_objects(combined_range,scevan_results,
                                      method_colname="scevan_cnv")

res<-get_all_metrics(combined_range,"scevan_cnv")

res_per_method<-rbind(res_per_method,res)

# ------------------------------------------------------------------------------
print("Evaluate Numbat (expr)")
# ------------------------------------------------------------------------------

#Get the directory name
input_numbat_obj<-dirname(input_numbat_obj)

numbat_results <- read_numbat_expr(input_numbat_obj)

combined_range<-combine_range_objects(wgs_results, numbat_results[[1]],
                                      method_colname="numbat")

res<-get_all_metrics(combined_range,"numbat")

res_per_method<-rbind(res_per_method,res)

# ------------------------------------------------------------------------------
print("Evaluate Numbat (CNV)")
# ------------------------------------------------------------------------------

numbat_results <- read_numbat_cnv(input_numbat_obj)

combined_range<-combine_range_objects(wgs_results, numbat_results,
                                      method_colname="numbat_cnv")

res<-get_all_metrics(combined_range,"numbat_cnv")

res_per_method<-rbind(res_per_method,res)

# ------------------------------------------------------------------------------
print("Save final results")
# ------------------------------------------------------------------------------

write.table(res_per_method,file=output_res,sep="\t",quote=FALSE,
            row.names=FALSE)

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()
