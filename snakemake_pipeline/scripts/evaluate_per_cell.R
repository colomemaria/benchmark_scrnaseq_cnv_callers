# ------------------------------------------------------------------------------
# Per cell comparison for the RNA benchmarking (pearson correlation, AUC gains &
# AUC losses) - exemplarily shown for SNU601
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
library(ggplot2)
library(ggpubr)
library(ROCR)
library(pROC)
library(viridis)
library(numbat)
library(crfsuite)

theme_set(theme_bw())

#Source script with help functions
source("scripts/lib.R")

# ------------------------------------------------------------------------------
# Input variables
# ------------------------------------------------------------------------------

#Input files
input_annot<-snakemake@input$annot
input_ref_groups<-snakemake@input$ref_groups
input_combined_res<-snakemake@input$combined_res

input_infercnv_cnv<-snakemake@input$infercnv_cnv
input_infercnv_expr<-snakemake@input$infercnv_expr
input_infercnv_gene_pos<-snakemake@input$infercnv_gene_pos

input_casper_grange<-snakemake@input$casper_grange
input_casper_cells<-snakemake@input$casper_cells

input_copykat<-snakemake@input$copykat

input_scevan_expr<-snakemake@input$scevan
input_scevan_gene_pos<-snakemake@input$scevan_gene_pos
input_scevan_clone1<-snakemake@input$scevan_clone1
input_scevan_annot<-snakemake@input$scevan_annot

input_numbat_obj<-snakemake@input$numbat_obj
input_numbat_obj<-dirname(input_numbat_obj)

input_CONICSmat_chrom_pos<-snakemake@input$CONICSmat_chr_pos
input_CONICSmat_cnv<-snakemake@input$CONICSmat_cnv

#Output files
output_total_res<-snakemake@output$total_res
output_plot_corr<-snakemake@output$plot_corr
output_plot_auc_gain<-snakemake@output$plot_auc_gain
output_plot_auc_loss<-snakemake@output$plot_auc_loss
output_plot_auc_trunc_gain<-snakemake@output$plot_auc_trunc_gain
output_plot_auc_trunc_loss<-snakemake@output$plot_auc_trunc_loss

# ------------------------------------------------------------------------------
print("Read each method and evaluate the results")
# ------------------------------------------------------------------------------

total_res_df<-NULL

#Look at the same regions as for the complete evaluation (=> only intersect per methods)
print("Get genomic ground truth (already filtered for methods intersect")
combined_res<-fread(input_combined_res)
wgs_results<-makeGRangesFromDataFrame(combined_res)
wgs_results$wgs_mean<-combined_res$wgs_mean
rm(combined_res)

print("Evaluate cells from InferCNV (CNV)")
infercnv_cells<-read_infercnv_6state_model_individual_cells(input_infercnv_cnv,
                                            input_infercnv_gene_pos)
infercnv_res_df<- evaluate_metrics_per_cell(wgs_results,infercnv_cells[[1]],
                                           infercnv_cells[[2]])
infercnv_res_df$method<-"InferCNV (CNV)"
total_res_df<-rbind(total_res_df,infercnv_res_df)

print("Evaluate cells from InferCNV (Expr)")
infercnv_cells<-read_infercnv_expr_individual_cells(input_infercnv_expr,
                                    input_infercnv_gene_pos)
infercnv_res_df<- evaluate_metrics_per_cell(wgs_results,infercnv_cells[[1]],
                                            infercnv_cells[[2]])
infercnv_res_df$method<-"InferCNV (Expr)"
total_res_df<-rbind(total_res_df,infercnv_res_df)

print("Evaluate cells from CaSpER")
casper_cells<-read_casper_individual_cells(input_casper_grange,input_casper_cells)
casper_res_df<- evaluate_metrics_per_cell(wgs_results,casper_cells[[1]],
                                                casper_cells[[2]])
casper_res_df$method<-"CaSpER"
total_res_df<-rbind(total_res_df,casper_res_df)

print("Evaluate cells from copyKat")
copykat_cells<-read_copykat_individual_cells(input_copykat,input_annot,input_ref_groups)
copykat_res_df<- evaluate_metrics_per_cell(wgs_results,copykat_cells[[1]],
                                           copykat_cells[[2]])
copykat_res_df$method<-"CopyKat"
total_res_df<-rbind(total_res_df,copykat_res_df)

print("Evaluate cells from SCEVAN (Expr)")
scevan_cells<-read_scevan_individual_cells(input_scevan_expr,input_scevan_gene_pos,
                                           input_annot,input_ref_groups)
scevan_res_df<- evaluate_metrics_per_cell(wgs_results,scevan_cells[[1]],
                                           scevan_cells[[2]])
scevan_res_df$method<-"SCEVAN (Expr)"
total_res_df<-rbind(total_res_df,scevan_res_df)

print("Evaluate cells from SCEVAN (CNV)")
scevan_cells<-read_scevan_cn_status_individual_cells(input_scevan_clone1, input_scevan_annot,
                                                     wgs_results)
scevan_res_df<- evaluate_metrics_per_cell(wgs_results,scevan_cells[[1]],
                                          scevan_cells[[2]])
scevan_res_df$method<-"SCEVAN (CNV)"
total_res_df<-rbind(total_res_df,scevan_res_df)

print("Evaluate cells from Numbat (CNV)")
numbat_cells<-read_numbat_cnv_individual_cells(input_numbat_obj)
numbat_res_df<- evaluate_metrics_per_cell(wgs_results,numbat_cells[[1]],
                                          numbat_cells[[2]])
numbat_res_df$method<-"Numbat (CNV)"
total_res_df<-rbind(total_res_df,numbat_res_df)

print("Evaluate cells from Numbat (Expr)")
numbat_cells <- read_numbat_expr_individual_cells(input_numbat_obj)
numbat_res_df<- evaluate_metrics_per_cell(wgs_results,numbat_cells[[1]],
                                          numbat_cells[[2]])
numbat_res_df$method<-"Numbat (Expr)"
total_res_df<-rbind(total_res_df,numbat_res_df)

print("Evaluate cells from CONICSmat")
CONICSmat_cells <- read_CONICSmat_individual_cells(input_CONICSmat_chrom_pos, 
                                                   input_CONICSmat_cnv,
                                  input_annot, input_ref_groups)

CONICSmat_res_df<- evaluate_metrics_per_cell(wgs_results,CONICSmat_cells[[1]],
                                          CONICSmat_cells[[2]])
CONICSmat_res_df$method<-"CONICSmat"
total_res_df<-rbind(total_res_df,CONICSmat_res_df)

# ------------------------------------------------------------------------------
print("Plot and save all the results")
# ------------------------------------------------------------------------------

write.table(total_res_df,file=output_total_res,
            quote = FALSE,row.names=FALSE,sep="\t")

#Plot overview
g<-ggplot(total_res_df,aes(x=cor,fill=method))+
  geom_histogram()+
  xlab("Pearson correlation")+ylab("Number of cells")+
  facet_wrap(~method)+
  theme(legend.position = "none")
ggsave(g,file=output_plot_corr,width=12,height=8)
  
g<-ggplot(total_res_df,aes(x=auc_gain,fill=method))+
  geom_histogram()+
  xlab("AUC (gains)")+ylab("Number of cells")+
  facet_wrap(~method)+
  theme(legend.position = "none")
ggsave(g,file=output_plot_auc_gain,width=12,height=8)

g<-ggplot(total_res_df,aes(x=auc_loss,fill=method))+
  geom_histogram()+
  xlab("AUC (losses)")+ylab("Number of cells")+
  facet_wrap(~method)+
  theme(legend.position = "none")
ggsave(g,file=output_plot_auc_loss,width=12,height=8)

g<-ggplot(total_res_df,aes(x=auc_trunc_gain,fill=method))+
  geom_histogram()+
  xlab("Truncated AUC (gains)")+ylab("Number of cells")+
  facet_wrap(~method)+
  theme(legend.position = "none")
ggsave(g,file=output_plot_auc_trunc_gain,width=12,height=8)

g<-ggplot(total_res_df,aes(x=auc_trunc_loss,fill=method))+
  geom_histogram()+
  xlab("Truncated AUC (losses)")+ylab("Number of cells")+
  facet_wrap(~method)+
  theme(legend.position = "none")
ggsave(g,file=output_plot_auc_trunc_loss,width=12,height=8)

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------

sessionInfo()





  

    
  
  

