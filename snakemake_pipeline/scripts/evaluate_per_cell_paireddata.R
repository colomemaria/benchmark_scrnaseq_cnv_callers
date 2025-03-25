# ------------------------------------------------------------------------------
# Per cell comparison for the RNA benchmarking (pearson correlation, AUC gains &
# AUC losses)
# Variant implemented for paired data (so matching cells)
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
input_wgs_cells<-snakemake@input$wgs_cells
input_pb_res<-snakemake@input$pb_res
  
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

input_CONICSmat_chrom_pos<-snakemake@input$CONICSmat_chr_pos
input_CONICSmat_cnv<-snakemake@input$CONICSmat_cnv

#Output files
output_total_res<-snakemake@output$total_res
output_plot_corr<-snakemake@output$plot_corr
output_plot_auc_gain<-snakemake@output$plot_auc_gain
output_plot_auc_loss<-snakemake@output$plot_auc_loss
output_plot_auc_trunc_gain<-snakemake@output$plot_auc_trunc_gain
output_plot_auc_trunc_loss<-snakemake@output$plot_auc_trunc_loss

output_plot_karyograms<-snakemake@output$plot_karyograms

# ------------------------------------------------------------------------------
print("Read each method and evaluate the results")
# ------------------------------------------------------------------------------

total_res_df<-NULL
karyoplots <- list()

print("Get genomic ground truth per cell")
combined_res<-fread(input_pb_res)
wgs_cell<-fread(input_wgs_cells)
  
#Look at the same regions as for the complete evaluation (=> only intersect per methods)
wgs_pb<-makeGRangesFromDataFrame(combined_res)
wgs_cell_range<-makeGRangesFromDataFrame(wgs_cell)
overlaps<-as.data.frame(findOverlaps(wgs_cell_range,wgs_pb))
wgs_cell_range<-wgs_cell_range[overlaps$queryHits,]
wgs_cell<-wgs_cell[overlaps$queryHits,]
rm(combined_res,wgs_pb)

#Convert into a matrix
wgs_cell<-wgs_cell[,-c("seqnames","start","end")]
wgs_cell<-as.matrix(wgs_cell)

print("Evaluate cells from InferCNV (CNV)")
infercnv_cells<-read_infercnv_6state_model_individual_cells(input_infercnv_cnv,
                                            input_infercnv_gene_pos)
infercnv_res_df<- evaluate_metrics_per_cell_paired(wgs_cell_range,wgs_cell,infercnv_cells[[1]],
                                           infercnv_cells[[2]])
infercnv_res_df$method<-"InferCNV (CNV)"
infercnv_res_df$cells<-rownames(infercnv_res_df)
total_res_df<-rbind(total_res_df,infercnv_res_df)

#Plot also the karyogram
g<-plot_karyogram_paired(wgs_cell_range,wgs_cell,infercnv_cells[[1]],
                         infercnv_cells[[2]],"InferCNV (CNV)")
karyoplots <- c(karyoplots,list(g))

print("Evaluate cells from InferCNV (Expr)")
infercnv_cells<-read_infercnv_expr_individual_cells(input_infercnv_expr,
                                    input_infercnv_gene_pos)
infercnv_res_df<- evaluate_metrics_per_cell_paired(wgs_cell_range,wgs_cell,infercnv_cells[[1]],
                                            infercnv_cells[[2]])
infercnv_res_df$method<-"InferCNV (Expr)"
infercnv_res_df$cells<-rownames(infercnv_res_df)
total_res_df<-rbind(total_res_df,infercnv_res_df)

#Plot also the karyogram
g<-plot_karyogram_paired(wgs_cell_range,wgs_cell,infercnv_cells[[1]],
                         infercnv_cells[[2]],"InferCNV (Expr)")
karyoplots <- c(karyoplots,list(g))

#Check if CaSpER output files are provided
if(! is.null(input_casper_grange)){
  
  print("Evaluate cells from CaSpER")
  casper_cells<-read_casper_individual_cells(input_casper_grange,input_casper_cells)
  
  if(! is.null(casper_cells)){
    casper_res_df<- evaluate_metrics_per_cell_paired(wgs_cell_range,wgs_cell,casper_cells[[1]],
                                                     casper_cells[[2]])
    casper_res_df$method<-"CaSpER"
    casper_res_df$cells<-rownames(casper_res_df)
    total_res_df<-rbind(total_res_df,casper_res_df)
    
    #Plot also the karyogram
    g<-plot_karyogram_paired(wgs_cell_range,wgs_cell,casper_cells[[1]],
                             casper_cells[[2]],"CaSpER")
    karyoplots <- c(karyoplots,list(g))
    
  } else {
    print("Skipping casper (no CNVs found).")
  }
} else {
  print("Skipping CaSpER (no output file provided).")  
}

print("Evaluate cells from copyKat")
copykat_cells<-read_copykat_individual_cells(input_copykat,input_annot,input_ref_groups)
copykat_res_df<- evaluate_metrics_per_cell_paired(wgs_cell_range,wgs_cell,copykat_cells[[1]],
                                           copykat_cells[[2]])
copykat_res_df$method<-"CopyKat"
copykat_res_df$cells<-rownames(copykat_res_df)
total_res_df<-rbind(total_res_df,copykat_res_df)

#Plot also the karyogram
g<-plot_karyogram_paired(wgs_cell_range,wgs_cell,copykat_cells[[1]],
                         copykat_cells[[2]],"CopyKat")
karyoplots <- c(karyoplots,list(g))

print("Evaluate cells from SCEVAN (Expr)")
scevan_cells<-read_scevan_individual_cells(input_scevan_expr,input_scevan_gene_pos,
                                           input_annot,input_ref_groups)
scevan_res_df<- evaluate_metrics_per_cell_paired(wgs_cell_range,wgs_cell,scevan_cells[[1]],
                                           scevan_cells[[2]])
scevan_res_df$method<-"SCEVAN (Expr)"
scevan_res_df$cells<-rownames(scevan_res_df)
total_res_df<-rbind(total_res_df,scevan_res_df)

#Plot also the karyogram
g<-plot_karyogram_paired(wgs_cell_range,wgs_cell,scevan_cells[[1]],
                         scevan_cells[[2]],"SCEVAN (Expr)")
karyoplots <- c(karyoplots,list(g))

#Plot also the karyogram
g<-plot_karyogram_paired(wgs_cell_range,wgs_cell,scevan_cells[[1]],
                         scevan_cells[[2]],"SCEVAN (Expr)")
karyoplots <- c(karyoplots,list(g))

print("Evaluate cells from SCEVAN (CNV)")
scevan_cells<-read_scevan_cn_status_individual_cells(input_scevan_clone1, input_scevan_annot,
                                                     input_annot, input_ref_groups,
                                                     wgs_cell_range)
scevan_res_df<- evaluate_metrics_per_cell_paired(wgs_cell_range,wgs_cell,scevan_cells[[1]],
                                                 scevan_cells[[2]])
scevan_res_df$method<-"SCEVAN (CNV)"
scevan_res_df$cells<-rownames(scevan_res_df)
total_res_df<-rbind(total_res_df,scevan_res_df)

#Plot also the karyogram
g<-plot_karyogram_paired(wgs_cell_range,wgs_cell,scevan_cells[[1]],
                         scevan_cells[[2]],"SCEVAN (CNV)")
karyoplots <- c(karyoplots,list(g))

#Check if Numbat output files are provided
if(! is.null(input_numbat_obj)){
  
  print("Evaluate cells from Numbat (CNV)")
  
  input_numbat_obj<-dirname(input_numbat_obj)
  numbat_cells<-read_numbat_cnv_individual_cells(input_numbat_obj)
  
  if(! is.null(numbat_cells)){
    numbat_res_df<- evaluate_metrics_per_cell_paired(wgs_cell_range,wgs_cell,numbat_cells[[1]],
                                                     numbat_cells[[2]])
    numbat_res_df$method<-"Numbat (CNV)"
    numbat_res_df$cells<-rownames(numbat_res_df)
    total_res_df<-rbind(total_res_df,numbat_res_df)
    
    #Plot also the karyogram
    g<-plot_karyogram_paired(wgs_cell_range,wgs_cell,numbat_cells[[1]],
                             numbat_cells[[2]],"Numbat (CNV)")
    karyoplots <- c(karyoplots,list(g))
    
  } else {
    print("Skipping Numbat (no CNVs found).")
  }

  print("Evaluate cells from Numbat (Expr)")
  numbat_cells <- read_numbat_expr_individual_cells(input_numbat_obj)
  numbat_res_df<- evaluate_metrics_per_cell_paired(wgs_cell_range,wgs_cell,numbat_cells[[1]],
                                            numbat_cells[[2]])
  numbat_res_df$method<-"Numbat (Expr)"
  numbat_res_df$cells<-rownames(numbat_res_df)
  total_res_df<-rbind(total_res_df,numbat_res_df)
  
  #Plot also the karyogram
  g<-plot_karyogram_paired(wgs_cell_range,wgs_cell,numbat_cells[[1]],
                           numbat_cells[[2]],"Numbat (Expr)")
  karyoplots <- c(karyoplots,list(g))
  

} else {
  print("Skipping Numbat (no output file provided).")  
}

print("Evaluate cells from CONICSmat")
CONICSmat_cells <- read_CONICSmat_individual_cells(input_CONICSmat_chrom_pos,
                                                   input_CONICSmat_cnv,
                                  input_annot, input_ref_groups)

CONICSmat_res_df<- evaluate_metrics_per_cell_paired(wgs_cell_range,wgs_cell,
                                                    CONICSmat_cells[[1]],CONICSmat_cells[[2]])
CONICSmat_res_df$method<-"CONICSmat"
CONICSmat_res_df$cells<-rownames(CONICSmat_res_df)
total_res_df<-rbind(total_res_df,CONICSmat_res_df)

#Plot also the karyogram
g<-plot_karyogram_paired(wgs_cell_range,wgs_cell,CONICSmat_cells[[1]],
                         CONICSmat_cells[[2]],"CONICSmat")
karyoplots <- c(karyoplots,list(g))

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

pdf(output_plot_karyograms, width = 12, height = 8)
for(kplot in karyoplots){
  grid::grid.draw(kplot)
}
dev.off()

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------

sessionInfo()





  

    
  
  

