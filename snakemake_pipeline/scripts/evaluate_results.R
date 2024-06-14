# ------------------------------------------------------------------------------
# Evaluate the results of the different CNV callers compared to a ground truth
# dataset (typically WGS)
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
  
input_hb_obj<-snakemake@input$hb_obj

input_numbat_obj<-snakemake@input$numbat_obj

input_CONICSmat_filtered_dims<-snakemake@input$CONICSmat_filtered_dims
input_CONICSmat_chrom_pos<-snakemake@input$CONICSmat_chr_pos
input_CONICSmat_cnv<-snakemake@input$CONICSmat_cnv

#Parameter how to interpret the WGS/WES data (default set to WGS)
param_genomic_format<-snakemake@params$genomic_format
if(is.null(param_genomic_format)){
  param_genomic_format<-"WGS"
}

#Parameter for plots - specifying type of genomic data
param_genomic_name<-snakemake@params$genomic_name
if(is.null(param_genomic_name)){
  param_genomic_name<-"WGS"
}

#Output files
output_merged_results<-snakemake@output$merged_results
output_text_filter<-snakemake@output$text_filter
output_plot_filter<-snakemake@output$plot_filter
output_plot_heatmap<-snakemake@output$plot_heatmap
output_text_scaling<-snakemake@output$text_scaling
output_text_corr<-snakemake@output$text_corr
output_plot_corr<-snakemake@output$plot_corr
output_text_auc<-snakemake@output$text_auc
output_plot_auc<-snakemake@output$plot_auc
output_text_f1<-snakemake@output$text_f1
output_plot_f1<-snakemake@output$plot_f1
output_plot_heatmap_binf1<-snakemake@output$plot_heatmap_binf1

# ------------------------------------------------------------------------------
print("Load the results of the different methods, calculate a pseudbulk version
      for each and combine them in one GRange object based on WGS bins")
# ------------------------------------------------------------------------------

# Save the number of evaluated genes per method
num_filtered_res<-NULL

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

print("Load inferCNV data (pseudobulk CNVs)")
infercnv_results<-read_infercnv_6state_model(input_infercnv_cnv,
                                             input_infercnv_gene_pos)

num_filtered_res<-rbind(num_filtered_res,
                        data.frame(method="inferCNV",
                                  annotated_genes=length(infercnv_results[[1]]),
                                  annotated_cells=infercnv_results[[2]]))

combined_range<-combine_range_objects(wgs_results,infercnv_results[[1]],
                                      method_colname="infercnv_cnv")

print("Load inferCNV data (normalized expression)")
infercnv_results<-read_infercnv_expr(input_infercnv_expr,
                                     input_infercnv_gene_pos)

combined_range<-combine_range_objects(combined_range,infercnv_results[[1]],
                                      method_colname="infercnv_expr")

#Check if Casper output files are provided
if(! is.null(input_casper)){
  print("Load casper results (filtered segments, aggregated per gene)")
  casper_results<-read_casper(input_casper)
  
  num_filtered_res<-rbind(num_filtered_res,
                          data.frame(method="CaSpER",
                                     annotated_genes=length(casper_results),
                                     annotated_cells=NA))
  
  combined_range<-combine_range_objects(combined_range,casper_results,
                                        method_colname="casper")
} else {
  print("Skipping casper (no output file provided).")
}

print("Load copykat data (normalized expression)")
copykat_results<-read_copykat(input_copykat,input_annot,input_ref_groups)

num_filtered_res<-rbind(num_filtered_res,
                        data.frame(method="copyKat",
                                   annotated_genes=length(copykat_results[[1]]),
                                   annotated_cells=copykat_results[[2]]))

combined_range<-combine_range_objects(combined_range,copykat_results[[1]],
                                      method_colname="copykat")

print("Load SCEVAN data (normalized expression)")
scevan_results<-read_scevan(input_scevan_expr,input_scevan_gene_pos,
                            input_annot,input_ref_groups)

num_filtered_res<-rbind(num_filtered_res,
                        data.frame(method="SCEVAN",
                                   annotated_genes=length(scevan_results[[1]]),
                                   annotated_cells=scevan_results[[2]]))

combined_range<-combine_range_objects(combined_range,scevan_results[[1]],
                                      method_colname="scevan")

print("Load SCEVAN data (CNV output)")
binned_genome<-combined_range
elementMetadata(binned_genome)<-NULL
scevan_results<-read_scevan_cn_status(input_scevan_clone1, input_scevan_annot,
                                      binned_genome)

combined_range<-combine_range_objects(combined_range,scevan_results,
                                      method_colname="scevan_cnv")

#Check if Numbat output files are provided
if(! is.null(input_numbat_obj)){
  
  #Get the directory name
  input_numbat_obj<-dirname(input_numbat_obj)
  
  print("Load Numbat results (expr)")
  
  numbat_results <- read_numbat_expr(input_numbat_obj)
  
  num_filtered_res<-rbind(num_filtered_res,
                          data.frame(method="Numbat",
                                     annotated_genes=length(numbat_results[[1]]),
                                     annotated_cells=numbat_results[[2]]))
  
  combined_range<-combine_range_objects(combined_range, numbat_results[[1]],
                                        method_colname="numbat")
  
  print("Load Numbat results (CNVs)")
  numbat_results <- read_numbat_cnv(input_numbat_obj)
  
  combined_range<-combine_range_objects(combined_range, numbat_results,
                                        method_colname="numbat_cnv")
} else {
  print("Skipping Numbat (no output file provided).")  
}

#Check if CONICSmat output files are provided
if(! is.null(input_CONICSmat_cnv)){
  print("Load CONICSmat results")
  CONICSmat_results <- read_CONICSmat(input_CONICSmat_chrom_pos, input_CONICSmat_cnv,
                                      input_annot, input_ref_groups)
    
  combined_range<-combine_range_objects(combined_range, CONICSmat_results,
                                        method_colname="CONICSmat")
  
  tmp<-fread(input_CONICSmat_filtered_dims)
  num_filtered_res<-rbind(num_filtered_res,
                          data.frame(method="CONICSmat",
                                     annotated_genes=tmp$ngenes[1],
                                     annotated_cells=tmp$ncells[1]))
  
} else {
  print("Skipping CONICSmat (no output file provided).")  
}


#Save the combined result data frame for later processing
write.table(combined_range,file=output_merged_results,
            quote = FALSE,row.names=FALSE,sep="\t")

# ------------------------------------------------------------------------------
print("Create filtering plot")
# ------------------------------------------------------------------------------

#Save the filtering results
write.table(num_filtered_res,file=output_text_filter,sep="\t",quote=FALSE,
            row.names=FALSE)

#Create plots of filtered cells and genes (or rather: how many are left)
g.1<-ggplot(num_filtered_res,aes(x=method,fill=method,y=annotated_cells))+
  geom_bar(stat="identity")+
  theme(legend.position="none",
        axis.text.x = element_text(angle = 45,vjust=1,hjust=1))+
  ylab("Analysed cells")+xlab("CNV Method")

g.2<-ggplot(num_filtered_res,aes(x=method,fill=method,y=annotated_genes))+
  geom_bar(stat="identity")+
  theme(legend.position="none",
        axis.text.x = element_text(angle = 45,vjust=1,hjust=1))+
  ylab("Annotated genes")+xlab("CNV Method")

g<-ggarrange(g.1,g.2,ncol=2)
ggsave(g,file=output_plot_filter,width=8,height=4)

# ------------------------------------------------------------------------------
print("Combined line plot / heatmap")
# ------------------------------------------------------------------------------

#Vector for renaming methods (official published names) and specifying
#the order in the plot
method_names<-setNames(c(param_genomic_name,"HoneyBADGER (CNV)","HoneyBADGER (Expr)",
                         "InferCNV (CNV)","InferCNV (Expr)","CaSpER",
                         "copyKat","SCEVAN (CNV)","SCEVAN (Expr)",
                         "Numbat (Expr)", "Numbat (CNV)", "CONICSmat"),
                       c("wgs_mean","hb_cnv","honeybadger",
                         "infercnv_cnv","infercnv_expr","casper",
                         "copykat","scevan_cnv","scevan",
                         "numbat","numbat_cnv", "CONICSmat"))

#Get results
combined_methods<-elementMetadata(combined_range)

#Reduce the method names to the ones that are really analyzed here
method_names<-method_names[names(method_names) %in% colnames(combined_methods)]

#Add position information
combined_methods$chr<-factor(combined_range@seqnames,
                             levels=combined_range@seqnames@values)
combined_methods$start_position<-combined_range@ranges@start
combined_methods<-as.data.frame(combined_methods)

#Add an artifical count through the whole genome and 
#get start positions for each new chromosome
combined_methods$counted_pos<-1:nrow(combined_methods)
chr_boundries<-combined_methods%>%
  group_by(chr)%>%
  summarize(start_chr=min(counted_pos),
            mean_chr=mean(counted_pos))

#Scale every dataset to have diploid values at 0 and a standard deviation of 1
scaling_factor<-NULL
scaled_methods<-combined_methods
for(method in names(method_names)){
  scaled_methods[,method]<-(scaled_methods[,method]-2) / 
    sd(scaled_methods[,method])
  
  #Save the standard deviation of each method to document the chosen normalization factor
  scaling_factor<-rbind(scaling_factor,
                        data.frame(method,
                                   sd=sd(combined_methods[,method]),
                                   mean=mean(combined_methods[,method]-2)))
}

write.table(scaling_factor,file=output_text_scaling,sep="\t",quote=FALSE,
            row.names=FALSE)

plot_data<-reshape2::melt(scaled_methods,
                          id.vars=c("chr","start_position","counted_pos"))

#Rename variable names
plot_data$variable<-method_names[as.character(plot_data$variable)]

#Order them respectively
plot_data$variable<-factor(plot_data$variable,
                           levels=method_names)


method_colors<-c('#1f77b4', '#ff7f0e', '#279e68', '#d62728', '#aa40fc', '#8c564b',
                 '#e377c2', '#b5bd61', '#17becf', '#aec7e8', '#ffbb78', '#98df8a',
                 '#ff9896', '#c5b0d5', '#c49c94', '#f7b6d2')

# #Combine heatmap and line plot
# g.1<-ggplot(plot_data,aes(x=counted_pos,y=value,color=variable))+geom_line()+
#   geom_vline(xintercept = chr_boundries$start_chr)+
#   ylab("Normalized score")+
#   scale_color_manual("Method",values=method_colors)+
#   scale_x_continuous(breaks=chr_boundries$mean_chr,
#                      labels=chr_boundries$chr)+
#   coord_cartesian(xlim=c(1,max(plot_data$counted_pos)),expand=FALSE)+
#   theme_bw()+
#   theme(legend.position="top",
#         #legend.key.size = unit(2,"line"),
#         axis.text.x=element_blank(),
#         axis.title.x=element_blank(),
#         text=element_text(size=21))

g.2<-ggplot(plot_data,aes(x=counted_pos,y=variable,fill=value))+geom_tile()+
  theme_bw()+
  scale_fill_gradient2("Score",low = "darkblue",
                        mid = "white",high = "darkred",midpoint = 0,
                       breaks=c(-5,0,5))+
  xlab("Chromosome position")+ylab("Method")+
  geom_vline(xintercept = chr_boundries$start_chr)+
  scale_x_continuous(breaks=chr_boundries$mean_chr,
                     labels=chr_boundries$chr)+
  scale_y_discrete(limits=rev)+
  coord_cartesian(xlim=c(1,max(plot_data$counted_pos)),expand=FALSE)+
  theme(legend.position="bottom",
        axis.text.x = element_text(angle=90,vjust=0.5,hjust=1),
        text=element_text(size=21))
#g<-ggarrange(g.1,g.2,ncol=1,align="v",heights=c(0.3,0.7))
ggsave(g.2, file=output_plot_heatmap,
       width=20,height=8)

# ------------------------------------------------------------------------------
print("Correlation calculation and plot")
# ------------------------------------------------------------------------------

#Create a heatmap which methods are how similar
compare_methods<-NULL
for(m1 in 1:(length(method_names)-1)){
  for(m2 in (m1+1):length(method_names)){
    res<-evaluate_correlation_colnames(combined_methods,
                                       names(method_names)[m1],
                                       names(method_names)[m2],
                                       printres=FALSE)
    compare_methods<-rbind(compare_methods,
                           data.frame(method1=method_names[m1],
                                      method2=method_names[m2],
                                      pearson=res[1],
                                      spearman=res[2],
                                      kendell=res[3],
                                      mse=res[4]))
  }
  
  #Add the diagonal for easier visualization
  compare_methods<-rbind(compare_methods,
                         data.frame(method1=method_names[m1],
                                    method2=method_names[m1],
                                    pearson=1,
                                    spearman=1,
                                    kendell=1,
                                    mse=0))
}

#Add the diagonal for the last element
compare_methods<-rbind(compare_methods,
                       data.frame(method1=method_names[m2],
                                  method2=method_names[m2],
                                  pearson=1,
                                  spearman=1,
                                  kendell=1,
                                  mse=0))

#Save result txt
write.table(compare_methods,file=output_text_corr,sep="\t",quote=FALSE,
            row.names=FALSE)

#Create plots
compare_methods$method1<-factor(compare_methods$method1,
                                levels=method_names)
compare_methods$method2<-factor(compare_methods$method2,
                                levels=method_names)

g<-ggplot(compare_methods,aes(x=method2,y=method1,fill=pearson))+
  geom_tile()+
  theme_bw()+
  geom_text(aes(label=round(pearson,3),
                color=ifelse(pearson<0.6,'white','black')),size=2)+
  scale_color_manual(values=c("black","white"))+
  scale_y_discrete(limits=rev)+
  xlab("Method")+
  ylab("Method")+
  scale_fill_viridis("Pearson\ncorrelation",limits=c(0,1))+
  guides(color="none")+
  theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1))

ggsave(g, file=output_plot_corr,
       width=6,height=4.5)

# ------------------------------------------------------------------------------
print("AUC evaluation and plots")
# ------------------------------------------------------------------------------

bins_res<-NULL
for(meth in setdiff(names(method_names),"wgs_mean")){
  
  #Get AUC and AUPRC values (looking at the complete curve)
  tmp<-evaluate_gains_loss_auc(combined_methods,meth,"wgs_mean")
  
  #Get truncated AUC values (only in the biological fair range)
  tmp2<-evaluate_gains_loss_auc_truncated(combined_methods,meth,"wgs_mean")
  
  bins_res<-rbind(bins_res,
                  data.frame(method=meth,
                             auc_gains=tmp[1],
                             auc_gains_trunc = tmp2[1],
                             aucpr_gains=tmp[2],
                             auc_losses=tmp[3],
                             auc_losses_trunc = tmp2[2],
                             aucpr_losses=tmp[4]))
}

#Save result txt
write.table(bins_res,file=output_text_auc,sep="\t",quote=FALSE,
            row.names=FALSE)

#Create barplot with all metrics
plot_bins_res<-melt(bins_res,id.vars="method")
plot_bins_res$method<-method_names[as.character(plot_bins_res$method)]
plot_bins_res$method<-factor(plot_bins_res$method,levels=method_names)

#Replace text values
metrics<-setNames(c("AUC (Gains)","AUC trunc (Gains)","AUCPR (Gains)",
                    "AUC (Losses)","AUC trunc (Losses)","AUCPR (Losses)"),
                  c("auc_gains","auc_gains_trunc","aucpr_gains",
                    "auc_losses","auc_losses_trunc", "aucpr_losses"))
plot_bins_res$variable<-metrics[plot_bins_res$variable]
plot_bins_res$variable<-factor(plot_bins_res$variable,
                               levels=metrics)

g<-ggplot(plot_bins_res,aes(x=method,y=value,fill=method))+
  geom_bar(stat="identity")+
  facet_wrap(~variable,ncol=3)+theme_bw()+
  xlab("CNV method")+ylab("Metric")+ylim(0,1)+
  theme(legend.position="none",
        axis.text.x = element_text(angle=65,vjust=1,hjust=1))

ggsave(g, file=output_plot_auc,
       width=12,height=6)

# ------------------------------------------------------------------------------
print("Sensitivity and specificity with F1 score optimal thresholds")
# ------------------------------------------------------------------------------

#Save the categorical CNV predictions (for plotting later)
binarized_methods<-combined_methods
binarized_methods$wgs_mean<-ifelse(binarized_methods$wgs_mean<1.5,1,
                                   ifelse(binarized_methods$wgs_mean>2.5,3,2))
bins_res<-NULL
for(meth in setdiff(names(method_names),"wgs_mean")){
  
  print(paste("Evaluating method -",meth,"- starting at",Sys.time()))
  
  tmp<-evaluate_optimal_f1_score(binarized_methods,meth,"wgs_mean")
  
  bins_res<-rbind(bins_res,
                  data.frame(method=meth,
                             max_f1=tmp[1],
                             cutoff_f1_gain=tmp[2],
                             sens_gains=tmp[3],
                             prec_gains=tmp[4],
                             cutoff_f1_loss=tmp[5],
                             sens_losses=tmp[6],
                             prec_losses=tmp[7]))
  
  binarized_methods[,meth]<-ifelse(binarized_methods[,meth]<tmp[5],1,
                                   ifelse(binarized_methods[,meth]>tmp[2],3,2))
}

#Save result txt
write.table(bins_res,file=output_text_f1,sep="\t",quote=FALSE,
            row.names=FALSE)

#Create barplot with all metrics
bins_res<-bins_res[,! (colnames(bins_res) %in% c("cutoff_f1_gain","cutoff_f1_loss"))]
plot_bins_res<-melt(bins_res,id.vars="method")
plot_bins_res$method<-method_names[as.character(plot_bins_res$method)]
plot_bins_res$method<-factor(plot_bins_res$method,levels=method_names)

#Replace text values
metrics<-setNames(c("F1 Score","Sensitivity (Gains)","Precision (Gains)",
                    "Sensitivity (Losses)","Precision (Losses)"),
                  c("max_f1","sens_gains","prec_gains",
                    "sens_losses","prec_losses"))
plot_bins_res$variable<-metrics[plot_bins_res$variable]
plot_bins_res$variable<-factor(plot_bins_res$variable,
                               levels=metrics)

g<-ggplot(plot_bins_res,aes(x=method,y=value,fill=method))+
  geom_bar(stat="identity")+
  facet_wrap(~variable)+theme_bw()+
  xlab("CNV method")+ylab("Metric")+ylim(0,1)+
  theme(legend.position="none",
        axis.text.x = element_text(angle=65,vjust=1,hjust=1))

ggsave(g, file=output_plot_f1,
       width=10,height=6)


# ------------------------------------------------------------------------------
print("Plot karyogram with F1 score optimal thresholds")
# ------------------------------------------------------------------------------

plot_data<-reshape2::melt(binarized_methods,
                          id.vars=c("chr","start_position","counted_pos"))

#Rename variable names
plot_data$variable<-method_names[as.character(plot_data$variable)]

#Order them respectively
plot_data$variable<-factor(plot_data$variable,
                           levels=method_names)

#Convert values to CNVs
plot_data$value<-ifelse(plot_data$value == 2,"normal",
                        ifelse(plot_data$value == 1,"loss","gain"))
plot_data$value<-factor(plot_data$value,levels=c("loss","normal","gain"))

g<-ggplot(plot_data,aes(x=counted_pos,y=variable,fill=value))+geom_tile()+
  theme_bw()+
  scale_fill_manual("CNV", values=c("darkblue","lightgreen","darkred"))+
  xlab("Chromosome position")+ylab("Method")+
  geom_vline(xintercept = chr_boundries$start_chr)+
  scale_x_continuous(breaks=chr_boundries$mean_chr,
                     labels=chr_boundries$chr)+
  scale_y_discrete(limits=rev)+
  coord_cartesian(xlim=c(1,max(plot_data$counted_pos)),expand=FALSE)+
  theme(legend.position="bottom",
        axis.text.x = element_text(angle=90,vjust=0.5,hjust=1),
        text=element_text(size=18))

ggsave(g, file=output_plot_heatmap_binf1,
       width=20,height=8)

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()
