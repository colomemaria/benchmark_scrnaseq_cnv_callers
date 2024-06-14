# ------------------------------------------------------------------------------
# Evaluate the results of the different CNV callers compared to a ground truth
# dataset - specific version of the script where the methods are compared that
# have the options to automatically identify the reference vs tumor cells
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
library(viridis)
library(numbat)

theme_set(theme_bw())

#Source script with help functions
source("scripts/lib.R")

# ------------------------------------------------------------------------------
print("Get input parameters from snakemake")
# ------------------------------------------------------------------------------

input_annot<-snakemake@input$annot
input_ref_groups<-snakemake@input$ref_groups
input_wes<-snakemake@input$wes

input_copykat_wref<-snakemake@input$copykat_wref
input_copykat_woref<-snakemake@input$copykat_woref

input_scevan_clone1_wref<-snakemake@input$scevan_clone1_wref
input_scevan_annot_wref<-snakemake@input$scevan_annot_wref
input_scevan_clone1_woref<-snakemake@input$scevan_clone1_woref
input_scevan_annot_woref<-snakemake@input$scevan_annot_woref

input_numbat_obj_wref<-dirname(snakemake@input$numbat_obj_wref)
input_numbat_obj_woref<-dirname(snakemake@input$numbat_obj_woref)

input_CONICSmat_chrom_pos<-snakemake@input$CONICSmat_chr_pos
input_CONICSmat_cnv_wref<-snakemake@input$CONICSmat_cnv_wref
input_CONICSmat_cnv_woref<-snakemake@input$CONICSmat_cnv_woref  

output_merged_results<-snakemake@output$merged_results 

output_plot_heatmap<-snakemake@output$plot_heatmap
output_text_corr<-snakemake@output$text_corr
output_plot_corr<-snakemake@output$plot_corr
output_text_auc<-snakemake@output$text_auc
output_plot_auc<-snakemake@output$plot_auc
output_text_f1<-snakemake@output$text_f1
output_plot_f1<-snakemake@output$plot_f1
output_plot_heatmap_binf1<-snakemake@output$plot_heatmap_binf1

# ------------------------------------------------------------------------------
print("Load the results of the different methods, calculate a pseudbulk version
      for each and combine them in one GRange object based on WES bins")
# ------------------------------------------------------------------------------

print("Load WES data")
wes_results<-read_wes_results_gatk(input_wes)

print("Load copykat data - with ref")
copykat_results<-read_copykat(input_copykat_wref,input_annot,input_ref_groups)
colnames(mcols(copykat_results[[1]]))<-c("copykat_wref")
  
combined_range<-combine_range_objects(wes_results,copykat_results[[1]],
                                      method_colname="copykat_wref")

print("Load copykat data - without ref")
copykat_results<-read_copykat(input_copykat_woref,input_annot,input_ref_groups)
colnames(mcols(copykat_results[[1]]))<-c("copykat_woref")

combined_range<-combine_range_objects(combined_range,copykat_results[[1]],
                                      method_colname="copykat_woref")


print("Load SCEVAN data - with ref")
binned_genome<-combined_range
elementMetadata(binned_genome)<-NULL
scevan_results<-read_scevan_cn_status(input_scevan_clone1_wref, input_scevan_annot_wref,
                                      binned_genome)
colnames(mcols(scevan_results))<-c("scevan_cnv_wref")
combined_range<-combine_range_objects(combined_range,scevan_results,
                                      method_colname="scevan_cnv_wref")

print("Load SCEVAN data - without ref")
binned_genome<-combined_range
elementMetadata(binned_genome)<-NULL
scevan_results<-read_scevan_cn_status(input_scevan_clone1_woref, input_scevan_annot_woref,
                                      binned_genome)
colnames(mcols(scevan_results))<-c("scevan_cnv_woref")
combined_range<-combine_range_objects(combined_range,scevan_results,
                                      method_colname="scevan_cnv_woref")


print("Load Numbat results - with ref")
numbat_results <- read_numbat_cnv(input_numbat_obj_wref)
colnames(mcols(numbat_results))<-c("numbat_wref")
combined_range<-combine_range_objects(combined_range, numbat_results,
                                      method_colname="numbat_wref")

print("Load Numbat results - without ref")
numbat_results <- read_numbat_cnv(input_numbat_obj_woref)
colnames(mcols(numbat_results))<-c("numbat_woref")
combined_range<-combine_range_objects(combined_range, numbat_results,
                                      method_colname="numbat_woref")

print("Load CONICSmat results - with ref")
CONICSmat_results <- read_CONICSmat(input_CONICSmat_chrom_pos, 
                                    input_CONICSmat_cnv_wref,
                                    input_annot, input_ref_groups)
colnames(mcols(CONICSmat_results))<-c("CONICSmat_wref")  
combined_range<-combine_range_objects(combined_range, CONICSmat_results,
                                      method_colname="CONICSmat_wref")

print("Load CONICSmat results - without ref")
CONICSmat_results <- read_CONICSmat(input_CONICSmat_chrom_pos, 
                                    input_CONICSmat_cnv_woref,
                                    input_annot, input_ref_groups)
colnames(mcols(CONICSmat_results))<-c("CONICSmat_woref")  
combined_range<-combine_range_objects(combined_range, CONICSmat_results,
                                      method_colname="CONICSmat_woref")

#Save the combined result data frame for later processing
write.table(combined_range,file=output_merged_results,
            quote = FALSE,row.names=FALSE,sep="\t")

# ------------------------------------------------------------------------------
print("Combined line plot / heatmap")
# ------------------------------------------------------------------------------

#Vector for renaming methods (official published names) and specifying
#the order in the plot
method_names<-setNames(c("WES","copyKat","copyKat (wo ref)",
                         "SCEVAN","SCEVAN (wo ref)",
                         "Numbat", "Numbat (wo ref)", 
                         "CONICSmat", "CONICSmat (wo ref)"),
                       c("wgs_mean","copykat_wref","copykat_woref",
                         "scevan_cnv_wref","scevan_cnv_woref",
                         "numbat_wref","numbat_woref", 
                         "CONICSmat_wref","CONICSmat_woref"))

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
scaled_methods<-combined_methods
for(method in names(method_names)){
  scaled_methods[,method]<-(scaled_methods[,method]-2) / 
    sd(scaled_methods[,method])
}

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
       width=20,height=10)

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
