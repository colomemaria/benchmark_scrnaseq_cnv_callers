# ------------------------------------------------------------------------------
# Evaluate the results for a diploid sample 
# (so in theory everything should be 0)
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

input_infercnv_cnv<-snakemake@input$infercnv_cnv
input_infercnv_expr<-snakemake@input$infercnv_expr
input_infercnv_gene_pos<-snakemake@input$infercnv_gene_pos

input_casper<-snakemake@input$casper
  
input_copykat<-snakemake@input$copykat

input_scevan_expr<-snakemake@input$scevan
input_scevan_gene_pos<-snakemake@input$scevan_gene_pos
input_scevan_clone1<-snakemake@input$scevan_clone1
input_scevan_annot<-snakemake@input$scevan_annot

input_numbat_obj<-dirname(snakemake@input$numbat_obj)

input_CONICSmat_chrom_pos<-snakemake@input$CONICSmat_chr_pos
input_CONICSmat_cnv<-snakemake@input$CONICSmat_cnv

input_scaling_ref<-snakemake@input$scaling_ref

#Output files
output_merged_results<-snakemake@output$merged_results

output_plot_heatmap<-snakemake@output$plot_heatmap
output_text_rmse<-snakemake@output$text_rmse
output_plot_rmse<-snakemake@output$plot_rmse

# ------------------------------------------------------------------------------
print("Load the results of the different methods, calculate a pseudbulk version
      for each and combine them in one GRange object")
# ------------------------------------------------------------------------------

print("Load inferCNV data (pseudobulk CNVs)")
infercnv_results<-read_infercnv_6state_model(input_infercnv_cnv,
                                             input_infercnv_gene_pos)

binned_genome<-create_binned_genome(infercnv_results[[1]])
combined_range<-combine_range_objects(binned_genome,infercnv_results[[1]],
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
  
  combined_range<-combine_range_objects(combined_range,casper_results,
                                        method_colname="casper")
} else {
  print("Skipping casper (no output file provided).")
}

print("Load copykat data (normalized expression)")
copykat_results<-read_copykat(input_copykat,input_annot,input_ref_groups)

combined_range<-combine_range_objects(combined_range,copykat_results[[1]],
                                      method_colname="copykat")

print("Load SCEVAN data (normalized expression)")
scevan_results<-read_scevan(input_scevan_expr,input_scevan_gene_pos,
                            input_annot,input_ref_groups)

combined_range<-combine_range_objects(combined_range,scevan_results[[1]],
                                      method_colname="scevan")

print("Load SCEVAN data (CNV output)")
binned_genome<-combined_range
elementMetadata(binned_genome)<-NULL
scevan_results<-read_scevan_cn_status(input_scevan_clone1, input_scevan_annot,
                                      binned_genome)

combined_range<-combine_range_objects(combined_range,scevan_results,
                                      method_colname="scevan_cnv")

print("Load Numbat results (expr)")
numbat_results <- read_numbat_expr(input_numbat_obj)

combined_range<-combine_range_objects(combined_range, numbat_results[[1]],
                                      method_colname="numbat")

print("Load Numbat results (CNVs)")
numbat_results <- read_numbat_cnv(input_numbat_obj)

#In case no CNVs were found, set everything to 2
if(is.null(numbat_results)){
  combined_range$numbat_cnv<-2
} else {
  combined_range<-combine_range_objects(combined_range, numbat_results,
                                        method_colname="numbat_cnv")
}

if(! is.null(input_CONICSmat_cnv)){
  print("Load CONICSmat results")
  CONICSmat_results <- read_CONICSmat(input_CONICSmat_chrom_pos, input_CONICSmat_cnv,
                                      input_annot, input_ref_groups)
    
  combined_range<-combine_range_objects(combined_range, CONICSmat_results,
                                        method_colname="CONICSmat")
} else {
  print("Skipping CONICSmat (no output file provided).")
}

#Save the combined result data frame for later processing
write.table(combined_range,file=output_merged_results,
            quote = FALSE,row.names=FALSE,sep="\t")

# ------------------------------------------------------------------------------
print("Combined line plot / heatmap")
# ------------------------------------------------------------------------------

#Vector for renaming methods (official published names) and specifying
#the order in the plot
method_names<-setNames(c("scWGS","InferCNV (CNV)","InferCNV (Expr)","CaSpER",
                         "copyKat","SCEVAN (CNV)","SCEVAN (Expr)",
                         "Numbat (Expr)", "Numbat (CNV)", "CONICSmat"),
                       c("wgs_mean","infercnv_cnv","infercnv_expr","casper",
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

#Load table with external scaling factors
scaling_table<-fread(input_scaling_ref)

#Scale every dataset (use an external reference to not falsify the current results)
scaled_methods<-combined_methods
for(method in names(method_names)){
  scaled_methods[,method]<-(scaled_methods[,method]-2) /
    scaling_table$sd[scaling_table$method==method]
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

#Combine heatmap and line plot
g.1<-ggplot(plot_data,aes(x=counted_pos,y=value,color=variable))+geom_line()+
  geom_vline(xintercept = chr_boundries$start_chr)+
  ylab("Normalized score")+
  scale_color_manual("Method",values=method_colors)+
  scale_x_continuous(breaks=chr_boundries$mean_chr,
                     labels=chr_boundries$chr)+
  coord_cartesian(xlim=c(1,max(plot_data$counted_pos)),expand=FALSE)+
  theme_bw()+
  theme(legend.position="top",
        #legend.key.size = unit(2,"line"),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        text=element_text(size=21))

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
g<-ggarrange(g.1,g.2,ncol=1,align="v",heights=c(0.3,0.7))
ggsave(g, file=output_plot_heatmap,
       width=20,height=10)

# ------------------------------------------------------------------------------
print("Calculate MSE (using data normalized with standard deviation from SNU601)")
# ------------------------------------------------------------------------------

metrics<-NULL
for(method in names(method_names)){
  
  metrics<-rbind(metrics,
                 data.frame(method,
                            rmse=sqrt(mean((scaled_methods[,method])^2))))
}

#Save result txt
write.table(metrics,file=output_text_rmse,sep="\t",quote=FALSE,
            row.names=FALSE)

g<-ggplot(metrics,aes(x=method,y=rmse,fill=method))+
  geom_bar(stat="identity")+
  theme_bw()+
  xlab("CNV method")+ylab("RMSE")+
  theme(legend.position="none",
        axis.text.x = element_text(angle=65,vjust=1,hjust=1))

ggsave(g, file=output_plot_rmse,
       width=8,height=6)

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()
