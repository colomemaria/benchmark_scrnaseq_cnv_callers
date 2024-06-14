# ------------------------------------------------------------------------------
# In depth evaluation of CNV predictions from InferCNV compared to ground truth
# - check how many genes are mapped to each bin
# - check if prediction performance is associated with gene expression level
# ------------------------------------------------------------------------------

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# ------------------------------------------------------------------------------
print("Load libraries")
# ------------------------------------------------------------------------------

library(Seurat)
library(data.table)
library(dplyr)
library(GenomicRanges)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(viridis)
library(clusterProfiler) #for GO enrichment
library(org.Hs.eg.db) #for GO enrichment

theme_set(theme_bw())

#Source script with help functions
source("scripts/lib.R")

# ------------------------------------------------------------------------------
print("Get input parameters from snakemake")
# ------------------------------------------------------------------------------

input_matrix<-snakemake@input$matrix
input_annot<-snakemake@input$annot
input_ref_groups<-snakemake@input$ref_groups

input_wgs<-snakemake@input$scwgs_truth
input_infercnv_cnv<-snakemake@input$infercnv_cnv
input_infercnv_gene_pos<-snakemake@input$infercnv_gene_pos

output_gene_dist<-snakemake@output$gene_dist
output_plot_count_dist<-snakemake@output$plot_count_dist
output_plot_count_norm<-snakemake@output$plot_count_norm
output_plot<-snakemake@output$overview_plot

output_boxplot_loss<-snakemake@output$boxplot_loss
output_boxplot_loss_lfc<-snakemake@output$boxplot_loss_lfc
output_boxplot_gain<-snakemake@output$boxplot_gain
output_boxplot_gain_lfc<-snakemake@output$boxplot_gain_lfc


# ------------------------------------------------------------------------------
print("Load the original count matrix and scWGS and RNA method results")
# ------------------------------------------------------------------------------

#Load scWGS data (already preprocessed to pseudobulk)
wgs_results<-read_wgs_results(input_wgs)

#Process inferCNV data (pseudobulk CNVs)
infercnv_results<-read_infercnv_6state_model(input_infercnv_cnv,
                                             input_infercnv_gene_pos)

# ------------------------------------------------------------------------------
print("Summary of genes per bin")
# ------------------------------------------------------------------------------

#Check overlaps (to get the bin sizes to a similar shape)
overlaps<-as.data.frame(findOverlaps(wgs_results,infercnv_results[[1]],type="any"))

genes_per_bin<-overlaps%>%
  group_by(queryHits)%>%
  summarize(n_genes=n())

print(summary(genes_per_bin$n_genes))

g<-ggplot(genes_per_bin,aes(x=n_genes))+
  geom_histogram(bins=30)+
  ylab("Number bins")+
  xlab("Number genes per bin")
ggsave(g,filename=output_gene_dist,width=4,height=3)

# ------------------------------------------------------------------------------
print("Extracting expression level of each gene")
# ------------------------------------------------------------------------------

#Load count matrix and
data_matrix<-fread(input_matrix)
#Format into a matrix
gene_names<-data_matrix$V1
data_matrix$V1<-NULL
data_matrix<-as.matrix(data_matrix)
rownames(data_matrix)<-gene_names

#Extract reference cells
annotation<-fread(input_annot, header=FALSE)
ref_groups<-read.table(input_ref_groups,header=TRUE)
ref_cells <- annotation$V1[annotation$V2 %in% ref_groups$ref_groups]

#Alternative approach: normalize counts before the analysis using Seurat
annotation<-as.data.frame(annotation)
rownames(annotation)<-annotation$V1
colnames(annotation)<-c("barcode","sample")

#Check raw UMI count distribution
#all(annotation$barcode==colnames(data_matrix))
annotation$raw_counts<-colSums(data_matrix)
annotation$group<-ifelse(annotation$sample %in% ref_groups$ref_groups,
                          "control","cancer")
g<-ggplot(annotation,aes(x=raw_counts,fill=group))+
  geom_histogram()+
  xlab("Raw UMI counts")
ggsave(g,filename=output_plot_count_dist,height=3,width=5)

#Normalize counts
seurat_obj <- CreateSeuratObject(counts = data_matrix, meta.data = annotation)
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize",
                            scale.factor = 10000)

#Check distribution after normalization
norm_counts<-seurat_obj@assays$RNA@data
#all(annotation$barcode==colnames(norm_counts))
annotation$norm_counts<-colSums(norm_counts)
  
g<-ggplot(annotation,aes(x=norm_counts,fill=group))+
  geom_histogram()+
  xlab("Normalized UMI counts")
ggsave(g,filename=output_plot_count_norm,height=3,width=5)

#Get mean expression in total (controls and cases)
mean_expr_total<-rowMeans(seurat_obj@assays$RNA@data)

#Get mean expression separately for cases and controls
seurat_control<-subset(seurat_obj,subset = sample %in% ref_groups$ref_groups)
mean_expr_controls<-rowMeans(seurat_control@assays$RNA@data)

ref_matrix <- as.matrix(seurat_control@assays$RNA@data)
seurat_cancer<-subset(seurat_obj,cells = 
                        annotation$barcode[!annotation$sample %in% ref_groups$ref_groups])
mean_expr_cancer <- rowMeans(seurat_cancer@assays$RNA@data)

#Combine all information in one df
expr_combined<-data.frame(genes=names(mean_expr_total),
                          mean_total=mean_expr_total,
                          mean_cancer=mean_expr_cancer,
                          mean_controls=mean_expr_controls)
rownames(expr_combined)<-expr_combined$genes

#Order the data frame corresponding to the inferCNV results
expr_combined<-expr_combined[infercnv_results[[1]]$gene_names,]

# ------------------------------------------------------------------------------
print("Combine all information on bin level")
# ------------------------------------------------------------------------------

#Filter first range (keep dimensions)
combined_range<-wgs_results[unique(overlaps$queryHits)]

#Average score of second method across bins
combined_range$infercnv<-sapply(unique(overlaps$queryHits),
        function(i) mean(infercnv_results[[1]]$infercnv_cnv[overlaps$subjectHits[overlaps$queryHits==i]]))                                 

#Add number of genes per bin
combined_range$number_genes<-genes_per_bin$n_genes

#Average expression per segment
combined_range$mean_total<-sapply(unique(overlaps$queryHits),
                                function(i) mean(expr_combined$mean_total[overlaps$subjectHits[overlaps$queryHits==i]]))                                 
combined_range$mean_controls<-sapply(unique(overlaps$queryHits),
                                  function(i) mean(expr_combined$mean_controls[overlaps$subjectHits[overlaps$queryHits==i]]))                                 
combined_range$mean_cancer<-sapply(unique(overlaps$queryHits),
                                  function(i) mean(expr_combined$mean_cancer[overlaps$subjectHits[overlaps$queryHits==i]]))                                 

# ------------------------------------------------------------------------------
print("Create overview plots")
# ------------------------------------------------------------------------------

#Get results
combined_methods<-elementMetadata(combined_range)

#Add position information
combined_methods$chr<-factor(combined_range@seqnames,levels=combined_range@seqnames@values)
combined_methods$start_position<-combined_range@ranges@start
combined_methods<-as.data.frame(combined_methods)

#Add an artifical count through the whole genome and 
#get start positions for each new chromosome
combined_methods$counted_pos<-1:nrow(combined_methods)
chr_boundries<-combined_methods%>%
  group_by(chr)%>%
  summarize(start_chr=min(counted_pos),
            mean_chr=mean(counted_pos))

#Lineplot with the number of genes
g.1<-ggplot(combined_methods,aes(x=counted_pos,y=number_genes))+
  geom_bar(stat="identity")+
  geom_vline(xintercept = chr_boundries$start_chr)+
  ylab("# genes\nper bin")+
  scale_x_continuous(breaks=chr_boundries$mean_chr,
                     labels=chr_boundries$chr)+
  coord_cartesian(xlim=c(1,max(combined_methods$counted_pos)),expand=FALSE)+
  theme_bw()+
  theme(axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        text=element_text(size=21))

#Lineplot with the expression levels
plot_data<-reshape2::melt(combined_methods[,c("chr","start_position",
                                              "counted_pos","mean_total",
                                              "mean_controls","mean_cancer")],
                          id.vars=c("chr","start_position","counted_pos"))
g.2<-ggplot(plot_data,aes(x=counted_pos,y=variable,fill=value))+
  geom_tile()+
  geom_vline(xintercept = chr_boundries$start_chr)+
  ylab("Expression")+
  scale_x_continuous(breaks=chr_boundries$mean_chr,
                     labels=chr_boundries$chr)+
  scale_fill_viridis("Expression")+
  coord_cartesian(xlim=c(1,max(combined_methods$counted_pos)),expand=FALSE)+
  theme_bw()+
  theme(legend.position="bottom",
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        text=element_text(size=21))

#Heatmap with the CNV results
plot_data<-reshape2::melt(combined_methods[,c("chr","start_position",
                                              "counted_pos","wgs_mean","infercnv")],
                          id.vars=c("chr","start_position","counted_pos"))
g.3<-ggplot(plot_data,aes(x=counted_pos,y=variable,fill=value))+geom_tile()+
  theme_bw()+
  scale_fill_gradient2("Score",low = "darkblue",
                       mid = "white",high = "darkred",midpoint = 2)+
  xlab("Chromosome position")+ylab("Method")+
  geom_vline(xintercept = chr_boundries$start_chr)+
  scale_x_continuous(breaks=chr_boundries$mean_chr,
                     labels=chr_boundries$chr)+
  scale_y_discrete(limits=rev)+
  coord_cartesian(xlim=c(1,max(plot_data$counted_pos)),expand=FALSE)+
  theme(legend.position="bottom",
        axis.text.x = element_text(angle=90,vjust=0.5,hjust=1),
        text=element_text(size=21))


g<-ggarrange(g.1,g.2,g.3,ncol=1,align="v",heights=c(0.15,0.45,0.4))
ggsave(g, file=output_plot,width=20,height=10)

print(paste("General correlation between the methods:",
      round(cor(combined_methods$wgs_mean,combined_methods$infercnv),3)))

# ------------------------------------------------------------------------------
print("Create boxplots for the different conditions - For loss")
# ------------------------------------------------------------------------------

#Discretize predictions into loss and gains
combined_methods$wgs_bin<-ifelse(combined_methods$wgs_mean<1.5,1,
                                 ifelse(combined_methods$wgs_mean>2.5,3,2))
combined_methods$infercnv_bin<-ifelse(combined_methods$infercnv<1.5,1,
                                 ifelse(combined_methods$infercnv>2.5,3,2))

#Define categories (TP/FP/FN/TN)
combined_methods$categories_loss<-with(combined_methods,
                                       case_when(wgs_bin==1 & infercnv_bin==1 ~ "TP",
                                                 wgs_bin!=1 & infercnv_bin==1 ~ "FP",
                                                 wgs_bin==1 & infercnv_bin!=1 ~ "FN",
                                                 wgs_bin!=1 & infercnv_bin!=1 ~ "TN"))

#Create boxplots  
plot_data<-reshape2::melt(combined_methods[,c("categories_loss","mean_total",
                                              "mean_controls","mean_cancer")],
                          id.vars=c("categories_loss"))

g<-ggplot(plot_data,aes(x=variable,y=value,fill=categories_loss))+
  geom_boxplot()+
  ylab("Expression")+xlab("")+
  theme(legend.position = "bottom")
ggsave(g,file=output_boxplot_loss,width=5,height=5)

#Check also foldchange between groups
combined_methods$lfc<-log2(combined_methods$mean_cancer/
                             combined_methods$mean_controls)

g<-ggplot(combined_methods,aes(x=categories_loss,y=lfc,fill=categories_loss))+
  geom_boxplot()+
  ylab("Log fold change bin")+xlab("Categorie")+
  theme(legend.position = "none")
ggsave(g,file=output_boxplot_loss_lfc,width=5,height=5)

# ------------------------------------------------------------------------------
print("Create boxplots for the different conditions - For gain")
# ------------------------------------------------------------------------------

#Define categories (TP/FP/FN/TN)
combined_methods$categories_gain<-with(combined_methods,
                                       case_when(wgs_bin==3 & infercnv_bin==3 ~ "TP",
                                                 wgs_bin!=3 & infercnv_bin==3 ~ "FP",
                                                 wgs_bin==3 & infercnv_bin!=3 ~ "FN",
                                                 wgs_bin!=3 & infercnv_bin!=3 ~ "TN"))

#Create boxplots  
plot_data<-reshape2::melt(combined_methods[,c("categories_gain","mean_total",
                                              "mean_controls","mean_cancer")],
                          id.vars=c("categories_gain"))

g<-ggplot(plot_data,aes(x=variable,y=value,fill=categories_gain))+
  geom_boxplot()+
  ylab("Expression")+xlab("")+
  theme(legend.position = "bottom")
ggsave(g,file=output_boxplot_gain,width=5,height=5)

#Check also foldchange between groups
g<-ggplot(combined_methods,aes(x=categories_gain,y=lfc,fill=categories_gain))+
  geom_boxplot()+
  ylab("Log fold change bin")+xlab("Categorie")+
  theme(legend.position = "none")
ggsave(g,file=output_boxplot_gain_lfc,width=5,height=5)

# ------------------------------------------------------------------------------
print("GO enrichment of subgroup")
# ------------------------------------------------------------------------------

#Map contingency categories back to overlaps object
combined_methods$bin_num_wgs<-unique(overlaps$queryHits)

all_genes<-infercnv_results[[1]]$gene_names  
for(cat in c("TP","TN","FN","FP")){
  
  #Get gene names associated with the bins of the respective category
  bin_numbers<-combined_methods$bin_num_wgs[combined_methods$categories_loss==cat]
  selected_genes<-infercnv_results[[1]]$gene_names[overlaps$subjectHits[overlaps$queryHits %in% bin_numbers]]
  
  enrich_out <-enrichGO(gene=unique(selected_genes),
                        OrgDb='org.Hs.eg.db',
                        keyType="SYMBOL",
                        qvalueCutoff = 0.2,
                        pAdjustMethod = "BH",
                        universe = all_genes,
                        ont="all",
                        minGSSize=5)
  
  print(paste("Enrichment",cat,"(losses)"))
  print(head(enrich_out@result))
  
  
  #Get gene names associated with the bins of the respective category
  bin_numbers<-combined_methods$bin_num_wgs[combined_methods$categories_gain==cat]
  selected_genes<-infercnv_results[[1]]$gene_names[overlaps$subjectHits[overlaps$queryHits %in% bin_numbers]]
  
  enrich_out <-enrichGO(gene=unique(selected_genes),
                        OrgDb='org.Hs.eg.db',
                        keyType="SYMBOL",
                        qvalueCutoff = 0.2,
                        pAdjustMethod = "BH",
                        universe = all_genes,
                        ont="all",
                        minGSSize=5)
  
  print(paste("Enrichment",cat,"(gains)"))
  print(head(enrich_out@result))
  
}
  
  
# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()
