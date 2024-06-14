# ------------------------------------------------------------------------------
# Visualize the performance of the methods for the MCF7, the BCC06 and the MM dataset
# with different reference datasets
# ------------------------------------------------------------------------------

library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(viridis)
library(GenomicRanges)
library(numbat)

library(Seurat) #for better count normalization
library(sctransform)

#Source script with help functions
source("scripts/lib.R")

theme_set(theme_bw())

#Vector for renaming methods (official published names)
method_names<-setNames(c("InferCNV (CNV)","InferCNV (Expr)","CaSpER",
                         "copyKat","SCEVAN (CNV)","SCEVAN (Expr)",
                         "Numbat (Expr)", "Numbat (CNV)", "CONICSmat"),
                       c("infercnv_cnv","infercnv_expr","casper",
                         "copykat","scevan_vega","scevan",
                         "numbat","numbat_cnv", "CONICSmat"))

# ------------------------------------------------------------------------------
# Dataset 1 - MCF7
# ------------------------------------------------------------------------------

#Load the different result files
dataset_names<-c("MCF7","MCF7_stromal","MCF7_endothelial","MCF7_immune","MCF7_SNU601")

#Load the files
all_evals<-NULL
for(dataset in dataset_names){
  
  corrs<-fread(paste0("results/output_",dataset,
                      "/evaluation/evaluation_cnv_prediction_corr.tsv"))
  
  #Filter for comparison with ground truth
  corrs<-corrs[corrs$method1 %in% c("scWGS","WGS","WES")]
  corrs<-corrs[,c("method2","pearson")]
  colnames(corrs)<-c("method","value")
  corrs$metric<-"pearson"
  corrs$dataset<-dataset
  
  #Set negative correlation values to 0
  corrs$norm_value<-corrs$value
  corrs$norm_value[corrs$norm_value<0]<-0
  
  #Remove row with scWGS/WES/WGS results
  corrs<-corrs[! corrs$method %in% c("scWGS","WGS","WES"),]
  all_evals<-rbind(all_evals,corrs)
  
  #Add also AUC values (separate values for gain and loss)
  auc<-fread(paste0("results/output_",dataset,
                    "/evaluation/evaluation_cnv_prediction_auc.txt"))
  auc$method<-method_names[auc$method]
  auc<-melt(auc[,c("method","auc_gains","auc_losses")],id.vars="method")
  auc<-auc[,c("method","value","variable")]
  colnames(auc)[3]<-"metric"
  auc$dataset<-dataset
  
  #Scale AUC to also range from 0-1 (set value < 0.5 to 0)
  auc$norm_value<-auc$value
  auc$norm_value<-(auc$norm_value-0.5)*2
  auc$norm_value[auc$norm_value<0]<-0
  
  all_evals<-rbind(all_evals,auc)
  
}

#Order methods based on overall score
method_ordering<-all_evals%>%
  filter(dataset == "MCF7")%>%
  filter(metric %in% c("pearson","auc_gains","auc_losses"))%>%
  group_by(method)%>%
  summarize(norm_value=mean(norm_value))%>%
  arrange(norm_value)

#Specify "MCF7" better: "MCF7_epithelial"
all_evals$dataset[all_evals$dataset=="MCF7"]<-"MCF7_epithelial"
dataset_names[dataset_names=="MCF7"]<-"MCF7_epithelial"

all_evals$method<-factor(all_evals$method,levels=method_ordering$method)
all_evals$dataset<-factor(all_evals$dataset,levels=dataset_names)

g<-ggplot(all_evals,aes(y=method,x=dataset,fill=norm_value))+
  geom_tile()+
  scale_fill_viridis("Normalized\nScore",limits=c(0,1))+
  geom_text(aes(label=round(value,2),
                color=ifelse(norm_value<0.6,'white','black')),
            size=2.8)+
  scale_color_manual(values=c("black","white"))+
  facet_grid(~metric,scales="free",space="free")+
  theme(axis.title.y = element_blank(),
        axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
        legend.position="none")+
  guides(color="none")
ggsave(g,file="~/Desktop/reference_dataset_MCF7_eval.png",width=12,height=5)


# ------------------------------------------------------------------------------
# Plot karyograms of the same method with different reference datasets to improve
# comparison - Version 1: Numbat CNV
# ------------------------------------------------------------------------------

dataset_names<-c("MCF7","MCF7_stromal","MCF7_endothelial","MCF7_immune","MCF7_SNU601")
method_names<-c("wgs_mean",dataset_names)

combined_range<-read_wgs_results("data/MCF7_WGS_groundtruth/wgs_results_formated.csv")

for(output_version in dataset_names){
  numbat_results <- read_numbat_cnv(paste0("results/output_",output_version,"/numbat/"))
  colnames(mcols(numbat_results))<-c(output_version) 
  
  combined_range<-combine_range_objects(combined_range, numbat_results,
                                        method_colname=output_version)
}

#Get results
combined_methods<-elementMetadata(combined_range)

#Replace MCF7 with MCF7_epithelial
method_names[method_names=="MCF7"]<-"MCF7_epithelial"
colnames(combined_methods)[colnames(combined_methods)=="MCF7"]<-"MCF7_epithelial"

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
for(method in method_names){
  scaled_methods[,method]<-(scaled_methods[,method]-2) / 
    sd(scaled_methods[,method])
  
  #Save the standard deviation of each method to document the chosen normalization factor
  scaling_factor<-rbind(scaling_factor,
                        data.frame(method,
                                   sd=sd(combined_methods[,method]),
                                   mean=mean(combined_methods[,method]-2)))
}

plot_data<-reshape2::melt(scaled_methods,
                          id.vars=c("chr","start_position","counted_pos"))

#Order them respectively
plot_data$variable<-factor(plot_data$variable,
                           levels=method_names)

g <- ggplot(plot_data,aes(x=counted_pos,y=variable,fill=value))+geom_tile()+
  theme_bw()+
  scale_fill_gradient2("Score",low = "darkblue",
                       mid = "white",high = "darkred",midpoint = 0,
                       breaks=c(-5,0,5))+
  xlab("Chromosome position")+ylab("Reference dataset")+
  geom_vline(xintercept = chr_boundries$start_chr)+
  scale_x_continuous(breaks=chr_boundries$mean_chr,
                     labels=chr_boundries$chr)+
  scale_y_discrete(limits=rev)+
  coord_cartesian(xlim=c(1,max(plot_data$counted_pos)),expand=FALSE)+
  theme(legend.position="bottom",
        axis.text.x = element_text(angle=90,vjust=0.5,hjust=1),
        text=element_text(size=21))
ggsave(g, file="~/Desktop/MCF_numbat_diffrefs.pdf",
       width=20,height=6)

# ------------------------------------------------------------------------------
# Version 2: InferCNV
# ------------------------------------------------------------------------------

dataset_names<-c("MCF7","MCF7_stromal","MCF7_endothelial","MCF7_immune","MCF7_SNU601")
method_names<-c("wgs_mean",dataset_names)

#Get CNV results for each method
combined_range<-read_wgs_results("data/MCF7_WGS_groundtruth/wgs_results_formated.csv")
for(output_version in dataset_names){
  
  infercnv_results<-read_infercnv_6state_model(paste0("results/output_",output_version,
                                                      "/infercnv/infercnv.20_HMM_predHMMi6.leiden.hmm_mode-subclusters.",
                                                      "Pnorm_0.5.repr_intensities.observations.txt"),
                                               "data/annotations/hg38_gencode_v27.txt")[[1]]
  infercnv_results$gene<-NULL
  colnames(mcols(infercnv_results))<-c(output_version) 
  
  combined_range<-combine_range_objects(combined_range, infercnv_results,
                                        method_colname=output_version)
}

#Gene positions
gene_pos<-fread("data/annotations/hg38_gencode_v27.txt")
colnames(gene_pos)<-c("genes","chr","start","end")

#Get mean expression of each reference dataset
gene_values<-combined_range
mcols(gene_values)<-NULL
for(output_version in dataset_names){
  
  #Load count matrix and annotations
  counts<-suppressWarnings(fread(paste0("data/input_",output_version,"/count_matrix.txt")))
  gene_names<-counts$V1
  counts$V1<-NULL
  counts<-as.matrix(counts)
  rownames(counts)<-gene_names
  
  # #Library size normalization to 10,000 counts per cell and logarithmize
  # counts<-t(t(counts)/colSums(counts)*10000)
  # norm_counts<-log(counts+1)
  
  #Size normalization with sctransform
  seurat <- CreateSeuratObject(counts = counts)
  seurat <- SCTransform(seurat, verbose = TRUE)
  norm_counts<-as.matrix(seurat@assays$SCT@data)
    
  annotation<-fread(paste0("data/input_",output_version,"/sample_annotation.txt"), header=FALSE)
  ref_groups<-fread(paste0("data/input_",output_version,"/ref_groups.txt"))
  ref_cells<-annotation$V1[annotation$V2 %in% ref_groups$ref_groups]
  
  #Get gene mean expression
  gene_means<-data.frame(genes=rownames(norm_counts),
                         mean_expr=rowMeans(norm_counts[,ref_cells]))
  gene_means<-merge(gene_pos,gene_means,by="genes")
  gene_means_gr<-makeGRangesFromDataFrame(gene_means,keep.extra.columns = TRUE)
  colnames(mcols(gene_means_gr))<-c("genes",output_version) 
  
  #Map mean gene expression to the bins
  gene_values<-combine_range_objects(gene_values, gene_means_gr,
                                        method_colname=output_version)
}

#Get a lineplot for the gene values
gene_expr<-as.data.frame(elementMetadata(gene_values))

#Replace MCF7 with MCF7_epithelial
colnames(gene_expr)[colnames(gene_expr)=="MCF7"]<-"MCF7_epithelial"
method_names[method_names=="MCF7"]<-"MCF7_epithelial"

gene_expr$counted_pos<-1:length(combined_range)
gene_expr<-reshape2::melt(gene_expr,id.vars=c("counted_pos"))

#Order them respectively
gene_expr$variable<-factor(gene_expr$variable,levels=method_names)

#Get CNV results
combined_methods<-elementMetadata(combined_range)

#Replace MCF7 with MCF7_epithelial
colnames(combined_methods)[colnames(combined_methods)=="MCF7"]<-"MCF7_epithelial"

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
for(method in method_names){
  scaled_methods[,method]<-(scaled_methods[,method]-2) / 
    sd(scaled_methods[,method])
  
  #Save the standard deviation of each method to document the chosen normalization factor
  scaling_factor<-rbind(scaling_factor,
                        data.frame(method,
                                   sd=sd(combined_methods[,method]),
                                   mean=mean(combined_methods[,method]-2)))
}

plot_data<-reshape2::melt(scaled_methods,
                          id.vars=c("chr","start_position","counted_pos"))

#Order them respectively
plot_data$variable<-factor(plot_data$variable,
                           levels=method_names)

g.1<-ggplot(gene_expr,aes(x=counted_pos,fill=value,y=variable))+
  geom_tile()+
  #scale_fill_gradient("Mean expression",low="white",high="darkred")+
  scale_fill_viridis("Mean expression",option="inferno",direction=-1,
                     breaks=c(0,1,2))+
  scale_y_discrete(limits=rev)+
  geom_vline(xintercept = chr_boundries$start_chr)+
  coord_cartesian(xlim=c(1,max(plot_data$counted_pos)),expand=FALSE)+
  ylab("Reference dataset")+
  theme(legend.position="bottom",
        axis.text.x = element_blank(),
        axis.title.x=element_blank(),
        text=element_text(size=21))

g.2 <- ggplot(plot_data,aes(x=counted_pos,y=variable,fill=value))+geom_tile()+
  theme_bw()+
  scale_fill_gradient2("Score",low = "darkblue",
                       mid = "white",high = "darkred",midpoint = 0,
                       breaks=c(-5,0,5))+
  xlab("Chromosome position")+ylab("Reference dataset")+
  geom_vline(xintercept = chr_boundries$start_chr)+
  scale_x_continuous(breaks=chr_boundries$mean_chr,
                     labels=chr_boundries$chr)+
  scale_y_discrete(limits=rev)+
  coord_cartesian(xlim=c(1,max(plot_data$counted_pos)),expand=FALSE)+
  theme(legend.position="bottom",
        axis.text.x = element_text(angle=90,vjust=0.5,hjust=1),
        text=element_text(size=21))

g<-ggarrange(g.1,g.2,nrow=2)
ggsave(g, file="~/Desktop/MCF_infercnv_diffrefs.pdf",
       width=20,height=12)


# ------------------------------------------------------------------------------
# Dataset 2 - BCC06
# ------------------------------------------------------------------------------

#Load the different result files
dataset_names<-c("BCC06","BCC06_immune","BCC06_TS","BCC06_SNU601")

#Load the files
all_evals<-NULL
for(dataset in dataset_names){
  
  corrs<-fread(paste0("results/output_",dataset,
                      "/evaluation/evaluation_cnv_prediction_corr_wes.tsv"))
  
  #Filter for comparison with ground truth
  corrs<-corrs[corrs$method1 %in% c("scWGS","WGS","WES","GATK")]
  corrs<-corrs[,c("method2","pearson")]
  colnames(corrs)<-c("method","value")
  corrs$metric<-"pearson"
  corrs$dataset<-dataset
  
  #Set negative correlation values to 0
  corrs$norm_value<-corrs$value
  corrs$norm_value[corrs$norm_value<0]<-0
  
  #Remove row with scWGS/WES/WGS results
  corrs<-corrs[! corrs$method %in% c("scWGS","WGS","WES","GATK"),]
  all_evals<-rbind(all_evals,corrs)
  
  #Add also AUC values (separate values for gain and loss)
  auc<-fread(paste0("results/output_",dataset,
                    "/evaluation/evaluation_cnv_prediction_auc_wes.txt"))
  auc$method<-method_names[auc$method]
  auc<-melt(auc[,c("method","auc_gains","auc_losses")],id.vars="method")
  auc<-auc[,c("method","value","variable")]
  colnames(auc)[3]<-"metric"
  auc$dataset<-dataset
  
  #Scale AUC to also range from 0-1 (set value < 0.5 to 0)
  auc$norm_value<-auc$value
  auc$norm_value<-(auc$norm_value-0.5)*2
  auc$norm_value[auc$norm_value<0]<-0
  
  all_evals<-rbind(all_evals,auc)
  
}

#Order methods based on overall score
method_ordering<-all_evals%>%
  filter(dataset == "BCC06")%>%
  filter(metric %in% c("pearson","auc_gains","auc_losses"))%>%
  group_by(method)%>%
  summarize(norm_value=mean(norm_value))%>%
  arrange(norm_value)


all_evals$method<-factor(all_evals$method,levels=method_ordering$method)
all_evals$dataset<-factor(all_evals$dataset,levels=dataset_names)

g<-ggplot(all_evals,aes(y=method,x=dataset,fill=norm_value))+
  geom_tile()+
  scale_fill_viridis("Normalized\nScore",limits=c(0,1))+
  geom_text(aes(label=round(value,2),
                color=ifelse(norm_value<0.6,'white','black')),
            size=2.8)+
  scale_color_manual(values=c("black","white"))+
  facet_grid(~metric,scales="free",space="free")+
  theme(axis.title.y = element_blank(),
        axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
        legend.position="none")+
  guides(color="none")
ggsave(g,file="~/Desktop/reference_dataset_BCC06_eval.png",width=9,height=5)

# ------------------------------------------------------------------------------
# Dataset 3 - MM
# ------------------------------------------------------------------------------

#Load the different result files
dataset_names<-c("MM","MM_tcell","MM_bcell","MM_mono","MM_SNU601")

#Load the files
all_evals<-NULL
for(dataset in dataset_names){
  
  corrs<-fread(paste0("results/output_",dataset,
                      "/evaluation/evaluation_cnv_prediction_corr_wes.tsv"))
  
  #Filter for comparison with ground truth
  corrs<-corrs[corrs$method1 %in% c("scWGS","WGS","WES","GATK")]
  corrs<-corrs[,c("method2","pearson")]
  colnames(corrs)<-c("method","value")
  corrs$metric<-"pearson"
  corrs$dataset<-dataset
  
  #Set negative correlation values to 0
  corrs$norm_value<-corrs$value
  corrs$norm_value[corrs$norm_value<0]<-0
  
  #Remove row with scWGS/WES/WGS results
  corrs<-corrs[! corrs$method %in% c("scWGS","WGS","WES","GATK"),]
  all_evals<-rbind(all_evals,corrs)
  
  #Add also AUC values (separate values for gain and loss)
  auc<-fread(paste0("results/output_",dataset,
                    "/evaluation/evaluation_cnv_prediction_auc_wes.txt"))
  auc$method<-method_names[auc$method]
  auc<-melt(auc[,c("method","auc_gains","auc_losses")],id.vars="method")
  auc<-auc[,c("method","value","variable")]
  colnames(auc)[3]<-"metric"
  auc$dataset<-dataset
  
  #Scale AUC to also range from 0-1 (set value < 0.5 to 0)
  auc$norm_value<-auc$value
  auc$norm_value<-(auc$norm_value-0.5)*2
  auc$norm_value[auc$norm_value<0]<-0
  
  all_evals<-rbind(all_evals,auc)
  
}

#Order methods based on overall score
method_ordering<-all_evals%>%
  filter(dataset == "MM")%>%
  filter(metric %in% c("pearson","auc_gains","auc_losses"))%>%
  group_by(method)%>%
  summarize(norm_value=mean(norm_value))%>%
  arrange(norm_value)


all_evals$method<-factor(all_evals$method,levels=method_ordering$method)
all_evals$dataset<-factor(all_evals$dataset,levels=dataset_names)

g<-ggplot(all_evals,aes(y=method,x=dataset,fill=norm_value))+
  geom_tile()+
  scale_fill_viridis("Normalized\nScore",limits=c(0,1))+
  geom_text(aes(label=round(value,2),
                color=ifelse(norm_value<0.6,'white','black')),
            size=2.8)+
  scale_color_manual(values=c("black","white"))+
  facet_grid(~metric,scales="free",space="free")+
  theme(axis.title.y = element_blank(),
        axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
        legend.position="none")+
  guides(color="none")
ggsave(g,file="~/Desktop/reference_dataset_MM_eval.png",width=9,height=5)

#Get the deviations between datasets
tmp<-dcast(all_evals,method+metric~dataset,value.var="value")
orginal_vals<-tmp$MM
for(dataset in dataset_names){
  tmp[,dataset]<-tmp[,..dataset]-orginal_vals
}
tmp<-melt(tmp,id.vars=c("method","metric"))
colnames(tmp)[3:4]<-c("dataset","deviation")
all_evals<-merge(all_evals,tmp,by=c("method","metric","dataset"))


g<-ggplot(all_evals,aes(y=method,x=dataset,fill=deviation))+
  geom_tile()+
  scale_fill_gradient2("Deviation",low="darkblue",mid="white",high="darkred",limits=c(-1,1))+
  geom_text(aes(label=round(value,2),
                color=ifelse(abs(deviation)>0.5,'white','black')),
            size=2.8)+
  scale_color_manual(values=c("black","white"))+
  facet_grid(~metric,scales="free",space="free")+
  theme(axis.title.y = element_blank(),
        axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
        legend.position="none")+
  guides(color="none")

ggsave(g,file="~/Desktop/reference_dataset_MM_eval_dev.png",width=9,height=5)
