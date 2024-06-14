# ------------------------------------------------------------------------------
# Check how different the CNVs are across the four individuals
# ------------------------------------------------------------------------------

library(data.table)
library(GenomicRanges)
library(dplyr)
library(ggplot2)
library(ggdendro)
library(cowplot)

source("scripts/lib.R")

sample_names<-c("BCC05","BCC06","BCC07","BCC08")

# Read input data from copykat
combined_input<-NULL
sample_annot<-NULL
for(sN in sample_names){
  
  matrix_path<-paste0("results/output_",sN,"/copykat/",sN,
                      "_copykat_CNA_raw_results_gene_by_cell.txt")
  annot_path<-paste0("data/input_",sN,"/sample_annotation.txt")
  ref_path<-paste0("data/input_",sN,"/ref_groups.txt")
  
  copykat_res<-fread(matrix_path)
  copykat_res$chromosome_name<-paste0("chr",copykat_res$chromosome_name)
  
  #Convert into a Grange
  copykat_grange<-makeGRangesFromDataFrame(copykat_res,
                                           seqnames.field="chromosome_name",
                                           start.field="start_position",
                                           end.field="end_position")
  
  #Get cell matrix
  copykat_cells<-as.matrix(copykat_res[,8:ncol(copykat_res)])
  
  #Remove reference cells
  annot<-fread(annot_path,header=FALSE)
  ref_groups<-fread(ref_path)
  cancer_cells<-annot$V1[!(annot$V2 %in% ref_groups$ref_groups)]
  copykat_cells<-copykat_cells[,colnames(copykat_cells) %in% cancer_cells]
  colnames(copykat_cells)<-paste0(sN,"_",colnames(copykat_cells))
  mcols(copykat_grange)<-copykat_cells
  
  #Annotation data frame with samples
  sample_annot<-rbind(sample_annot,
                      data.frame(cell=colnames(copykat_cells),
                                 sample=sN))
  if(is.null(combined_input)){
    binned_genome<-create_binned_genome(copykat_grange)
    combined_input<-combine_range_objects_allcols(binned_genome,copykat_grange)
  } else {
    combined_input<-combine_range_objects_allcols(combined_input,copykat_grange)
  }
}


# ------------------------------------------------------------------------------
# Create the heatmap
# ------------------------------------------------------------------------------

combined_methods<-elementMetadata(combined_input)

#Add position information
combined_methods$chr<-factor(combined_input@seqnames,
                             levels=combined_input@seqnames@values)
combined_methods$start_position<-combined_input@ranges@start
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
for(method in setdiff(colnames(scaled_methods),c("chr","start_position","counted_pos"))){
  scaled_methods[,method]<-scaled_methods[,method] / 
    sd(scaled_methods[,method],na.rm=TRUE)
}

#Cluster methods on how similar the clonal structure is
scaled_data<-scaled_methods[,setdiff(colnames(scaled_methods),
                                     c("chr","start_position","counted_pos"))]
cluster_clones<-hclust(dist(t(scaled_data)), method="average")

# Build dendrogram object from hclust results
dend <- as.dendrogram(cluster_clones)
dend_data <- dendro_data(dend) #, type = "rectangle")

# Plot line segments and add labels
g.1 <- ggplot(dend_data$segments) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend))+
  scale_y_reverse()+
  scale_x_reverse()+
  coord_flip(expand=FALSE)+
  ylab("")+
  theme(panel.background=element_blank(),
        axis.ticks=element_blank(),
        axis.text=element_blank(),
        axis.line=element_blank(),
        axis.title.y=element_blank(),
        axis.title.x=element_text(size=18),
        panel.grid.major = element_blank(),
        panel.grid.minor=element_blank(),
        panel.border = element_blank())

plot_data<-reshape2::melt(scaled_methods,
                          id.vars=c("chr","start_position","counted_pos"))

#Order them respectively
plot_data$variable<-factor(plot_data$variable,
                           levels=cluster_clones$labels[cluster_clones$order])

g.2<-ggplot(plot_data,aes(x=counted_pos,y=variable,fill=value))+geom_tile()+
  theme_bw()+
  scale_fill_gradient2("Score",low = "darkblue",
                       mid = "white",high = "darkred",midpoint = 0,
                       breaks=c(-5,0,5))+
  xlab("Chromosome position")+
  geom_vline(xintercept = chr_boundries$start_chr)+
  scale_x_continuous(breaks=chr_boundries$mean_chr,
                     labels=chr_boundries$chr)+
  scale_y_discrete(limits=rev)+
  coord_cartesian(xlim=c(1,max(plot_data$counted_pos)),expand=FALSE)+
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x=element_text(size=18),
        legend.position = 'none',
        axis.title.y=element_blank(),
        axis.text.y = element_blank())

#Add a barplot with the proportions per donor
sample_annot$cell<-gsub("-",".",as.character(sample_annot$cell))
sample_annot$cell<-factor(sample_annot$cell,
                          levels=rev(cluster_clones$labels[cluster_clones$order]))
g.3<-ggplot(sample_annot,aes(x=cell,y=1,fill=sample))+
  geom_tile()+
  coord_flip()+
  labs(fill="Donor annotation")+
  theme(legend.title=element_text(size=16),
        legend.text=element_text(size=16))+
  guides(fill=guide_legend(nrow=1,byrow=TRUE))

gglegend<-get_legend(g.3)

g.3<-g.3+
  ylab("")+
  theme(axis.title.y=element_blank(),
        axis.title.x=element_text(size=18),
        axis.ticks=element_blank(),
        axis.text=element_blank(),
        legend.position="none")

combiplot <- cowplot::plot_grid(g.1, g.2, g.3, ncol = 3, 
                                rel_widths = c(0.1,1,0.03), 
                                axis='b', align = 'hw')
combiplot <- cowplot::plot_grid(combiplot,gglegend,ncol=1,
                                rel_heights=c(1,0.05))
ggsave("../figure_plots/supp_figure_BCC_karyogram.png",combiplot, width = 38, height=20, units = "in")

