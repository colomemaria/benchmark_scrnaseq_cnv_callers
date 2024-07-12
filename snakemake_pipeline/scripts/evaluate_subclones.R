# ------------------------------------------------------------------------------
# Evaluate how similar the clones are identified by the different methods
# Requires an (artifical) ground truth (first test: pooling samples 
# from different donors together)
# ------------------------------------------------------------------------------

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# ------------------------------------------------------------------------------
print("Load libraries")
# ------------------------------------------------------------------------------

library(data.table)
library(mclust) #for the adjusted Rand Index
library(clevr) # to get the homogeneity between clusters
library(ggplot2)
library(viridis)
library(infercnv)
library(CONICSmat)
library(numbat)
library(GenomicRanges)
library(dplyr)
library(ggdendro)
library(ggpubr)

theme_set(theme_bw())

source("scripts/lib.R")

# ------------------------------------------------------------------------------
print("Get input parameters from snakemake")
# ------------------------------------------------------------------------------

input_expr_mat<-snakemake@input$expr_mat
input_annot<-snakemake@input$annot
input_ref_groups<-snakemake@input$ref_groups
input_annot_truth<-snakemake@input$annot_truth

input_clusters_copykat<-snakemake@input$clusters_copykat
input_matrix_copykat<-snakemake@input$matrix_copykat

input_clones_numbat<-snakemake@input$clones_numbat

input_infercnv_cnv<-snakemake@input$infercnv_cnv
input_infercnv_gene_pos<-snakemake@input$infercnv_gene_pos
input_clusters_infercnv<-snakemake@input$clusters_infercnv

input_clones_scevan<-snakemake@input$clones_scevan
input_scevan_subclones<-snakemake@input$scevan_subclones

input_l_CONICSmat<-snakemake@input$l_CONICSmat
input_lr_CONICSmat<-snakemake@input$lr_CONICSmat
input_genes_CONICSmat<-snakemake@input$genes_CONICSmat
input_CONICSmat_chrom_pos<-snakemake@input$CONICSmat_chrom_pos
input_CONICSmat_cnv<-snakemake@input$CONICSmat_cnv

input_casper_hclust<-snakemake@input$casper_hclust
input_casper_cellmatrix<-snakemake@input$casper_cellmatrix
input_casper_gene_annot<-snakemake@input$casper_gene_annot

param_dataset<-snakemake@params$dataset

output_merged_clusters<-snakemake@output$merged_clusters
output_merged_results<-snakemake@output$merged_results
output_conicsmat_histo<-snakemake@output$conicsmat_histo
output_text_ari<-snakemake@output$text_ari
output_plot_ari<-snakemake@output$plot_ari
output_plot_clust_dist<-snakemake@output$plot_clust_dist
output_plot_clones_clustered<-snakemake@output$plot_clones_clustered

# ------------------------------------------------------------------------------
# Load files
# ------------------------------------------------------------------------------

#Load ground truth
ground_truth<-fread(input_annot_truth)
colnames(ground_truth)[2]<-"true_clusters"

#Define the number of clusters dependent on how many clusters are in the ground truth
param_splitted_clusters<-length(unique(ground_truth$true_clusters))

#Load dendogram from copykat
clusters_copykat<-readRDS(input_clusters_copykat)

#Split copykat into the same number of groups (number defined as parameter)
groups_copykat<-cutree(clusters_copykat,k=param_splitted_clusters)
clustering_copykat<-data.frame(cell=names(groups_copykat),
                               copykat=groups_copykat)
#Replace trailing .1 (not formated correctly)
clustering_copykat$cell<-gsub("\\.1","-1",clustering_copykat$cell)
rownames(clustering_copykat)<-clustering_copykat$cell
all_clusters<-merge(ground_truth,clustering_copykat,by="cell",all=TRUE)

#Load clones from numbat
clones_numbat<-fread(input_clones_numbat)
clones_numbat<-clones_numbat[,c("cell","clone_opt")]
colnames(clones_numbat)<-c("cell","numbat")
all_clusters<-merge(all_clusters,clones_numbat,by="cell",all=TRUE)

#Load subclusters from InferCNV
infercnv_obj <- readRDS(input_clusters_infercnv)
subclusters_tumor<-infercnv_obj@tumor_subclusters$subclusters[[param_dataset]]
cell_names<-colnames(infercnv_obj@expr.data)
clustering_infercnv<-NULL
for(subc in names(subclusters_tumor)){
  clustering_infercnv<-rbind(clustering_infercnv,
            data.frame(cell=cell_names[subclusters_tumor[[subc]]],
                       infercnv=subc))
                              
}
rm(infercnv_obj, subclusters_tumor)
gc()
clustering_infercnv$infercnv<-gsub(".*_s","",clustering_infercnv$infercnv)
all_clusters<-merge(all_clusters,clustering_infercnv,by="cell",all=TRUE)

#CONICSmat (need to calculate the clusters on the fly)
l<-as.matrix(read.table(input_l_CONICSmat,sep="\t",header=T,row.names=1,check.names=F))
lrbic<-read.table(input_lr_CONICSmat,
                  sep="\t",header=T,row.names=1,check.names=F)
candRegions<-rownames(lrbic)[which(lrbic[,"BIC difference"]>200 &
                                     lrbic[,"LRT adj. p-val"]<0.01)]

#Requires additionally filtered expression matrix (using same code as in CONICSmat)
expr_df <- as.matrix(read.table(input_expr_mat,sep="\t",header=T,row.names=1,check.names=F))
gene_pos <- read.table(file = input_genes_CONICSmat, header = TRUE, sep = "\t")
expr_df[which(is.na(expr_df))] <- 0
expr_df <- filterMatrix(expr_df,gene_pos[,"hgnc_symbol"],minCells=5)

pdf(output_conicsmat_histo)
hi<-plotHistogram(l[,candRegions],expr_df,
                  clusters=param_splitted_clusters,zscoreThreshold=4)
dev.off()

clusters_conicsmat<-data.frame(cell=names(hi),
                               conicsmat=hi)
all_clusters<-merge(all_clusters,clusters_conicsmat,by="cell",all=TRUE)

# Get CaSpER clones
clusters_casper<-readRDS(input_casper_hclust)
groups_casper<-cutree(clusters_casper,k=param_splitted_clusters)
clustering_casper<-data.frame(cell=names(groups_casper),
                              casper=groups_casper)
all_clusters<-merge(all_clusters,clustering_casper,by="cell",all=TRUE)

#SCEVAN (need to run specific mode to find subclusters)
clustering_scevan<-fread(input_clones_scevan)
#Remove filtered and normal cells (no subclone defined for those)
clustering_scevan<-clustering_scevan[clustering_scevan$class == "tumor",]
clustering_scevan<-clustering_scevan[,c("V1","subclone")]
colnames(clustering_scevan)<-c("cell","scevan")
#Replace all . by _ because SCEVAN causes problems here
all_clusters$cell<-gsub("\\.","_",all_clusters$cell)
all_clusters<-merge(all_clusters,clustering_scevan,by="cell",all=TRUE)

#Remove all healthy cells
annotation<-fread(input_annot, header=FALSE)
ref_groups<-read.table(input_ref_groups,header=TRUE)
cancer_cells <- annotation$V1[! (annotation$V2 %in% ref_groups$ref_groups)]
#Replace all . by _ because SCEVAN causes problems here
cancer_cells<-gsub("\\.","_",cancer_cells)
all_clusters<-all_clusters[all_clusters$cell %in% cancer_cells,]

#Save the combined result data frame for later processing
write.table(all_clusters,file=output_merged_clusters,
            quote = FALSE,row.names=FALSE,sep="\t")

# ------------------------------------------------------------------------------
print("Pairwise comparisons using adjusted rand index")
# ------------------------------------------------------------------------------

#Convert to call columns
all_clusters<-as.data.frame(all_clusters)

methods<-c("true_clusters","copykat","infercnv","scevan",
           "conicsmat","numbat","casper")
compare_methods<-NULL
for(m1 in 1:length(methods)){
  for(m2 in 1:length(methods)){
    compare_methods<-rbind(compare_methods,
                           data.frame(method1=methods[m1],
                                method2=methods[m2],
                                adj_rand_index=mclust::adjustedRandIndex(
                                  all_clusters[,methods[m1]],
                                  all_clusters[,methods[m2]]),
                                homogen = clevr::homogeneity(
                                  all_clusters[,methods[m1]],
                                  all_clusters[,methods[m2]]),
                                merged_ari=recursive_ari(
                                  all_clusters[,methods[m1]],
                                  all_clusters[,methods[m2]])))
  }
}

#Save result txt
write.table(compare_methods,file=output_text_ari,sep="\t",quote=FALSE,
            row.names=FALSE)

#Plot results
compare_methods$method1<-factor(compare_methods$method1,levels=methods)
compare_methods$method2<-factor(compare_methods$method2,levels=methods)

g.1<-ggplot(compare_methods,aes(x=method2,y=method1,fill=adj_rand_index))+
  geom_tile()+
  theme_bw()+
  geom_text(aes(label=round(adj_rand_index,2),
                color=ifelse(adj_rand_index<0.6,'white','black')),size=2.8)+
  scale_color_manual(values=c("black","white"))+
  scale_y_discrete(limits=rev)+
  xlab("Method")+
  ylab("Method")+
  scale_fill_viridis("Adj Rand Index",limits=c(0,1))+
  guides(color="none")+
  theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1),
        legend.position = "bottom")

g.2<-ggplot(compare_methods,aes(x=method2,y=method1,fill=merged_ari))+
  geom_tile()+
  theme_bw()+
  geom_text(aes(label=round(merged_ari,2),
                color=ifelse(merged_ari<0.6,'white','black')),size=2.8)+
  scale_color_manual(values=c("black","white"))+
  scale_y_discrete(limits=rev)+
  xlab("Method")+
  ylab("Method")+
  scale_fill_viridis("Merged ARI",limits=c(0,1))+
  guides(color="none")+
  theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1),
        legend.position = "bottom")

g.3<-ggplot(compare_methods,aes(x=method2,y=method1,fill=homogen),
            legend.position = "bottom")+
  geom_tile()+
  theme_bw()+
  geom_text(aes(label=round(homogen,2),
                color=ifelse(homogen<0.6,'white','black')),size=2.8)+
  scale_color_manual(values=c("black","white"))+
  scale_y_discrete(limits=rev)+
  xlab("Method")+
  ylab("Method")+
  scale_fill_viridis("Homogeneity",limits=c(0,1))+
  guides(color="none")+
  theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1),
        legend.position = "bottom")

g<-ggarrange(g.1,g.2,g.3,nrow=1)  
ggsave(g, file=output_plot_ari,
       width=11,height=4.5)

# ------------------------------------------------------------------------------
print("Create a stacked barplot how the donors are distributed for each method")
# ------------------------------------------------------------------------------

comp_clusters<-melt(all_clusters,id.vars=c("cell","true_clusters"))
comp_clusters<-comp_clusters%>%
  group_by(variable,value,true_clusters)%>%
  summarize(n_cells=n())

comp_clusters$cluster_name<-paste(comp_clusters$variable,comp_clusters$value,sep="_")

g<-ggplot(comp_clusters,aes(x=cluster_name,y=n_cells,fill=true_clusters))+
  geom_bar(stat="identity")+xlab("Cluster")+ylab("Clone")+
  scale_fill_discrete("Donor")+
  theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5))
ggsave(g,file=output_plot_clust_dist,width=6,height=4)

# ------------------------------------------------------------------------------
print("Get per subclone and method pseudobulk results")
# ------------------------------------------------------------------------------

print("Load copykat data (normalized expression)")
copykat_results<-read_copykat(input_matrix_copykat,input_annot,input_ref_groups,
                              clones=clustering_copykat)[[1]]
colnames(mcols(copykat_results))<-gsub("clone","copykat",
                                       colnames(mcols(copykat_results)))

#Create genomic range object
binned_genome<-create_binned_genome(copykat_results)
combined_clones<-combine_range_objects_allcols(binned_genome,copykat_results)

print("Load SCEVAN data (segments per clone)")
scevan_results<-read_scevan_subclones(input_scevan_subclones,binned_genome)
colnames(mcols(scevan_results))<-gsub("CNV_subclone","scevan_",
                                      colnames(mcols(scevan_results)))
combined_clones<-combine_range_objects_allcols(combined_clones,scevan_results)

print("Load inferCNV data (pseudobulk CNVs)")
infercnv_results<-read_infercnv_6state_model(input_infercnv_cnv,
                                             input_infercnv_gene_pos,
                                             clones=clustering_infercnv)[[1]]
infercnv_results$gene<-NULL
colnames(mcols(infercnv_results))<-paste0("infercnv_",
                                        colnames(mcols(infercnv_results)))
combined_clones<-combine_range_objects_allcols(combined_clones,infercnv_results)

print("Load CONICSmat data")
CONICSmat_results <- read_CONICSmat(input_CONICSmat_chrom_pos, input_CONICSmat_cnv,
                                    input_annot, input_ref_groups,
                                    clones=clusters_conicsmat)
colnames(mcols(CONICSmat_results))<-gsub("clone","conicsmat_",
                                          colnames(mcols(CONICSmat_results)))
combined_clones<-combine_range_objects_allcols(combined_clones,CONICSmat_results)

print("Load Numbat data")
numbat_results <- read_numbat_cnv(dirname(input_clones_numbat),clones=TRUE)
colnames(mcols(numbat_results))<-gsub("cnv_","",colnames(mcols(numbat_results)))
combined_clones<-combine_range_objects_allcols(combined_clones,numbat_results)

print("Load Casper data")
casper_results <- read_casper_cnv(input_casper_cellmatrix,input_casper_gene_annot,
                                  clustering_casper)
combined_clones<-combine_range_objects_allcols(combined_clones,casper_results)

#Save the combined result data frame for later processing
write.table(combined_clones,file=output_merged_results,
            quote = FALSE,row.names=FALSE,sep="\t")

# ------------------------------------------------------------------------------
print(paste0("Generate a heatmap per subclone and method ",
             "combined with a barplot of cells per donor"))
# ------------------------------------------------------------------------------

combined_methods<-elementMetadata(combined_clones)

#Add position information
combined_methods$chr<-factor(combined_clones@seqnames,
                             levels=combined_clones@seqnames@values)
combined_methods$start_position<-combined_clones@ranges@start
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
  scaled_methods[,method]<-(scaled_methods[,method]-2) / 
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
  theme(panel.background=element_blank(), 
        axis.ticks=element_blank(), 
        axis.text=element_blank(), 
        axis.line=element_blank(), 
        axis.title=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor=element_blank(),
        panel.border = element_blank(),
        plot.title = element_text(size=18),
        plot.subtitle=element_text(size=12))
  
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
  theme(legend.position="bottom",
        axis.text.x = element_text(angle=90,vjust=0.5,hjust=1),
        axis.title.y=element_blank(),
        text=element_text(size=20))

#Add a barplot with the proportions per donor
comp_clusters$cluster_name<-factor(comp_clusters$cluster_name,
                           levels=rev(cluster_clones$labels[cluster_clones$order]))
comp_clusters<-comp_clusters[!is.na(comp_clusters$cluster_name),]
g.3<-ggplot(comp_clusters,aes(x=cluster_name,y=n_cells,fill=true_clusters))+
  geom_bar(stat="identity")+ylab("# cells per donor")+
  scale_fill_discrete("Donor")+
  coord_flip(expand=FALSE)+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position="bottom",
        text=element_text(size=16))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))

g<-cowplot::plot_grid(g.1,g.2,g.3,ncol=3,rel_widths=c(0.1,1,0.2),axis="b",align="h")
ggsave(g, file=output_plot_clones_clustered,width=20,height=8)

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()
