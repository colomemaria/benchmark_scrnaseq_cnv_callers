# ------------------------------------------------------------------------------
# Evaluate how similar the clones are identified by the different methods
# ------------------------------------------------------------------------------

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# ------------------------------------------------------------------------------
print("Load libraries")
# ------------------------------------------------------------------------------

library(data.table)
library(mclust) #for the adjusted Rand Index
library(ggplot2)
library(viridis)
library(infercnv)
library(CONICSmat)

source("scripts/lib.R")

# ------------------------------------------------------------------------------
print("Get input parameters from snakemake")
# ------------------------------------------------------------------------------

path_expr_mat<-"data/input_SNU601/count_matrix.txt"
path_annot<-"data/input_SNU601/sample_annotation.txt"
path_ref_groups<-"data/input_SNU601/ref_groups.txt"

path_clusters_copykat<-"results/output_SNU601/copykat/SNU601_copykat_clustering_results.rds"
input_matrix_copykat<-"results/output_SNU601/copykat/SNU601_copykat_CNA_raw_results_gene_by_cell.txt"

path_clones_numbat<-"results/output_SNU601/numbat/clone_post_2.tsv"

input_infercnv_cnv<-"results/output_SNU601/infercnv/infercnv.20_HMM_predHMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.repr_intensities.observations.txt"
input_infercnv_expr<-"results/output_SNU601/infercnv/infercnv.observations.txt"
input_infercnv_gene_pos<-"data/annotations/hg38_gencode_v27.txt"
path_clusters_infercnv<-"results/output_SNU601/infercnv/run.final.infercnv_obj"

path_clones_scevan<-"results/output_SNU601/scevan_find_clones/output/SNU601_scevan_prediction.txt"
input_scevan_expr<-"results/output_SNU601/scevan_find_clones//output/SNU601_CNAmtx.RData"
input_scevan_gene_pos<-"results/output_SNU601/scevan_find_clones//output/SNU601_count_mtx_annot.RData"
input_scevan_subclones<-"results/output_SNU601/scevan_find_clones/output/SNU601_subclone1_CN.seg"

path_l_CONICSmat<-"results/output_SNU601/CONICSmat/cnv_posterior_probs.tsv"
path_lr_CONICSmat<-"results/output_SNU601/CONICSmat/CONICSmat_CNV_BIC_LR.txt"
path_genes_CONICSmat<-"results/output_SNU601/CONICSmat/gene_annotation_conicsmat_SNU601.tsv"
input_CONICSmat_chrom_pos<-"data/annotations/chromosome_arm_positions_grch38.txt"
input_CONICSmat_cnv<-"results/output_SNU601/CONICSmat/cnv_types.tsv"

path_clusters_casper<-"results/output_SNU601/casper/SNU601_casper_hclust.RDS"
param_splitted_clusters<-6

plot_conicsmat_histo<-"results/output_SNU601/CONICSmat/histogram_clones.pdf"
output_text_ari
output_plot_ari

# ------------------------------------------------------------------------------
# Load files
# ------------------------------------------------------------------------------

#Load dendogram from copykat
clusters_copykat<-readRDS(path_clusters_copykat)
clusters_copykat$labels<-gsub("\\.","-",clusters_copykat$labels)

#Split copykat into the same number of groups (number defined as parameter)
groups_copykat<-cutree(clusters_copykat,k=param_splitted_clusters)
clustering_copykat<-data.frame(cell=names(groups_copykat),
                               copykat=groups_copykat)
all_clusters<-clustering_copykat

#Load clones from numbat
clones_numbat<-fread(path_clones_numbat)
clones_numbat<-clones_numbat[,c("cell","clone_opt")]
colnames(clones_numbat)<-c("cell","numbat")
all_clusters<-merge(all_clusters,clones_numbat,by="cell",all=TRUE)

#Load subclusters from InferCNV
infercnv_obj <- readRDS(path_clusters_infercnv)
subclusters_tumor<-infercnv_obj@tumor_subclusters$subclusters$SNU601
cell_names<-colnames(infercnv_obj@expr.data)
clustering_infercnv<-NULL
for(subc in names(subclusters_tumor)){
  clustering_infercnv<-rbind(clustering_infercnv,
            data.frame(cell=cell_names[subclusters_tumor[[subc]]],
                       infercnv=subc))
                              
}
rm(infercnv_obj, subclusters_tumor)
gc()
all_clusters<-merge(all_clusters,clustering_infercnv,by="cell",all=TRUE)

#SCEVAN (need to run specific mode to find subclusters)
clustering_scevan<-fread(path_clones_scevan)
#Remove filtered and normal cells (no subclone defined for those)
clustering_scevan<-clustering_scevan[clustering_scevan$class == "tumor",]
clustering_scevan<-clustering_scevan[,c("V1","subclone")]
colnames(clustering_scevan)<-c("cell","scevan")
all_clusters<-merge(all_clusters,clustering_scevan,by="cell",all=TRUE)

#CONICSmat (need to calculate the clusters on the fly)
l<-as.matrix(read.table(path_l_CONICSmat,sep="\t",header=T,row.names=1,check.names=F))
lrbic<-read.table(path_lr_CONICSmat,
                  sep="\t",header=T,row.names=1,check.names=F)
candRegions<-rownames(lrbic)[which(lrbic[,"BIC difference"]>200 &
                                     lrbic[,"LRT adj. p-val"]<0.01)]

#Requires additionally filtered expression matrix (using same code as in CONICSmat)
expr_df <- as.matrix(read.table(path_expr_mat,sep="\t",header=T,row.names=1,check.names=F))
gene_pos <- read.table(file = path_genes_CONICSmat, header = TRUE, sep = "\t")
expr_df[which(is.na(expr_df))] <- 0
expr_df <- filterMatrix(expr_df,gene_pos[,"hgnc_symbol"],minCells=5)

pdf(plot_conicsmat_histo)
hi<-plotHistogram(l[,candRegions],expr_df,
                  clusters=param_splitted_clusters,zscoreThreshold=4)
dev.off()

clusters_conicsmat<-data.frame(cell=names(hi),
                               conicsmat=hi)
all_clusters<-merge(all_clusters,clusters_conicsmat,by="cell",all=TRUE)

#Add casper results
hclust_casper<-readRDS(path_clusters_casper)
groups_casper<-cutree(hclust_casper,k=param_splitted_clusters)
clustering_casper<-data.frame(cell=names(groups_casper),
                               casper=groups_casper)
all_clusters<-merge(all_clusters,clustering_casper,by="cell",all=TRUE)

#Remove all healthy cells
annotation<-fread(path_annot, header=FALSE)
ref_groups<-read.table(path_ref_groups,header=TRUE)
cancer_cells <- annotation$V1[! (annotation$V2 %in% ref_groups$ref_groups)]

all_clusters<-all_clusters[all_clusters$cell %in% cancer_cells,]

# ------------------------------------------------------------------------------
print("Pairwise comparisons using adjusted rand index")
# ------------------------------------------------------------------------------

methods<-c("copykat","numbat","infercnv","scevan","conicsmat","casper")
compare_methods<-NULL
for(m1 in 1:(length(methods)-1)){
  for(m2 in (m1+1):length(methods)){
    compare_methods<-rbind(compare_methods,
                           data.frame(method1=methods[m1],
                                method2=methods[m2],
                                adj_rand_index=mclust::adjustedRandIndex(
                                  all_clusters[,methods[m1]],
                                  all_clusters[,methods[m2]])))
  }
  compare_methods<-rbind(compare_methods,
                         data.frame(method1=methods[m1],
                              method2=methods[m1],
                              adj_rand_index=1))
}
compare_methods<-rbind(compare_methods,
                       data.frame(method1=methods[m2],
                            method2=methods[m2],
                            adj_rand_index=1))

#Save result txt
write.table(compare_methods,file=output_text_ari,sep="\t",quote=FALSE,
            row.names=FALSE)

compare_methods$method1<-factor(compare_methods$method1,levels=methods)
compare_methods$method2<-factor(compare_methods$method2,levels=methods)

g<-ggplot(compare_methods,aes(x=method2,y=method1,fill=adj_rand_index))+
  geom_tile()+
  theme_bw()+
  geom_text(aes(label=round(adj_rand_index,3),
                color=ifelse(adj_rand_index<0.6,'white','black')),size=2)+
  scale_color_manual(values=c("black","white"))+
  scale_y_discrete(limits=rev)+
  xlab("Method")+
  ylab("Method")+
  scale_fill_viridis("Adj Rand\nIndex",limits=c(0,1))+
  guides(color="none")+
  theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1))

ggsave(g, file=output_plot_ari,
       width=6,height=4.5)

# ------------------------------------------------------------------------------
print("Get per subclone and method pseudobulk results")
# ------------------------------------------------------------------------------

print("Load copykat data (normalized expression)")
copykat_results<-read_copykat(input_matrix_copykat,path_annot,path_ref_groups,
                              clones=clustering_copykat)[[1]]
colnames(mcols(copykat_results))<-paste0("copykat_",
                                         colnames(mcols(copykat_results)))

#Create genomic range object
binned_genome<-create_binned_genome(copykat_results)
combined_clones<-combine_range_objects_allcols(binned_genome,copykat_results)

# print("Load SCEVAN data (normalized expression)")
# scevan_results<-read_scevan(input_scevan_expr,input_scevan_gene_pos,
#                             path_annot,path_ref_groups,
#                             clones=clustering_scevan)[[1]]
# colnames(mcols(scevan_results))<-paste0("scevan_expr_",
#                                          colnames(mcols(scevan_results)))
# combined_clones<-combine_range_objects_allcols(combined_clones,scevan_results)

print("Load SCEVAN data (segments per clone)")
scevan_results<-read_scevan_subclones(input_scevan_subclones,binned_genome)
colnames(mcols(scevan_results))<-paste0("scevan_",
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

print("Load Numbat data")
numbat_results <- read_numbat_cnv(input_numbat_obj)

print("Load CONICSmat data")
CONICSmat_results <- read_CONICSmat(input_CONICSmat_chrom_pos, input_CONICSmat_cnv,
                                    path_annot, path_ref_groups,
                                    clones=clusters_conicsmat)
colnames(mcols(CONICSmat_results))<-paste0("conicsmat_",
                                          colnames(mcols(CONICSmat_results)))
combined_clones<-combine_range_objects_allcols(combined_clones,CONICSmat_results)

# ------------------------------------------------------------------------------
print("Generate a heatmap per subclone and method")
# ------------------------------------------------------------------------------

#Get results
combined_methods<-elementMetadata(combined_clones)

pheatmap::pheatmap(t(combined_methods),scale="row",cluster_cols=FALSE, cluster_rows=TRUE,
                   clustering_method = "average",file="~/Desktop/test_heatmap.pdf")

#ordering_clones<-hclust(dist(combined_methods,method="euclidean"),method="average")
# #Add position information
# combined_methods$chr<-factor(combined_clones@seqnames,
#                              levels=combined_clones@seqnames@values)
# combined_methods$start_position<-combined_clones@ranges@start
# combined_methods<-as.data.frame(combined_methods)
# 
# #Add an artifical count through the whole genome and 
# #get start positions for each new chromosome
# combined_methods$counted_pos<-1:nrow(combined_methods)
# chr_boundries<-combined_methods%>%
#   group_by(chr)%>%
#   summarize(start_chr=min(counted_pos),
#             mean_chr=mean(counted_pos))
# 
# scaled_methods<-combined_methods
# 
# plot_data<-reshape2::melt(scaled_methods,
#                           id.vars=c("chr","start_position","counted_pos"))
# 
# g.2<-ggplot(plot_data,aes(x=counted_pos,y=variable,fill=value))+geom_tile()+
#   theme_bw()+
#   scale_fill_gradient2("Score",low = "darkblue",
#                        mid = "white",high = "darkred",midpoint = 0,
#                        breaks=c(-5,0,5))+
#   xlab("Chromosome position")+ylab("Method")+
#   geom_vline(xintercept = chr_boundries$start_chr)+
#   scale_x_continuous(breaks=chr_boundries$mean_chr,
#                      labels=chr_boundries$chr)+
#   scale_y_discrete(limits=rev)+
#   coord_cartesian(xlim=c(1,max(plot_data$counted_pos)),expand=FALSE)+
#   theme(legend.position="bottom",
#         axis.text.x = element_text(angle=90,vjust=0.5,hjust=1),
#         text=element_text(size=21))
# g<-ggarrange(g.1,g.2,ncol=1,align="v",heights=c(0.3,0.7))

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()
