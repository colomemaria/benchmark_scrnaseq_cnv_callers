# ------------------------------------------------------------------------------
# Check exemplarily for the MM dataset how similar the gene annotations are
# in order to decide how to best further process them
# Remark: CONICSmat gives only gene-level predictions
# ------------------------------------------------------------------------------

library(data.table)
library(ggplot2)
library(viridis)
library(numbat)

# Load annotation InferCNV
infercnv<-fread("snakemake_pipeline/data/annotations/hg38_gencode_v27.txt")
colnames(infercnv) <- c("gene","chr_infercnv","start_infercnv","end_infercnv")

#Load annotation CopyKat
copykat<-fread("snakemake_pipeline/results/output_MM/copykat/MM_copykat_CNA_raw_results_gene_by_cell.txt")
copykat<-copykat[,c(2:6)]

#Check that the hgnc_symbol is unique
length(copykat$hgnc_symbol) == length(unique(copykat$hgnc_symbol))

copykat$ensembl_gene_id<-NULL
colnames(copykat)<-c("chr_copykat","start_copykat","end_copykat","gene")
#Replace chr23 with chrX
copykat$chr_copykat<-paste0("chr",copykat$chr_copykat)
copykat$chr_copykat[copykat$chr_copykat=="chr23"]<-"chrX"

#Combine and compare the two annotations
combined_res<-merge(infercnv,copykat,by="gene")
mean(combined_res$chr_infercnv==combined_res$chr_copykat)
mean(combined_res$start_infercnv==combined_res$start_copykat)
mean(combined_res$end_infercnv==combined_res$end_copykat)

#Load annotation SCEVAN
load("snakemake_pipeline/results/output_MM/scevan/output/MM_count_mtx_annot.RData") #count_mtx_annot
scevan<-count_mtx_annot
rm(count_mtx_annot)

#Check that the gene name is unique
length(scevan$gene_name) == length(unique(scevan$gene_name))

scevan$gene_id<-NULL
colnames(scevan)<-c("chr_scevan","start_scevan","end_scevan","gene")
scevan$chr_scevan<-paste0("chr",scevan$chr_scevan)

#Combine and compare the two annotations
combined_res<-merge(combined_res,scevan,by="gene")

#Check how many genes are still there
nrow(combined_res)

mean(combined_res$chr_infercnv==combined_res$chr_scevan)
mean(combined_res$start_infercnv==combined_res$start_scevan)
mean(combined_res$end_infercnv==combined_res$end_scevan)
mean(combined_res$chr_copykat==combined_res$chr_scevan)
mean(combined_res$start_copykat==combined_res$start_scevan)
mean(combined_res$end_copykat==combined_res$end_scevan)

#Load annotation Casper
casper<-fread("snakemake_pipeline/results/output_MM/casper/gene_annotation_casper_MM.tsv")
casper<-casper[,c(2:5)]
colnames(casper)<-c("gene","chr_casper","start_casper","end_casper")
casper$chr_casper<-paste0("chr",casper$chr_casper)

#Combine and compare the two annotations
combined_res<-merge(combined_res,casper,by="gene")

#Check how many genes are still there
nrow(combined_res)

mean(combined_res$chr_infercnv==combined_res$chr_casper)
mean(combined_res$start_infercnv==combined_res$start_casper)
mean(combined_res$end_infercnv==combined_res$end_casper)

#Check how distant they are on average
combined_res$dist_annot<-combined_res$start_infercnv-combined_res$start_casper
hist(log(abs(combined_res$dist_annot)))

#Load annotation numbat
numbat_obj<-Numbat$new("snakemake_pipeline/results/output_MM/numbat")
numbat<-as.data.frame(numbat_obj$gtf)
numbat<-numbat[,1:4]
colnames(numbat)<-c("gene","start_numbat","end_numbat","chr_numbat")
numbat$chr_numbat<-paste0("chr",numbat$chr_numbat)

#Combine and compare the two annotations
combined_res<-merge(combined_res,numbat,by="gene")

#Check how many genes are still there
nrow(combined_res)

mean(combined_res$chr_infercnv==combined_res$chr_numbat)
mean(combined_res$start_infercnv==combined_res$start_numbat)
mean(combined_res$end_infercnv==combined_res$end_numbat)

mean(combined_res$chr_casper==combined_res$chr_numbat)
mean(combined_res$start_casper==combined_res$start_numbat)
mean(combined_res$end_casper==combined_res$end_numbat)

#To visualize it somehow: plot to what extend the start position matches
start_comp<-NULL
start_positions<-paste0("start_",c("infercnv","copykat","scevan","casper","numbat"))
for(m1 in 1:(length(start_positions)-1)){
  for(m2 in (m1+1):length(start_positions)){
    start_comp<-rbind(start_comp,
                      data.frame(method1=gsub("start_","",start_positions[m1]),
                           method2=gsub("start_","",start_positions[m2]),
                           overlap=mean(combined_res[[start_positions[m1]]]==
                                          combined_res[[start_positions[m2]]]),
                           overlap_100=mean(abs(combined_res[[start_positions[m1]]]-
                                              combined_res[[start_positions[m2]]])<100)))
  }
  
  #Add the diagonal for easier visualization
  start_comp<-rbind(start_comp,
                    data.frame(method1=gsub("start_","",start_positions[m1]),
                               method2=gsub("start_","",start_positions[m1]),
                               overlap=1,
                               overlap_100=1))
}

#Add the diagonal for easier visualization
start_comp<-rbind(start_comp,
                  data.frame(method1=gsub("start_","",start_positions[m2]),
                             method2=gsub("start_","",start_positions[m2]),
                             overlap=1,
                             overlap_100=1))

#Create plots
start_comp$method1<-factor(start_comp$method1,
                           levels=c("infercnv","copykat","scevan","casper","numbat"))
start_comp$method2<-factor(start_comp$method2,
                           levels=c("infercnv","copykat","scevan","casper","numbat"))

g<-ggplot(start_comp,aes(x=method2,y=method1,fill=overlap))+
  geom_tile()+
  theme_bw()+
  geom_text(aes(label=round(overlap,3),
                color=ifelse(overlap<0.6,'white','black')),size=2)+
  scale_color_manual(values=c("black","white"))+
  scale_y_discrete(limits=rev)+
  xlab("Method")+
  ylab("Method")+
  scale_fill_viridis("Overlap\nstart pos",limits=c(0,1))+
  guides(color="none")+
  theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1))

ggsave(g, file="snakemake_pipeline/results/output_MM/evaluation/overlap_gene_annots_MM.png",
       width=4,height=3)

g<-ggplot(start_comp,aes(x=method2,y=method1,fill=overlap_100))+
  geom_tile()+
  theme_bw()+
  geom_text(aes(label=round(overlap_100,3),
                color=ifelse(overlap_100<0.6,'white','black')),size=2)+
  scale_color_manual(values=c("black","white"))+
  scale_y_discrete(limits=rev)+
  xlab("Method")+
  ylab("Method")+
  scale_fill_viridis("Overlap\nstart pos",limits=c(0,1))+
  guides(color="none")+
  theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1))

ggsave(g, file="snakemake_pipeline/results/output_MM/evaluation/overlap_gene_100bp_annots_MM.png",
       width=4,height=3)
