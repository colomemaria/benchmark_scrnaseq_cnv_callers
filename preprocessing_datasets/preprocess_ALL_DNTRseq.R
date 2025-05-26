# ------------------------------------------------------------------------------
# Preprocess ALL samples ALL1 and ALL2 from DNTRseq paper
# Explanations from the first author about the shared data:
# Only cell identifiers for which DNA was sequenced are in the clone file. 
# Clone = NA means it was filtered out (low counts for example). 
# Clone = 0 is that it could not be reliably assigned. 
# ALL1 has 3 copy states: 1_1 (diploid), 2_3 (tumor) and 3_1 (tumor)
# ALL2 has 2 copy states: 2_1m (diploid) and 1_1m (tumor)
#
# ------------------------------------------------------------------------------

library(data.table)
library(ggplot2)
library(GenomicRanges)
library(dplyr)

#General gene annotation file
gene_annot<-fread("geneInfo.tab",skip=1,header=FALSE)

# ------------------------------------------------------------------------------
# Processing sample ALL 1
# ------------------------------------------------------------------------------

cnv<-fread("ALL1-cn_segments.txt")

clones<-fread("ALL1-clones_final.txt")
clones$clone[is.na(clones$clone)]<-"filtered"

ggplot(clones,aes(x=dna_reads,y=rna_counts,color=clone))+
  geom_point()+scale_x_log10()+scale_y_log10()

counts<-fread("ALL1-rna_counts.tsv.gz")
gene_names<-counts$gene_id
counts$gene_id<-NULL
counts<-as.matrix(counts)
rownames(counts)<-gene_names

#DNA library is labeled with D, RNA with R, so need to update the names
clones$rna_index<-gsub("D","R",clones$dna_library_id)
clones$rna_index<-gsub("_","-",clones$rna_index)

#Remove cells from RNA matrix without annotation and vice versa
clones<-clones[clones$rna_index %in% colnames(counts),]
counts<-counts[,clones$rna_index]

#Check whether RNA counts are matching
clones$rna_counts_new<-colSums(counts)

#Not identical but very close
ggplot(clones,aes(rna_counts,rna_counts_new))+
  geom_point()

#Remove filtered cells, cells with clone 0 and cells with less than 20,000 counts
clones<-clones[! clones$clone %in% c("0","filtered"),]
clones<-clones[clones$rna_counts_new > 20000,]

#Remove these also from the count matrix
counts<-counts[,clones$rna_index]

#Remove empty genes from matrix
counts<-counts[rowSums(counts)>0,]

#Shift count matrix from ENSG to gene symbols (to fit the pipeline)
gene_annot_all1<-data.frame(gene_annot)
rownames(gene_annot_all1)<-gene_annot_all1$V1
gene_annot_all1<-gene_annot_all1[rownames(counts),]

#Merge counts for ENSG IDs with same gene symbol
counts<-apply(counts, 2, tapply, as.factor(gene_annot_all1$V2),sum, na.rm=T)

#Save sample annotation
sample_annot<-data.frame(barcode=clones$rna_index,
                         sample=ifelse(clones$clone=="1_1","diploid","ALL1"))

write.table(sample_annot, file="input_ALL1/sample_annotation.txt",
            sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)


df_refs<-data.frame(ref_groups=setdiff(unique(sample_annot$sample),"ALL1"))
write.table(df_refs, file="input_ALL1/ref_groups.txt",
            sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)

#Save the reduced count matrix
write.table(counts,
            file="input_ALL1/count_matrix.txt",
            sep="\t",quote=FALSE,row.names=TRUE)

#Reformat the CNVs into 100kb bins
cnv_grange<-makeGRangesFromDataFrame(cnv, start.field="start.pos",
                                     end.field="end.pos")

#Add CNV information (weighted mean between both clones)
perc_clone1<-sum(clones$clone=="2_3") / sum(clones$clone %in% c("2_3","3_1"))
perc_clone2<-sum(clones$clone=="3_1") / sum(clones$clone %in% c("2_3","3_1"))
cnv_grange$wgs_gain<-perc_clone1 * as.numeric(cnv[,"2_3"]>2) + perc_clone2 * as.numeric(cnv[,"3_1"]>2)
cnv_grange$wgs_base<-perc_clone1 * as.numeric(cnv[,"2_3"]==2) + perc_clone2 * as.numeric(cnv[,"3_1"]==2)
cnv_grange$wgs_loss<-perc_clone1 * as.numeric(cnv[,"2_3"]<2) + perc_clone2 * as.numeric(cnv[,"3_1"]<2)


#Determine the annotation region of each chromsome
cnv_limits<-cnv%>%group_by(chr)%>%
  summarize(global_start=min(start.pos),global_end=max(end.pos))

#Round it to 100,000 borders (to have bins starting at 0)
cnv_limits$global_start_rounded<-floor(cnv_limits$global_start/100000)*100000
cnv_limits$global_end_rounded<-ceiling(cnv_limits$global_end/100000)*100000-1

#Use Genomic Ranges to create 100,000 kB windows
cnv_limits_gr<-makeGRangesFromDataFrame(cnv_limits,
                                        seqnames.field = "chr",
                                        start.field="global_start_rounded",
                                        end.field="global_end_rounded")
cnv_limits_gr<-cnv_limits_gr[order(cnv_limits_gr)]

cnv_bins<-tile(cnv_limits_gr,width=100000)
cnv_bins<-unlist(cnv_bins)

#Map it back to the genomic segments
overlaps<-as.data.frame(findOverlaps(cnv_bins,cnv_grange,type="any"))
cnv_bins<-cnv_bins[unique(overlaps$queryHits)]
cnv_bins$gain_wgs <- sapply(unique(overlaps$queryHits),
                            function(i) mean(cnv_grange$wgs_gain[overlaps$subjectHits[overlaps$queryHits==i]]))
cnv_bins$loss_wgs <- sapply(unique(overlaps$queryHits),
                            function(i) mean(cnv_grange$wgs_loss[overlaps$subjectHits[overlaps$queryHits==i]]))
cnv_bins$base_wgs <- sapply(unique(overlaps$queryHits),
                            function(i) mean(cnv_grange$wgs_base[overlaps$subjectHits[overlaps$queryHits==i]]))

#Create also a clone specific matrix (part 1)
cnv_bins_per_clone<-cnv_bins
elementMetadata(cnv_bins_per_clone)<-NULL

#Remove sex chromosomes
cnv_bins<-cnv_bins[! cnv_bins@seqnames %in% c("chrX","chrY"),]

#Save the results
cnv_bins<-as.data.frame(cnv_bins)
cnv_bins$strand<-NULL
cnv_bins$width<-NULL
colnames(cnv_bins)[1]<-"chr"

#Remove the "chr" to match the Aneufinder chromosome encoding
cnv_bins$chr<-gsub("chr","",cnv_bins$chr)

write.table(cnv_bins,file="ALL1_WGS_groundtruth/wgs_results_formated.csv",
            sep=",",quote=FALSE)

#Create also a clone specific matrix (part 2)
cnv<-as.data.frame(cnv)
cnv_bins_per_clone$clone2_3<-sapply(unique(overlaps$queryHits),
                                    function(i) mean(cnv[overlaps$subjectHits[overlaps$queryHits==i],"2_3"]))
cnv_bins_per_clone$clone3_1<-sapply(unique(overlaps$queryHits),
                                    function(i) mean(cnv[overlaps$subjectHits[overlaps$queryHits==i],"3_1"]))

clone_matrix<-matrix(c(rep(cnv_bins_per_clone$clone2_3,sum(clones$clone=="2_3")),
                       rep(cnv_bins_per_clone$clone3_1,sum(clones$clone=="3_1"))),
                     ncol=sum(clones$clone %in% c("2_3","3_1")))
colnames(clone_matrix)<-c(clones$rna_index[clones$clone == "2_3"],
                          clones$rna_index[clones$clone == "3_1"])
clone_res<-as.data.frame(cnv_bins_per_clone)
clone_res<-cbind(clone_res[,c("seqnames","start","end")],clone_matrix)

#Filter sex chromsomes and NA values (NA in all cells in these cases)
clone_res<-clone_res[! clone_res$seqnames %in% c("chrX","chrY"),]
clone_res<-clone_res[! is.na(clone_res$`VZA01001R-A01`),]

write.table(clone_res,file="ALL1_WGS_groundtruth/cnvs_per_cell_filtered.tsv",
            quote=FALSE,sep="\t",row.names=FALSE)

# ------------------------------------------------------------------------------
# Processing sample ALL 2
# ------------------------------------------------------------------------------

cnv<-fread("ALL2-cn_segments.txt")

clones<-fread("ALL2-clones_final.txt")
clones$clone[is.na(clones$clone)]<-"filtered"

ggplot(clones,aes(x=dna_reads,y=rna_counts,color=clone))+
  geom_point()+scale_x_log10()+scale_y_log10()


counts<-fread("ALL2-rna_counts.tsv.gz")
gene_names<-counts$gene_id
counts$gene_id<-NULL
counts<-as.matrix(counts)
rownames(counts)<-gene_names

#DNA library is labeled with D, RNA with R, so need to update the names
clones$rna_index<-gsub("D","R",clones$dna_library_id)
clones$rna_index<-gsub("_","-",clones$rna_index)

#Remove cells from RNA matrix without annotation and vice versa
clones<-clones[clones$rna_index %in% colnames(counts),]
counts<-counts[,clones$rna_index]

#Check whether RNA counts are matching
clones$rna_counts_new<-colSums(counts)

#Not identical but very close
ggplot(clones,aes(rna_counts,rna_counts_new))+
  geom_point()

#Remove filtered cells, cells with clone 0 and cells with less than 20,000 counts
clones<-clones[! clones$clone %in% c("0","filtered"),]
clones<-clones[clones$rna_counts_new > 20000,]

#Remove these also from the count matrix
counts<-counts[,clones$rna_index]

#Remove empty genes from matrix
counts<-counts[rowSums(counts)>0,]

#Shift count matrix from ENSG to gene symbols (to fit the pipeline)
gene_annot_all2<-data.frame(gene_annot)
rownames(gene_annot_all2)<-gene_annot_all2$V1
gene_annot_all2<-gene_annot_all2[rownames(counts),]

#Merge counts for ENSG IDs with same gene symbol
counts<-apply(counts, 2, tapply, as.factor(gene_annot_all2$V2),sum, na.rm=T)

#Save sample annotation
sample_annot<-data.frame(barcode=clones$rna_index,
                         sample=ifelse(clones$clone=="2_1m","diploid","ALL2"))

write.table(sample_annot, file="input_ALL2/sample_annotation.txt",
            sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)


df_refs<-data.frame(ref_groups=setdiff(unique(sample_annot$sample),"ALL2"))
write.table(df_refs, file="input_ALL2/ref_groups.txt",
            sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)

#Save the reduced count matrix
write.table(counts,
            file="input_ALL2/count_matrix.txt",
            sep="\t",quote=FALSE,row.names=TRUE)

#Reformat the CNVs into 100kb bins
cnv_grange<-makeGRangesFromDataFrame(cnv, start.field="start.pos",
                                     end.field="end.pos")

#Add CNV information
cnv_grange$wgs_gain<-as.numeric(cnv[,"1_1m"]>2)
cnv_grange$wgs_base<-as.numeric(cnv[,"1_1m"]==2)
cnv_grange$wgs_loss<-as.numeric(cnv[,"1_1m"]<2)


#Determine the annotation region of each chromsome
cnv_limits<-cnv%>%group_by(chr)%>%
  summarize(global_start=min(start.pos),global_end=max(end.pos))

#Round it to 100,000 borders (to have bins starting at 0)
cnv_limits$global_start_rounded<-floor(cnv_limits$global_start/100000)*100000
cnv_limits$global_end_rounded<-ceiling(cnv_limits$global_end/100000)*100000-1

#Use Genomic Ranges to create 100,000 kB windows
cnv_limits_gr<-makeGRangesFromDataFrame(cnv_limits,
                                        seqnames.field = "chr",
                                        start.field="global_start_rounded",
                                        end.field="global_end_rounded")
cnv_limits_gr<-cnv_limits_gr[order(cnv_limits_gr)]

cnv_bins<-tile(cnv_limits_gr,width=100000)
cnv_bins<-unlist(cnv_bins)


#Map it back to the genomic segments
overlaps<-as.data.frame(findOverlaps(cnv_bins,cnv_grange,type="any"))
cnv_bins<-cnv_bins[unique(overlaps$queryHits)]
cnv_bins$gain_wgs <- sapply(unique(overlaps$queryHits),
                            function(i) mean(cnv_grange$wgs_gain[overlaps$subjectHits[overlaps$queryHits==i]]))
cnv_bins$loss_wgs <- sapply(unique(overlaps$queryHits),
                            function(i) mean(cnv_grange$wgs_loss[overlaps$subjectHits[overlaps$queryHits==i]]))
cnv_bins$base_wgs <- sapply(unique(overlaps$queryHits),
                            function(i) mean(cnv_grange$wgs_base[overlaps$subjectHits[overlaps$queryHits==i]]))


#Remove sex chromosomes
cnv_bins<-cnv_bins[! cnv_bins@seqnames %in% c("chrX","chrY"),]

#Save the results
cnv_bins<-as.data.frame(cnv_bins)
cnv_bins$strand<-NULL
cnv_bins$width<-NULL
colnames(cnv_bins)[1]<-"chr"

#Remove the "chr" to match the Aneufinder chromosome encoding
cnv_bins$chr<-gsub("chr","",cnv_bins$chr)

write.table(cnv_bins,file="ALL2_WGS_groundtruth/wgs_results_formated.csv",
            sep=",",quote=FALSE)


#Combine loss, base and gain to one number
cnv_bins$wgs_mean<-with(cnv_bins,loss_wgs+(base_wgs*2)+gain_wgs*3)

#Save the per cell clone matrix (here completely redundant of course)
clone_matrix<-matrix(rep(cnv_bins$wgs_mean,sum(clones$clone=="1_1m")),
                     ncol=sum(clones$clone=="1_1m"))
colnames(clone_matrix)<-clones$rna_index[clones$clone == "1_1m"]
clone_res<-as.data.frame(cnv_bins)
clone_res<-cbind(clone_res[,c("chr","start","end")],clone_matrix)
colnames(clone_res)[1]<-"seqnames"

#Add the chr again
clone_res$seqnames<-paste0("chr",clone_res$seqnames)

#Filter NA values (NA in all cells in these cases)
clone_res<-clone_res[! is.na(clone_res$`VZA01101R-A05`),]

write.table(clone_res,file="ALL2_WGS_groundtruth/cnvs_per_cell_filtered.tsv",
            quote=FALSE,sep="\t",row.names=FALSE)
