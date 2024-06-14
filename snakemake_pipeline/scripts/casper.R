# ------------------------------------------------------------------------------
# Run CaSpER
# ------------------------------------------------------------------------------

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# ------------------------------------------------------------------------------
print("Load libraries")
# ------------------------------------------------------------------------------

library(CaSpER)
library(data.table)
library(Seurat)

# ------------------------------------------------------------------------------
print("Get input parameters from snakemake")
# ------------------------------------------------------------------------------

input_file<-snakemake@input$matrix
input_annotations<-snakemake@input$annot
input_ref_groups<-snakemake@input$ref_groups
input_gene_annot<-snakemake@input$gene_annot
input_af<-snakemake@input$af

input_segment_gamma<-as.numeric(snakemake@params$segment_gamma)

#In order to try different cutoffs, as default value is very strict
input_expr_cutoff<-snakemake@params$expr_cutoff
if(is.null(input_expr_cutoff)){
  #input_expr_cutoff<-4.5 # the default threshold
  input_expr_cutoff<-0.1 # from the 10X tutorial
} else {
  input_expr_cutoff<-as.numeric(input_expr_cutoff)
}
print(paste("Running CaSpEr with the following expression cutoff",input_expr_cutoff))

#output_casper_object<-snakemake@output$casper_object
#output_casper_segments<-snakemake@output$casper_segments
output_casper_cellmatrix<-snakemake@output$casper_cellmatrix
output_casper_pseudobulk<-snakemake@output$casper_pseudobulk
output_casper_hclust<-snakemake@output$casper_clusters

output_plot_density<-snakemake@output$casper_plot_density 
output_plot_large_events<-snakemake@output$casper_plot_large_events
output_plot_heatmap<-snakemake@output$casper_plot_heatmap
output_plot_baf<-snakemake@output$casper_plot_baf

# ------------------------------------------------------------------------------
print("Execute Casper")
# ------------------------------------------------------------------------------

#Load count matrix
count_matrix<-fread(input_file)

#Transform it into a numeric matrix
gene_names<-count_matrix$V1
count_matrix$V1<-NULL
count_matrix<-as.matrix(count_matrix)
rownames(count_matrix)<-gene_names

#Get cytoband information
data(hg38_cytoband)

#Load metadata
meta_data <- fread(input_annotations,header=FALSE)
colnames(meta_data)<-c("cell","sample")

#Load count matrix into seurat for filtering and normalization 
#(following the respective single cell tutorial)
count_seurat <- CreateSeuratObject(counts = count_matrix, project = "bcc", 
                                   min.cells = 3, min.features = 200)

count_seurat <- NormalizeData(count_seurat , scale.factor = 1e6, 
                              normalization.method = "RC")

log.ge <- as.matrix(count_seurat@assays$RNA@data)
log.ge <- log2(log.ge +1)

rm(count_seurat)
gc()

#Filter gene matrix to contain only annotated cells in the right order
annotation<-read.table(input_gene_annot,
                       header=TRUE,stringsAsFactors = FALSE)
annotation<-annotation[annotation$Gene %in% rownames(log.ge),]
log.ge <- log.ge[match(annotation$Gene,rownames(log.ge)) , ]
all(rownames(log.ge)==annotation$Gene)

#Read reference groups (saved in one tsv file)
ref_groups<-read.table(input_ref_groups,header=TRUE)

#Extract barcodes of healthy cells as controls
control_cells <- meta_data$cell[meta_data$sample %in% ref_groups$ref_groups]

loh <- readBAFExtractOutput(path=dirname(input_af), 
                            sequencing.type="single-cell",
                            suffix="af")

#Rename files to match sample name
names(loh)<-gsub("_BAFExtract.af","",names(loh))

#Create sample - barcode mapping for LOH annotation
loh_name_mapping <- data.frame(loh.name= meta_data$sample, 
                               sample.name=meta_data$cell)

#Initialize CaSpER object
object <- CreateCasperObject(raw.data=log.ge,
                             loh.name.mapping=loh_name_mapping, 
                             sequencing.type="single-cell", 
                             cnv.scale=3, loh.scale=3, 
                             expr.cutoff=input_expr_cutoff, 
                             filter="median", 
                             matrix.type="normalized",
                             annotation=annotation, method="iterative", 
                             loh=loh, 
                             control.sample.ids=control_cells, 
                             cytoband=cytoband_hg38)

# Print some general statistics
print("Dimensions of the raw count matrix:")
print(dim(object@raw.data))
print("Dimensions of the filtered count matrix:")
print(dim(object@data))

# Run CaSpER
pdf(file=output_plot_density)
final.objects <- runCaSpER(object, removeCentromere=T, cytoband=cytoband_hg38, 
                           method="iterative")
dev.off()

# #Save results
# saveRDS(final.objects,file=output_casper_object)

## plot large scale events
finalChrMat <- extractLargeScaleEvents(final.objects, thr=0.75)

plot.data <- reshape2::melt(finalChrMat)
plot.data$value2 <- "neutral"
plot.data$value2[plot.data$value > 0] <- "amplification"
plot.data$value2[plot.data$value < 0] <- "deletion"
plot.data$value2 <- factor(plot.data$value2, levels = c("amplification",
                                                        "deletion", "neutral"))
plot.data$Var2 <- factor(plot.data$Var2, levels = colnames(finalChrMat))
p <- ggplot(aes(x = Var2, y = Var1, fill = value2), data = plot.data) +
  geom_tile(colour = "white", size = 0.01)+
  labs(x = "",
       y = "") + scale_fill_manual(values = c(amplification = muted("red"),
                                              deletion = muted("blue"),
                                              neutral = "white"))+
  theme_grey(base_size = 6) +
  theme(legend.position = "right", legend.direction = "vertical",
        legend.title = element_blank(), strip.text.x = element_blank(),
        legend.text = element_text(colour = "black", size = 7,
                                   face = "bold"),
        legend.key.height = grid::unit(0.8, "cm"), legend.key.width = grid::unit(0.5, "cm"),
        axis.text.x = element_text(size = 5, colour = "black",
                                   angle = -45, hjust = 0),
        axis.text.y = element_text(size = 6, vjust = 0.2, colour = "black"),
        axis.ticks = element_line(size = 0.4),
        plot.title = element_text(colour = "black", hjust = 0,
                                  size = 6, face = "bold"))
ggsave(p,file=output_plot_large_events,
       width=12,height=10)

# ------------------------------------------------------------------------------
print("Get clustering")
# ------------------------------------------------------------------------------

#Heatmap with normalized expression
obj <- final.objects[[9]]
plotHeatmap10x(object=obj, fileName=output_plot_heatmap,
             cnv.scale= 3, cluster_cols = F,
             cluster_rows = T, show_rownames = F, only_soi = T)
print("Heatmap saved")

#Redo clustering as implemetned in plotHeatmap10x to be save directly the dendogram
cnv.scale <- 3
data <- obj@control.normalized.noiseRemoved[[cnv.scale]]
x.center <- mean(data)
quantiles = quantile(data[data != x.center], c(0.01, 0.99))
delta = max(abs(c(x.center - quantiles[1], quantiles[2] - 
                    x.center)))
low_threshold = x.center - delta
high_threshold = x.center + delta
x.range = c(low_threshold, high_threshold)
data[data < low_threshold] <- low_threshold
data[data > high_threshold] <- high_threshold
data <- data[, !(colnames(data) %in% obj@control.sample.ids)]
clusters<-hclust(dist(t(data), "euclidean"), method="complete")
saveRDS(clusters,file=output_casper_hclust)

# ------------------------------------------------------------------------------
print("plot BAF deviation")
# ------------------------------------------------------------------------------

plotBAFAllSamples (loh = final.objects[[9]]@loh.median.filtered.data,
                   fileName=output_plot_baf)

# ------------------------------------------------------------------------------
print("Extract CNV segments with score")
# ------------------------------------------------------------------------------

print(paste("Chosen thresold to define segments (larger or equal):",
            input_segment_gamma))
segment.summary <- extractSegmentSummary(final.objects)
loss <- segment.summary$all.summary.loss
gain <- segment.summary$all.summary.gain
#loh <- segment.summary$all.summary.loh #ignore loss of heterozygosity for the moment

loss.final <- loss[loss$count>=input_segment_gamma, ]
gain.final <- gain[gain$count>=input_segment_gamma, ]
segments<-rbind(loss.final,gain.final)
# write.table(segments,quote=FALSE, sep="\t",file=output_casper_segments)
# print("Saved segmentation results detailed!")

# ------------------------------------------------------------------------------
print("Create pseudobulk aggregate over all cells")
# ------------------------------------------------------------------------------

#Filter the segments for the cancer cells and create a Grange
cancer_cells <- meta_data$cell[!(meta_data$sample %in% ref_groups$ref_groups)]
segments <- segments[segments$ID %in% cancer_cells,]
seg_grange <-  GRanges(seqnames = Rle(paste0("chr",gsub("p|q", "", segments$seqnames))), 
                       IRanges(segments$start, segments$end)) 

#Create a Grange from the gene annotation 
#(use the filtered annotation not the original one!)
annot_filtered<-final.objects[[1]]@annotation.filt
annot_filtered$Chr<-paste0("chr",annot_filtered$Chr)
ann_gr <- makeGRangesFromDataFrame(annot_filtered, 
                                   keep.extra.columns = TRUE, seqnames.field="Chr")

#Get gene-wise CNV annotations following the code from the CaSpER tutorial
genes <- splitByOverlap(ann_gr, seg_grange, "GeneSymbol")
genes_ann <- lapply(genes, function(x) x[!(x=="")])
rna_matrix <- gene.matrix(seg=segments, all.genes=unique(annot_filtered$GeneSymbol), 
                          all.samples=cancer_cells, 
                          genes.ann=genes_ann)

#Save also the RNA matrix (per cell results)
write.table(rna_matrix,file=output_casper_cellmatrix,
            sep="\t",quote=FALSE)

#Get average results across all cells
ann_gr$mean_loss<-rowMeans(rna_matrix == -1)
ann_gr$mean_base<-rowMeans(rna_matrix == 0)
ann_gr$mean_gain<-rowMeans(rna_matrix == 1)

#Delete unused positions
ann_gr@elementMetadata[,c("Gene","band","cytoband","isCentromer",
                          "Position","new_positions")]<-NULL

saveRDS(ann_gr,file=output_casper_pseudobulk)

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()
