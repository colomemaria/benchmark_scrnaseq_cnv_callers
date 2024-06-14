# ------------------------------------------------------------------------------
# Processing ccRCC data to get cell type annotations
# 
# Code downloaded from github 
# https://github.com/lessonskit/Single-cell-multi-omics-profiling-of-ccRCC
# associated with the publication
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9978887/
#
# Slightly adapted to get it running ....
# ------------------------------------------------------------------------------

###scRNA-seq

library(Seurat)
library(magrittr)
library(cowplot)
library(harmony)
library(dplyr)

data_path<-"~/Documents/CNV_RNAseq_benchmark/data/ccRCC_2023/scRNA-seq/"
  
K81.data <- Read10X(data.dir = paste0(data_path,"RCC81"))
kid81 <- CreateSeuratObject(counts = K81.data, project = "mRCC81", min.cells = 8, min.features = 500)
kid81[["percent.mt"]] <- PercentageFeatureSet(kid81, pattern = "^MT-")
#VlnPlot(kid81, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
kid81 <- subset(kid81, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 10)

K84.data <- Read10X(data.dir = paste0(data_path,"RCC84"))
kid84 <- CreateSeuratObject(counts = K84.data, project = "mRCC84", min.cells = 7, min.features = 500)
kid84[["percent.mt"]] <- PercentageFeatureSet(kid84, pattern = "^MT-")
#VlnPlot(kid84, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
kid84 <- subset(kid84, subset = nFeature_RNA > 500 & nFeature_RNA < 4700 & percent.mt < 10)

K86.data <- Read10X(data.dir = paste0(data_path,"RCC86"))
kid86 <- CreateSeuratObject(counts = K86.data, project = "mRCC86", min.cells = 9, min.features = 500)
kid86[["percent.mt"]] <- PercentageFeatureSet(kid86, pattern = "^MT-")
#VlnPlot(kid86, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
kid86 <- subset(kid86, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10)

K87.data <- Read10X(data.dir = paste0(data_path,"RCC87"))
kid87 <- CreateSeuratObject(counts = K87.data, project = "mRCC87", min.cells = 9, min.features = 500)
kid87[["percent.mt"]] <- PercentageFeatureSet(kid87, pattern = "^MT-")
#VlnPlot(kid87, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
kid87 <- subset(kid87, subset = nFeature_RNA > 500 & nFeature_RNA < 2900 & percent.mt < 10)

K94.data <- Read10X(data.dir = paste0(data_path,"RCC94"))
kid94 <- CreateSeuratObject(counts = K94.data, project = "mRCC94", min.cells = 8, min.features = 500)
kid94[["percent.mt"]] <- PercentageFeatureSet(kid94, pattern = "^MT-")
#VlnPlot(kid94, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
kid94 <- subset(kid94, subset = nFeature_RNA > 500 & nFeature_RNA < 3700 & percent.mt < 10)

K96.data <- Read10X(data.dir = paste0(data_path,"RCC96"))
kid96 <- CreateSeuratObject(counts = K96.data, project = "mRCC96", min.cells = 8, min.features = 500)
kid96[["percent.mt"]] <- PercentageFeatureSet(kid96, pattern = "^MT-")
#VlnPlot(kid96, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
kid96 <- subset(kid96, subset = nFeature_RNA > 500 & nFeature_RNA < 3300 & percent.mt < 10)

K99.data <- Read10X(data.dir = paste0(data_path,"RCC99"))
kid99 <- CreateSeuratObject(counts = K99.data, project = "mRCC99", min.cells = 9, min.features = 500)
kid99[["percent.mt"]] <- PercentageFeatureSet(kid99, pattern = "^MT-")
#VlnPlot(kid99, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
kid99 <- subset(kid99, subset = nFeature_RNA > 500 & nFeature_RNA < 2800 & percent.mt < 10)

K100.data <- Read10X(data.dir = paste0(data_path,"RCC100"))
kid100 <- CreateSeuratObject(counts = K100.data, project = "mRCC100", min.cells = 10, min.features = 500)
kid100[["percent.mt"]] <- PercentageFeatureSet(kid100, pattern = "^MT-")
#VlnPlot(kid100, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
kid100 <- subset(kid100, subset = nFeature_RNA > 500 & nFeature_RNA < 3800 & percent.mt < 10)

K101.data <- Read10X(data.dir = paste0(data_path,"RCC101"))
kid101 <- CreateSeuratObject(counts = K101.data, project = "mRCC101", min.cells = 10, min.features = 500)
kid101[["percent.mt"]] <- PercentageFeatureSet(kid101, pattern = "^MT-")
#VlnPlot(kid101, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
kid101 <- subset(kid101, subset = nFeature_RNA > 500 & nFeature_RNA < 2400 & percent.mt < 10)

K103.data <- Read10X(data.dir = paste0(data_path,"RCC103"))
kid103 <- CreateSeuratObject(counts = K103.data, project = "mRCC103", min.cells = 10, min.features = 500)
kid103[["percent.mt"]] <- PercentageFeatureSet(kid103, pattern = "^MT-")
#VlnPlot(kid103, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
kid103 <- subset(kid103, subset = nFeature_RNA > 500 & nFeature_RNA < 3300 & percent.mt < 10)

K104.data <- Read10X(data.dir = paste0(data_path,"RCC104"))
kid104 <- CreateSeuratObject(counts = K104.data, project = "mRCC104", min.cells = 10, min.features = 500)
kid104[["percent.mt"]] <- PercentageFeatureSet(kid104, pattern = "^MT-")
#VlnPlot(kid104, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
kid104 <- subset(kid104, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10)

K106.data <- Read10X(data.dir = paste0(data_path,"RCC106"))
kid106 <- CreateSeuratObject(counts = K106.data, project = "mRCC106", min.cells = 10, min.features = 500)
kid106[["percent.mt"]] <- PercentageFeatureSet(kid106, pattern = "^MT-")
#VlnPlot(kid106, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
kid106 <- subset(kid106, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & percent.mt < 10)

K112.data <- Read10X(data.dir = paste0(data_path,"RCC112"))
kid112 <- CreateSeuratObject(counts = K112.data, project = "mRCC112", min.cells = 10, min.features = 500)
kid112[["percent.mt"]] <- PercentageFeatureSet(kid112, pattern = "^MT-")
#VlnPlot(kid112, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
kid112 <- subset(kid112, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & percent.mt < 10)

K113.data <- Read10X(data.dir = paste0(data_path,"RCC113"))
kid113 <- CreateSeuratObject(counts = K113.data, project = "mRCC113", min.cells = 9, min.features = 500)
kid113[["percent.mt"]] <- PercentageFeatureSet(kid113, pattern = "^MT-")
#VlnPlot(kid113, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
kid113 <- subset(kid113, subset = nFeature_RNA > 500 & nFeature_RNA < 4500 & percent.mt < 10)

K114.data <- Read10X(data.dir = paste0(data_path,"RCC114"))
kid114 <- CreateSeuratObject(counts = K114.data, project = "mRCC114", min.cells = 7, min.features = 500)
kid114[["percent.mt"]] <- PercentageFeatureSet(kid114, pattern = "^MT-")
#VlnPlot(kid114, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
kid114 <- subset(kid114, subset = nFeature_RNA > 500 & nFeature_RNA < 3200 & percent.mt < 10)

K115.data <- Read10X(data.dir = paste0(data_path,"RCC115"))
kid115 <- CreateSeuratObject(counts = K115.data, project = "mRCC115", min.cells = 9, min.features = 500)
kid115[["percent.mt"]] <- PercentageFeatureSet(kid115, pattern = "^MT-")
#VlnPlot(kid115, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
kid115 <- subset(kid115, subset = nFeature_RNA > 500 & nFeature_RNA < 2400 & percent.mt < 10)

K116.data <- Read10X(data.dir = paste0(data_path,"RCC116"))
kid116 <- CreateSeuratObject(counts = K116.data, project = "mRCC116", min.cells = 6, min.features = 500)
kid116[["percent.mt"]] <- PercentageFeatureSet(kid116, pattern = "^MT-")
#VlnPlot(kid116, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
kid116 <- subset(kid116, subset = nFeature_RNA > 500 & nFeature_RNA < 4500 & percent.mt < 10)

K119.data <- Read10X(data.dir = paste0(data_path,"RCC119"))
kid119 <- CreateSeuratObject(counts = K119.data, project = "mRCC119", min.cells = 10, min.features = 500)
kid119[["percent.mt"]] <- PercentageFeatureSet(kid119, pattern = "^MT-")
#VlnPlot(kid119, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
kid119 <- subset(kid119, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & percent.mt < 10)

K120.data <- Read10X(data.dir = paste0(data_path,"RCC120"))
kid120 <- CreateSeuratObject(counts = K120.data, project = "mRCC120", min.cells = 10, min.features = 500)
kid120[["percent.mt"]] <- PercentageFeatureSet(kid120, pattern = "^MT-")
#VlnPlot(kid120, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
kid120 <- subset(kid120, subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 10)

#Remove unused data
rm(K81.data,K84.data,K86.data,K87.data,K94.data,K96.data,K99.data,
   K100.data,K101.data,K103.data,K104.data,K106.data,K112.data,K113.data,
   K114.data,K115.data,K116.data,K119.data,K120.data)
gc()

#Merge the individual plots
mRCC <- merge(x = kid81, y = list(kid84, kid86, kid87, kid94, kid96, kid99, 
                                  kid100, kid101, kid103, kid104, kid106, kid112, kid113, 
                                  kid114, kid115, kid116, kid119, kid120))

#Remove unused data
rm(kid81,kid84, kid86, kid87, kid94, kid96, kid99, 
    kid100, kid101, kid103, kid104, kid106, kid112, kid113, 
    kid114, kid115, kid116, kid119, kid120)
gc()

mRCC[["percent.mt"]] <- PercentageFeatureSet(mRCC, pattern = "^MT-")

#Analysis plots (currently not saved, could be saved with ggsave)
# VlnPlot(mRCC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size= 0)
# plot1 <- FeatureScatter(mRCC, feature1 = "nCount_RNA", feature2 = "percent.mt")
# plot2 <- FeatureScatter(mRCC, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# CombinePlots(plots = list(plot1, plot2))

#Strange that the normalization is done twice .. 
#mRCC <- NormalizeData(mRCC, normalization.method = "LogNormalize", scale.factor = 10000)
mRCC <- NormalizeData(mRCC)

mRCC <- FindVariableFeatures(mRCC, selection.method = "vst", nfeatures = 2000)

# top10 <- head(VariableFeatures(mRCC), 10)
# plot1 <- VariableFeaturePlot(mRCC)
# plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
# CombinePlots(plots = list(plot1, plot2))

#Get cell cycle genes saved in Seurat and correct for cell cycle score
s.genes <-cc.genes$s.genes
g2m.genes<-cc.genes$g2m.genes
mRCC <- CellCycleScoring(mRCC, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
mRCC <- ScaleData(mRCC, vars.to.regress = c("S.Score", "G2M.Score"))

mRCC <- RunPCA(mRCC, pc.genes = mRCC@var.genes, npcs = 30, verbose = FALSE)
options(repr.plot.height = 2.5, repr.plot.width = 6)
mRCC <- mRCC %>% 
    RunHarmony("orig.ident", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(mRCC, 'harmony')
#harmony_embeddings[1:5, 1:5]
mRCC <- mRCC %>% 
    RunUMAP(reduction = "harmony", dims = 1:30) %>% 
    FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
    FindClusters(resolution = 0.6) %>% 
    identity()

DimPlot(mRCC, reduction = "umap", label = TRUE, pt.size = .5)

#Find all the differences between clusters
#mRCC.markers <- FindAllMarkers(mRCC, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
#write.table(mRCC.markers,sep="\t",file="/data_8t/file/teacher/yuzhenyuan/ATAC_mRNA/metadata/mRNA19_PC30.xls")

#Renaming clusters doesn't work for me
#new.cluster.ids <- c(1, 2, 3, 4, 5, 6, 7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22)
#names(new.cluster.ids) <- levels(mRCC)
#mRCC <- RenameIdents(mRCC, new.cluster.ids)

#DimPlot(mRCC, reduction = "umap", label = TRUE, pt.size = .5)
#DimPlot(mRCC, reduction = "umap", split.by = "orig.ident",label = TRUE, pt.size = .5, ncol = 5)
##save object
save(mRCC, file= paste0(data_path,"mRNA19_merged_seurat_PC30.RData"))

#Marker genes
features1 <- c("NDUFA4L2","CA9","KRT18","KRT8","KLRD1","KLRB1","GNLY",
               "NKG7","CD3D","CD3E","CD8A","CD8B","IL7R","CD68","CD163","GPNMB",
               "SLC40A1","MSR1","ACTA2","PDGFRB","COL1A2","PECAM1","KDR","CDH5",
               "S100A8","S100A9","LYZ","KCNQ1OT1","CP","CTLA4","FOXP3","CD1C","CD1E",
               "ACKR1","VWF","MKI67","TOP2A","CD79A","CD79B","IGKC","IGLC2","MS4A1",
               "TPSB2","TPSAB1","KIT","KRT19","WFDC2")
g<-DotPlot(mRCC, features = features1, cols = c("blue", "red"), dot.scale = 8) + RotatedAxis()
ggsave(g,file=paste0(data_path,"marker_genes.png"),width=14,height=6)

#Add assign cell types based on marker genes
new.cluster.ids <- c("0 NK cells","1 ccRCC 1", 
                     "2 Exhau CD8+ T", "3 CD4+ T cells", 
                     "4 TAM 1", "5 CAF", "6 Endo cell 1", 
                     "7 Monocytes","8 ccRCC 2","9 Treg cells",
                     "10 TAM 2","11 DCs",
                     "12 ccRCC 4","13 Pro CD8+ T","14 Endo cell 2","15 ccRCC 3","16 Plasma cells",
                     "17 B cells", "18 Mast cells",
                     "19 Plasma cells","20 E-P CD8+ T")
names(new.cluster.ids) <- levels(mRCC)
mRCC <- RenameIdents(mRCC, new.cluster.ids)
g<-DimPlot(mRCC, reduction = "umap", label = TRUE, pt.size = .5)
g<-g+theme(legend.position = "none")
ggsave(g,file=paste0(data_path,"umap_annotated_ccRCC.png"),
       width=10,height=6)

#Create a data frame with cell type annotations per cell and sample
mRCC@meta.data$cell_type<-Idents(mRCC)
rename_cts<-setNames(
              c("NK cells","ccRCC", 
                "Exhau CD8+ T", "CD4+ T cells", 
                "TAM", "CAF", "Endo cell", 
                "Monocytes","ccRCC","Treg cells",
                "TAM","DCs",
                "ccRCC","Pro CD8+ T","Endo cell","ccRCC","Plasma cells",
                "B cells", "Mast cells",
                "Plasma cells","E-P CD8+ T"),
              c("0 NK cells","1 ccRCC 1", 
                "2 Exhau CD8+ T", "3 CD4+ T cells", 
                "4 TAM 1", "5 CAF", "6 Endo cell 1", 
                "7 Monocytes","8 ccRCC 2","9 Treg cells",
                "10 TAM 2","11 DCs",
                "12 ccRCC 4","13 Pro CD8+ T","14 Endo cell 2","15 ccRCC 3","16 Plasma cells",
                "17 B cells", "18 Mast cells",
                "19 Plasma cells","20 E-P CD8+ T"))
mRCC@meta.data$cell_type_renamed<-rename_cts[as.character(mRCC@meta.data$cell_type)]
DimPlot(mRCC, reduction = "umap", group.by="cell_type_renamed",label = TRUE, pt.size = .5)

save(mRCC, file= paste0(data_path,"mRNA19_merged_seurat_annotated.RData"))

# Save cell type annotation (all samples together)
meta_data<-mRCC@meta.data
meta_data<-meta_data[,c("orig.ident","cell_type_renamed")]
meta_data$barcode<-gsub("_.*","",rownames(meta_data))
colnames(meta_data)<-c("sample","cell_type","barcode")
write.table(meta_data,file=paste0(data_path,"cell_type_annotation.tsv"),
            sep="\t",quote=FALSE,row.names=FALSE)

# ------------------------------------------------------------------------------
# Prepare input files for RNA benchmarking - start with sample 86
# ------------------------------------------------------------------------------

output_dir<-"~/Documents/CNV_RNAseq_benchmark/data/ccRCC_2023/input_RCC86/"
  
RCC86<-mRCC[,mRCC$orig.ident=="mRCC86"]
raw_counts<-as.matrix(RCC86@assays$RNA@counts)

#Rename barcodes (remove _3)
colnames(raw_counts)<-gsub("_.*","",colnames(raw_counts))

#Get annotation file
meta_data_sample<-meta_data[meta_data$sample=="mRCC86",]
#all(meta_data_sample$barcode == colnames(raw_counts))
meta_data_sample<-meta_data_sample[,c("barcode","cell_type")]
meta_data_sample$cell_type[meta_data_sample$cell_type=="ccRCC"]<-"RCC86"

#Replace all white spaces with "_"
meta_data_sample$cell_type<-gsub(" ","_",meta_data_sample$cell_type)

#Save result data frames
write.table(raw_counts,
            file=paste0(output_dir,"count_matrix.txt"),
            sep="\t",quote=FALSE,row.names=TRUE)

write.table(meta_data_sample, 
            file=paste0(output_dir,"sample_annotation.txt"),
            sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)

#Save the reference groupfs
df_refs<-data.frame(ref_groups=setdiff(unique(meta_data_sample$cell_type),"RCC86"))
write.table(df_refs, 
            file=paste0(output_dir,"ref_groups.txt"),
            sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)

