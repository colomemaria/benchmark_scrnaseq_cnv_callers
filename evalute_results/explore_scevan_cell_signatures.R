# ------------------------------------------------------------------------------
# Explore how the confidential normal is working for SCEVAN
# => because it doesn't seem to work correctly for the iAMP21 dataset
# Test it for different datasets
# ------------------------------------------------------------------------------

library(SCEVAN)
library(data.table)

#Gene set:
#names(SCEVAN:::geneSet)

for(dataset in c("MM","BCC06","BCC06post")){ #iAMP21
  
  print("Getting confidentally normal genes for dataset:")
  print(dataset)
  
  #Load iAMP21 data
  input_file<-paste0("data/input_",dataset,"/count_matrix.txt")
  
  data_matrix<-fread(input_file)
  #Format into a matrix
  gene_names<-data_matrix$V1
  data_matrix$V1<-NULL
  data_matrix<-as.matrix(data_matrix)
  rownames(data_matrix)<-gene_names
  
  count_mtx=data_matrix
  sample=""
  
  #Following pipeline of SCEVAN::pipeline step by step
  
  #dir.create(file.path("./output"), showWarnings = FALSE)
  setwd(paste0("results/output_",dataset,"/scevan"))
  
  res_proc <- SCEVAN:::preprocessingMtx(count_mtx, sample, par_cores = 1,
                               findConfident = TRUE,
                               AdditionalGeneSets = NULL,
                               SCEVANsignatures = TRUE, organism = "human")
  saveRDS(res_proc, 
          file="output/first_healthy_cell_annotation.RDS")

  setwd("../../..")
  
}

# #Inside SCEVAN:::preprocessingMtx
# perc_genes = 0.1
# organism="human"
# 
# print(paste(" raw data - genes: ", nrow(count_mtx), " cells: ", 
#             ncol(count_mtx), sep = ""))
# print("1) Filter: cells > 200 genes")
# genes.raw <- apply(count_mtx, 2, function(x) (sum(x > 0)))
# if (sum(genes.raw > 200) == 0) 
#   stop("none cells have more than 200 genes")
# if (sum(genes.raw < 100) > 1) {
#   count_mtx <- count_mtx[, -which(genes.raw < 200)]
#   print(paste("filtered out ", sum(genes.raw <= 200), " cells past filtering ", 
#               ncol(count_mtx), " cells", sep = ""))
# }
# der <- apply(count_mtx, 1, function(x) (sum(x > 0)))/ncol(count_mtx)
# if (sum(der > perc_genes) > 7000) {
#   print(paste0("2) Filter: genes > ", perc_genes * 100, 
#                "% of cells"))
#   count_mtx <- count_mtx[which(der > perc_genes), ]
# } else {
#   perc_genes <- perc_genes - 0.05
#   print("low data quality")
#   print(paste0("2) Filter: genes > ", perc_genes * 100, 
#                "% of cells"))
#   count_mtx <- count_mtx[which(der > perc_genes), ]
# }
# print(paste(nrow(count_mtx), " genes past filtering", sep = ""))
# print("3) Annotations gene coordinates")
# count_mtx_annot <- annotateGenes(count_mtx, organism)
# count_mtx <- count_mtx_annot[, -c(1:5)]
# rownames(count_mtx) <- count_mtx_annot$gene_name
# 
# AdditionalGeneSets=NULL
# SCEVANsignatures=TRUE
# par_cores=1
# norm_cell <- getConfidentNormalCells(count_mtx, sample, 
#                par_cores = par_cores, AdditionalGeneSets = AdditionalGeneSets, 
#                SCEVANsignatures = SCEVANsignatures, organism = organism)
# 
# #Remark: this is based on the yaGST mwwGST package to run a competitive test to highlight
# #whether a gene set is highly ranked in a sequence of gene values on the genes outside the gene-set
# SCEVAN:::ssMwwGst(geData = count_mtx, geneSet = SCEVAN:::geneSet, ncore = par_cores, 
#          minLenGeneSet = 5, filename = paste0("./output/", sample, 
#                                               "_confidentNormal"), standardize = TRUE)
# 
# load(paste0("./output/", sample, "_confidentNormal_MWW.RData")) #FDR, NES, pValue

#SCEVAN:::getConfidentNormalCells
#SCEVAN:::ssMwwGst

