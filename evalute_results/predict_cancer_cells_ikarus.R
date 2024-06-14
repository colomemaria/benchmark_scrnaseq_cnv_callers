# ------------------------------------------------------------------------------
# Get cancer cell prediction from ICARUS as a ground truth 
# to validate annotations with the different CNV tools
# ------------------------------------------------------------------------------

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# ------------------------------------------------------------------------------
print("Load libraries")
# ------------------------------------------------------------------------------

library(data.table)
library(SingleCellExperiment)
library(basilisk)
library(zellkonverter)
library(reticulate)

# ------------------------------------------------------------------------------
print("Get input parameters from snakemake")
# ------------------------------------------------------------------------------

input_file <- "data/input_MM/count_matrix.txt"
conda_path <- "~/miniconda3/envs/ikarus/"

signatures_path <- "~/Documents/CNV_RNAseq_benchmark/results/test_ikarus/signatures.gmt"
model_path <- "~/Documents/CNV_RNAseq_benchmark/results/test_ikarus/core_model.joblib"
ikarus_output_path <- "~/Desktop/ikarus_output"

# ------------------------------------------------------------------------------
print("Run ikarus")
# ------------------------------------------------------------------------------

#Activate python version of conda for reticulate
use_condaenv(conda_path)

#Load dataset matrix
data_matrix<-fread(input_file)
#Format into a matrix
gene_names<-data_matrix$V1
data_matrix$V1<-NULL
data_matrix<-as.matrix(data_matrix)
rownames(data_matrix)<-gene_names

#Check whether the gene names are matching
genSig<-fread(signatures_path,header=FALSE)
tumor_genes<-unlist(genSig[genSig$V1=="Tumor",3:ncol(genSig)])
tumor_genes<-tumor_genes[tumor_genes != ""]
mean(tumor_genes %in% gene_names)

normal_genes<-unlist(genSig[genSig$V1=="Normal",3:ncol(genSig)])
normal_genes<-normal_genes[normal_genes != ""]
mean(normal_genes %in% gene_names)

#Create a single cell experiment object
sce <- SingleCellExperiment(assays = list(counts = data_matrix),
                            colData = data.frame(cell=colnames(data_matrix)))

#Covert sce object to an anndata object (using zellkonverter)
adata <- basiliskRun(fun = function(sce) {
  SCE2AnnData(sce)}, env = conda_path, sce = sce)
rm(sce)

#Load the trained ikarus model (problem with load_core_model function and paths?)
setwd(dirname(model_path))
ikarus <- import("ikarus")
model  <- ikarus$classifier$Ikarus(signatures_gmt = basename(signatures_path), 
                                   out_dir=ikarus_output_path)
model$load_core_model(basename(model_path))

#Preprocess adata object
adata <- ikarus$data$preprocess_adata(adata)

#Predict cancer cells
class_predict <- model$predict(adata, "test")


# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()