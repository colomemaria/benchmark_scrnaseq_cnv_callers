# -----------------------------------------------------------------------------
# Extract dataset characteristics
# -----------------------------------------------------------------------------

library(data.table)
library(Seurat)

dataset_names<-c("SNU601","NCIN87","MKN45","MCF7","KATOIII","NUGC4","SNU638",
                 "HGC27","SNU16","SNU668",
                 "MM","BCC06","BCC06post","COLO320",
                 "SNU601_sample20","SNU601_sample40","SNU601_sample60",
                 "SNU601_sample80")

dataset_charac<-NULL
for(dataset in dataset_names){
  
  print(dataset)
  
  #Load dataset matrix
  data_matrix<-suppressWarnings(fread(paste0("data/input_",dataset,"/count_matrix.txt")))
  #Format into a matrix
  gene_names<-data_matrix$V1
  data_matrix$V1<-NULL
  data_matrix<-as.matrix(data_matrix)
  rownames(data_matrix)<-gene_names
  
  #Get cancer cell annotation
  annotation<-fread(paste0("data/input_",dataset,"/sample_annotation.txt"), 
                    header=FALSE)
  ref_groups<-read.table(paste0("data/input_",dataset,"/ref_groups.txt"),
                         header=TRUE)
  ref_cells <- annotation$V1[annotation$V2 %in% ref_groups$ref_groups]
  cancer_cells <- annotation$V1[! (annotation$V2 %in% ref_groups$ref_groups)]
  
  #Get coefficient of variation for cancer cells (first log-normalize with seurat)
  seurat_obj <- CreateSeuratObject(counts = data_matrix[,cancer_cells],
                                   min.cells = 3)
  seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize")
  #norm_matrix <- seurat_obj@assays$RNA@data
  norm_matrix <- seurat_obj@assays$RNA$data
  means<-rowMeans(norm_matrix)
  sds<-apply(norm_matrix,1,sd)
  coef_var<-sqrt(exp(sds^2)-1)

  
  #Save number cells and mean UMI counts per cell (separately for cancer and ref)
  dataset_charac<-rbind(dataset_charac,
            data.frame(num_cancer_cells=length(cancer_cells),
                       num_ref_cells=length(ref_cells),
                       num_genes=nrow(data_matrix),
                        mean_umi_cancer=mean(colSums(data_matrix[,cancer_cells])),
                        mean_umi_ref=mean(colSums(data_matrix[,ref_cells])),
                        mean_umi_total=mean(colSums(data_matrix)),
                        mean_drop_cell=mean(colMeans(data_matrix==0)),
                        mean_nonZeroGenes_cell = mean(colSums(data_matrix!=0)),
                        mean_drop_gene=mean(rowMeans(data_matrix==0)),
                        mean_coef_var=mean(coef_var),
                        dataset))
  
}

fwrite(dataset_charac,file="data/general_dataset_characteristics.tsv",
       quote=FALSE,sep="\t")
