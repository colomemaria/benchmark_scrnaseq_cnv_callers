# ------------------------------------------------------------------------------
# Run copykat (both modes possible, specifying reference cells or not)
# ------------------------------------------------------------------------------

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# ------------------------------------------------------------------------------
print("Load libraries")
# ------------------------------------------------------------------------------

library(copykat)
library(data.table)

# ------------------------------------------------------------------------------
print("Get input parameters from snakemake")
# ------------------------------------------------------------------------------

input_file<-snakemake@input$matrix

#Potentially getting reference annotations (not required)
input_annotations<-snakemake@input$annot
input_ref_groups<-snakemake@input$ref_groups

output_file<-snakemake@output$cnv_file
  
# ------------------------------------------------------------------------------
print("Execute the copyKat in the chosen settings.")
# ------------------------------------------------------------------------------

#Load dataset matrix (annotation file not required in this modus)
data_matrix<-fread(input_file)
#Format into a matrix
gene_names<-data_matrix$V1
data_matrix$V1<-NULL
data_matrix<-as.matrix(data_matrix)
rownames(data_matrix)<-gene_names

#Define reference cells in case they are specified
if(is.null(input_annotations)){
  print("No reference cells defined, copyKat will identify cancer cells automatically")
  
  ref_cells <- ""
  
} else {
  print("Reference cells for copyKat defined")
  
  #Extract reference cells
  annotation<-fread(input_annotations, header=FALSE)
  
  #Read reference groups (saved in one tsv file)
  ref_groups<-read.table(input_ref_groups,header=TRUE)
  
  ref_cells <- annotation$V1[annotation$V2 %in% ref_groups$ref_groups]
  
}

#Create the output directory and define it as a working directory
output_dir<-dirname(output_file)
dir.create(output_dir,recursive=TRUE) # gives a warning if the directory exists already
setwd(output_dir)

#Extract name of the dataset
dataset_name<-gsub("output_","",unlist(strsplit(output_dir,split="/"))[2])

#Run copyKat
copykat.test <- copykat(rawmat=data_matrix, #2d matrix with gene expression counts
                        id.type="S", #gene id type (symbol or ensemble)
                        cell.line="no", #if data is from pure cell line
                        ngene.chr=5,
                        LOW.DR = 0.05,
                        UP.DR = 0.1,
                        win.size = 25, #minimal window size for segmentation
                        norm.cell.names = ref_cells, #not specifying the reference cells
                        sam.name=dataset_name, #sample name used for output files
                        distance="euclidean",
                        output.seg="FALSE",
                        plot.genes="TRUE",
                        genome = "hg20",
                        n.cores=1)

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()
