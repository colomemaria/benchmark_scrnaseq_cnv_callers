# ------------------------------------------------------------------------------
# Run CONICSmat
# ------------------------------------------------------------------------------

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# ------------------------------------------------------------------------------
print("Load libraries")
# ------------------------------------------------------------------------------

library(CONICSmat)
library(data.table)

# ------------------------------------------------------------------------------
print("Get input parameters from snakemake")
# ------------------------------------------------------------------------------

input_expression_df <- snakemake@input$expression_df
input_chrom_pos_df <- snakemake@input$chrom_pos_df
input_gene_pos_df <- snakemake@input$gene_pos_df
input_annotations <- snakemake@input$annot
input_ref_groups <- snakemake@input$ref_groups

output_filtering_dims <- snakemake@output$filtering_dims
out_cnv_types <- snakemake@output$cnv_types
out_tumor_pred <- snakemake@output$tumor_pred
out_cnv_binary_mat <- snakemake@output$cnv_binary_mat
out_p_val <- snakemake@output$p_val
out_post_prob<-snakemake@output$post_prob
  
output_boxplot_path <-snakemake@output$boxplot_path
output_heatmap_expr_path <- snakemake@output$heatmap_expr_path
output_histo_path <- snakemake@output$histo_path

output_path <- dirname(out_cnv_types)
  
# ------------------------------------------------------------------------------
print("Execute CONICSmat")
# ------------------------------------------------------------------------------

#load inputs
expr_df <- as.matrix(read.table(input_expression_df,sep="\t",
                                header=T,row.names=1,check.names=F))
chrom_regions <- read.table(input_chrom_pos_df,sep="\t",row.names = 1,header = T)
gene_pos <- read.table(file = input_gene_pos_df, header = TRUE, sep = "\t")

#Normalize expr_df to counts per million
expr_df<-t(t(expr_df)/colSums(expr_df)*1e6)
#Followed by normalization to log2(CPM/10+1) (following CONICSmat tutorial)
expr_df<-log2(expr_df/10+1)
  
#replace missing values with zeroes
expr_df[which(is.na(expr_df))] <- 0 #convert NA values to 0's

#filter uninformative genes (expressed in few cells)
expr_df <- filterMatrix(expr_df,gene_pos[,"hgnc_symbol"],minCells=5)

#save the number of genes after filtering
filtered_res<-data.frame(ncells=ncol(expr_df),ngenes=nrow(expr_df))
write.table(filtered_res, file = output_filtering_dims, sep = "\t",row.names = FALSE)

#calculate a normalization factor for each cell
normFactor <- calcNormFactors(expr_df)

#Define reference cells in case they are specified
if(is.null(input_annotations)){
  print("No reference cells defined, CONICSmat will identify cancer cells automatically")
  
  #determine if avg gene expression in any region shows bimodal distrib. across cells
  l <- plotAll(expr_df,normFactor,chrom_regions,gene_pos, 
               fname = paste(output_path, "CONICSmat_CNV", sep = "//"))
  
} else {
  print("Reference cells for CONICSmat defined")
  
  #Extract reference cells
  annotation<-fread(input_annotations,header=FALSE)
  
  #Read reference groups (saved in one tsv file)
  ref_groups<-read.table(input_ref_groups,header=TRUE)
  
  ref_cells <- annotation$V1[annotation$V2 %in% ref_groups$ref_groups]
  ref_indices <- which(colnames(expr_df) %in% ref_cells)
  
  #determine if avg gene expression in any region shows bimodal distrib. across cells
  l <- plotAll(expr_df,normFactor,chrom_regions,gene_pos, 
               fname = paste(output_path, "CONICSmat_CNV", sep = "//"),
               normal = ref_indices)
  
}

#create a heatmap of posterior prob. of cells for component2 of each region
#plotHistogram doesn't return any data, so below is the calculation part of it.
#hi=plotHistogram(l,suva_expr,clusters=2,zscoreThreshold=4,patients)
t = 4
l_1 = scale(l)
if (max(l_1) > t) {
  l_1[which(l_1 > t)] = t
  l_1[which(l_1 < (-t))] = (-t)
} else {
  mx = min(max(l_1), abs(min(l)))
  sc = t/mx
  l_1 = l_1 * sc
  l_1[which(l_1 > t)] = t
  l_1[which(l_1 < (-t))] = (-t)
}

#Determine which cells are malignant
#filter out uninformative regions based on results of likelihood ratio test and BIC
lrbic=read.table(paste(output_path, "CONICSmat_CNV_BIC_LR.txt", sep = "//"),
                 sep="\t",header=T,row.names=1,check.names=F)
colnames(lrbic)
candRegions=rownames(lrbic)[which(lrbic[,"BIC difference"]>200 & 
                                    lrbic[,"LRT adj. p-val"]<0.01)]

#Check if any significant candidate regions were found 
#(otherwise second part can not be run)
if(length(candRegions)){
  
  #generate a heatmap of posterior probabilities for candidate Regions
  png(output_histo_path)
  hi=plotHistogram(l[,candRegions],expr_df,clusters=2,zscoreThreshold=4)
  dev.off()
  
  #assign labels to subsets
  normal <- which(hi==1)
  tumor <- which(hi!=1)
  
  is_tumor <- data.frame( cell = names(hi), class = ifelse(hi == 1, 0, 1))
  
  #plot posterior probabilities with statistic for normal and tumor cells
  redu <- plotAll(expr_df,normFactor,chrom_regions[candRegions,],gene_pos, 
                  fname = paste(output_path, "CONICSmat_CNV_with_info", sep = "//"),
                  normal=normal,tumor=tumor)
  
  #Generate the binary matrix, inferring CNV presence through 
  #thresholding posterior probabilities
  bin_mat <- binarizeMatrix(redu,normal,tumor,0.8)
  
  #obtain p-values for each cell and cnv
  png(output_boxplot_path)
  r <- generatePvalMat(expr_df,chrom_regions[candRegions,],
                       normFactor,normal,tumor,gene_pos,threshold=0.8)
  binr <- ifelse(r>0.1,0,1)
  boxplot(r)
  dev.off()
  
  #visualize chromosomal alterations
  png(output_heatmap_expr_path)
  plotChromosomeHeatmap(expr_df,normal = normal, plotcells = c(1:ncol(expr_df)), 
                        gene_pos = gene_pos, windowsize = 121, chr=T, expThresh=0.2, 
                        thresh = 1)
  dev.off()
  
#Create empty files so that snakemake doesn't complain
} else{
  
  is_tumor<-data.frame()
  r<-data.frame()
  binr<-data.frame()
  
  png(output_histo_path)
  plot.new()
  dev.off()
  
  png(output_boxplot_path)
  plot.new()
  dev.off()
  
  png(output_heatmap_expr_path)
  plot.new()
  dev.off()
}

#Save posterior probabilities
write.table(l,file=out_post_prob,sep = "\t",quote=FALSE)

#output relevant data
write.table(l_1, file = out_cnv_types, sep = "\t")
write.table(is_tumor, file = out_tumor_pred, sep = "\t")
write.table(binr, file = out_cnv_binary_mat, sep = "\t")
write.table(r, file = out_p_val, sep = "\t")

# ------------------------------------------------------------------------------
print("SessionInfo:")
# ------------------------------------------------------------------------------
sessionInfo()
