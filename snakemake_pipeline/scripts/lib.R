# ------------------------------------------------------------------------------
# Collection of small help scripts:
# - Read scWGS results and output of the different CNV methods
# ------------------------------------------------------------------------------

#Read WGS results (need to be preprocessed to pseudobulk)
read_wgs_results<-function(path,add_chr_symbol=TRUE){
  
  wgs<-fread(path)
  wgs$V1<-NULL
  
  #Convert chromosome names to match EpiAneufinder and inferCNV results
  if(add_chr_symbol){
    wgs$chr<-paste0("chr",wgs$chr)
  }
  
  #Combine loss, base and gain to one number
  wgs$wgs_mean<-with(wgs,loss_wgs+(base_wgs*2)+gain_wgs*3)
  
  #Create Grange results
  wgs_grange<-makeGRangesFromDataFrame(wgs, keep.extra.columns=FALSE)
  wgs_grange$wgs_mean<-wgs$wgs_mean
  return(wgs_grange)
}

#Read WES results (using per exon log fold change from CNVkit)
read_wes_results_lfc<-function(wes_path,gene_annot){
  
  #Create a gene grange as ground truth
  gene_annot<-fread(gene_annot)
  colnames(gene_annot)<-c("gene","chr","start","end")
  gene_grange<-makeGRangesFromDataFrame(gene_annot, keep.extra.columns=TRUE)
  
  #Mapp the exons to the genes
  wes<-fread(wes_path)
  #Remove antitarget regions
  wes<-wes[wes$gene!="Antitarget",]
  wes_grange<-makeGRangesFromDataFrame(wes,keep.extra.columns = FALSE)
  
  #Map exons to genes
  overlaps<-as.data.frame(findOverlaps(gene_grange,wes_grange,type="any"))
  
  print(paste("From",length(wes_grange),"exons (wo Antitargets) mapping",
              length(unique(overlaps$subjectHits)),"to genes."))
  
  #Filter first gene range (keep dimensions)
  gene_grange<-gene_grange[unique(overlaps$queryHits)]
  
  #Calculate exon length to average the score
  wes$length<-wes$end-wes$start
  
  # #Average log fold change
  # gene_grange$wes_lfc <- sapply(unique(overlaps$queryHits),
  #                                        function(i) mean(wes$log2[overlaps$subjectHits[overlaps$queryHits==i]]))   
  
  #Average log fold change (weighted by exon length)
  gene_grange$wes_lfc <- sapply(unique(overlaps$queryHits),
              function(i) sum(wes$log2[overlaps$subjectHits[overlaps$queryHits==i]]*
                                wes$length[overlaps$subjectHits[overlaps$queryHits==i]])/
                sum(wes$length[overlaps$subjectHits[overlaps$queryHits==i]]))               
  
  #Shift it again by 2 to set the normal level at 2
  gene_grange$wes_lfc<-gene_grange$wes_lfc+2
  return(gene_grange)
}


#Read WGS/WES results processed with CNVkit (long segments)
read_results_cnvkit_segs<-function(wes_path){
  wes<-fread(wes_path)
  
  wes$cn_binarized<-ifelse(wes$cn<1,1,
                           ifelse(wes$cn>3,3,wes$cn))
  
  wes_gr<-makeGRangesFromDataFrame(wes)
  wes_gr$call<-wes$cn_binarized
  
  #Determine the annotation region of each chromsome
  wes_limits<-wes%>%group_by(chromosome)%>%
    summarize(global_start=min(start),global_end=max(end))
  
  #Round it to 100,000 borders (to have bins starting at 0)
  wes_limits$global_start_rounded<-floor(wes_limits$global_start/100000)*100000
  wes_limits$global_end_rounded<-ceiling(wes_limits$global_end/100000)*100000-1
  
  #Use Genomic Ranges to create 100,000 kB windows
  wes_limits_gr<-makeGRangesFromDataFrame(wes_limits,
                                          start.field="global_start_rounded",
                                          end.field="global_end_rounded")
  wes_limits_gr<-wes_limits_gr[order(wes_limits_gr)]
  
  wes_bins<-tile(wes_limits_gr,width=100000)
  wes_bins<-unlist(wes_bins)
  
  
  #Map it back to the genomic segments
  overlaps<-as.data.frame(findOverlaps(wes_bins,wes_gr,type="any"))
  wes_bins<-wes_bins[unique(overlaps$queryHits)]
  wes_bins$wgs_mean <- sapply(unique(overlaps$queryHits),
                             function(i) mean(wes_gr$call[overlaps$subjectHits[overlaps$queryHits==i]]))
  
  return(wes_bins)
}


#Read WES results from GATK (file: <filename>_clean.called.seg)
read_wes_results_gatk<-function(wes_path){
  
  #Start reading at the line starting with CONTIG (everything before is header)
  wes<-fread(wes_path,skip="CONTIG")
  
  wes$CALL_num<-ifelse(wes$CALL=="-",1,
                       ifelse(wes$CALL=="+",3,2))
  
  wes_gr<-makeGRangesFromDataFrame(wes,seqnames.field="CONTIG",
                                   start.field="START",end.field="END")
  wes_gr$call<-wes$CALL_num
  
  #Determine the annotation region of each chromsome
  wes_limits<-wes%>%group_by(CONTIG)%>%
    summarize(global_start=min(START),global_end=max(END))
  
  #Round it to 100,000 borders (to have bins starting at 0)
  wes_limits$global_start_rounded<-floor(wes_limits$global_start/100000)*100000
  wes_limits$global_end_rounded<-ceiling(wes_limits$global_end/100000)*100000-1
  
  #Use Genomic Ranges to create 100,000 kB windows
  wes_limits_gr<-makeGRangesFromDataFrame(wes_limits,
                                          seqnames.field = "CONTIG",
                                          start.field="global_start_rounded",
                                          end.field="global_end_rounded")
  wes_limits_gr<-wes_limits_gr[order(wes_limits_gr)]
  
  wes_bins<-tile(wes_limits_gr,width=100000)
  wes_bins<-unlist(wes_bins)
  
  
  #Map it back to the genomic segments
  overlaps<-as.data.frame(findOverlaps(wes_bins,wes_gr,type="any"))
  wes_bins<-wes_bins[unique(overlaps$queryHits)]
  wes_bins$wgs_mean <- sapply(unique(overlaps$queryHits),
                             function(i) mean(wes_gr$call[overlaps$subjectHits[overlaps$queryHits==i]]))
  
  return(wes_bins)
}

#Read results from inferCNV (6 state HMM model)
read_infercnv_6state_model<-function(matrix_path,gene_path,clones=NULL){
  
  suppressWarnings(res_inferCNV<-fread(matrix_path))
  
  #Transform into a numeric matrix
  gene_names<-res_inferCNV$V1
  res_inferCNV$V1<-NULL
  res_inferCNV<-as.matrix(res_inferCNV)
  rownames(res_inferCNV)<-gene_names
  
  #Combine the different loss states and gain states
  res_inferCNV[res_inferCNV==0.5]<-0
  res_inferCNV[res_inferCNV == 1.5 | res_inferCNV == 2 | res_inferCNV == 3]<-2
  
  #Create a grange object from the gene position information
  gene_pos <- fread(gene_path)
  colnames(gene_pos)<-c("gene","chr","start","end")
  gene_pos<-as.data.frame(gene_pos)
  rownames(gene_pos)<-gene_pos$gene
  gene_pos<-gene_pos[rownames(res_inferCNV),]
  
  gene_pos_range<-makeGRangesFromDataFrame(gene_pos,keep.extra.columns=TRUE)
  
  #Calculate the pseudobulk combined over all cells unless the clonal structure is given
  #(then pseudobulk per subclone)
  if(is.null(clones)){
    gene_pos_range$infercnv_cnv<-rowMeans(res_inferCNV)+1
  } else {
    rownames(clones)<-clones$cell
    clones<-clones[colnames(res_inferCNV),]
    clones$infercnv<-as.factor(clones$infercnv)
    pbulk_clones<-t(apply(res_inferCNV,1,tapply,clones$infercnv,mean))+1
    mcols(gene_pos_range)<-cbind(mcols(gene_pos_range),pbulk_clones)
  }

  return(list(gene_pos_range,ncol(res_inferCNV)))
}

#Read results from inferCNV (6 state HMM model) - per cell
read_infercnv_6state_model_individual_cells<-function(matrix_path,gene_path){
  
  suppressWarnings(res_inferCNV<-fread(matrix_path))
  
  #Transform into a numeric matrix
  gene_names<-res_inferCNV$V1
  res_inferCNV$V1<-NULL
  res_inferCNV<-as.matrix(res_inferCNV)
  rownames(res_inferCNV)<-gene_names
  
  #Combine the different loss states and gain states
  res_inferCNV[res_inferCNV==0.5]<-0
  res_inferCNV[res_inferCNV == 1.5 | res_inferCNV == 2 | res_inferCNV == 3]<-2
  
  #Shift by one to have baseline at 0
  res_inferCNV<-res_inferCNV+1
  
  #Create a grange object from the gene position information
  gene_pos <- fread(gene_path)
  colnames(gene_pos)<-c("gene","chr","start","end")
  gene_pos<-as.data.frame(gene_pos)
  rownames(gene_pos)<-gene_pos$gene
  gene_pos<-gene_pos[rownames(res_inferCNV),]
  
  gene_pos_range<-makeGRangesFromDataFrame(gene_pos,keep.extra.columns=TRUE)
  
  return(list(gene_pos_range,res_inferCNV))
}

#Read results from inferCNV (normalized expression)
read_infercnv_expr<-function(matrix_path,gene_path){
  
  suppressWarnings(res_inferCNV<-fread(matrix_path))
  
  #Transform into a numeric matrix
  gene_names<-res_inferCNV$V1
  res_inferCNV$V1<-NULL
  res_inferCNV<-as.matrix(res_inferCNV)
  rownames(res_inferCNV)<-gene_names
  
  #Get average results across all cells
  pseudobulk_res_inferCNV<-data.frame(gene_names,
                                      mean_expr=rowMeans(res_inferCNV))
  
  #Read gene position information
  gene_pos <- fread(gene_path)
  colnames(gene_pos)<-c("gene","chr","start","end")
  
  #Filter for genes that were not filtered out
  pseudobulk_res_inferCNV<-merge(pseudobulk_res_inferCNV,gene_pos,
                                 by.x="gene_names",by.y="gene",sort=FALSE)
  
  gene_pos_range<-makeGRangesFromDataFrame(pseudobulk_res_inferCNV, keep.extra.columns=TRUE)
  
  
  #Combine values to a signal score (shift so that baseline is at 2)
  gene_pos_range$infercnv_expr<-gene_pos_range$mean_expr+1
  
  
  return(list(gene_pos_range,ncol(res_inferCNV)))
}

#Read results from inferCNV (normalized expression) - per cell
read_infercnv_expr_individual_cells<-function(matrix_path,gene_path){
  
  suppressWarnings(res_inferCNV<-fread(matrix_path))
  
  #Transform into a numeric matrix
  gene_names<-res_inferCNV$V1
  res_inferCNV$V1<-NULL
  res_inferCNV<-as.matrix(res_inferCNV)
  rownames(res_inferCNV)<-gene_names
  
  #Shift it so that the baseline is 2
  res_inferCNV<-res_inferCNV+1
  
  #Read gene position information
  gene_pos <- fread(gene_path)
  colnames(gene_pos)<-c("gene","chr","start","end")
  gene_pos<-gene_pos[gene_pos$gene %in% rownames(res_inferCNV),]
  
  gene_pos_range<-makeGRangesFromDataFrame(gene_pos)
  
  return(list(gene_pos_range,res_inferCNV))
}

#Read results from casper (normalized expression)
read_casper<-function(input_file){
  
  casper_complete<-readRDS(input_file)
  casper_score_merged<-with(casper_complete,mean_loss+(mean_base*2)+mean_gain*3)
  casper_complete$casper<-casper_score_merged
    
  return(casper_complete)
}

read_casper_individual_cells<-function(input_casper_grange,input_casper_cells){
  
  #Get position information
  casper_pseudobulk<-readRDS(input_casper_grange)
  elementMetadata(casper_pseudobulk)<-NULL
  
  #Convert cells to a numeric matrix
  suppressWarnings(casper_cells<-fread(input_casper_cells))
  genes<-casper_cells$V1
  casper_cells$V1<-NULL
  casper_cells<-as.matrix(casper_cells)
  rownames(casper_cells)<-genes
  
  #Convert them to have a mean of 2
  casper_cells<-casper_cells+2
  
  return(list(casper_pseudobulk,casper_cells))
}

#Read results from copyKat (normalized expression)
read_copykat<-function(matrix_path,annot_path,ref_path,
                       clones=NULL){
  
  suppressWarnings(copykat_res<-fread(matrix_path))
  
  #Generate a grange object with pseudobulk results
  copykat_res$chromosome_name<-paste0("chr",copykat_res$chromosome_name)
  
  #Replace chr23 by chrX (to align with other methods)
  copykat_res$chromosome_name[copykat_res$chromosome_name=="chr23"]<-"chrX"
    
  copykat_grange<-makeGRangesFromDataFrame(copykat_res,
                                           seqnames.field="chromosome_name",
                                           start.field="start_position",
                                           end.field="end_position")
  
  #Get cell matrix
  copykat_cells<-as.matrix(copykat_res[,8:ncol(copykat_res)])
  
  #Remove reference cells
  annot<-fread(annot_path,header=FALSE)
  ref_groups<-fread(ref_path)
  cancer_cells<-annot$V1[!(annot$V2 %in% ref_groups$ref_groups)]
  copykat_cells<-copykat_cells[,colnames(copykat_cells) %in% cancer_cells]
  
  #Calculate the pseudobulk combined over all cells unless the clonal structure is given
  #(then pseudobulk per subclone)
  if(is.null(clones)){
    copykat_grange$copykat<-rowMeans(copykat_cells)+2
  } else {
    clones<-clones[colnames(copykat_cells),]
    clones$copykat<-as.factor(clones$copykat)
    pbulk_clones<-t(apply(copykat_cells,1,tapply,clones$copykat,mean))+2
    colnames(pbulk_clones)<-paste0("clone_",colnames(pbulk_clones))
    mcols(copykat_grange)<-cbind(mcols(copykat_grange),pbulk_clones)
  }
  
  return(list(copykat_grange,ncol(copykat_cells)))
}

#Read copykat results without aggregating it to pseudobulk
read_copykat_individual_cells<-function(matrix_path,annot_path,ref_path){
  
  suppressWarnings(copykat_res<-fread(matrix_path))
  
  #Generate a grange object with pseudobulk results
  copykat_res$chromosome_name<-paste0("chr",copykat_res$chromosome_name)
  
  #Replace chr23 by chrX (to align with other methods)
  copykat_res$chromosome_name[copykat_res$chromosome_name=="chr23"]<-"chrX"
  
  copykat_grange<-makeGRangesFromDataFrame(copykat_res,
                                           seqnames.field="chromosome_name",
                                           start.field="start_position",
                                           end.field="end_position")
  
  #Get cell matrix
  copykat_cells<-as.matrix(copykat_res[,8:ncol(copykat_res)])
  
  #Move mean to 2
  copykat_cells<-copykat_cells+2
  
  #Remove reference cells
  annot<-fread(annot_path,header=FALSE)
  ref_groups<-fread(ref_path)
  cancer_cells<-annot$V1[!(annot$V2 %in% ref_groups$ref_groups)]
  copykat_cells<-copykat_cells[,colnames(copykat_cells) %in% cancer_cells]
  
  return(list(copykat_grange,copykat_cells))
}

#Read results from SCEVAN (normalized expression)
read_scevan<-function(matrix_path,gene_annot,annot_path,ref_path,
                      clones=NULL){

  #Load normalized count matrix and gene positions
  load(matrix_path) #results.com (when looking at the subclonal matrix)
  CNA_mtx_relat<-results.com
  load(gene_annot) #count_mtx_annot
  
  #Generate a grange object with pseudobulk results
  count_mtx_annot$seqnames<-paste0("chr",count_mtx_annot$seqnames)
  scevan_grange<-makeGRangesFromDataFrame(count_mtx_annot,
                                           seqnames.field="seqnames",
                                           start.field="start",
                                           end.field="end")
  
  #Remove reference cells
  annot<-fread(annot_path,header=FALSE)
  ref_groups<-fread(ref_path)
  cancer_cells<-annot$V1[!(annot$V2 %in% ref_groups$ref_groups)]
  CNA_mtx_relat<-CNA_mtx_relat[,colnames(CNA_mtx_relat) %in% cancer_cells]
  
  #Calculate the pseudobulk combined over all cells unless the clonal structure is given
  #(then pseudobulk per subclone)
  if(is.null(clones)){
    #Calculate mean over all cancer cells (and shift baseline to 2)
    scevan_grange$scevan<-rowMeans(CNA_mtx_relat)+2
  } else {
    
    #Filter cells
    common_cells<-intersect(colnames(CNA_mtx_relat),clones$cell)
    CNA_mtx_relat<-CNA_mtx_relat[,common_cells]
    clones<-as.data.frame(clones)
    rownames(clones)<-clones$cell
    clones<-clones[common_cells,]
    
    clones$scevan<-as.factor(clones$scevan)
    pbulk_clones<-t(apply(CNA_mtx_relat,1,tapply,clones$scevan,mean))+2
    colnames(pbulk_clones)<-paste0("clone_",colnames(pbulk_clones))
    mcols(scevan_grange)<-cbind(mcols(scevan_grange),pbulk_clones)
  }
  
  return(list(scevan_grange,ncol(CNA_mtx_relat)))
}


#Read results from SCEVAN per cell (normalized expression)
read_scevan_individual_cells<-function(matrix_path,gene_annot,annot_path,ref_path,
                      clones=NULL){
  
  #Load normalized count matrix and gene positions
  load(matrix_path) #results.com (when looking at the subclonal matrix)
  CNA_mtx_relat<-results.com
  load(gene_annot) #count_mtx_annot
  
  #Generate a grange object with pseudobulk results
  count_mtx_annot$seqnames<-paste0("chr",count_mtx_annot$seqnames)
  scevan_grange<-makeGRangesFromDataFrame(count_mtx_annot,
                                          seqnames.field="seqnames",
                                          start.field="start",
                                          end.field="end")
  
  #Remove reference cells
  annot<-fread(annot_path,header=FALSE)
  ref_groups<-fread(ref_path)
  cancer_cells<-annot$V1[!(annot$V2 %in% ref_groups$ref_groups)]
  CNA_mtx_relat<-CNA_mtx_relat[,colnames(CNA_mtx_relat) %in% cancer_cells]
  
  #Move so that the mean is 2
  CNA_mtx_relat<-CNA_mtx_relat+2
  
  return(list(scevan_grange,CNA_mtx_relat))
}


#Read results from SCEVAN (CNV predictions (VEGA output))
read_scevan_vega<-function(file_path){
  
  #Read input file
  suppressWarnings(scevan_res<-fread(file_path))
                           
  #Create grange object
  scevan_res$Chr<-paste0("chr",scevan_res$Chr)
  scevan_grange<-makeGRangesFromDataFrame(scevan_res,keep.extra.columns = FALSE)
  
  scevan_grange$scevan_vega<-scevan_res$Mean+2
  return(scevan_grange)
}

#Read SCEVAN CN status
read_scevan_cn_status<-function (clone1_path, clone_annot_path,binned_genome=NULL){
  
  require(data.table)
  require(GenomicRanges)
  
  #Load how many cells there are per subclone (to combine them to a weighted sum)
  suppressWarnings(clone_annot<-fread(clone_annot_path))
  
  #Check if significant differently subclones were found 
  #(otherwise the column named subclone is missing)
  if("subclone" %in% colnames(clone_annot)){
    
    all_files<-list.files(dirname(clone1_path),
                          pattern="*_subclone[1-9]_CN.seg",full.names=TRUE)
    
    #Read the first subclone
    suppressWarnings(clone<-fread(clone1_path))
    clone$Chr<-paste0("chr",clone$Chr)
    clone_gr<-makeGRangesFromDataFrame(clone,
                                       seqnames.field="Chr",start.field="Pos",
                                       end.field="End")
    clone_gr$CNV_subclone1<-ifelse(clone$CN<2,1,ifelse(clone$CN>2,3,2))
    
    #Separate the first clone into 100kb bins to faciliate combination with other clones
    if(is.null(binned_genome)){
      binned_genome<-create_binned_genome(clone_gr)
    }
    combined_clones<-combine_range_objects(binned_genome,clone_gr,
                                           "CNV_subclone1")
    
    #Read all other clones and combine add it to the Grange object of the first clone
    for(fl in setdiff(all_files,clone1_path)){
      clonename<-paste0("CNV_subclone",unlist(strsplit(gsub(".*subclone","",fl),split="_"))[[1]])
      suppressWarnings(clone<-fread(fl))
      clone$Chr<-paste0("chr",clone$Chr)
      clone_gr<-makeGRangesFromDataFrame(clone,
                                         seqnames.field="Chr",start.field="Pos",
                                         end.field="End")
      mcols(clone_gr)[,clonename]<-ifelse(clone$CN<2,1,ifelse(clone$CN>2,3,2))
      combined_clones<-combine_range_objects(combined_clones,clone_gr,
                                             clonename)
    }
    
  
    clone_freq<-as.data.frame(table(clone_annot$subclone))
    clone_freq$perc<-clone_freq$Freq/sum(clone_freq$Freq)
    
    cnv_matrix<-as.matrix(elementMetadata(combined_clones))
    merged_cnv_vector<-colSums(t(cnv_matrix) * clone_freq$perc)
    elementMetadata(combined_clones)<-NULL
    combined_clones$scevan_cnv<-merged_cnv_vector
    
  } else {
    
    filename<-list.files(dirname(clone1_path),
               pattern="*_Clonal_CN.seg",full.names=TRUE)
    
    #Read the general clonal profile
    suppressWarnings(clone<-fread(filename[1]))
    clone$Chr<-paste0("chr",clone$Chr)
    clone_gr<-makeGRangesFromDataFrame(clone,
                                       seqnames.field="Chr",start.field="Pos",
                                       end.field="End")
    clone_gr$scevan_cnv<-ifelse(clone$CN<2,1,ifelse(clone$CN>2,3,2))
    
    #Separate the first clone into 100kb bins to faciliate combination with other clones
    if(is.null(binned_genome)){
      binned_genome<-create_binned_genome(clone_gr)
    }
    combined_clones<-combine_range_objects(binned_genome,clone_gr,
                                           "scevan_cnv")
  }
  
  return(combined_clones)
  
}

#Read SCEVAN CN status - get results per cell
read_scevan_cn_status_individual_cells<-function (clone1_path,clone_annot_path,binned_genome=NULL){
  
  require(data.table)
  require(GenomicRanges)
  
  all_files<-list.files(dirname(clone1_path),
                        pattern="*_subclone[1-9]_CN.seg",full.names=TRUE)
  
  #Read the first subclone
  suppressWarnings(clone<-fread(clone1_path))
  clone$Chr<-paste0("chr",clone$Chr)
  clone_gr<-makeGRangesFromDataFrame(clone,
                                     seqnames.field="Chr",start.field="Pos",
                                     end.field="End")
  clone_gr$CNV_subclone1<-ifelse(clone$CN<2,1,ifelse(clone$CN>2,3,2))
  
  #Separate the first clone into 100kb bins to faciliate combination with other clones
  if(is.null(binned_genome)){
    binned_genome<-create_binned_genome(clone_gr)
  }
  combined_clones<-combine_range_objects(binned_genome,clone_gr,
                                         "CNV_subclone1")
  
  #Read all other clones and combine add it to the Grange object of the first clone
  for(fl in setdiff(all_files,clone1_path)){
    clonename<-paste0("CNV_subclone",unlist(strsplit(gsub(".*subclone","",fl),split="_"))[[1]])
    suppressWarnings(clone<-fread(fl))
    clone$Chr<-paste0("chr",clone$Chr)
    clone_gr<-makeGRangesFromDataFrame(clone,
                                       seqnames.field="Chr",start.field="Pos",
                                       end.field="End")
    mcols(clone_gr)[,clonename]<-ifelse(clone$CN<2,1,ifelse(clone$CN>2,3,2))
    combined_clones<-combine_range_objects(combined_clones,clone_gr,
                                           clonename)
  }
  
  #Load how many cells there are per subclone (to combine them to a weighted sum)
  suppressWarnings(clone_annot<-fread(clone_annot_path))
  clone_freq<-as.data.frame(table(clone_annot$subclone))
  
  scevan_cells<-NULL
  for(clone_number in as.numeric(clone_freq$Var1)){
    vector<-elementMetadata(combined_clones)[,paste0("CNV_subclone",clone_number)]
    ncells<-clone_freq$Freq[clone_number]
    clone_matrix<-matrix(rep(vector,ncells),ncol=ncells)
    scevan_cells<-cbind(scevan_cells,clone_matrix)
  }
  
  elementMetadata(combined_clones)<-NULL
  return(list(combined_clones,scevan_cells))
}

#Read CNV predictions for each genome
read_scevan_subclones<-function(clone1_path,binned_genome=NULL){
  all_files<-list.files(dirname(clone1_path),
                        pattern="*_subclone[1-9]_CN.seg",full.names=TRUE)
  
  #Read the first subclone
  suppressWarnings(clone<-fread(clone1_path))
  clone$Chr<-paste0("chr",clone$Chr)
  clone_gr<-makeGRangesFromDataFrame(clone,
                                     seqnames.field="Chr",start.field="Pos",
                                     end.field="End")
  clone_gr$CNV_subclone1<-ifelse(clone$CN<2,1,ifelse(clone$CN>2,3,2))
  
  if(is.null(binned_genome)){
    binned_genome<-create_binned_genome(clone_gr)
  }
  
  combined_clones<-combine_range_objects(binned_genome,clone_gr,
                                         "CNV_subclone1")
  
  for(fl in setdiff(all_files,clone1_path)){
    clonename<-paste0("CNV_",unlist(strsplit(basename(fl),split="_"))[2])
    suppressWarnings(clone<-fread(fl))
    clone$Chr<-paste0("chr",clone$Chr)
    clone_gr<-makeGRangesFromDataFrame(clone,
                                       seqnames.field="Chr",start.field="Pos",
                                       end.field="End")
    mcols(clone_gr)[,clonename]<-ifelse(clone$CN<2,1,ifelse(clone$CN>2,3,2))
    combined_clones<-combine_range_objects(combined_clones,clone_gr,
                                           clonename)
  }
  
  return(combined_clones)
}

#Read results from HoneyBadger
read_honeybadger_expr<-function(file_path){
  
  #Extract normalized counts from HB object
  hb<-readRDS(file_path)
  norm_counts<-hb$gexp.norm
  
  #Get normalized expression
  hb_range<-hb$genes
  hb_range$honeybadger<-rowMeans(norm_counts)+2
  
  return(list(hb_range,ncol(norm_counts)))
}

read_honeybadger_cnv<-function(file_path){
  
  #Get cnv predictions
  hb<-readRDS(file_path)
  cnvs<-hb$cnvs$`gene-based`
  
  # Assign all amplicifications with 3 and deletions with 1
  amp<-cnvs$amp
  amp$hb_cnv<-3
  
  del<-cnvs$del
  if(length(del)>0){
    del$hb_cnv<-1
  }
  
  cnvs_combined<-c(amp,del)
  
  # #Posterior probability
  # post_prob<-hb$summary$`gene-based`
  
  return(cnvs_combined)
}

#Important remark: similar to "combine_range_objects" but as we have only
#amplifications (3) and deletions (1) predictions and no information about
#normal(2) assuming that everything else is normal that is missing in range 2
add_honeybadger_to_range<-function(range_1,range_2,method_colname){
  
  if(! method_colname %in% colnames(elementMetadata(range_2))){
    stop(paste("Specified column name",method_colname,"not found in range 2!"))
  }
  
  #Check overlaps (to get the bin sizes to a similar shape)
  overlaps<-as.data.frame(findOverlaps(range_1,range_2,type="any"))
  
  #Filter first range (keep dimensions)
  elementMetadata(range_1)[method_colname] <- 2
  elementMetadata(range_1)[overlaps$queryHits,method_colname]<-
    elementMetadata(range_2)[overlaps$subjectHits,method_colname]
  
  #Return range 1 with the new additional column from range2                                                   
  return(range_1)
}

read_numbat_expr<-function(file_path){
  
  #Extract gene positions
  gene_position<-as.data.frame(numbat::gtf_hg38)
  rownames(gene_position)<-gene_position$gene
  
  #Load expression matrix (broken in current numbat version (1.2.1), 
  #so currently loaded manually)
  #expr_matrix<-numbat_obj$gexp_roll_wide
  expr_matrix<-fread(paste0(file_path,"/gexp_roll_wide.tsv.gz"))
  cells<-expr_matrix$cell
  expr_matrix$cell<-NULL
  expr_matrix<-as.matrix(expr_matrix)
  rownames(expr_matrix)<-cells
  
  #Filter for genes that are part of the expression matrix
  gene_position<-gene_position[colnames(expr_matrix),]
  
  #Match chromosome naming
  gene_position$CHROM<-paste0("chr",gene_position$CHROM)
  
  #Create a Grange object from genomic positions
  numbat_grange<-makeGRangesFromDataFrame(gene_position,
                                              start.field="gene_start",
                                              end.field="gene_end",
                                              seqnames.field="CHROM")
  
  #calculate mean expression levels and write down cell count
  numbat_grange$numbat<-colMeans(expr_matrix)+2
  numbat_expr<-list()
  numbat_expr[[1]] <- numbat_grange
  numbat_expr[[2]] <- nrow(expr_matrix)
  
  return(numbat_expr)
}

read_numbat_expr_individual_cells<-function(file_path){
  
  #Extract gene positions
  gene_position<-as.data.frame(numbat::gtf_hg38)
  rownames(gene_position)<-gene_position$gene
  
  #Load expression matrix (broken in current numbat version (1.2.1), 
  #so currently loaded manually)
  #expr_matrix<-numbat_obj$gexp_roll_wide
  expr_matrix<-fread(paste0(file_path,"/gexp_roll_wide.tsv.gz"))
  cells<-expr_matrix$cell
  expr_matrix$cell<-NULL
  expr_matrix<-as.matrix(expr_matrix)
  rownames(expr_matrix)<-cells
  
  #Move mean to 2
  expr_matrix<-expr_matrix+2
  
  #Filter for genes that are part of the expression matrix
  gene_position<-gene_position[colnames(expr_matrix),]
  
  #Match chromosome naming
  gene_position$CHROM<-paste0("chr",gene_position$CHROM)
  
  #Create a Grange object from genomic positions
  numbat_grange<-makeGRangesFromDataFrame(gene_position,
                                          start.field="gene_start",
                                          end.field="gene_end",
                                          seqnames.field="CHROM")
  
  return(list(numbat_grange,t(expr_matrix)))
}

read_numbat_cnv<-function(file_path,clones=FALSE){

  #Check if numbat found any CNVs (could be empty in diploid tests)
  if(! file.exists(paste(file_path,"bulk_clones_final.tsv.gz",sep="/"))){
    print("No CNVS found for Numbat!")
   return(NULL) 
  }

  #create a Numbat object
  numbat_obj<-Numbat$new(file_path)
  
  #remove SNPs without genes assigned
  numbat_pseudobulk<-numbat_obj$bulk_clones[numbat_obj$bulk_clones$gene!='',]
  
  #convert cnv states to index form
  numbat_pseudobulk$cnv_state_int <- -1
  #filter by deletion cnv state
  numbat_pseudobulk[numbat_pseudobulk$cnv_state == "del",]$cnv_state_int <- 1
  #filter by neutral cnv state (and loh since they are not implemented in all methods)
  numbat_pseudobulk[numbat_pseudobulk$cnv_state %in% c("loh", "neu"),]$cnv_state_int <- 2
  #filter by amplification cnv states
  numbat_pseudobulk[numbat_pseudobulk$cnv_state %in% c("amp", "amp|bamp", "bamp"),]$cnv_state_int <- 3
  
  #Create one pseudobulk over all cells unless specified differently
  if(! clones){
    #get weighted indices for multiple occurrences of one gene
    numbat_cnv<-data.frame()
    
    gene_list<-unique(numbat_pseudobulk$gene)
    for (i in 1:length(gene_list)){
      gene_df<-numbat_pseudobulk[numbat_pseudobulk$gene==gene_list[i],]
      
      numbat_cnv[i,1]<-gene_list[i]
      numbat_cnv[i,2]<-sum(gene_df$cnv_state_int*gene_df$n_cells)/sum(gene_df$n_cells)
    }
    
    colnames(numbat_cnv)<-c("gene","numbat_cnv")
  } else {
    
    numbat_cnv<-NULL
    for(clone in unique(numbat_pseudobulk$sample)){
      #Get per clone annotations
      per_clone<-numbat_pseudobulk[numbat_pseudobulk$sample == clone,]
      
      #Merge in cases with multiple annotations per gene 
      #(when there are multiple SNPs)
      cnv_clone<-per_clone%>%
                  group_by(gene)%>%
                  summarize(numbat_cnv = mean(cnv_state_int))%>%
                  as.data.frame()
      colnames(cnv_clone)<-c("gene",paste0("numbat_cnv_",clone))
      if(is.null(numbat_cnv)){
        numbat_cnv<-cnv_clone
      } else {
        numbat_cnv<-merge(numbat_cnv,cnv_clone,by="gene",all=TRUE)
      }
    }
  }
  
  #Extract gene positions
  gene_position<-as.data.frame(numbat_obj$gtf)
  
  #Match chromosome naming
  gene_position$CHROM<-paste0("chr",gene_position$CHROM)
  
  #Add cnv states
  gene_position<-merge(gene_position, numbat_cnv)
  
  #Delete columns that should not be part of Grange object
  gene_position[,c("gene","gene_length")]<-NULL
  
  #Create a Grange object from genomic positions
  numbat_grange<-makeGRangesFromDataFrame(gene_position,
                                          start.field="gene_start",
                                          end.field="gene_end",
                                          seqnames.field="CHROM",
                                          keep.extra.columns=TRUE)
  
  numbat_grange<-numbat_grange[order(numbat_grange)]

  return(numbat_grange)
}

read_numbat_cnv_individual_cells<-function(file_path){
  
  #Check if numbat found any CNVs (could be empty in diploid tests)
  if(! file.exists(paste(file_path,"bulk_clones_final.tsv.gz",sep="/"))){
    print("No CNVS found for Numbat!")
    return(NULL) 
  }
  
  #create a Numbat object
  numbat_obj<-Numbat$new(file_path)
  
  #remove SNPs without genes assigned
  numbat_pseudobulk<-numbat_obj$bulk_clones[numbat_obj$bulk_clones$gene!='',]
  
  #convert cnv states to index form
  numbat_pseudobulk$cnv_state_int <- -1
  #filter by deletion cnv state
  numbat_pseudobulk[numbat_pseudobulk$cnv_state == "del",]$cnv_state_int <- 1
  #filter by neutral cnv state (and loh since they are not implemented in all methods)
  numbat_pseudobulk[numbat_pseudobulk$cnv_state %in% c("loh", "neu"),]$cnv_state_int <- 2
  #filter by amplification cnv states
  numbat_pseudobulk[numbat_pseudobulk$cnv_state %in% c("amp", "amp|bamp", "bamp"),]$cnv_state_int <- 3
  
  #Get the number of cells per clone
  ncells_clone<-data.frame(table(numbat_pseudobulk$sample,
                                 numbat_pseudobulk$n_cells))
  ncells_clone<-ncells_clone[ncells_clone$Freq!=0,]
  colnames(ncells_clone)<-c("clone","ncells","freq")
  ncells_clone$clone<-as.numeric(as.character(ncells_clone$clone))
  ncells_clone$ncells<-as.numeric(as.character(ncells_clone$ncells))
  
  #Get a per cell CNV state (only possible from clone level)
  numbat_cnv<-NULL
  for(clone in unique(numbat_pseudobulk$sample)){
    
    #Get per clone annotations
    per_clone<-numbat_pseudobulk[numbat_pseudobulk$sample == clone,]
    
    #Merge in cases with multiple annotations per gene 
    #(when there are multiple SNPs)
    cnv_clone<-per_clone%>%
      group_by(gene)%>%
      summarize(numbat_cnv = mean(cnv_state_int))%>%
      as.data.frame()
    
    #Transfer it to a matrix
    ncells<-ncells_clone$ncells[ncells_clone$clone==clone]
    clone_matrix<-matrix(rep(cnv_clone$numbat_cnv,ncells),ncol=ncells)
    rownames(clone_matrix)<-cnv_clone$gene
    
    #Combine the two matrices (need to use merge as rownames are not identical)
    if(is.null(numbat_cnv)){
      numbat_cnv<-clone_matrix
    } else {
      combined_matrix<-merge(numbat_cnv,clone_matrix,by="row.names",all=T)
      genes<-combined_matrix[,1]
      combined_matrix<-as.matrix(combined_matrix[,-1])
      rownames(combined_matrix)<-genes
      colnames(combined_matrix)<-paste("V",1:ncol(combined_matrix))
      numbat_cnv<-combined_matrix
    }
  }
  
  #Extract gene positions and filter them
  gene_position<-as.data.frame(numbat_obj$gtf)
  gene_position$CHROM<-paste0("chr",gene_position$CHROM)
  gene_position<-gene_position[gene_position$gene %in% rownames(numbat_cnv),]
  
  #Create a Grange object from genomic positions
  numbat_grange<-makeGRangesFromDataFrame(gene_position,
                                          start.field="gene_start",
                                          end.field="gene_end",
                                          seqnames.field="CHROM",
                                          keep.extra.columns=TRUE)
  
  numbat_cnv<-numbat_cnv[gene_position$gene,]
  
  return(list(numbat_grange,numbat_cnv))
}

# TODO: Potentially using directly the larger segment from numbat 
# (but missing subclonal information here ...)
# read_numbat_seg<-function(file_path){
#   
#   #create a Numbat object
#   numbat_obj<-Numbat$new(file_path)
#   
#   per_cell_annot<-numbat_obj$joint_post[,c("cell","seg","cnv_state",
#                                            "CHROM","seg_start","seg_end")]
#   
#   per_cell_annot$cnv_state_int<-ifelse(per_cell_annot$cnv_state=="amp",3,
#                                       ifelse(per_cell_annot$cnv_state=="del",1,2))
#   seg_annot<-per_cell_annot%>%
#     group_by(seg,CHROM,seg_start,seg_end)%>%
#     summarize(numbat_seg=mean(cnv_state_int))
# }

read_CONICSmat<-function(input_CONICSmat_chrom_pos, input_CONICSmat_cnv, 
                         input_annot, input_ref_groups,
                         clones=NULL){
  #read data from files
  chrom_pos<-read.table(input_CONICSmat_chrom_pos, sep="\t", header = T)
  cnv_states<-read.table(input_CONICSmat_cnv, sep="\t", row.names = 1, header = T)
  colnames(cnv_states)<-gsub("X", "",colnames(cnv_states))
  annot<-read.table(input_annot, sep="\t", header = F)
  colnames(annot)<-c("barcode", "label")
  ref_groups<-read.table(input_ref_groups, header = T)
  
  #filter out reference samples
  query_samples<-annot[!(annot$label %in% ref_groups[,1]),]
  cnv_states<-cnv_states[(rownames(cnv_states) %in% query_samples[,1]),]
  
  #get a mean value of cnv states for each chromosome arm
  if(is.null(clones)){
    mean_cnv_states<-colMeans(cnv_states)+2
    transposed_cnv<-data.frame(Idf = colnames(cnv_states), 
                               CONICSmat = mean_cnv_states)
  } else {
    #Filter cells
    common_cells<-intersect(rownames(cnv_states),clones$cell)
    cnv_states<-cnv_states[common_cells,]
    clones<-clones[common_cells,]
    
    clones$conicsmat<-as.factor(clones$conicsmat)
    transposed_cnv<-t(apply(cnv_states,2,tapply,clones$conicsmat,mean))+2
    colnames(transposed_cnv)<-paste0("clone",colnames(transposed_cnv))
    transposed_cnv<-as.data.frame(transposed_cnv)
    transposed_cnv$Idf<-rownames(transposed_cnv)
  }
 
  #assemble a dataframe with all relevant data
  CONICSmat_data<-full_join(chrom_pos, transposed_cnv, by = "Idf", multiple = "all")
  
  #substitute NA with 0's
  CONICSmat_data[is.na(CONICSmat_data)] <- 0
  
  #Match chromosome naming
  CONICSmat_data$Chrom<-paste0("chr", CONICSmat_data$Chrom)
  
  #Create a Grange object from chromosome arm positions
  CONICSmat_grange<-makeGRangesFromDataFrame(CONICSmat_data,
                                             start.field="Start",
                                             end.field="End",
                                             seqnames.field="Chrom",
                                             keep.extra.columns=TRUE)
  
  #Delete unnecessary columns
  CONICSmat_grange$Idf<-NULL
  CONICSmat_grange$Length<-NULL
  
  return(CONICSmat_grange)
}

read_CONICSmat_individual_cells<-function(input_CONICSmat_chrom_pos, input_CONICSmat_cnv, 
                         input_annot, input_ref_groups){
  
  #read data from files
  chrom_pos<-read.table(input_CONICSmat_chrom_pos, sep="\t", header = T)
  cnv_states<-read.table(input_CONICSmat_cnv, sep="\t", row.names = 1, header = T)
  colnames(cnv_states)<-gsub("X", "",colnames(cnv_states))
  annot<-read.table(input_annot, sep="\t", header = F)
  colnames(annot)<-c("barcode", "label")
  ref_groups<-read.table(input_ref_groups, header = T)
  
  #filter out reference samples
  query_samples<-annot[!(annot$label %in% ref_groups[,1]),]
  cnv_states<-cnv_states[(rownames(cnv_states) %in% query_samples[,1]),]
  
  #substitute NA with 0's
  cnv_states[is.na(cnv_states)] <- 0
  
  #Convert to numeric
  cnv_states<-cnv_states+2
  
  #Filter chromosome arms to the ones with cnv values
  chrom_pos<-chrom_pos[chrom_pos$Idf %in% colnames(cnv_states),]
  
  #Match chromosome naming
  chrom_pos$Chrom<-paste0("chr", chrom_pos$Chrom)
  
  #Create a Grange object from chromosome arm positions
  CONICSmat_grange<-makeGRangesFromDataFrame(chrom_pos,
                                             start.field="Start",
                                             end.field="End",
                                             seqnames.field="Chrom",
                                             keep.extra.columns=TRUE)
  
  #Delete unnecessary columns
  CONICSmat_grange$Idf<-NULL
  CONICSmat_grange$Length<-NULL
  
  return(list(CONICSmat_grange,t(cnv_states)))
}

#Function to combine output from two different functions 
#bin size of the first range is taken and method column of the second range is transferred
combine_range_objects<-function(range_1,range_2,method_colname){
  
  if(! method_colname %in% colnames(elementMetadata(range_2))){
    stop(paste("Specified column name",method_colname,"not found in range 2!"))
  }
  #Check overlaps (to get the bin sizes to a similar shape)
  overlaps<-as.data.frame(findOverlaps(range_1,range_2,type="any"))
  
  #Filter first range (keep dimensions)
  range_1<-range_1[unique(overlaps$queryHits)]
  
  #Average score of second method across bins
  elementMetadata(range_1)[method_colname] <- sapply(unique(overlaps$queryHits),
                                                     function(i) mean(elementMetadata(range_2)[overlaps$subjectHits[overlaps$queryHits==i], method_colname]))                                 
  
  #Return range 1 with the new additional column from range2                                                   
  return(range_1)
}

#Version of combine_range_objects when all columns of range2 should be added 
#to the grange object (columns need to be numeric)
combine_range_objects_allcols<-function(range_1,range_2){
  
  #Check overlaps (to get the bin sizes to a similar shape)
  overlaps<-as.data.frame(findOverlaps(range_1,range_2,type="any"))
  
  #Filter first range (keep dimensions)
  range_1<-range_1[unique(overlaps$queryHits)]
  
  #Average score of second method across bins
  meta_data_range2<-as.matrix(mcols(range_2))
  combined_scores <- sapply(unique(overlaps$queryHits),
          function(i) {
            indices<-overlaps$subjectHits[overlaps$queryHits==i]
            if(length(indices)==1){
              meta_data_range2[indices,]
            } else {
              colMeans(meta_data_range2[indices,])
            }}) 
  mcols(range_1)<-cbind(mcols(range_1),as.data.frame(t(combined_scores)))
  
  return(range_1)
}


#Evaluate correlation between two methods in the same data frame 
#(specified by column names)
evaluate_correlation_colnames<-function(bin_df,var_name_one,var_name_two,printres=TRUE){
  
  #Calculate correlation values
  pearson_corr<-round(cor(bin_df[[var_name_one]],bin_df[[var_name_two]],
                          method="pearson"),3)
  spearman_corr<-round(cor(bin_df[[var_name_one]],bin_df[[var_name_two]],
                           method="spearman"),3)
  kendall_corr<-round(cor(bin_df[[var_name_one]],bin_df[[var_name_two]],
                          method="kendall"),3)
  
  #Calculate mean square error
  mse <- round(mean((bin_df[[var_name_one]]-bin_df[[var_name_two]])^2),3)
  
  #Print results
  if(printres){
    print(paste("Comparison between",var_name_one,"and",var_name_two))
    print(paste("Pearson correlation:",pearson_corr))
    print(paste("Spearman correlation:",spearman_corr))
    print(paste("Kendell correlation:",kendall_corr))
    print(paste("Mean square error:",mse))
  }
  
  return(c(pearson_corr,spearman_corr,kendall_corr,mse))
  
}

#Evaluate sensitivity and precision for the prediction of gains and losses (separate values)
evaluate_gains_loss<-function(binarized_methods,method_1,method_2="wgs"){
  
  #Evaluate gains (encoded as 3)
  tp_gains <- sum(binarized_methods[,method_1]==3 & binarized_methods[,method_2]==3)
  sens_gains <- tp_gains / sum(binarized_methods[,method_2]==3)
  prec_gains <- tp_gains / sum(binarized_methods[,method_1]==3)
  
  #Evaluate losses (encoded as 1)
  tp_losses <- sum(binarized_methods[,method_1]==1 & binarized_methods[,method_2]==1)
  sens_losses <- tp_losses / sum(binarized_methods[,method_2]==1)
  prec_losses <- tp_losses / sum(binarized_methods[,method_1]==1)
  
  return(c(sens_gains,prec_gains,sens_losses,prec_losses))
}

#Evaluate AUC and AUPRC (again separately for gains and losses)
evaluate_gains_loss_auc<-function(combined_methods,method_1,method_2="wgs"){
  
  library(ROCR)
  
  #Evaluate gains (currentyl defined as >2.5)
  truth <- ifelse(combined_methods[,method_2]>2.5,1,0)
  
  #Check if both classes are present (esp. if there are gain values)
  if(length(table(truth))>1){
    #Allow that AUC thresholds only above 2 are tested
    gain_scores<-combined_methods[,method_1]
    
    pr<-prediction(gain_scores, truth)
    prf<-performance(pr,'auc')
    auc_gains <- unlist(prf@y.values)
    
    prf<-performance(pr,'aucpr')
    aucpr_gains <- unlist(prf@y.values)
  } else {
    auc_gains<-NA
    aucpr_gains<-NA
  }
  
  #Evaluate losses (encoded as 1)
  truth <- ifelse(combined_methods[,method_2]<1.5,0,1)
  
  #Check if both classes are present (esp. if there are loss values)
  if(length(table(truth))>1){
    pr<-prediction(combined_methods[,method_1], truth)
    prf<-performance(pr,'auc')
    auc_losses <- unlist(prf@y.values)
    
    prf<-performance(pr,'aucpr')
    aucpr_losses <- unlist(prf@y.values)
  } else {
    auc_losses<-NA
    aucpr_losses<-NA
  }
  
  return(c(auc_gains,aucpr_gains,auc_losses,aucpr_losses))
}


#Evaluate AUC truncated:evaluate the thresholds only in the logical region 
#(for gain, scores above 2, for loss, scores below 2)
evaluate_gains_loss_auc_truncated<-function(combined_methods,method_1,method_2="wgs"){
  
  library(pROC)
  
  #Evaluate gains (defined as >2.5)
  gain_scores<-combined_methods[,method_1]
  truth <- ifelse(combined_methods[,method_2]>2.5,1,0)
  
  #Check if both classes are present (esp. if there are gain values)
  if(length(table(truth))>1){
    #Calculate sensitivity at the gain threshold
    gain_scores_threshold<-ifelse(gain_scores>2,1,0)
    sens_threshold<-sum(gain_scores_threshold==1 & truth==1)/sum(truth==1)
    
    #Using the pROC package to get the partial AUC curve
    roc_gain <- roc(response=truth, predictor=gain_scores)
    auc_gain <- auc(roc_gain, partial.auc=c(0,sens_threshold), partial.auc.focus="sens")
  } else {
    auc_gain<-NA
  }
    
  #Evaluate loss (defined as <1.5)
  loss_scores <- gain_scores * -1 #invert to get loss as positive class
  truth <- ifelse(combined_methods[,method_2]<1.5,1,0)
  
  #Check if both classes are present (esp. if there are loss values)
  if(length(table(truth))>1){
    #Calculate sensitivity at the loss threshold
    loss_scores_threshold<-ifelse(loss_scores>-2,1,0)
    sens_threshold<-sum(loss_scores_threshold==1 & truth==1)/sum(truth==1)
    
    #Using the pROC package to get the partial AUC curve
    roc_loss <- roc(response=truth, predictor=loss_scores)
    auc_loss <- auc(roc_loss, partial.auc=c(0,sens_threshold), partial.auc.focus="sens")
  } else {
    auc_loss<-NA
  }
  
  return(c(auc_gain,auc_loss))
}


#Calculate the F1 score for the given thresholds
# @params threshold_1 Loss threshold
# @params threshold_2 Gain threshold
get_f1_score_unweighted<-function(continuous_scores,threshold_1,
                                  threshold_2,ground_truth){
  
  require(crfsuite)
  
  #Discretize the score dependent on the thresholds
  discrete_scores<-ifelse(continuous_scores<threshold_1,1,
                          ifelse(continuous_scores>threshold_2,3,2))

  #Multi-class metrics from the crf package
  multi_metrics<-crf_evaluation(pred=discrete_scores,
                                obs=ground_truth)
  
  #Set bylabel F1 NA values to 0 (happens when it predicts no gain/loss at all)
  by_label<-multi_metrics$bylabel
  by_label$f1[is.na(by_label$f1)]<-0
  
  #Return the accuracy, unweighted F1 scores and weighted F1 scores
  return(mean(by_label$f1))
}

# Identify the optimal threshold based on the F1 score per method
# Important remark: only the biological meaning thresholds for gain and loss
# are explored (i.e. gain scores > 2 and loss scores < 2)
evaluate_optimal_f1_score<-function(combined_methods,method_1,method_2="wgs"){
  
  #Identify all thresholds that need to be evaluated
  scores<-sort(unique(combined_methods[,method_1]))
  
  #Only test loss scores < 2
  loss_scores<-scores[scores<2]
  #Subsample scores in case there are too many - otherwise not feasible runtimewise
  if(length(loss_scores)>200){
    loss_scores<-loss_scores[round(seq(1,length(loss_scores),length.out=200))]
  #Set one default loss scores in case no loss score exist (otherwise the function crashes)
  } else if(length(loss_scores)==0){
    loss_scores<-c(1)
  }
  
  #Only test gain scores > 2
  gain_scores<-scores[scores>2]
  #Subsample scores in case there are too many - otherwise not feasible runtimewise
  if(length(gain_scores)>200){
    gain_scores<-gain_scores[round(seq(1,length(gain_scores),length.out=200))]
  #Set one default loss scores in case no loss score exist (otherwise the function crashes)
  } else if(length(gain_scores)==0){
    gain_scores<-c(3)
  }

  #Get all pairwise threshold combinations that should be tested
  combis<-expand.grid(loss_scores,gain_scores)
  
  combis$f1<-sapply(1:nrow(combis),
                             function(i) get_f1_score_unweighted(combined_methods[,method_1],
                                                   combis$Var1[i],
                                                   combis$Var2[i],
                                                   combined_methods[,method_2]))
  
  max_combi<-combis[which.max(combis$f1),]
  
  #Get sens and spec for gains
  tp_gains <- sum(combined_methods[,method_2]==3 & combined_methods[,method_1]>max_combi$Var2)
  sens_gains <- tp_gains / sum(combined_methods[,method_2]==3)
  prec_gains <- tp_gains / sum(combined_methods[,method_1]>max_combi$Var2)
  
  #Get sens and spec for losses
  tp_losses <- sum(combined_methods[,method_1]<max_combi$Var1 & combined_methods[,method_2]==1)
  sens_losses <- tp_losses / sum(combined_methods[,method_2]==1)
  prec_losses <- tp_losses / sum(combined_methods[,method_1]<max_combi$Var1)
  
  return(c(max_combi$f1,max_combi$Var2,sens_gains,prec_gains,max_combi$Var1,sens_losses,prec_losses))
  
}
  
#Create a binned version of the whole genome for the respective Grange object
create_binned_genome <- function(grange_obj,bin_size=100000){
  
  require(GenomicRanges)
  require(tidyverse)
  
  #Determine the covered region of each chromsome
  grange_df<-as.data.frame(grange_obj)
  gr_limits<-grange_df%>%group_by(seqnames)%>%
    summarize(global_start=min(start),global_end=max(end))
  
  #Round the borders to fit the bin_size (to have bins starting at 0)
  gr_limits$global_start_rounded<-floor(gr_limits$global_start/bin_size)*bin_size
  gr_limits$global_end_rounded<-ceiling(gr_limits$global_end/bin_size)*bin_size-1
  
  #Use Genomic Ranges to create bin_size windows
  limits_gr<-makeGRangesFromDataFrame(gr_limits,
                                          start.field="global_start_rounded",
                                          end.field="global_end_rounded")
  limits_gr<-limits_gr[order(limits_gr)]
  
  gr_bins<-tile(limits_gr,width=bin_size)
  gr_bins<-unlist(gr_bins)
  
  return(gr_bins)
}

#Recursive adjusted rand index - merge cluster until the "correct number of clones"
#Group 1 = ground truth with correct number of clones
recursive_ari<-function(group1,group2){
  
  #Remove NA values
  group1<-group1[! is.na(group2)]
  group2<-group2[! is.na(group2)]
  
  #Correct number of clones
  n_clones<-length(unique(group1))
  
  #Merge subclones from group2 until the same number of clones as in group1 is found
  while(length(unique(group2))>n_clones){
    clusters<-unique(group2)
    
    max_ari<-NULL
    for(i in 1:(length(clusters)-1)){
      for(j in (i+1):length(clusters)){
        tmp<-ifelse(group2==clusters[j],clusters[i],group2)
        max_ari<-rbind(max_ari,
                       data.frame(i,j,ari=mclust::adjustedRandIndex(group1,tmp)))
        
      }
    }
    
    ari_pos<-max_ari[which.max(max_ari$ari),]
    group2<-ifelse(group2==clusters[ari_pos$i] | group2==clusters[ari_pos$j],
                   paste0(clusters[ari_pos$i],clusters[ari_pos$j]),group2)
  }
  
  return(mclust::adjustedRandIndex(group1,group2))
  
}

# Get evaluation scores for each individual cell compared to a WGS ground truth
# Calculate Pearson correlation, AUC gain and loss as well as AUC trunc gain and loss
evaluate_metrics_per_cell<-function(wgs_grange,method_grange,method_matrix){
  
  require(pROC)
  
  #Combine WGS and RNA results
  overlaps<-as.data.frame(findOverlaps(wgs_grange,method_grange,type="any"))
  
  #Filter first range (keep dimensions)
  wgs_filtered<-wgs_grange[unique(overlaps$queryHits)]
  
  #Map per gene counts for the methods to the WGS bins
  combined_scores <- sapply(unique(overlaps$queryHits),
                            function(i) {
                              indices<-overlaps$subjectHits[overlaps$queryHits==i]
                              if(length(indices)==1){
                                method_matrix[indices,]
                              } else {
                                colMeans(method_matrix[indices,])
                              }}) 
  combined_scores<-t(combined_scores)
  
  #Get the pearson correlation estimates per cell
  cor_per_cell<-apply(combined_scores,2,function(i)cor(i,wgs_filtered$wgs_mean,
                                                       use="complete.obs"))
  
  #Get the AUC gain values per cell
  wgs_filtered$bin_gain<-as.numeric(wgs_filtered$wgs_mean>2.5)
  auc_gain_per_cell<-apply(combined_scores,2,function(i){
    
    #Remove NA values
    df<-data.frame(rna=i,
                   wgs=wgs_filtered$bin_gain)
    df<-df[!is.na(df$rna),]
    
    #Calculate AUC
    pr<-prediction(df$rna, df$wgs)
    prf<-performance(pr,'auc')
    return(unlist(prf@y.values))})
  
  #Get the AUC loss values per cell
  wgs_filtered$bin_loss<-as.numeric(wgs_filtered$wgs_mean>=1.5)
  auc_loss_per_cell<-apply(combined_scores,2,function(i){
    
    #Remove NA values
    df<-data.frame(rna=i,
                   wgs=wgs_filtered$bin_loss)
    df<-df[!is.na(df$rna),]
    
    pr<-prediction(df$rna, df$wgs)
    prf<-performance(pr,'auc')
    return(unlist(prf@y.values))})
  
  #Get truncated AUC gain value per cell
  auc_trunc_gain_per_cell<-apply(combined_scores,2,function(i){
    
    #Remove NA values
    df<-data.frame(rna=i,
                   wgs=wgs_filtered$bin_gain)
    df<-df[!is.na(df$rna),]
    
    #Calculate sensitivity at the gain threshold
    gain_scores_threshold<-ifelse(df$rna>2,1,0)
    sens_threshold<-sum(gain_scores_threshold==1 & df$wgs==1)/sum(df$wgs==1)
    
    #Using the pROC package to get the partial AUC curve
    suppressMessages(roc_gain <- roc(response=df$wgs, predictor=df$rna))
    auc_gain <- auc(roc_gain, partial.auc=c(0,sens_threshold), partial.auc.focus="sens")
    return(as.numeric(auc_gain))})
  
  #Get truncated AUC gain value per cell
  auc_trunc_loss_per_cell<-apply(combined_scores,2,function(i){
    
    #Remove NA values (invert so that loss is the positive class)
    df<-data.frame(rna=i*(-1),
                   wgs=wgs_filtered$bin_loss*(-1))
    df<-df[!is.na(df$rna),]
    
    #Calculate sensitivity at the loss threshold
    loss_scores_threshold<-ifelse(df$rna>-2,1,0)
    sens_threshold<-sum(loss_scores_threshold==1 & df$wgs==0)/sum(df$wgs==0)
    
    #Using the pROC package to get the partial AUC curve
    suppressMessages(roc_loss <- roc(response=df$wgs, predictor=df$rna))
    auc_loss <- auc(roc_loss, partial.auc=c(0,sens_threshold), partial.auc.focus="sens")
    return(as.numeric(auc_loss))})
  
  res_df<-data.frame(cor=cor_per_cell,
             auc_gain=auc_gain_per_cell,
             auc_loss=auc_loss_per_cell,
             auc_trunc_gain=auc_trunc_gain_per_cell,
             auc_trunc_loss=auc_trunc_loss_per_cell)
  return(res_df)
}

#Combine the different metrics and run them all for one specific methods
get_all_metrics<-function(combined_range, method_name){
  
  combined_methods<-elementMetadata(combined_range)
  
  #Get correlation results
  corr<-evaluate_correlation_colnames(combined_methods,"wgs_mean",method_name,printres=FALSE)
  
  #Get AUC and AUPRC values (looking at the complete curve)
  auc<-evaluate_gains_loss_auc(combined_methods,method_name,"wgs_mean")
  
  #Get truncated AUC values (only in the biological fair range)
  auc_trunc<-evaluate_gains_loss_auc_truncated(combined_methods,method_name,"wgs_mean")
  
  #Binarize WGS score
  binarized_methods<-combined_methods
  binarized_methods$wgs_mean<-ifelse(binarized_methods$wgs_mean<1.5,1,
                                     ifelse(binarized_methods$wgs_mean>2.5,3,2))
  
  #Optimal F1 score
  f1<-evaluate_optimal_f1_score(binarized_methods,method_name,"wgs_mean")
  
  res<-data.frame(method=method_name,
                  num_bins=length(combined_range),
                  pearson=corr[1],
                  auc_gains=auc[1],
                  auc_gains_trunc = auc_trunc[1],
                  aucpr_gains=auc[2],
                  auc_losses=auc[3],
                  auc_losses_trunc = auc_trunc[2],
                  aucpr_losses=auc[4],
                  max_f1=f1[1],
                  cutoff_f1_gain=f1[2],
                  sens_gains=f1[3],
                  prec_gains=f1[4],
                  cutoff_f1_loss=f1[5],
                  sens_losses=f1[6],
                  prec_losses=f1[7])
  
  return(res)
}

