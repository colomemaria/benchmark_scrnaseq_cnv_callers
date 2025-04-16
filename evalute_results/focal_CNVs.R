# ------------------------------------------------------------------------------
# Explore the CNV types for each dataset with scWGS data
# ------------------------------------------------------------------------------

library(data.table)
library(GenomicRanges)
#library(ggplot2)

source("scripts/lib.R")

dataset_names<-c("SNU601","NCIN87","MKN45","KATOIII",
                 "NUGC4","SNU638","SNU668","HGC27", "SNU16","iAMP21")

#Get chromosome arm definition (from CONICSmat, but why not)
chr_arm<-fread("data/annotations/chromosome_arm_positions_grch38.txt")
chr_arm$Chrom<-paste0("chr",chr_arm$Chrom)

chr_arm_gr<-makeGRangesFromDataFrame(chr_arm,keep.extra.columns = TRUE)

focal_dists<-NULL
summary_eval_focal<-NULL
summary_eval_focal2<-NULL
summary_eval_focal3<-NULL
for(dataset in dataset_names){
  
  print(dataset)
  
  pb_wgs<-fread(paste0("data/",dataset,"_WGS_groundtruth/wgs_results_formated.csv"))
  
  #Reformat chromosome column
  pb_wgs$chr<-paste0("chr",pb_wgs$chr)
  
  #Discretize the WGS status into 1,2,3
  pb_wgs$wgs_mean<-with(pb_wgs,loss_wgs+(base_wgs*2)+gain_wgs*3)
  pb_wgs$wgs_discrete<-ifelse(pb_wgs$wgs_mean<1.5,1,
                              ifelse(pb_wgs$wgs_mean>2.5,3,2))
  pb_wgs<-pb_wgs[,c("chr","start","end","wgs_discrete")]
  
  pb_grange<-makeGRangesFromDataFrame(pb_wgs,keep.extra.columns = TRUE)
  
  #Creates a Grange list with three entries for each wgs discrete class
  pb_splitted<-split(pb_grange,pb_grange$wgs_discrete)
  
  #Merge entries next to each other
  pb_reduced <- lapply(pb_splitted, reduce)
  
  #Get the losses and gains
  cnv_lengths<-rbind(data.frame(cnv_width=width(pb_reduced[["1"]]), class="loss"),
                     data.frame(cnv_width=width(pb_reduced[["3"]]), class="gain"))
  
  # ----------------------------------------------------------------------------
  # Definition 1: focal meaning smaller than 3 MB
  # ----------------------------------------------------------------------------
  
  #Check how many are less than 3Mbp long
  cnv_lengths$cnv_focal<-(cnv_lengths$cnv_width<=3e6)
  
  # 1) all CNVs
  pb_all_cnvs <- lapply(names(pb_reduced), function(status) {
    gr <- pb_reduced[[status]]
    mcols(gr)$focal_cnvs2 <- status
    gr
  })
  
  #Save for later all the regions
  pb_all_combined<-c(pb_all_cnvs[[1]],pb_all_cnvs[[2]],pb_all_cnvs[[3]])
  pb_all_combined<-sort(pb_all_combined)
  
  #Then combine back to one Grange (only the losses and gains)
  pb_all_cnvs<-c(pb_all_cnvs[[1]],pb_all_cnvs[[3]])
  pb_all_cnvs<-sort(pb_all_cnvs)
  
  #Convert the status to numeric
  pb_all_cnvs$focal_cnvs2<-as.numeric(pb_all_cnvs$focal_cnvs2)
  

  # 2) Focal CNVS:
  # Add 1 for focal losses, 3 for focal gains and 2 otherwise for combining it again to one grange
  pb_focal <- lapply(names(pb_reduced), function(status) {
    gr <- pb_reduced[[status]]
    mcols(gr)$focal_cnvs <- ifelse(width(gr)<=3e6,status,2)
    gr
  })
  
  #Combine back to one Grange
  pb_focal_combined<-c(pb_focal[[1]],pb_focal[[2]],pb_focal[[3]])
  
  #Convert the status to numeric
  pb_focal_combined$focal_cnvs<-as.numeric(pb_focal_combined$focal_cnvs)
  
  #Combine with existing Grange containing RNA data
  eval_res<-fread(paste0("results/output_",dataset,"/evaluation/outputs_allmethods_combined.tsv"))
  eval_gr<-makeGRangesFromDataFrame(eval_res,keep.extra.columns = TRUE)
  
  #Check how many of the focal CNVs overlap the genomic region
  pb_focal_filtered<-pb_focal_combined[pb_focal_combined$focal_cnvs!=2,]
  focal_overlaps<-as.data.frame(findOverlaps(pb_focal_filtered, eval_gr))
  pb_focal_filtered_exp<-pb_focal_filtered[unique(focal_overlaps$queryHits)]
  
  #Check the same for the large CNVs
  pb_all_filtered<-pb_all_cnvs[pb_all_cnvs$focal_cnvs2!=2,]
  all_overlaps<-as.data.frame(findOverlaps(pb_all_filtered, eval_gr))
  pb_all_filtered_exp<-pb_all_filtered[unique(all_overlaps$queryHits)]

  eval_gr<-combine_range_objects(eval_gr,pb_focal_combined,method_colname="focal_cnvs")
  
  #Evaluate each method for sensitivity and specificity using the F1 score determined before
  f1_scores<-fread(paste0("results/output_",dataset,"/evaluation/evaluation_cnv_prediction_f1.txt"))
  
  eval_focal<-NULL
  for(meth in f1_scores$method){
    gain_cutoff <- f1_scores$cutoff_f1_gain[f1_scores$method==meth]
    loss_cutoff<- f1_scores$cutoff_f1_loss[f1_scores$method==meth]
    
    #Get sens and spec for gains
    tp_gains <- sum(eval_gr$focal_cnvs == 3 & mcols(eval_gr)[,meth]>gain_cutoff)
    sens_gains <- tp_gains / sum(eval_gr$focal_cnvs==3)
    #prec_gains <- tp_gains / sum(mcols(eval_gr)[,meth]>gain_cutoff)
    
    #Get sens and spec for losses
    tp_losses <- sum(mcols(eval_gr)[,meth]<loss_cutoff & eval_gr$focal_cnvs==1)
    sens_losses <- tp_losses / sum(eval_gr$focal_cnvs==1)
    #prec_losses <- tp_losses / sum(mcols(eval_gr)[,meth]<loss_cutoff)
    
    eval_focal<-rbind(eval_focal,
                      data.frame(dataset,
                                 method=meth,
                                 sens_gains,
                                 sens_losses))
  }
  
  summary_eval_focal<-rbind(summary_eval_focal,eval_focal)
  
  
  # ----------------------------------------------------------------------------
  # Definition 2: focal CNVs are smaller than the chromosome arm
  # ----------------------------------------------------------------------------
  
  #Annotate this is the data frame again
  overlaps<-as.data.frame(findOverlaps(eval_gr,pb_all_cnvs))
  
  eval_gr$merged_cnvs<-2
  eval_gr$merged_cnvs[overlaps$queryHits]<-pb_all_cnvs$focal_cnvs2[overlaps$subjectHits]
  
  inside_chrarm<-findOverlaps(pb_all_cnvs,chr_arm_gr,type="within")
  inside_chrarm<-as.data.frame(inside_chrarm)
  
  #Filter for length (using two different thresholds, 50% and 95%)
  inside_smaller<-inside_chrarm[width(pb_all_cnvs)[inside_chrarm$queryHits] <
                                  (0.50*width(chr_arm_gr)[inside_chrarm$subjectHits]),]
  
  pb_focal_2<-pb_all_cnvs[inside_smaller$queryHits]
  
  inside_smaller<-inside_chrarm[width(pb_all_cnvs)[inside_chrarm$queryHits] <
                                  (0.95*width(chr_arm_gr)[inside_chrarm$subjectHits]),]
  
  pb_focal_3<-pb_all_cnvs[inside_smaller$queryHits]
  
  focal_dists<-rbind(focal_dists,
                     data.frame(dataset,
                                all_gains_unfiltered=sum(cnv_lengths$class=="gain"),
                                all_losses_unfiltered=sum(cnv_lengths$class=="loss"),
                                focal1_gains_unfiltered=sum(cnv_lengths$class=="gain" & cnv_lengths$cnv_focal),
                                focal1_losses_unfiltered=sum(cnv_lengths$class=="loss" & cnv_lengths$cnv_focal),
                                all_gains=sum(pb_all_filtered_exp$focal_cnvs2==3),
                                all_losses=sum(pb_all_filtered_exp$focal_cnvs2==1),
                                focal1_gains=sum(pb_focal_filtered_exp$focal_cnvs==3),
                                focal1_losses=sum(pb_focal_filtered_exp$focal_cnvs==1),
                                focal2_gains=sum(pb_focal_2$focal_cnvs2==3),
                                focal2_losses=sum(pb_focal_2$focal_cnvs2==1),
                                focal3_gains=sum(pb_focal_3$focal_cnvs2==3),
                                focal3_losses=sum(pb_focal_3$focal_cnvs2==1)))
  
  
  #Find the non-unique intersect again
  overlaps<-as.data.frame(findOverlaps(pb_all_combined, pb_focal_2))
  pb_all_combined$focal_cnvs2<-2  
  pb_all_combined$focal_cnvs2[overlaps$queryHits]<-pb_focal_2$focal_cnvs2[overlaps$subjectHits]
  
  eval_gr<-combine_range_objects(eval_gr,pb_all_combined,method_colname="focal_cnvs2")
  
  eval_focal<-NULL
  for(meth in f1_scores$method){
    gain_cutoff <- f1_scores$cutoff_f1_gain[f1_scores$method==meth]
    loss_cutoff<- f1_scores$cutoff_f1_loss[f1_scores$method==meth]
    
    #Get sens and spec for gains
    tp_gains <- sum(eval_gr$focal_cnvs2 == 3 & mcols(eval_gr)[,meth]>gain_cutoff)
    sens_gains <- tp_gains / sum(eval_gr$focal_cnvs2==3)
    #prec_gains <- tp_gains / sum(mcols(eval_gr)[,meth]>gain_cutoff)
    
    #Get sens and spec for losses
    tp_losses <- sum(mcols(eval_gr)[,meth]<loss_cutoff & eval_gr$focal_cnvs2==1)
    sens_losses <- tp_losses / sum(eval_gr$focal_cnvs2==1)
    #prec_losses <- tp_losses / sum(mcols(eval_gr)[,meth]<loss_cutoff)
    
    eval_focal<-rbind(eval_focal,
                      data.frame(dataset,
                                 method=meth,
                                 sens_gains,
                                 sens_losses))
  }
  
  summary_eval_focal2<-rbind(summary_eval_focal2,eval_focal)
  
  
  #Find the non-unique intersect again
  overlaps<-as.data.frame(findOverlaps(pb_all_combined, pb_focal_3))
  pb_all_combined$focal_cnvs3<-2  
  pb_all_combined$focal_cnvs3[overlaps$queryHits]<-pb_focal_3$focal_cnvs2[overlaps$subjectHits]
  
  eval_gr<-combine_range_objects(eval_gr,pb_all_combined,method_colname="focal_cnvs3")
  
  eval_focal<-NULL
  for(meth in f1_scores$method){
    gain_cutoff <- f1_scores$cutoff_f1_gain[f1_scores$method==meth]
    loss_cutoff<- f1_scores$cutoff_f1_loss[f1_scores$method==meth]
    
    #Get sens and spec for gains
    tp_gains <- sum(eval_gr$focal_cnvs3 == 3 & mcols(eval_gr)[,meth]>gain_cutoff)
    sens_gains <- tp_gains / sum(eval_gr$focal_cnvs3==3)
    #prec_gains <- tp_gains / sum(mcols(eval_gr)[,meth]>gain_cutoff)
    
    #Get sens and spec for losses
    tp_losses <- sum(mcols(eval_gr)[,meth]<loss_cutoff & eval_gr$focal_cnvs3==1)
    sens_losses <- tp_losses / sum(eval_gr$focal_cnvs3==1)
    #prec_losses <- tp_losses / sum(mcols(eval_gr)[,meth]<loss_cutoff)
    
    eval_focal<-rbind(eval_focal,
                      data.frame(dataset,
                                 method=meth,
                                 sens_gains,
                                 sens_losses))
  }
  
  summary_eval_focal3<-rbind(summary_eval_focal3,eval_focal)
  
}

#Save the two result files
write.table(focal_dists,file="results/focal_cnvs_distribution.tsv",
            quote=FALSE,sep="\t",row.names = FALSE)
write.table(summary_eval_focal,file="results/focal_definition1_cnvs_evaluation.tsv",
            quote=FALSE,sep="\t",row.names = FALSE)
write.table(summary_eval_focal2,file="results/focal_definition2_cnvs_evaluation.tsv",
            quote=FALSE,sep="\t",row.names = FALSE)
write.table(summary_eval_focal3,file="results/focal_definition3_cnvs_evaluation.tsv",
            quote=FALSE,sep="\t",row.names = FALSE)

