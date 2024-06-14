# ------------------------------------------------------------------------------
# Test different versions for multi-class performance metrics for one example
# and decide afterwards which ones to include in the general workflow
# ------------------------------------------------------------------------------

library(data.table)
library(dplyr)
library(GenomicRanges)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(ROCR)
library(viridis)
library(numbat)

library(crfsuite) #new package to calculate multi-class metrics

theme_set(theme_bw())

#Source script with help functions
source("scripts/lib.R")

#Load test data: as usual SNU601 and inferCNV
input_wgs<-"../../../results/SNU601_scWGS_results/wgs_results_aneufinder_pseudobulk.csv"
input_infercnv_cnv<-"results/output_SNU601/infercnv/infercnv.20_HMM_predHMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.repr_intensities.observations.txt"
input_infercnv_gene_pos<-"data/annotations/hg38_gencode_v27.txt"
  
wgs_results<-read_wgs_results(input_wgs)

#Distribution WGS scores
plot_wgs<-elementMetadata(wgs_results)
plot_wgs<-as.data.frame(plot_wgs)
g<-ggplot(plot_wgs,aes(x=wgs_mean))+
  xlab("Pseudobulk CNV score for scWGS")+ylab("Number regions")+
  geom_histogram(bins=50)+
  geom_vline(xintercept=c(1.5,2.5))
ggsave(g,file="~/Desktop/wgs_cnv_scores.png",height=3,width=4)
  
infercnv_results<-read_infercnv_6state_model(input_infercnv_cnv,
                                             input_infercnv_gene_pos)

#Distribution InferCNV scores
plot_wgs<-elementMetadata(infercnv_results[[1]])
plot_wgs<-as.data.frame(plot_wgs)
g<-ggplot(plot_wgs,aes(x=infercnv_cnv))+
  xlab("Pseudobulk CNV score for InferCNV")+ylab("Number regions")+
  geom_histogram(bins=50)+
  geom_vline(xintercept=c(1.05,2.6))
ggsave(g,file="~/Desktop/infercnv_scores.png",height=3,width=4)

combined_range<-combine_range_objects(wgs_results,infercnv_results[[1]],
                                      method_colname="infercnv_cnv")
combined_methods<-elementMetadata(combined_range)

#Iterate over each possible score combinations
combined_methods$wgs_mean<-ifelse(combined_methods$wgs_mean<1.5,1,
                                   ifelse(combined_methods$wgs_mean>2.5,3,2))
scores<-sort(unique(combined_methods$infercnv_cnv))

#Get all pairwise combinations that should be tested
combis<-expand.grid(1:(length(scores)-1),2:length(scores))
combis<-combis[combis$Var1<combis$Var2,]

#Function to calculate the evaluation metrics
get_score<-function(index){
  combined_methods$infercnv_discrete<-ifelse(combined_methods$infercnv_cnv<scores[combis$Var1[index]],1,
                              ifelse(combined_methods$infercnv_cnv>scores[combis$Var2[index]],3,2))
  
  #Get the accuracy
  acc<-mean(combined_methods$infercnv_discrete==combined_methods$wgs_mean)
  
  #Multi-class metrics from the crf package
  multi_metrics<-crf_evaluation(pred=combined_methods$infercnv_discrete,
                                obs=combined_methods$wgs_mean)
  
  #Set bylabel F1 NA values to 0 (happens when it predicts no gain/loss at all)
  by_label<-multi_metrics$bylabel
  by_label$f1[is.na(by_label$f1)]<-0
  
  #Return the accuracy, unweighted F1 scores and weighted F1 scores
  return(c(acc,mean(by_label$f1),
           sum(by_label$f1*by_label$support)/sum(by_label$support)))
}

multi_class_scores<-sapply(1:nrow(combis),get_score)
multi_class_scores<-t(multi_class_scores)
multi_class_scores<-as.data.frame(multi_class_scores)
colnames(multi_class_scores)<-c("acc","f1","f1_mean")

multi_class_scores$score_loss<-scores[combis$Var1]
multi_class_scores$score_gain<-scores[combis$Var2]

multi_class_scores[which.max(multi_class_scores$acc),]
multi_class_scores[which.max(multi_class_scores$f1),]
multi_class_scores[which.max(multi_class_scores$f1_mean),]

g<-ggplot(multi_class_scores,aes(x=as.factor(score_loss),y=as.factor(score_gain),fill=acc))+
  geom_tile()+xlab("Loss score")+ylab("Gain score")+
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank())
ggsave(g,file="~/Desktop/accuracy_2d.png",height=3,width=4)

g<-ggplot(multi_class_scores,aes(x=as.factor(score_loss),y=as.factor(score_gain),fill=f1))+
  geom_tile()+xlab("Loss score")+ylab("Gain score")+
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank())
ggsave(g,file="~/Desktop/f1_unweighted_2d.png",height=3,width=4)

g<-ggplot(multi_class_scores,aes(x=as.factor(score_loss),y=as.factor(score_gain),fill=f1_mean))+
  geom_tile()+xlab("Loss score")+ylab("Gain score")+
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank())
