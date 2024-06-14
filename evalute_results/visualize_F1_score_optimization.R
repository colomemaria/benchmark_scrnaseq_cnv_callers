# ------------------------------------------------------------------------------
# Create a example visualization how the F1 score is calculated
# ------------------------------------------------------------------------------

library(fread)
library(ggplot2)

source("scripts/lib.R")

#Example dataset
res<-fread("results/output_SNU601/evaluation/outputs_allmethods_combined.tsv")

#Binarize WGS results
res$wgs_mean<-ifelse(res$wgs_mean<1.5,1,ifelse(res$wgs_mean>2.5,3,2))

# ------------------------------------------------------------------------------
# Code from evaluate_optimal_f1_score
# ------------------------------------------------------------------------------

combined_methods<-as.data.frame(res)
method_1<-"infercnv_cnv"
method_2<-"wgs_mean"

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

#Create a plot with results
g<-ggplot(combis,aes(x=Var1,y=Var2,color=f1^10))+
  geom_point()+
  scale_color_viridis("F1 score",option="turbo",breaks=c(0.5^10,0.65^10),
                      labels=c(0.5,0.65))+
  xlab("Loss score")+ylab("Gain score")
ggsave(g,filename="~/Desktop/optimize_f1.png",width=5,height=4)

#Create a histogram showing the chosen cutoffs
opt_combi<-combis[which.max(combis$f1),]

g<-ggplot(combined_methods,aes(x=infercnv_cnv))+
  geom_histogram()+
  geom_vline(xintercept=c(opt_combi$Var1+0.05,opt_combi$Var2))+
  ylab("Number of regions")+
  xlab("Pseudobulk CNV score for InferCNV")
ggsave(g,filename="~/Desktop/optimize_f1_cutoffs.png",width=5,height=4)
