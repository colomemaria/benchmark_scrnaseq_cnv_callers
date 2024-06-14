# ------------------------------------------------------------------------------
# Compare whether the scaling factors for each method are consistent
# across datasets
# ------------------------------------------------------------------------------

library(data.table)
library(ggplot2)

theme_set(theme_bw())

dataset_ending<-setNames(c("","_wohb","_WES","_WES","_cnvkit"),
                         c("SNU601","MCF7","MM","BCC08","COLO320"))

#Load the files
scaling_factors<-NULL
for(dataset in names(dataset_ending)){
  
  tmp<-fread(paste0("results/output_",dataset,
                      "/evaluation/scaling_factors",
                      dataset_ending[dataset],".tsv"))
  tmp$dataset<-dataset
  scaling_factors<-rbind(scaling_factors,tmp)
}

#Group WGS and WES results together
scaling_factors$method[scaling_factors$method %in% c("wgs_mean","wes_lfc")]<-"ground_truth"

#Order methods by SNU601 scaling factor size
tmp<-scaling_factors[scaling_factors$dataset=="SNU601",]
method_order<-tmp$method[order(tmp$sd)]
scaling_factors$method<-factor(scaling_factors$method,
                               levels=method_order)

g<-ggplot(scaling_factors,aes(x=method,y=sd,color=dataset))+
  geom_point()+ylab("Standard deviation (unnormalized)")+xlab("Method")+
  theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))

ggsave(g,file="~/Desktop/scaling_factors.png",width=6,height=5)


