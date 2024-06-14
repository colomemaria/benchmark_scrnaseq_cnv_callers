# ------------------------------------------------------------------------------
# Visualize the performance of the methods across the different datasets
# ------------------------------------------------------------------------------

library(data.table)
library(dplyr)
library(ggplot2)
library(viridis)

theme_set(theme_bw())

#Load the different result files
dataset_names<-c("SNU601","MCF7","MM","BCC06","BCC06post","COLO320")
dataset_ending<-setNames(c("","","_WES","_WES","_WES","_cnvkit"),
                         dataset_names)

#Vector for renaming methods (official published names)
method_names<-setNames(c("InferCNV (CNV)","InferCNV (Expr)","CaSpER",
                         "copyKat","SCEVAN (CNV)","SCEVAN (Expr)",
                         "Numbat (Expr)", "Numbat (CNV)", "CONICSmat"),
                       c("infercnv_cnv","infercnv_expr","casper",
                         "copykat","scevan_vega","scevan",
                         "numbat","numbat_cnv", "CONICSmat"))

#Match resources (sometimes multiple variants for one method)
mapping_resources<-data.frame(methods_full=c("InferCNV (CNV)","InferCNV (Expr)","CaSpER",
                                         "copyKat","SCEVAN (CNV)","SCEVAN (Expr)",
                                         "Numbat (Expr)", "Numbat (CNV)", "CONICSmat"),
                              methods_reduced=c("InferCNV","InferCNV","CaSpER",
                                         "copyKat","SCEVAN","SCEVAN",
                                         "Numbat","Numbat", "CONICSmat"))

#Load the files
all_evals<-NULL
f1_optimal_cutoffs<-NULL
for(dataset in dataset_names){
  
  corrs<-fread(paste0("results/output_",dataset,
                      "/evaluation/evaluation_cnv_prediction_corr",
                      dataset_ending[dataset],".tsv"))
  
  #Filter for comparison with ground truth
  corrs<-corrs[corrs$method1 %in% c("scWGS","WGS","WES","GATK")]
  corrs<-corrs[,c("method2","pearson")]
  colnames(corrs)<-c("method","value")
  corrs$metric<-"pearson"
  corrs$dataset<-dataset
  
  #Set negative correlation values to 0
  corrs$norm_value<-corrs$value
  corrs$norm_value[corrs$norm_value<0]<-0
  
  #Remove row with scWGS/WES/WGS results
  corrs<-corrs[! corrs$method %in% c("scWGS","WGS","WES","GATK"),]
  all_evals<-rbind(all_evals,corrs)
  
  #Add also AUC values (separate values for gain and loss)
  auc<-fread(paste0("results/output_",dataset,
                    "/evaluation/evaluation_cnv_prediction_auc",
                    dataset_ending[dataset],".txt"))
  auc$method<-method_names[auc$method]
  auc<-melt(auc[,c("method","auc_gains","auc_losses")],id.vars="method")
  auc<-auc[,c("method","value","variable")]
  colnames(auc)[3]<-"metric"
  auc$dataset<-dataset
  
  #Scale AUC to also range from 0-1 (set value < 0.5 to 0)
  auc$norm_value<-auc$value
  auc$norm_value<-(auc$norm_value-0.5)*2
  auc$norm_value[auc$norm_value<0]<-0
  
  all_evals<-rbind(all_evals,auc)
  
  #Add sens and spec values
  sensspec<-fread(paste0("results/output_",dataset,
                      "/evaluation/evaluation_cnv_prediction_senspec",
                      dataset_ending[dataset],".txt"))
  sensspec<-melt(sensspec,id.vars="method")
  sensspec<-sensspec[,c("method","value","variable")]
  colnames(sensspec)[3]<-"metric"
  sensspec$dataset<-dataset
  sensspec$norm_value<-sensspec$value
  all_evals<-rbind(all_evals,sensspec)
  
  #Add sens and spec values for F1 scores
  filename<-paste0("results/output_",dataset,
                   "/evaluation/evaluation_cnv_prediction_f1",
                   dataset_ending[dataset],".txt")
  if(file.exists(filename)){
    sensspecf1<-fread(filename)
    sensspecf1$method<-method_names[sensspecf1$method]
    
    #Save separately the optimal cutoffs
    f1_cutoff<-sensspecf1[,c("method","cutoff_f1_gain","cutoff_f1_loss")]
    f1_cutoff$dataset<-dataset
    f1_optimal_cutoffs<-rbind(f1_optimal_cutoffs,f1_cutoff)
    
    colnames(sensspecf1)[c(4,5,7,8)]<-paste0(colnames(sensspecf1)[c(4,5,7,8)],
                                             "_f1")
    sensspecf1<-melt(sensspecf1[,!c("cutoff_f1_gain","cutoff_f1_loss")],id.vars="method")
    sensspecf1<-sensspecf1[,c("method","value","variable")]
    colnames(sensspecf1)[3]<-"metric"
    sensspecf1$dataset<-dataset
    sensspecf1$norm_value<-sensspecf1$value
    all_evals<-rbind(all_evals,sensspecf1)
  }
  
  #Get memory and runtime
  filename<-paste0("results/output_",dataset,
                   "/evaluation/",dataset,"_resources_required.txt")
  if(file.exists(filename)){

    resources<-fread(filename)
    resources<-resources[,c("method","runtime_h","max_vms")]
    
    #Remove methods run without ground truth
    resources<-resources[! endsWith(resources$method,"(wo ref)"),]
    
    #Convert max_vms from MB to GB (for better readability)
    resources$max_vms<-resources$max_vms/1024
    
    #Scale runtime and memory both between 0 and 1
    resources$runtime_h_norm<-(resources$runtime_h-min(resources$runtime_h))/
      (max(resources$runtime_h)-min(resources$runtime_h))
    resources$max_vms_norm<-(resources$max_vms-min(resources$max_vms))/
      (max(resources$max_vms)-min(resources$max_vms))
    
    #Revert the scale so that efficient methods get a score of 1
    resources$runtime_h_norm<-1-resources$runtime_h_norm
    resources$max_vms_norm<-1-resources$max_vms_norm
    
    #Some methods are evaluated multiple times (duplicate these entries)
    resources<-merge(resources,mapping_resources,
                     by.x="method",by.y="methods_reduced")
    resources$method<-NULL
    
    #Reformat dataframe to match other entries
    resources_time<-resources[,c("methods_full","runtime_h","runtime_h_norm")]
    resources_time$metric<-"runtime_h"
    resources_time$dataset<-dataset
    colnames(resources_time)<-c("method","value","norm_value","metric","dataset")
    resources_time<-resources_time[,c("method","value","metric","dataset","norm_value")]
    all_evals<-rbind(all_evals,resources_time)
    
    resources_mem<-resources[,c("methods_full","max_vms","max_vms_norm")]
    resources_mem$metric<-"max_vms_gb"
    resources_mem$dataset<-dataset
    colnames(resources_mem)<-c("method","value","norm_value","metric","dataset")
    resources_mem<-resources_mem[,c("method","value","metric","dataset","norm_value")]
    all_evals<-rbind(all_evals,resources_mem)
    
    
  }
}

per_categorie_metric<-all_evals%>%
  group_by(method,metric)%>%
  summarize(value=mean(value,na.rm=TRUE),norm_value=mean(norm_value,na.rm=TRUE))
per_categorie_metric$dataset<-"combined"
per_categorie_metric<-per_categorie_metric[,c("method","value",
                                              "metric","dataset","norm_value")]
#all_evals<-rbind(all_evals,per_categorie_metric)

#Order methods based on overall score
method_ordering<-per_categorie_metric%>%
  filter(metric %in% c("pearson","auc_gains","auc_losses","max_f1"))%>%
  group_by(method)%>%
  summarize(norm_value=mean(norm_value))%>%
  arrange(norm_value)

# Create one plot only for the metrics
colnames(method_ordering)<-c("method","value")
method_ordering$metric<-"overall"
method_ordering$dataset<-"combined"
method_ordering$norm_value<-method_ordering$value
scores<-all_evals[all_evals$metric %in% 
                    c("pearson","max_f1","auc_gains","auc_losses"),]
scores<-rbind(scores,method_ordering)
scores$method<-factor(scores$method,levels=c(method_ordering$method,"combined"))
g<-ggplot(scores,aes(y=method,x=dataset,fill=norm_value))+
  geom_tile()+xlab("Data set")+
  scale_fill_viridis("Normalized\nScore",limits=c(0,1))+
  geom_text(aes(label=round(value,2),
                color=ifelse(norm_value<0.6,'white','black')),
            size=2.8)+
  scale_color_manual(values=c("black","white"))+
  facet_grid(~metric,scales = "free",space="free")+
  theme(axis.title.y = element_blank(),
        axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
        legend.position="none")+
  guides(color="none")
ggsave(g,file="~/Desktop/scoring_mean.png",width=10,height=4)

all_evals$method<-factor(all_evals$method,levels=method_ordering$method)
all_evals$dataset<-factor(all_evals$dataset,levels=dataset_names)

# #Dotplot version
# g<-ggplot(all_evals,aes(y=method,x=dataset,color=value,size=value))+
#   geom_point()+
#   scale_color_viridis("Normalized\nScore")+
#   facet_grid(~metric,scales="free")+
#   theme(axis.title.y = element_blank(),
#         axis.text.x=element_text(angle=90))+
#   guides(size="none")

#Split the large overview into two parts (does not fit into one plot)
all_evals_part1<-all_evals[all_evals$metric %in% 
                             c("pearson","max_f1","auc_gains","auc_losses",
                               "runtime_h","max_vms_gb"),]
all_evals_part1$metric<-factor(all_evals_part1$metric,
                               levels=c("pearson","auc_gains","auc_losses","max_f1",
                                        "runtime_h","max_vms_gb"))
g<-ggplot(all_evals_part1,aes(y=method,x=dataset,fill=norm_value))+
  geom_tile()+xlab("Data set")+
  scale_fill_viridis("Normalized\nScore",limits=c(0,1))+
  geom_text(aes(label=round(value,2),
                color=ifelse(norm_value<0.6,'white','black')),
            size=2.8)+
  scale_color_manual(values=c("black","white"))+
  facet_wrap(~metric,ncol=3)+
  theme(axis.title.y = element_blank(),
        axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
        legend.position="none")+
  guides(color="none")
ggsave(g,file="~/Desktop/scoring_overview.png",width=10,height=6)

all_evals_part2<-all_evals[all_evals$metric %in% 
                             c("sens_gains","sens_gains_f1",
                               "prec_gains","prec_gains_f1",
                               "sens_losses","sens_losses_f1",
                               "prec_losses","prec_losses_f1"),]

g<-ggplot(all_evals_part2,aes(y=method,x=dataset,fill=norm_value))+
  geom_tile()+
  scale_fill_viridis("Normalized\nScore",limits=c(0,1))+
  geom_text(aes(label=round(value,2),
                color=ifelse(norm_value<0.6,'white','black')),
            size=2.8)+
  scale_color_manual(values=c("black","white"))+
  facet_wrap(~metric,ncol=4)+
  theme(axis.title.y = element_blank(),
        axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
        legend.position="none")+
  guides(color="none")
ggsave(g,file="~/Desktop/scoring_overview_sens_spec.png",width=12,height=8)

all_evals_part3<-all_evals[all_evals$metric %in% 
                             c("auc_gains","sens_gains_f1","prec_gains_f1",
                               "auc_losses","sens_losses_f1","prec_losses_f1"),]
all_evals_part3$metric<-factor(all_evals_part3$metric,
                               levels=c("auc_gains","sens_gains_f1","prec_gains_f1",
                                        "auc_losses","sens_losses_f1","prec_losses_f1"))
g<-ggplot(all_evals_part3,aes(y=method,x=dataset,fill=norm_value))+
  geom_tile()+
  scale_fill_viridis("Normalized\nScore",limits=c(0,1))+
  geom_text(aes(label=round(value,2),
                color=ifelse(norm_value<0.6,'white','black')),
            size=2.8)+
  scale_color_manual(values=c("black","white"))+
  facet_wrap(~metric,ncol=3)+
  theme(axis.title.y = element_blank(),
        axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
        legend.position="none")+
  guides(color="none")
ggsave(g,file="~/Desktop/scoring_overview_f1.png",width=10,height=8)

#Overall summary score


#Plot chosen optimal thresholds
cutoffs<-melt(f1_optimal_cutoffs,id.vars=c("method","dataset"))
g<-ggplot(cutoffs,aes(x=method,y=value))+
  geom_boxplot()+geom_point(aes(color=dataset))+
  xlab("Method")+ylab("Cutoff")+
  facet_wrap(~variable)+
  theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),
        axis.title.x=element_blank(),
        legend.position="bottom")
ggsave(g,file="~/Desktop/cutoffs_f1.png",width=8,height=5)

