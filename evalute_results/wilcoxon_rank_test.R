# ------------------------------------------------------------------------------
# Estimate significant differences for the tested methods across datasets
# According to googling: using Wilcoxon signed rank test (because the data 
# is paired) rather than Wilcoxon rank sum test
# ------------------------------------------------------------------------------

library(data.table)
library(dplyr)
library(ggplot2)

theme_set(theme_bw())

# ------------------------------------------------------------------------------
#Load performance values
# ------------------------------------------------------------------------------

#Load the different result files
dataset_names<-c("SNU601","NCIN87","MKN45","KATOIII",
                 "NUGC4","SNU638","SNU668","HGC27", "SNU16",
                 "MCF7","iAMP21",
                 "COLO320","MM","BCC06","BCC06post")
dataset_ending<-setNames(c("","","","",
                           "","","","", "",
                           "","",
                           "_cnvkit","_wes","_wes","_wes"),
                         dataset_names)


#Vector for renaming methods (official published names)
method_names<-setNames(c("InferCNV (CNV)","InferCNV (Expr)","CaSpER",
                         "copyKat","SCEVAN (CNV)","SCEVAN (CNV)","SCEVAN (Expr)",
                         "Numbat (Expr)", "Numbat (CNV)", "CONICSmat"),
                       c("infercnv_cnv","infercnv_expr","casper",
                         "copykat","scevan_vega","scevan_cnv","scevan",
                         "numbat","numbat_cnv", "CONICSmat"))

#Load the files
all_evals<-NULL
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
  
  #Get standard AUC values
  auc_basic<-melt(auc[,c("method","auc_gains","auc_losses")],id.vars="method")
  auc_basic<-auc_basic[,c("method","value","variable")]
  colnames(auc_basic)[3]<-"metric"
  auc_basic$dataset<-dataset
  
  #Scale AUC to also range from 0-1 (set value < 0.5 to 0)
  auc_basic$norm_value<-auc_basic$value
  auc_basic$norm_value<-(auc_basic$norm_value-0.5)*2
  auc_basic$norm_value[auc_basic$norm_value<0]<-0
  
  all_evals<-rbind(all_evals,auc_basic)
  
  #Get AUPRC and truncated AUC values
  auc_others<-melt(auc[,c("method","auc_gains_trunc","aucpr_gains",
                          "auc_losses_trunc","aucpr_losses")],id.vars="method")
  auc_others<-auc_others[,c("method","value","variable")]
  colnames(auc_others)[3]<-"metric"
  auc_others$dataset<-dataset
  auc_others$norm_value<-auc_others$value
  all_evals<-rbind(all_evals,auc_others) 
  
  #Add sens and spec values for F1 scores
  filename<-paste0("results/output_",dataset,
                   "/evaluation/evaluation_cnv_prediction_f1",
                   dataset_ending[dataset],".txt")
  sensspecf1<-fread(filename)
  sensspecf1$method<-method_names[sensspecf1$method]
  
  colnames(sensspecf1)[c(4,5,7,8)]<-paste0(colnames(sensspecf1)[c(4,5,7,8)],
                                           "_f1")
  sensspecf1<-melt(sensspecf1[,!c("cutoff_f1_gain","cutoff_f1_loss")],id.vars="method")
  sensspecf1<-sensspecf1[,c("method","value","variable")]
  colnames(sensspecf1)[3]<-"metric"
  sensspecf1$dataset<-dataset
  sensspecf1$norm_value<-sensspecf1$value
  all_evals<-rbind(all_evals,sensspecf1)
  
}

# ------------------------------------------------------------------------------
# Test for statistical significance
# ------------------------------------------------------------------------------

#Get the ordering as in Figure 2a
method_order<-all_evals%>%
  filter(metric=="max_f1")%>%
  group_by(method)%>%
  summarize(avg_value=mean(value))

method_order<-method_order[order(method_order$avg_value,decreasing=TRUE),]

#Test each method against each other
res_combined<-NULL
for(metric_type in c("max_f1","pearson","auc_gains_trunc","auc_losses_trunc")){
  
  res_frame<-NULL
  for(m1 in 1:(nrow(method_order))){
    for(m2 in (m1):nrow(method_order)){
      
        vals1<-all_evals$value[all_evals$metric==metric_type & all_evals$method==method_order$method[m1]]
        vals2<-all_evals$value[all_evals$metric==metric_type & all_evals$method==method_order$method[m2]]
        wtest<-wilcox.test(vals1,vals2,alternative="two.sided",paired=TRUE)
        
        res_frame<-rbind(res_frame,data.frame(method1=method_order$method[m1],
                                              method2=method_order$method[m2],
                                              metric = metric_type,
                                              statistics=wtest$statistic,
                                              p_val=wtest$p.value))
    }
  }
  
  #Multiple testing adjustment
  res_frame$p_adjust_BH<-p.adjust(res_frame$p_val,method="BH")
  res_frame$signif<-res_frame$p_adjust_BH<0.05
  
  res_combined<-rbind(res_combined, res_frame)
}

res_combined$method1<-factor(res_combined$method1,levels=method_order$method)
res_combined$method2<-factor(res_combined$method2,levels=rev(method_order$method))


#Rename the metrics
rename_metric<-setNames(c("Maximal F1 Score","Correlation",
                          "Partial AUC (gain)","Partial AUC (loss)"),
                        c("max_f1","pearson","auc_gains_trunc","auc_losses_trunc"))
res_combined$metric<-rename_metric[res_combined$metric]
res_combined$metric<-factor(res_combined$metric,
                            levels=rename_metric)
  
#Set NA values to 1
res_combined$signif[is.na(res_combined$signif)]<-FALSE

g<-ggplot(res_combined,aes(x=method1,y=method2,fill=signif))+
  geom_tile()+
  scale_fill_manual("FDR < 0.05",values=c("orange","lightgreen"))+
  geom_text(aes(label=round(p_adjust_BH,3)),
           size=2.8)+
  facet_wrap(~metric,ncol=2)+
  theme(axis.text.x=element_text(angle=60,vjust=1,hjust=1),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position = "bottom")
  
ggsave(g,file=paste0("../figure_plots/supp_figure_significant_differences_combined.png"),
       width=8,height=6)
