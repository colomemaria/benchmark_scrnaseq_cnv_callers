# ------------------------------------------------------------------------------
# For better explanation, visualize the difference between the AUC
# and the truncated AUC
#
# Take as an example MCF7, Numbat_CNV
# ------------------------------------------------------------------------------

library(data.table)
library(pROC)
library(ggplot2)

theme_set(theme_bw())

#Load the results
combined_methods<-fread("results/output_MCF7/evaluation/outputs_allmethods_combined.tsv")

# ------------------------------------------------------------------------------
# One plot for gain truncated AUC values
# ------------------------------------------------------------------------------

gain_scores<-combined_methods$copykat
gain_scores<-gain_scores-2
truth <- ifelse(combined_methods$wgs_mean>2.5,1,0)

gain_scores_threshold<-ifelse(gain_scores>0,1,0)
sens_threshold<-sum(gain_scores_threshold==1 & truth==1)/sum(truth==1)
  
roc_gain <- roc(response=truth, predictor=gain_scores)

plot_roc<-data.frame(sens=roc_gain$sensitivities,
                     fpr=1-roc_gain$specificities,
                     threshold=roc_gain$thresholds)

#Reorder them so that the line plot is generated correctly
plot_roc<-plot_roc[order(plot_roc$fpr, plot_roc$sens), ]


g<-ggplot(plot_roc,aes(x=fpr,y=sens,color=threshold))+
  geom_abline(linetype = "dashed", color = "gray")+
  geom_line(linewidth=1.2)+
  geom_hline(yintercept=sens_threshold)+
  scale_color_gradient("Copykat\ngradient",low = "blue", high = "red") +
  xlab("False positive rate")+
  ylab("Sensitvity")

ggsave(g,file="roc_gains_copykat_MCF7.pdf",width=5,height=5)

