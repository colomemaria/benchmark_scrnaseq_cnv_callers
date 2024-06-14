# ------------------------------------------------------------------------------
# Evaluate the prediction of cancer cells when not specifying the reference
# explicetly (only done for a subset of the methods)
# ------------------------------------------------------------------------------

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

# ------------------------------------------------------------------------------
print("Load libraries")
# ------------------------------------------------------------------------------

library(data.table)
library(ggplot2)
library(ggpubr)
library(numbat)
theme_set(theme_bw())

# ------------------------------------------------------------------------------
print("Get input parameters from snakemake")
# ------------------------------------------------------------------------------

input_annot<-snakemake@input$annot_truth
input_ref_groups<-snakemake@input$ref_groups
  
input_copykat<-snakemake@input$copykat_pred
input_scevan<-snakemake@input$scevan_pred
input_numbat<-snakemake@input$numbat_pred
  
output_file<-snakemake@output$text
output_plot<-snakemake@output$plot
  
# ------------------------------------------------------------------------------
print("Evaluation of the cancer cell prediction of the different methods.")
# ------------------------------------------------------------------------------

#Read annotation input files
annot<-fread(input_annot,header=FALSE)
ref_groups<-fread(input_ref_groups)
cancer_cells<-annot$V1[! annot$V2 %in% ref_groups$ref_groups]

general_performance<-NULL

#Evaluate copykat input
pred_copykat <-fread(input_copykat)
pred_copykat$truth<-ifelse(pred_copykat$cell.names %in% cancer_cells,"aneuploid","diploid")

print("Evaluation for copykat")
print(table(pred_copykat$copykat.pred,pred_copykat$truth))

general_performance<-rbind(general_performance,
                           data.frame(
                             method="copyKat",
                             acc_general=mean(pred_copykat$copykat.pred==pred_copykat$truth),
                             acc_woNA=mean((pred_copykat$copykat.pred==pred_copykat$truth)
                                           [!pred_copykat$copykat.pred=="not.defined"])
                           ))
                           
#Evaluate scevan input
suppressWarnings(pred_scevan <-fread(input_scevan))
pred_scevan$truth<-ifelse(pred_scevan$V1 %in% cancer_cells,"tumor","normal")

print("Evaluation for SCEVAN")
print(table(pred_scevan$class,pred_scevan$truth))

general_performance<-rbind(general_performance,
                           data.frame(
                             method="SCEVAN",
                             acc_general=mean(pred_scevan$class==pred_scevan$truth),
                             acc_woNA=mean((pred_scevan$class==pred_scevan$truth)
                                           [!pred_scevan$class=="filtered"])
                           ))

#Evaluate numbat input 
numbat_obj<-Numbat$new(dirname(input_numbat))
pred_numbat<-numbat_obj$clone_post[,c("cell","compartment_opt")]
pred_numbat$truth<-ifelse(pred_numbat$cell %in% cancer_cells,"tumor","normal")

print("Evaluation for numbat")
print(table(pred_numbat$compartment_opt,pred_numbat$truth))

general_performance<-rbind(general_performance,
                           data.frame(
                             method="Numbat",
                             acc_general=mean(pred_numbat$compartment_opt==pred_numbat$truth),
                             acc_woNA=mean(pred_numbat$compartment_opt==pred_numbat$truth)
                           ))

#Save text file with accuracy results
write.table(general_performance,file=output_file,sep="\t",
            row.names=FALSE,quote=FALSE)

#Save plot with accuracy
g.1<-ggplot(general_performance,aes(x=method,fill=method,y=acc_general))+
  geom_bar(stat="identity")+ylim(0,1)+
  theme(legend.position="none",
        axis.text.x = element_text(angle=45,vjust=1,hjust=1))+
  ylab("Accuracy (including filtered cells)")+xlab("CNV Method")

g.2<-ggplot(general_performance,aes(x=method,fill=method,y=acc_woNA))+
  geom_bar(stat="identity")+ylim(0,1)+
  theme(legend.position="none",
        axis.text.x = element_text(angle=45,vjust=1,hjust=1))+
  ylab("Accuracy (wo filtered cells)")+xlab("CNV Method")

g<-ggarrange(g.1,g.2,ncol=2)
ggsave(g,file=output_plot,width=6,height=4)

