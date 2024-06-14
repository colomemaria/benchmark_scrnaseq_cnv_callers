# ------------------------------------------------------------------------------
# Evaluate Casper results with stricter filtering
# ------------------------------------------------------------------------------

library(data.table)
library(GenomicRanges)

source("scripts/lib.R")

#Evaluation results from all other methods
combined_output<-fread("results/output_SNU601/evaluation/outputs_allmethods_combined.tsv")
combined_output<-makeGRangesFromDataFrame(combined_output,keep.extra.columns = TRUE)

#Add filtered CaSpER results
casper_filtered<-readRDS("results/output_SNU601/casper_strict/SNU601_casper_pseudobulk_aggregate.RDS")
length(casper_filtered)
casper_filtered$casper_filtered<-with(casper_filtered,mean_loss+(mean_base*2)+mean_gain*3)

#Add it to the combined file (check file before and afterwards)
length(combined_output)
combined_output<-combine_range_objects(combined_output,casper_filtered,
                                       method_colname="casper_filtered")
length(combined_output)

cor(combined_output$wgs_mean,combined_output$casper)
cor(combined_output$wgs_mean,combined_output$casper_filtered)

#Evaluate Casper and Casper (filtered)
eval_1<-get_all_metrics(combined_output,"casper_filtered")
eval_2<-get_all_metrics(combined_output,"casper")
all_res<-rbind(eval_1,eval_2)
write.table(all_res,quote=FALSE,sep="\t",
            file="results/output_SNU601/evaluation/eval_results_filtered_casper.tsv")

