# ------------------------------------------------------------------------------
# Evaluate for each method the runtime and memory requirements
# by plotting two bar plots
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
theme_set(theme_bw())

# ------------------------------------------------------------------------------
print("Get input parameters from snakemake")
# ------------------------------------------------------------------------------

all_files<-snakemake@input

output_plot<-snakemake@output$plot
output_txt<-snakemake@output$file

# ------------------------------------------------------------------------------
print("Evaluate runtime and memory of all specified methods")
# ------------------------------------------------------------------------------

#Iterate over all files in the benchmarking directory
benchmark_combined<-NULL
for(filename in all_files){
  res_benchmark<-fread(filename)
  res_benchmark$method<-gsub("_[a-zA-Z0-9]+_benchmark.txt","",basename(filename))
  
  benchmark_combined<-rbind(benchmark_combined,
                            res_benchmark)
}

#Specific case for SNU601 due to MCF7_SNU601 being also a test dataset,
#both version are mixed for SNU601 and the MCF7 files need to be filtered out
benchmark_combined<-benchmark_combined[! endsWith(benchmark_combined$method,"_MCF7") &
                                         ! endsWith(benchmark_combined$method,"_MM") &
                                         ! endsWith(benchmark_combined$method,"_BCC06") &
                                         ! endsWith(benchmark_combined$method,"_BCC08"),]

#Rename the second column name (double point makes trouble)
colnames(benchmark_combined)[2]<-"h_m_s"

#Combine Honeybadger results (if they exist both)
h_1<-benchmark_combined[benchmark_combined$method == "honeybadger_expr_part1",]
h_2<-benchmark_combined[benchmark_combined$method == "honeybadger_expr_part2",]
benchmark_combined<-benchmark_combined[! benchmark_combined$method %in%
                                         c("honeybadger_expr_part1",
                                           "honeybadger_expr_part2"),]
if(nrow(h_1) > 0 & nrow(h_2) > 0){
  
  benchmark_combined<-rbind(benchmark_combined,
                            data.frame(s=(h_1$s + h_2$s),
                                       h_m_s="",
                                       max_rss=max(h_1$max_rss,h_2$max_rss),
                                       max_vms=max(h_1$max_vms,h_2$max_vms),
                                       max_uss=max(h_1$max_uss,h_2$max_uss),
                                       max_pss=max(h_1$max_pss,h_2$max_pss),
                                       io_in=NA,
                                       io_out=NA,
                                       mean_load=NA,
                                       cpu_time=NA,
                                       method="honeybadger_expr"))
}

#Combine Numbat results (if they exist both)
h_1<-benchmark_combined[benchmark_combined$method == "numbat_preprocessing",]
h_2<-benchmark_combined[benchmark_combined$method == "numbat",]
h_3<-benchmark_combined[benchmark_combined$method == "numbat_predict",]
benchmark_combined<-benchmark_combined[! benchmark_combined$method %in%
                                         c("numbat_preprocessing",
                                           "numbat","numbat_predict"),]
if(nrow(h_1) > 0 & nrow(h_2) > 0){
  
  benchmark_combined<-rbind(benchmark_combined,
                            data.frame(s=(h_1$s + h_2$s),
                                       h_m_s="",
                                       max_rss=max(h_1$max_rss,h_2$max_rss),
                                       max_vms=max(h_1$max_vms,h_2$max_vms),
                                       max_uss=max(h_1$max_uss,h_2$max_uss),
                                       max_pss=max(h_1$max_pss,h_2$max_pss),
                                       io_in=NA,
                                       io_out=NA,
                                       mean_load=NA,
                                       cpu_time=NA,
                                       method="numbat"))
  
  #In case Numbat was run also with the "predict cancer cell" option
  if(nrow(h_3)>0){
    
    benchmark_combined<-rbind(benchmark_combined,
                              data.frame(s=(h_1$s + h_3$s),
                                         h_m_s="",
                                         max_rss=max(h_1$max_rss,h_3$max_rss),
                                         max_vms=max(h_1$max_vms,h_3$max_vms),
                                         max_uss=max(h_1$max_uss,h_3$max_uss),
                                         max_pss=max(h_1$max_pss,h_3$max_pss),
                                         io_in=NA,
                                         io_out=NA,
                                         mean_load=NA,
                                         cpu_time=NA,
                                         method="numbat_predict"))
    
  }
}


#Combine casper results (if they exist all three)
h_1<-benchmark_combined[benchmark_combined$method == "casper_preprocess",]
h_2<-benchmark_combined[benchmark_combined$method == "casper_gene_annot",]
h_3<-benchmark_combined[benchmark_combined$method == "casper",]
benchmark_combined<-benchmark_combined[! benchmark_combined$method %in%
                                         c("casper_preprocess",
                                           "casper_gene_annot",
                                           "casper"),]

if(nrow(h_1) > 0 & nrow(h_2) > 0 & nrow(h_3) > 0){
  
  benchmark_combined<-rbind(benchmark_combined,
                            data.frame(s=(h_1$s + h_2$s + h_3$s),
                                       h_m_s="",
                                       max_rss=max(h_1$max_rss,h_2$max_rss,h_3$max_rss),
                                       max_vms=max(h_1$max_vms,h_2$max_vms,h_3$max_vms),
                                       max_uss=max(h_1$max_uss,h_2$max_uss,h_3$max_uss),
                                       max_pss=max(h_1$max_pss,h_2$max_pss,h_3$max_pss),
                                       io_in=NA,
                                       io_out=NA,
                                       mean_load=NA,
                                       cpu_time=NA,
                                       method="casper"))
}

#Combine CONICSmat results (if they exist both)
h_1<-benchmark_combined[benchmark_combined$method == "CONICSmat_import",]
h_2<-benchmark_combined[benchmark_combined$method == "CONICSmat",]
h_3<-benchmark_combined[benchmark_combined$method == "CONICSmat_predict",]
benchmark_combined<-benchmark_combined[! benchmark_combined$method %in%
                                         c("CONICSmat_import",
                                           "CONICSmat"),]
if(nrow(h_1) > 0 & nrow(h_2) > 0){
  
  benchmark_combined<-rbind(benchmark_combined,
                            data.frame(s=(h_1$s + h_2$s),
                                       h_m_s="",
                                       max_rss=max(h_1$max_rss,h_2$max_rss),
                                       max_vms=max(h_1$max_vms,h_2$max_vms),
                                       max_uss=max(h_1$max_uss,h_2$max_uss),
                                       max_pss=max(h_1$max_pss,h_2$max_pss),
                                       io_in=NA,
                                       io_out=NA,
                                       mean_load=NA,
                                       cpu_time=NA,
                                       method="CONICSmat"))
  
  #In case CONICSmat was run also with the "predict cancer cell" option
  if(nrow(h_3)>0){
    
    benchmark_combined<-rbind(benchmark_combined,
                              data.frame(s=(h_1$s + h_3$s),
                                         h_m_s="",
                                         max_rss=max(h_1$max_rss,h_3$max_rss),
                                         max_vms=max(h_1$max_vms,h_3$max_vms),
                                         max_uss=max(h_1$max_uss,h_3$max_uss),
                                         max_pss=max(h_1$max_pss,h_3$max_pss),
                                         io_in=NA,
                                         io_out=NA,
                                         mean_load=NA,
                                         cpu_time=NA,
                                         method="CONICSmat_predict"))
    
  }
}

#Convert runtime in seconds to hours
benchmark_combined$runtime_h<-benchmark_combined$s / 3600

#Vector for renaming methods (official published names) and specifying
#the order in the plot
method_names<-setNames(c("HoneyBADGER","CaSpER",
                         "Numbat","Numbat (wo ref)","InferCNV",
                         "CONICSmat", "CONICSmat (wo ref)",
                         "copyKat","copyKat (wo ref)",
                         "SCEVAN","SCEVAN (wo ref)",
                         "SCEVAN (subclones)"),
                       c("honeybadger_expr","casper",
                         "numbat","numbat_predict","infercnv",
                         "CONICSmat", "CONICSmat_predict",
                         "copykat","copykat_predict",
                         "scevan","scevan_predict",
                         "scevan_clones"))

benchmark_combined$method<-method_names[benchmark_combined$method]
benchmark_combined$method<-factor(benchmark_combined$method,
                                  levels=method_names)
  
#Create plots of runtime and memory requirements
g.1<-ggplot(benchmark_combined,aes(x=method,fill=method,y=runtime_h))+
  geom_bar(stat="identity")+
  theme(legend.position="none",
        axis.text.x = element_text(angle=45,vjust=1,hjust=1))+
  ylab("Runtime (in h)")+xlab("CNV Method")

g.2<-ggplot(benchmark_combined,aes(x=method,fill=method,y=max_vms))+
  geom_bar(stat="identity")+
  theme(legend.position="none",
        axis.text.x = element_text(angle=45,vjust=1,hjust=1))+
  ylab("Memory (max vms) (in MB)")+xlab("CNV Method")

g<-ggarrange(g.1,g.2,ncol=2)
ggsave(g,file=output_plot,width=10,height=4)

#Save the combined runtime and benchmarking also as a combined tsv file
write.table(benchmark_combined,file=output_txt,
            sep="\t",quote=FALSE,row.names=FALSE)

