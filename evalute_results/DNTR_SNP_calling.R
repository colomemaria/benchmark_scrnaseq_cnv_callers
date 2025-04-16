# ------------------------------------------------------------------------------
# Check how many SNPs for HCT116 and A375
# ------------------------------------------------------------------------------

library(data.table)
library(ggplot2)

theme_set(theme_bw())

#Get previous SNP results
snp_calls<-fread("results/snp_calling_summary.tsv")
snp_calls<-snp_calls[! snp_calls$dataset %in% c("SNU601_sample20","SNU601_sample40",
                                                "SNU601_sample60","SNU601_sample80"),]

# ------------------------------------------------------------------------------
# Add Numbat and CaSpER results
# ------------------------------------------------------------------------------

#Get SNPs for HCT116
ds<-"HCT116"
numbat_afs<-fread(paste0("results/output_",ds,"/numbat_smartseq/",ds,"_allele_counts.tsv.gz"))


casper_af<-fread(paste0("results/output_",ds,"/casper_smartseq/",ds,
                        "_BAFExtract/",ds,"_BAFExtract.af"))

snp_calls<-rbind(snp_calls,
                   data.frame(dataset=ds,
                              numbat_nsnps=nrow(numbat_afs),
                              numbat_cells_wsnp=length(unique(numbat_afs$cell)),
                              numbat_mean_DP=mean(numbat_afs$DP),
                              casper_nsnps = nrow(casper_af),
                              casper_mean_DP=mean(casper_af$V6)))

#Get SNPs for A375
ds<-"A375"
numbat_afs<-fread(paste0("results/output_",ds,"/numbat_smartseq/",ds,"_allele_counts.tsv.gz"))

casper_af<-fread(paste0("results/output_",ds,"/casper_smartseq/",ds,
                        "_BAFExtract/",ds,"_BAFExtract.af"))

snp_calls<-rbind(snp_calls,
                 data.frame(dataset=ds,
                            numbat_nsnps=nrow(numbat_afs),
                            numbat_cells_wsnp=length(unique(numbat_afs$cell)),
                            numbat_mean_DP=mean(numbat_afs$DP),
                            casper_nsnps = nrow(casper_af),
                            casper_mean_DP=mean(casper_af$V6)))


snp_calls$type<-ifelse(snp_calls$dataset %in% c("HCT116","A375"),
                       snp_calls$dataset,"10X")

g<-ggplot(snp_calls,aes(x=numbat_nsnps,y=casper_nsnps,color=type))+
  geom_point()+scale_x_log10()+scale_y_log10()+
  xlab("Number SNPs Numbat")+ylab("Number SNPs CaSpER")
ggsave(g,file="number_SNP_Smartseq2.png",
       width=6,height=5)

g<-ggplot(snp_calls,aes(x=numbat_mean_DP,y=casper_mean_DP,color=type))+
  geom_point()+scale_x_log10()+scale_y_log10()+
  xlab("Mean coverage Numbat (per cell)")+ylab("Mean coverage CaSpER (pseudobulk)")
ggsave(g,file="coverage_SNP_Smartseq2.png",
       width=6,height=5)
