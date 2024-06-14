# ------------------------------------------------------------------------------
# Check how many SNPs were found in each dataset for Numbat and Casper
# ------------------------------------------------------------------------------

library(data.table)

dataset_names<-c("SNU601","SNU601_sample20","SNU601_sample40","SNU601_sample60",
                 "SNU601_sample80","NCIN87","MKN45","KATOIII",
                 "NUGC4","SNU638","SNU16","SNU668","HGC27",
                 "MCF7","COLO320","MM","BCC06","BCC06post")


#Get Numbat SNPs
snp_calling<-NULL
for(ds in dataset_names){
  
  #Get Numbat SNPs
  numbat_afs<-fread(paste0("results/output_",ds,"/numbat/",ds,"_allele_counts.tsv.gz"))
  
  #Get Casper SNPs
  casper_af<-fread(paste0("results/output_",ds,"/casper/",ds,
                          "_BAFExtract/",ds,"_BAFExtract.af"))
  
  snp_calling<-rbind(snp_calling,
                     data.frame(dataset=ds,
                                numbat_nsnps=nrow(numbat_afs),
                                numbat_cells_wsnp=length(unique(numbat_afs$cell)),
                                numbat_mean_DP=mean(numbat_afs$DP),
                                casper_nsnps = nrow(casper_af),
                                casper_mean_DP=mean(casper_af$V6)))
}

write.table(snp_calling,sep="\t",quote=FALSE,row.names=FALSE,
            file="snp_calling_summary.tsv")



