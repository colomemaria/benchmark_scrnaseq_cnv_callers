# ------------------------------------------------------------------------------
# Extract unfiltered distribution of CNVs in the genome
# ------------------------------------------------------------------------------

library(GenomicRanges)
library(data.table)
library(dplyr)

source("scripts/lib.R")

dataset_names<-c("SNU601","NCIN87","MKN45","KATOIII",
                 "NUGC4","SNU638","SNU668","HGC27", "SNU16",
                 "MCF7","iAMP21",
                 "COLO320","MM","BCC06","BCC06post",
                  "mouse","HCT116","A375","ALL1","ALL2")
wgs_type<-setNames(c("WGS","WGS","WGS","WGS",
                           "WGS","WGS","WGS","WGS", "WGS",
                           "WGS","WGS",
                           "CNVKIT","GATK","GATK","GATK",
                     "WGS","WGS","WGS","WGS", "WGS"),
                         dataset_names)


dist_genomic<-NULL
for(ds in dataset_names){
  
  #Process WGS data
  if(wgs_type[ds]=="WGS"){
    file_path<-paste0("data/",ds,"_WGS_groundtruth/wgs_results_formated.csv")
    wgs<-read_wgs_results(file_path)
    
    #Discretize values
    wgs$wgs_discrete<-ifelse(wgs$wgs_mean>2.5,3,
                             ifelse(wgs$wgs_mean<1.5,1,2))
    
    dist_genomic<-rbind(dist_genomic,
                        data.frame(dataset=ds,
                                   type="WGS",
                                   gains=sum(wgs$wgs_discrete==3, na.rm=TRUE),
                                   losses=sum(wgs$wgs_discrete==1, na.rm=TRUE),
                                   total=length(wgs)))
    
  } else if (wgs_type[ds]=="GATK"){
    file_path<-paste0("data/",ds,"_WES_GATK/",ds,"_tumor_clean.called.seg")
    wgs<-read_wes_results_gatk(file_path)
    
    #Discretize values
    wgs$wgs_discrete<-ifelse(wgs$wgs_mean>2.5,3,
                             ifelse(wgs$wgs_mean<1.5,1,2))
    
    dist_genomic<-rbind(dist_genomic,
                        data.frame(dataset=ds,
                                   type="WES",
                                   gains=sum(wgs$wgs_discrete==3),
                                   losses=sum(wgs$wgs_discrete==1),
                                   total=length(wgs)))
  } else if (wgs_type[ds]=="CNVKIT"){
    file_path<-"data/COLO320_WGS/COLO320_HSR_WGS.aln.sort.rmdup.call.cns"
    wgs<-read_results_cnvkit_segs(file_path)
    
    #Discretize values
    wgs$wgs_discrete<-ifelse(wgs$wgs_mean>2.5,3,
                             ifelse(wgs$wgs_mean<1.5,1,2))
    
    dist_genomic<-rbind(dist_genomic,
                        data.frame(dataset=ds,
                                   type="WGS",
                                   gains=sum(wgs$wgs_discrete==3),
                                   losses=sum(wgs$wgs_discrete==1),
                                   total=length(wgs)))
  }
}

#Save the complete results somewhere
write.table(dist_genomic,file="data/genomice_CNV_distribution_unfiltered.tsv",
            sep="\t",row.names=FALSE, quote=FALSE)


