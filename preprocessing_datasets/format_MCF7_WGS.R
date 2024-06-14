# ------------------------------------------------------------------------------
# Process MCF7 WGS results (from ichorCNA) to match input requirements for 
# snakemake evaluation pipeline
# ------------------------------------------------------------------------------

library(fread)
library(GenomicRanges)

#Read CNV predictions per bin (1MB bins)
segs<-fread("../../../results/MCF7_AA_WGS/MCF7_AA.cna.seg")

#Create the file in the same format as the scWGS file
segs_formated<-segs[,c("chr","start","end")]
segs_formated$gain_wgs<-ifelse(segs$MCF7_AA.copy.number>2,1,0)
segs_formated$loss_wgs<-ifelse(segs$MCF7_AA.copy.number<2,1,0)
segs_formated$base_wgs<-ifelse(segs$MCF7_AA.copy.number==2,1,0)

#Split into a binsize of 100,000 instead of 1,000,000
segs_grange<-makeGRangesFromDataFrame(segs_formated,keep.extra.columns = TRUE)
segs_grange_bins<-tile(segs_grange,width=100000)
segs_grange_bins<-unlist(segs_grange_bins)

#Add meta data (=CNV status) again
ovs<-as.data.frame(findOverlaps(segs_grange_bins,segs_grange))
elementMetadata(segs_grange_bins)<-elementMetadata(segs_grange)[ovs$subjectHits,]

#Format it again as a data frame
segs_grange_bins<-as.data.frame(segs_grange_bins)
segs_grange_bins$width<-NULL
segs_grange_bins$strand<-NULL
colnames(segs_grange_bins)[1]<-c("chr")

write.table(segs_grange_bins,file="../../../results/MCF7_AA_WGS/wgs_results_formated.csv",
            quote=FALSE,sep=",")
