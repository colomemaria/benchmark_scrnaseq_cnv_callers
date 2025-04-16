# ------------------------------------------------------------------------------
# Check Aneufinder results
# For HCT116 and A375 filter on cells with sufficient coverage
# ------------------------------------------------------------------------------

library(AneuFinder)
library(epiAneufinder)

# ------------------------------------------------------------------------------
# A375
# ------------------------------------------------------------------------------

rdat_files<-list.files("A375_DNTRseq/aneufinder_results/MODELS/method-edivisive/",full.names=TRUE)

model_list<-list()
qual_df<-NULL
for(rfile in rdat_files){
  load(rfile) #model
  model_list<-c(model_list,list(model))
  
  #Save the quality info separate
  qual<-as.data.frame(model$qualityInfo)
  qual$barcode<-model$ID
  
  #Get the total copy number
  qual$mean_copy_number<-mean(model$bins$copy.number)
  
  qual_df<-rbind(qual_df,qual)
}

#Identified 
qual_df$strange_profile<-(qual_df$barcode %in% c("HCA00102.E14_dedup.bam","HCA00102.E16_dedup.bam",
                                              "HCA00102.C13_dedup.bam","HCA00102.M13_dedup.bam"))

ggplot(qual_df,aes(x=total.read.count,y=mean_copy_number,color=strange_profile))+
  geom_point()

#Remove cells with very few reads or odd profiles (huge number of gains, probably wrong baseline)
qual_df<-qual_df[qual_df$total.read.count>1e5,]
qual_df<-qual_df[qual_df$mean_copy_number<4.5,]

#Filter for the cells where we have also RNA data
rna_samples<-fread("A375_DNTRseq/input_A375/matching_sample_names.txt",header=FALSE)

qual_df$sample_name<-gsub("_dedup.bam","",qual_df$barcode)
qual_df<-qual_df[qual_df$sample_name %in% rna_samples$V1,]


#Merge all filtered cells to one grange
complete_grange<-NULL
for(i in 1:length(model_list)){
  
  #Check whether the cell is in the filtered set
  if(model_list[[i]]$ID %in% qual_df$barcode){
    
    #For the first cell, save the coordinates
    if(is.null(complete_grange)){
      complete_grange<-model_list[[i]]$bins
      elementMetadata(complete_grange)<-NULL
      
    }
    
    #Simplify CNV states to 1,2,3 (loss,base,gain)
    cnvs<-model_list[[i]]$bins$copy.number
    elementMetadata(complete_grange)[gsub("_dedup.bam","",model_list[[i]]$ID)]<-ifelse(cnvs<=1,1,ifelse(cnvs>3,3,cnvs))
  }
}


#Filter for the regions chosen in the pseudobulk CNV
pb_res<-fread("A375_DNTRseq/aneufinder_results/wgs_results_formated.csv")
pb_res$chr<-paste0("chr",pb_res$chr)
pb_grange<-makeGRangesFromDataFrame(pb_res)
ovs<-as.data.frame(findOverlaps(complete_grange,pb_grange))

complete_grange<-complete_grange[ovs$queryHits,]

#Save the results
complete_df<-as.data.frame(complete_grange)
complete_df$width<-NULL
complete_df$strand<-NULL

write.table(complete_df,file="A375_DNTRseq/aneufinder_results/cnvs_per_cell_filtered.tsv",
            quote=FALSE,sep="\t",row.names=FALSE)

#Tweak the results to use epiAneufinder for plotting
colnames(complete_df)[1]<-"seq"
colnames(complete_df)[4:ncol(complete_df)]<-paste0("cell-",colnames(complete_df)[4:ncol(complete_df)])
complete_df[,4:ncol(complete_df)]<-complete_df[,4:ncol(complete_df)]-1
epiAneufinder::plot_karyo_annotated(complete_df,plot_path="A375_DNTRseq/aneufinder_results/karyogram_filtered_cells.png")

# ------------------------------------------------------------------------------
# A375
# ------------------------------------------------------------------------------

rdat_files<-list.files("HCT116_DNTRseq/aneufinder_results/MODELS/method-edivisive/",full.names=TRUE)

model_list<-list()
qual_df<-NULL
for(rfile in rdat_files){
  load(rfile) #model
  model_list<-c(model_list,list(model))
  
  #Save the quality info separate
  qual<-as.data.frame(model$qualityInfo)
  qual$barcode<-model$ID
  
  #Get the total copy number
  qual$mean_copy_number<-mean(model$bins$copy.number)
  
  qual_df<-rbind(qual_df,qual)
}


#Remove one cell with 0 coverage
qual_df<-qual_df[qual_df$total.read.count>0,]

ggplot(qual_df,aes(x=total.read.count,y=mean_copy_number))+
  geom_point()

#Remove cells with very few reads or odd profiles (huge number of gains, probably wrong baseline)
qual_df<-qual_df[qual_df$total.read.count>1e5,]
qual_df<-qual_df[qual_df$mean_copy_number<2.2 & qual_df$mean_copy_number>1.7,]

#Filter for the cells where we have also RNA data
rna_samples<-fread("HCT116_DNTRseq/input_HCT116/matching_sample_names.txt",header=FALSE)

qual_df$sample_name<-gsub("_dedup.bam","",qual_df$barcode)
qual_df<-qual_df[qual_df$sample_name %in% rna_samples$V1,]

#Merge all filtered cells to one grange
complete_grange<-NULL
for(i in 1:length(model_list)){
  
  #Check whether the cell is in the filtered set
  if(model_list[[i]]$ID %in% qual_df$barcode){
    
    #For the first cell, save the coordinates
    if(is.null(complete_grange)){
      complete_grange<-model_list[[i]]$bins
      elementMetadata(complete_grange)<-NULL
      
    }
    
    #Simplify CNV states to 1,2,3 (loss,base,gain)
    cnvs<-model_list[[i]]$bins$copy.number
    elementMetadata(complete_grange)[gsub("_dedup.bam","",model_list[[i]]$ID)]<-ifelse(cnvs<=1,1,ifelse(cnvs>3,3,cnvs))
  }
}


#Filter for the regions chosen in the pseudobulk CNV
pb_res<-fread("A375_DNTRseq/aneufinder_results/wgs_results_formated.csv")
pb_res$chr<-paste0("chr",pb_res$chr)
pb_grange<-makeGRangesFromDataFrame(pb_res)
ovs<-as.data.frame(findOverlaps(complete_grange,pb_grange))

complete_grange<-complete_grange[ovs$queryHits,]

#Save the results
complete_df<-as.data.frame(complete_grange)
complete_df$width<-NULL
complete_df$strand<-NULL

write.table(complete_df,file="HCT116_DNTRseq/aneufinder_results/cnvs_per_cell_filtered.tsv",
            quote=FALSE,sep="\t",row.names=FALSE)

#Tweak the results to use epiAneufinder for plotting
colnames(complete_df)[1]<-"seq"
colnames(complete_df)[4:ncol(complete_df)]<-paste0("cell-",colnames(complete_df)[4:ncol(complete_df)])
complete_df[,4:ncol(complete_df)]<-complete_df[,4:ncol(complete_df)]-1
epiAneufinder::plot_karyo_annotated(complete_df,plot_path="HCT116_DNTRseq/aneufinder_results/karyogram_filtered_cells.png")





