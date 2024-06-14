for(comp in c("epithelial","endothelial","immune","stromal","SNU601")){
  
  #Define output directory for the files
  if(comp=="epithelial"){
    sample_name<-"MCF7"
  } else {
    sample_name<-paste0("MCF7_",comp)
  }
  
  annot<-fread(paste0("data/input_",sample_name,"/sample_annotation.txt"),header=FALSE)
  
  print(sample_name)
  print(sum(annot$V2==sample_name))
  print(sum(annot$V2!=sample_name))
}


for(comp in c("","_immune","_TS","_SNU601")){
  
  #Define output directory for the files
  sample_name<-paste0("BCC08",comp)
  
  annot<-fread(paste0("data/input_",sample_name,"/sample_annotation.txt"),header=FALSE)
  
  print(sample_name)
  print(sum(annot$V2==sample_name))
  print(sum(annot$V2!=sample_name))
}
