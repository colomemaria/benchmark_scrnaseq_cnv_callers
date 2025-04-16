# ------------------------------------------------------------------------------
# Put the karyograms of all datasets together in one document
# ------------------------------------------------------------------------------

library(data.table)
library(ggplot2)
library(dplyr)

dataset_names<-c("SNU601","NCIN87","MKN45","KATOIII",
                 "NUGC4","SNU638","SNU668","HGC27", "SNU16",
                 "MCF7","COLO320","MM","BCC06","BCC06post","iAMP21")
dataset_ending<-setNames(c("","","","",
                           "","","","", "",
                           "","_cnvkit","_wes","_wes","_wes",""),
                         dataset_names)

ground_truth_type<-setNames(c("scWGS","scWGS","scWGS","scWGS",
                           "scWGS","scWGS","scWGS","scWGS", "scWGS",
                           "WGS","WGS","WES","WES","WES","scWGS"),
                         dataset_names)

#Vector for renaming methods (official published names)
method_names<-setNames(c("InferCNV (CNV)","InferCNV (Expr)","CaSpER",
                         "copyKat","SCEVAN (CNV)","SCEVAN (CNV)","SCEVAN (Expr)",
                         "Numbat (Expr)", "Numbat (CNV)", "CONICSmat"),
                       c("infercnv_cnv","infercnv_expr","casper",
                         "copykat","scevan_vega","scevan_cnv","scevan",
                         "numbat","numbat_cnv", "CONICSmat"))

karyoplots <- list()
for(dataset in dataset_names){
  
  combined_methods<-fread(paste0("results/output_",dataset,
                                 "/evaluation/outputs_allmethods_combined",
                                 dataset_ending[dataset],".tsv"))
  method_names_extended<-c(setNames(ground_truth_type[dataset],"wgs_mean"),method_names)
  
  #Remove unused columns
  combined_methods<-combined_methods[,-c("end","width","strand")]
  #Remove first column to chr (instead of seqnames)
  colnames(combined_methods)[1:2]<-c("chr","start_position")
  
  #Sort chromosomes correctly
  combined_methods$chr<-factor(combined_methods$chr,
                               levels=unique(combined_methods$chr))
  combined_methods<-as.data.frame(combined_methods)
  
  #Add an artifical count through the whole genome and 
  #get start positions for each new chromosome
  combined_methods$counted_pos<-1:nrow(combined_methods)
  chr_boundries<-combined_methods%>%
    group_by(chr)%>%
    summarize(start_chr=min(counted_pos),
              mean_chr=mean(counted_pos))
  
  #Scale every dataset to have diploid values at 0 and a standard deviation of 1
  scaling_factor<-NULL
  scaled_methods<-combined_methods
  method_names_extended<-method_names_extended[names(method_names_extended) %in% colnames(scaled_methods)]
  for(method in names(method_names_extended)){
    scaled_methods[,method]<-(scaled_methods[,method]-2) / 
      sd(scaled_methods[,method])
    
    #Save the standard deviation of each method to document the chosen normalization factor
    scaling_factor<-rbind(scaling_factor,
                          data.frame(method,
                                     sd=sd(combined_methods[,method]),
                                     mean=mean(combined_methods[,method]-2)))
  }
  
  plot_data<-reshape2::melt(scaled_methods,
                            id.vars=c("chr","start_position","counted_pos"))
  
  #Rename variable names
  plot_data$variable<-method_names_extended[as.character(plot_data$variable)]
  
  #Order them respectively
  plot_data$variable<-factor(plot_data$variable,
                             levels=method_names_extended)
  
  g<-ggplot(plot_data,aes(x=counted_pos,y=variable,fill=value))+geom_tile()+
    theme_bw()+
    ggtitle(paste0(LETTERS[which(dataset_names==dataset)],". Karyogram of ",dataset))+
    scale_fill_gradient2("Score",low = "darkblue",
                         mid = "white",high = "darkred",midpoint = 0,
                         breaks=c(-5,0,5),limits = c(-5,5),
                         labels=c("loss","base","gain"))+
    xlab("Chromosome position")+ylab("Method")+
    geom_vline(xintercept = chr_boundries$start_chr)+
    scale_x_continuous(breaks=chr_boundries$mean_chr,
                       labels=chr_boundries$chr)+
    scale_y_discrete(limits=rev)+
    coord_cartesian(xlim=c(1,max(plot_data$counted_pos)),expand=FALSE)+
    theme(legend.position="bottom",
          axis.text.x = element_text(angle=90,vjust=0.5,hjust=1),
          text=element_text(size=21))
  
  karyoplots <- c(karyoplots,list(g))

}

pdf("../figure_plots/all_karyograms.pdf", width = 17, height = 8)
for(kplot in karyoplots){
  grid::grid.draw(kplot)
}
dev.off()
