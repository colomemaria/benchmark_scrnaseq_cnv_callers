# ------------------------------------------------------------------------
# Reformat Aneufinder output (create a table format of bins times cells)
# ------------------------------------------------------------------------

from collections import defaultdict
import pandas as pd
import gzip
from natsort import natsort_keygen

#Function for creating a dictionary from the Aneufinder input
#Argument chr_col_num: whether the first column with chromosome is an integer (1) or "chr1"
def createDictionaryFromBed(bedfile, binsize=100000):
    bed_dict=defaultdict(list)
    lines=bedfile.readlines()
    l=[]
    for line in lines:
        line = line.decode()
        if line.strip()[0] != "t":
            chrom = line.strip().split("\t")[0]
            start = int(line.strip().split("\t")[1])
            end = int(line.strip().split("\t")[2])
            somy = int(line.strip().split("\t")[3].split("-")[0])
            if somy > 3:
                somy=3
            
            #Split the large bins into equally sized regions
            start = start+1
            while start < end:
                new_end=start+binsize-1
                l.append([chrom, start, new_end, somy])
                start=start+binsize
        
    for chrom, start, end,somy in l:
        bed_dict[chrom,start,end].append(somy)
    return(bed_dict)

#Function for filtering the Aneufinder dictionary in the same way as the epiAneufinder results.
#Bins having more than 85% of the cells with zero status are removed
def filterDictionary(dict):
    new_dict = {k: v for k, v in bed_dict.items() if (bed_dict[k].count(0)/len(bed_dict[k]))<0.85}
    return(new_dict)

#Aggregate pseudbobulk for scWGS (adapted from function "calculatePopulationSomies")
def createPseudobulkWGS(wgs_dict):
    gain_wgs = []
    loss_wgs = []
    base_wgs = [] 
    for k in list(wgs_dict):
            #Calculating pseudobulk representation for the scWGS. 1 is loss, 2 is disomic and 3 is gain
            loss_wgs.append((wgs_dict[k].count(1)+wgs_dict[k].count(0))/len(wgs_dict[k]))
            base_wgs.append(wgs_dict[k].count(2) / len(wgs_dict[k]))
            gain_wgs.append(wgs_dict[k].count(3) / len(wgs_dict[k]))
            
    #Format results into a pandas data frame
    data_pdf = {"chr" : [k[0].split("r")[1] for k in list(wgs_dict)], #Remove the chr to have the same format as other files
                "start" : [k[1] for k in list(wgs_dict)],
                "end" : [k[2] for k in list(wgs_dict)],
                "gain_wgs" : gain_wgs,
                "loss_wgs" : loss_wgs,
                "base_wgs" : base_wgs}
    
    df = pd.DataFrame(data=data_pdf)

    #Order the dataframe after chromsomes and position
    df = df.sort_values(by=["chr","start"],key=natsort_keygen())
    return(df)

##Process the Aneufinder results for the HCT116 dataset
#print("Processing HCT116 DNTRseq dataset")
#fin=gzip.open("data/HCT116_DNTRseq/aneufinder_results/binsize_1e+05_stepsize_1e+05_CNV.bed.gz")
#bed_dict=createDictionaryFromBed(fin)
#filtered_dict=filterDictionary(bed_dict)
##Delete object not longer used
#del bed_dict
#wgs_df = createPseudobulkWGS(filtered_dict)
#wgs_df.to_csv("data/HCT116_DNTRseq/aneufinder_results/wgs_results_formated.csv")

##Process the Aneufinder results for the A375 dataset
#print("Processing A375 DNTRseq dataset")
#fin=gzip.open("data/A375_DNTRseq/aneufinder_results/binsize_1e+05_stepsize_1e+05_CNV.bed.gz")
#bed_dict=createDictionaryFromBed(fin)
#filtered_dict=filterDictionary(bed_dict)
##Delete object not longer used
#del bed_dict
#wgs_df = createPseudobulkWGS(filtered_dict)
#wgs_df.to_csv("data/A375_DNTRseq/aneufinder_results/wgs_results_formated.csv")

#Process the Aneufinder results for the A375 dataset
print("Processing mouse dataset")
fin=gzip.open("data/mouse_UMCG/aneufinder_results/binsize_1e+05_stepsize_1e+05_CNV.bed.gz")
bed_dict=createDictionaryFromBed(fin)
filtered_dict=filterDictionary(bed_dict)
#Delete object not longer used
del bed_dict
wgs_df = createPseudobulkWGS(filtered_dict)
wgs_df.to_csv("data/mouse_UMCG/aneufinder_results/wgs_results_formated.csv")

##Process the Aneufinder results for the iAMP21 dataset
#print("Processing iAMP dataset (ALL)")
#fin=gzip.open("data/iAMP21/aneufinder_results/binsize_1e+05_stepsize_1e+05_CNV.bed.gz")
#bed_dict=createDictionaryFromBed(fin)
#filtered_dict=filterDictionary(bed_dict)
##Delete object not longer used
#del bed_dict
#wgs_df = createPseudobulkWGS(filtered_dict)
#wgs_df.to_csv("data/iAMP21/aneufinder_results/wgs_results_formated.csv")

##Process the Aneufinder results for the different gastric cell lines
#for cell_line in ["NCIN87","MKN45","NUGC4","SNU638"]:
#for cell_line in ["KATOIII","SNU16","SNU668","HGC27"]:
#    print("Processing "+cell_line)
#    fin=gzip.open("data/gastric_cell_lines/"+cell_line+"/aneufinder_results/binsize_1e+05_stepsize_1e+05_CNV.bed.gz")
#    bed_dict=createDictionaryFromBed(fin)
#    filtered_dict=filterDictionary(bed_dict)
#    #Delete object not longer used
#    del bed_dict
#    wgs_df = createPseudobulkWGS(filtered_dict)
#    wgs_df.to_csv("data/gastric_cell_lines/"+cell_line+"/aneufinder_results/wgs_results_formated.csv")