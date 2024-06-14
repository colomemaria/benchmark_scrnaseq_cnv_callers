# Benchmarking scRNA-seq copy number variant callers

Katharina T. Schmid 1, Aikaterini Symeonidi 1,2, Dmytro Hlushchenko 1, Maria L. Richter 1, Maria Colomé-Tatché 1,2

1 Biomedical Center (BMC), Physiological Chemistry, Faculty of Medicine, LMU Munich, Planegg-Martinsried, Germany
2 Institute of Computational Biology, Computational Health Center, Helmholtz Zentrum München, German Research Center for Environmental Health, Neuherberg, Germany

Work in progress, publication link will follow here as soon as available.

## Project description

We evaluated six popular computational methods in their ability to detect CNVs in 14 scRNA-seq datasets, comprising cancer cell lines, primary cancer samples and one diploid PBMC dataset. We assessed the methods according to their ability to recover the ground truth CNVs, estimaed with (sc)WGS or WES, using a large set of performance metrics. Additionally, we explored whether they could correctly identify euploid cells, especially also in fully diploid samples, and subclonal structures in heterogeneous tumor samples. We assessed also the scalability of each method.

![Project outline](pictures/figure1.png)

Assessed methods:

* CaSpER: https://github.com/akdess/CaSpER

* CONICSmat: https://github.com/diazlab/CONICS

* copyKat: https://github.com/navinlabcode/copykat 

* InferCNV: https://github.com/broadinstitute/infercnv

* Numbat: https://github.com/kharchenkolab/numbat/

* SCEVAN: https://github.com/AntonioDeFalco/SCEVAN 


## Running the pipeline



## Main results

For the CNV analysis of aneuploid samples, the evaluated methods showed a very variable performance which was in large part dataset dependent (Figure 2A). We identified certain dataset characteristics that impacted the performance, e.g. total coverage, number of measured cells, and the relative frequency of gains and losses.