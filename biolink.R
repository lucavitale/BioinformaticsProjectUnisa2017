library(TCGAbiolinks)

RNAq <- GDCquery(project = "TCGA-BRCA", data.category = "Transcriptome Profiling",data.type = "Gene Expression Quantification", workflow.type = "HTSeq - FPKM")
GDCdownload(RNAq)


miRNAq <- GDCquery(project = "TCGA-BRCA", data.category = "Transcriptome Profiling",data.type = "miRNA Expression Quantification", workflow.type = "BCGSC miRNA Profiling")
GDCdownload(miRNAq)
