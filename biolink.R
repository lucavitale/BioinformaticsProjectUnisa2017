library(TCGAbiolinks)
library(SummarizedExperiment)

RNAq <- GDCquery(project = "TCGA-BRCA", data.category = "Transcriptome Profiling",data.type = "Gene Expression Quantification", workflow.type = "HTSeq - FPKM-UQ")
GDCdownload(RNAq,chunks.per.download= 150)

miRNAq <- GDCquery(project = "TCGA-BRCA", data.category = "Transcriptome Profiling",data.type = "miRNA Expression Quantification", workflow.type = "BCGSC miRNA Profiling")
GDCdownload(miRNAq, chunks.per.download= 250)

RNATable <-  GDCprepare(RNAq,mut.pipeline = NULL)
miRNATable <- GDCprepare(miRNAq,mut.pipeline = NULL)