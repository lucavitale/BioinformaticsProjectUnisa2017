library(DFP)

readData <- function() {
  fileExprsGenes <- file.path(getwd(), "all_genes.txt"); fileExprsGenes
  fileExprsMirnas <- file.path(getwd(), "all_mirnas.txt"); fileExprsMirnas
  filePhenodata <- file.path(getwd(), "patient_classes.txt"); filePhenodata
  
  mirnaTable <- as.matrix(read.table(fileExprsMirnas))
  phenoData <- read.table(filePhenodata)
  mirnaTable <- mirnaTable[,which(colnames(mirnaTable) %in% rownames(phenoData))]
  phenoData <- AnnotatedDataFrame(phenoData)
  
  return(ExpressionSet(mirnaTable, phenoData = phenoData))
}

mirnaSet <- readData()