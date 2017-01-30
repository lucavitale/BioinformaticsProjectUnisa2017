library(DFP)

readData <- function() {
  fileExprsGenes <- file.path(getwd(), "all_genes.txt"); fileExprsGenes
  fileExprsMirnas <- file.path(getwd(), "all_mirnas.txt"); fileExprsMirnas
  filePhenodata <- file.path(getwd(), "patient_classes.txt"); filePhenodata
  mirnaData <- readExpressionSet(fileExprsMirnas, phenoDataFile = filePhenodata, exprsArgs = list(), phenoDataArgs = list(quote="\""))
}

readData()

readData2 <- function() {
  fileExprsGenes <- file.path(getwd(), "all_genes.txt"); fileExprsGenes
  fileExprsMirnas <- file.path(getwd(), "all_mirnas.txt"); fileExprsMirnas
  filePhenodata <- file.path(getwd(), "patient_classes.txt"); filePhenodata
  
  mirnaTable <- as.matrix(read.table(fileExprsMirnas))
  phenoData <- read.table(filePhenodata)
}