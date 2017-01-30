library(DFP)

readData <- function() {
  fileExprsGenes <- file.path(getwd(), "all_genes.txt"); fileExprsGenes
  fileExprsMirnas <- file.path(getwd(), "all_mirnas.txt"); fileExprsMirnas
  filePhenodata <- file.path(getwd(), "patient_classes.txt"); filePhenodata
  browser()
  geneData <- readExpressionSet(fileExprsMirnas, phenoDataFile = filePhenodata, exprsArgs = list(), phenoDataArgs = list(quote="\""))
}
readData()
