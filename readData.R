library(DFP)

readData <- function() {
  fileExprsGenes <- file.path(getwd(), "all_genes.txt"); fileExprsGenes
  fileExprsMirnas <- file.path(getwd(), "all_mirnas.txt"); fileExprsMirnas
  filePhenodata <- file.path(getwd(), "patient_classes.txt"); filePhenodata
  
  mirnaTable <- as.matrix(read.table(fileExprsMirnas))
  phenoData <- read.table(filePhenodata)
  names(phenoData) <- c("class", "pat_id")
  mirnaTable <- mirnaTable[,which(colnames(mirnaTable) %in% rownames(phenoData))]
  phenoData <- AnnotatedDataFrame(phenoData)
  
  return(ExpressionSet(mirnaTable, phenoData = phenoData))
}

mirnaSet <- readData()
res <- discriminantFuzzyPattern(mirnaSet, piVal = 0.6)
plotDiscriminantFuzzyPattern(res$discriminant.fuzzy.pattern)
