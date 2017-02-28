library(entropy)
library(parallel)
library(parmigene)

load("tcga_trainingIdx.RData")
load("tcgaDRPResults/paramList96.Rdata")
miRNA <- read.table("tcga_breast/miRNASeq.txt", header = T, sep = "\t")
miRNATraining <- miRNA[,trainingIdx]
RNA <- read.table("tcga_breast/RNASeq_full.txt", header = T, sep = "\t")
RNATraining <- RNA[rownames(paramList$dfps),trainingIdx]

cl <- makeCluster(8)

# empirical entropy (ML)
nbins <- floor(ncol(miRNATraining)/10)

nmi <- function(RNArow, miRNArow) {
  X <- miRNArow
  Y <- as.numeric(RNArow)
  freqTable <- discretize2d(X, Y, nbins, nbins)
  #mi <- mi.plugin(freqTable)
  mi <- mi.empirical(freqTable)
  res <- mi/(entropy(X) + entropy(Y))
  return(res)
  #knnmi(miRNArow,as.numeric(RNArow),3)
}

onAllGenes <- function(miRNArow) {
  apply(RNATraining, MARGIN = 1, FUN = nmi, as.numeric(miRNArow))
}

clusterExport(cl, varlist = c("nmi","RNATraining","miRNATraining","discretize2d","mi.plugin","nbins","entropy","knnmi","mi.empirical"))
nmi_miRNA_RNA <- parApply(cl, miRNATraining, MARGIN = 1, FUN = onAllGenes)

stopCluster(cl)
