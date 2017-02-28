library(entropy)
library(parallel)

load("tcga_trainingIdx.RData")
load("tcgaDRPResults/paramList96.Rdata")
miRNA <- read.table("tcga_breast/miRNASeq.txt", header = T, sep = "\t")
miRNATraining <- miRNA[,trainingIdx]
RNA <- read.table("tcga_breast/RNASeq_full.txt", header = T, sep = "\t")
RNATraining <- RNA[rownames(paramList$dfps),trainingIdx]

cl <- makeCluster(8)

# empirical entropy (ML)
nbins <- floor(ncol(miRNATraining)/30)

nmi <- function(RNArow, miRNArow) {
  X <- miRNArow
  Y <- as.numeric(RNArow)
  freqTable <- discretize2d(X, Y, nbins, nbins)
  mi <- mi.plugin(freqTable)
  res <- mi/(entropy(X) + entropy(Y))
  return(res)
}

onAllGenes <- function(miRNArow) {
  apply(RNATraining, MARGIN = 1, FUN = nmi, as.numeric(miRNArow))
}

clusterExport(cl, varlist = c("nmi","RNATraining","miRNATraining","discretize2d","mi.plugin","nbins","entropy"))
parApply(cl, miRNATraining, MARGIN = 1, FUN = onAllGenes)

stopCluster(cl)
