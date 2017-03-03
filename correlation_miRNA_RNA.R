#library(entropy)
library(parallel)
#library(parmigene)

load("perm_test_results.Rdata")
#MM <- MM[-which(MM$Pathway == "GO_REGULATION_OF_CELL_CYCLE_PHASE_TRANSITION"),]
pathwayPath <- paste("tcgaPathways/", MM$Pathway, ".txt", sep="")

getGenes <- function(path) {
  nm <- names(read.table(path, header = T, sep = " ", check.names=FALSE))
}

genes <- unlist(sapply(pathwayPath, getGenes))
genes <- unique(genes)

load("tcga_trainingIdx.RData")
#load("tcgaDRPResults/paramList96.Rdata")
miRNA <- read.table("tcga_breast/miRNASeq.txt", header = T, sep = "\t")
miRNATraining <- miRNA[,trainingIdx]
RNA <- read.table("tcga_breast/RNASeq_full.txt", header = T, sep = "\t")
RNATraining <- RNA[genes,trainingIdx]

cl <- makeCluster(8)

# empirical entropy (ML)
#nbins <- floor(ncol(miRNATraining)/10)

correlate <- function(RNArow, miRNArow) {
  X <- miRNArow
  Y <- as.numeric(RNArow)
  #freqTable <- discretize2d(X, Y, nbins, nbins)
  #mi <- mi.plugin(freqTable)
  #mi <- mi.empirical(freqTable)
  #res <- mi/(entropy(X) + entropy(Y))
  res <- cor(X, Y)
  return(res)
  #knnmi(miRNArow,as.numeric(RNArow),3)
}

onAllGenes <- function(miRNArow) {
  apply(RNATraining, MARGIN = 1, FUN = correlate, as.numeric(miRNArow))
}

#clusterExport(cl, varlist = c("correlate","RNATraining","miRNATraining","discretize2d","mi.plugin","nbins","entropy","knnmi","mi.empirical"))
clusterExport(cl, varlist = c("correlate","RNATraining","miRNATraining"))
cor_miRNA_RNA <- parApply(cl, miRNATraining, MARGIN = 1, FUN = onAllGenes)

stopCluster(cl)

save(cor_miRNA_RNA, file = "cor_miRNA_RNA.Rdata")
