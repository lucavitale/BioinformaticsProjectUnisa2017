source("DFP.R")
oxfordGenes <- read.table("tcga_breast/RNASeq_full.txt",sep="\t",header = T)
oxfordClass <- read.table("tcga_breast/patient_classes.txt",sep="\t",header = T)

testFraction <- 0.3

idx1 <- which(oxfordClass$x == "Basal")
idx2 <- which(oxfordClass$x == "LumA")
idx3 <- which(oxfordClass$x == "LumB")
idx4 <- which(oxfordClass$x == "Her2")

tr1 <- sample(idx1, length(idx1)*(1-testFraction))
tr2 <- sample(idx2, length(idx2)*(1-testFraction))
tr3 <- sample(idx3, length(idx3)*(1-testFraction))
tr4 <- sample(idx4, length(idx4)*(1-testFraction))

trainingIdx <- sort(c(tr1,tr2,tr3,tr4))
save(trainingIdx, file = "tcga_trainingIdx.RData")

#randomForest(x=t(oxfordGenes),y=oxfordClass$x)

multiDFP(RNAFinal = oxfordGenes[,trainingIdx],
         RNAPatientsFinal = oxfordClass[trainingIdx,],
         piVal = seq(0.4, 0.8, 0.05),
         skipFactor = c(0, 1, 2, 3),
         zeta = c(0.35, 0.4, 0.45, 0.5),
         overlapping = c(1, 2),
         core = 8,
         filterGenes = F,
         savePath = "tcgaDRPResults",
         datasetName = "tcga training dataset")

load("tcga_trainingIdx.RData")

getDFP(RNAFinal = oxfordGenes[,trainingIdx],
       RNAPatientsFinal = oxfordClass[trainingIdx,],
       piVal = seq(0.45, 0.7, 0.05),
       skipFactor = 2,
       zeta = 0.5,
       overlapping = 2,
       core = 7,
       filterGenes = F,
       savePath = "oxfordDFPResults",
       datasetName = "oxford training dataset")
