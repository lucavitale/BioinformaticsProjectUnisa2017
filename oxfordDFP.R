source("DFP.R")
oxfordGenes <- read.table("oxford/genes.txt",sep="\t",header = T)
oxfordClass <- read.table("oxford/patient_classes.txt",sep="\t",header = T)

testFraction <- 0.3

idx1 <- which(oxfordClass$x == "Grave_1")
idx2 <- which(oxfordClass$x == "Grave_2")
idx3 <- which(oxfordClass$x == "Grave_3")
idx4 <- which(oxfordClass$x == "Grave_4")

tr1 <- sample(idx1, length(idx1)*(1-testFraction))
tr2 <- sample(idx2, length(idx2)*(1-testFraction))
tr3 <- sample(idx3, length(idx3)*(1-testFraction))
tr4 <- sample(idx4, length(idx4)*(1-testFraction))
  
trainingIdx <- sort(c(tr1,tr2,tr3,tr4))
save(trainingIdx, file = "trainingIdx.RData")

#randomForest(x=t(oxfordGenes),y=oxfordClass$x)

multiDFP(RNAFinal = oxfordGenes[,trainingIdx],
         RNAPatientsFinal = oxfordClass[trainingIdx,],
         piVal = seq(0.4, 0.7, 0.05),
         skipFactor = c(0, 1, 2, 3),
         zeta = c(0.35, 0.4, 0.45, 0.5),
         overlapping = c(1, 2),
         core = 7,
         filterGenes = F,
         savePath = "oxfordDFPResults",
         datasetName = "oxford training dataset")