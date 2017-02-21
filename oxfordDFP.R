oxfordGenes <- read.table("oxford/genes.txt",sep="\t",header = T)

oxfordClass <- read.table("oxford/patient_classes.txt",sep="\t",header = T)


randomForest(x=t(oxfordGenes),y=oxfordClass$x)

getDFP(RNAFinal = oxfordGenes,RNAPatientsFinal = oxfordClass,piVal = c(0.9,0.8,0.7),skipFactor = 3,zeta = 0.5,overlapping = 1,core = 2,filterGenes = F)
