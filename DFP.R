library(DFP)

load("RNAFinal.Rdata")
load("RNAPatientsFinal.Rdata")

filterGenes = TRUE
numberOfGenes = nrow(RNAFinal)
saveData = TRUE

phenoData <- RNAPatientsFinal
names(phenoData)[length(names(phenoData))] <- "class"
phenoData <- new("AnnotatedDataFrame",data = phenoData)

if (filterGenes) {
  vars <- apply(RNAFinal, 1, var)
  RNAReduced <- RNAFinal[order(vars, decreasing = TRUE)[1:numberOfGenes],]
  esRNA <- ExpressionSet(as.matrix(RNAReduced),phenoData = phenoData)
} else {
  esRNA <- ExpressionSet(as.matrix(RNAFinal),phenoData = phenoData)
}

skipFactor <- 3 # Factor to skip odd values
zeta <- 0.5 # Threshold used in the membership functions to label the float values with a discrete value
piVal <- 0.5 # Percentage of values of a class to determine the fuzzy patterns
overlapping <- 2 # Determines the number of discrete labels


mfs <- calculateMembershipFunctions(esRNA, skipFactor); mfs[[1]]
if (filterGenes) {
  toremove <- c()
  for(i in 1:length(mfs)){
    if (is.na(mfs[[i]]$lel@center))
      toremove <- c(toremove, i)
  }
  mfs <- mfs[-toremove]
  RNAReduced <- RNAReduced[-toremove,]
  esRNA <- ExpressionSet(as.matrix(RNAReduced),phenoData = phenoData)
}

#plotMembershipFunctions(esRNA, mfs, featureNames(esRNA)[1:2])


dvs <- discretizeExpressionValues(esRNA, mfs, zeta, overlapping); dvs[1:4,1:10]
showDiscreteValues(dvs, featureNames(esRNA)[1:10])

fps <- calculateFuzzyPatterns(esRNA, dvs, piVal, overlapping); fps[1:30,]
showFuzzyPatterns(fps, "stage i")[21:50]


dfps <- calculateDiscriminantFuzzyPattern(esRNA, fps); dfps[1:5,]
plotDiscriminantFuzzyPattern(dfps, overlapping)

if (saveData) {
  paramList <- list(
    skipFactor = skipFactor,
    zeta = zeta,
    piVal = piVal,
    overlapping = overlapping,
    filterGenes = filterGenes,
    numberOfGenes = numberOfGenes,
    mfs = mfs,
    dvs = dvs,
    fps = fps,
    dfps = dfps
  )
  
  exFiles <- list.files(".", "paramList*")
  nFile <- length(exFiles) + 1
  save(paramList, file = paste("paramList", nFile, ".Rdata", sep = ""))
}