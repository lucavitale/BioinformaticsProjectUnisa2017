library(DFP)

#load("RNAFinal.Rdata")
#load("RNAPatientsFinal.Rdata")

getDFP <- function(RNAFinal, RNAPatientsFinal, datasetName, customFileName = NA, skipFactor = 3, zeta = 0.5, piVal = 0.5, overlapping = 1, filterGenes = TRUE, saveData = TRUE) {
  numberOfGenes = nrow(RNAFinal)
  
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
  shdaowDiscreteValues(dvs, featureNames(esRNA)[1:10])
  
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
      dfps = dfps,
      dataset = datasetName
    )
    
    if (is.na(customFileName)) {
      exFiles <- list.files(".", "paramList*")
      nFile <- length(exFiles) + 1
      save(paramList, file = paste("paramList", nFile, ".Rdata", sep = ""))
    } else {
      save(paramList, file = customFileName)
    }
  }

}