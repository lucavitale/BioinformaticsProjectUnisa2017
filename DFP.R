library(DFP)
library(parallel)

#load("RNAFinal.Rdata")
#load("RNAPatientsFinal.Rdata")

parDiscretizeExpressionValues <- function(cluster, rmadataset, mfs, zeta = 0.5, overlapping = 2) 
{
  rmam <- exprs(rmadataset)
  rmam[c(1:8), c(1:4)]
  rmav <- as.vector(pData(phenoData(rmadataset))$class)
  rmav
  names(rmav) <- sampleNames(rmadataset)
  rmav
  gene.names <- rownames(rmam)
  gene.names
  dvs <- NULL
  #for (ig in gene.names) {
  #  values <- rmam[ig, ]
  #  values
  #  disc.values <- .fuzzyDiscretization(mfs[[ig]]$lel, mfs[[ig]]$mel, 
  #                                      mfs[[ig]]$hel, values, zeta, overlapping)
  #  disc.values
  #  dvs <- rbind(dvs, disc.values)
  #  dvs
  #}
  
  doit <- function(ig) {
    #browser()
    values <- rmam[ig, ]
    disc.values <- DFP:::.fuzzyDiscretization(mfs[[ig]]$lel, mfs[[ig]]$mel, 
                                              mfs[[ig]]$hel, values, zeta, overlapping)
    #return(c(ig, disc.values))
    print("done")
    return(disc.values)
  }
  
  dvs <- t(parSapply(cluster, gene.names, doit))
  
  rownames(dvs) <- gene.names
  head(dvs)
  attr(dvs, "types") <- rmav
  dvs
  return(dvs)
}

skipOddValues <- function (values, skipFactor = 3) 
{
  if (skipFactor > 0) {
    orderv <- order(values)
    orderv
    vals <- values[orderv]
    vals
    first <- trunc(length(vals)/4)
    first
    third <- trunc(length(vals)/4) * 3
    third
    firstValue <- vals[first + 1]
    firstValue
    thirdValue <- vals[third + 1]
    thirdValue
    RIC <- thirdValue - firstValue
    RIC
    lowBarrier <- firstValue - (skipFactor * RIC)
    lowBarrier
    highBarrier <- thirdValue + (skipFactor * RIC)
    highBarrier
    isOutlier <- values < lowBarrier | values > highBarrier
    isOutlier
  } else if (skipFactor == 0) {
    isOutlier <- rep(FALSE, length(values))
  }
  return(isOutlier)
}

assignInNamespace(".skipOddValues", skipOddValues, "DFP")

getDFP <- function(RNAFinal, RNAPatientsFinal, datasetName, customFileName = NA, restoreFromDvs = NA, skipFactor = 3, zeta = 0.5, piVal = 0.5, overlapping = 1, filterGenes = TRUE, saveData = TRUE, core = 1) {
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
  
  print(paste("Number of selected genes:", nrow(RNAReduced)))
  #plotMembershipFunctions(esRNA, mfs, featureNames(esRNA)[1:2])
  
  if (class(restoreFromDvs) != "matrix" && is.na(restoreFromDvs)) {
    if (core > 1) {
      cl <- makeCluster(core, outfile = "progress.log")
      dvs <- parDiscretizeExpressionValues(cl, esRNA, mfs, zeta, overlapping); #dvs[1:4,1:10]
      stopCluster(cl)
    } else {
      dvs <- discretizeExpressionValues(esRNA, mfs, zeta, overlapping); #dvs[1:4,1:10]
    }
    assign("dvs", dvs, envir = .GlobalEnv)
    #showDiscreteValues(dvs, featureNames(esRNA)[1:10])
  } else {
    dvs <- restoreFromDvs
  }
  
  it <- 0
  for (val in piVal) {
    it <- it + 1
    print(paste("Calculating dfps", it, "of", length(piVal), "..."))
    
    fps <- calculateFuzzyPatterns(esRNA, dvs, val, overlapping); #fps[1:30,]
    #showFuzzyPatterns(fps, "stage i")[21:50]
    
    
    dfps <- calculateDiscriminantFuzzyPattern(esRNA, fps); #dfps[1:5,]
    #plotDiscriminantFuzzyPattern(dfps, overlapping)
    
    if (saveData) {
      paramList <- list(
        skipFactor = skipFactor,
        zeta = zeta,
        piVal = val,
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
        nFile <- 1
        while (file.exists(paste("paramList", nFile, ".Rdata", sep = ""))) {
          nFile <- nFile + 1
        }
        save(paramList, file = paste("paramList", nFile, ".Rdata", sep = ""))
      } else {
        save(paramList, file = customFileName)
      }
    }
  }
  file.remove("progress.log")
  
}
