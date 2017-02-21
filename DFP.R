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

parCalculateDiscriminantFuzzyPattern <- function (cluster, rmadataset, fps) 
{
  discriminants <- NULL
  #for (ig in featureNames(rmadataset)) {
  #  table.facFP <- table(factor(fps[ig, ]))
  #  max.facFP <- ifelse(sum(table.facFP) == 0, 0, max(table.facFP))
  #  max.facFP
  #  if (max.facFP > 0 & max.facFP < sum(table.facFP)) 
  #    discriminants <- c(discriminants, ig)
  #}
  
  doit <- function(ig) {
    table.facFP <- table(factor(fps[ig, ]))
    max.facFP <- ifelse(sum(table.facFP) == 0, 0, max(table.facFP))
    max.facFP
    res <- NA
    if (max.facFP > 0 & max.facFP < sum(table.facFP)) {
      #print(paste("max:", max.facFP, "| sum:", sum(table.facFP), "max>0?", max.facFP > 0, "| max < sum?", max.facFP < sum(table.facFP)))
      res <- ig
    }
    ##print("done")
    return(res)
  }
  
  discriminants <- parSapply(cluster, featureNames(rmadataset), doit)
  discriminants <- discriminants[!is.na(discriminants)]
  #browser()
  dfp <- fps[discriminants, ]
  dfp
  attr(dfp, "ifs") <- attr(fps, "ifs")[discriminants, ]
  dfp
  return(dfp)
}

parCalculateFuzzyPatterns <- function (cluster, rmadataset, dvs, piVal = 0.9, overlapping = 2) 
{
  if (overlapping == 1) {
    disc.alphab <- c("Low", "Medium", "High")
  }
  else if (overlapping == 2) {
    disc.alphab <- c("Low", "Low-Medium", "Medium", "Medium-High", 
                     "High")
  }
  else {
    disc.alphab <- c("Low", "Low-Medium", "Low-Medium-High", 
                     "Medium", "Medium-High", "High")
  }
  fps <- NULL
  ifs <- NULL
  #for (ig in featureNames(rmadataset)) {
  #  disc.values <- dvs[ig, ]
  #  disc.values
  #  attr(disc.values, "types") <- attr(dvs, "types")
  #  disc.values
  #  fuzzypat <- .fuzzyPatterns(disc.values, disc.alphab, 
  #                             piVal)
  #  fuzzypat
  #  fps <- rbind(fps, fuzzypat)
  #  fps
  #  ifs <- rbind(ifs, attr(fuzzypat, "ifs"))
  #  ifs
  #}
  
  doit <- function(ig) {
    disc.values <- dvs[ig, ]
    disc.values
    attr(disc.values, "types") <- attr(dvs, "types")
    disc.values
    fuzzypat <- DFP:::.fuzzyPatterns(disc.values, disc.alphab, 
                               piVal)
    fuzzypat
    #print(paste(fuzzypat, "|", attr(fuzzypat, "ifs")))
    return(list(fuzzypat, attr(fuzzypat, "ifs")))
  }
  
  result <- parLapply(cluster, featureNames(rmadataset), doit)
  
  fps <- parLapply(cluster, result, "[[", 1)
  ifs <- parLapply(cluster, result, "[[", 2)
  fps <- t(simplify2array(fps))
  ifs <- t(simplify2array(ifs))
  
  #if (class(fps == "character")) { # if it's a vector, then it has only 1 row
  #  fps <- matrix(fps, nrow = 1)
  #  ifs <- matrix(attr(fps, "ifs"), nrow = 1)
  #}
  
  rownames(fps) <- featureNames(rmadataset)
  head(fps)
  rownames(ifs) <- featureNames(rmadataset)
  head(ifs)
  attr(fps, "ifs") <- ifs
  head(fps)
  return(fps)
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
    print(paste("Calculating dfps ", it, " of ", length(piVal), "...", sep = ""))
    
    if (core > 1) {
      #cl <- makeCluster(core, outfile = "progressFP.log")
      cl <- makeCluster(core)
      fps <- parCalculateFuzzyPatterns(cl, esRNA, dvs, val, overlapping);
      stopCluster(cl)
    } else {
      fps <- calculateFuzzyPatterns(esRNA, dvs, val, overlapping); #fps[1:30,]
      #showFuzzyPatterns(fps, "stage i")[21:50]
    }
    
    if (core > 1) {
      #cl <- makeCluster(core, outfile = "progressDFP.log")
      cl <- makeCluster(core)
      dfps <- parCalculateDiscriminantFuzzyPattern(cl, esRNA, fps);
      stopCluster(cl)
    } else {
      dfps <- calculateDiscriminantFuzzyPattern(esRNA, fps); #dfps[1:5,]
      #plotDiscriminantFuzzyPattern(dfps, overlapping)
    }
    
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
    
    gc() # execute garbage collection between calls
  }
  #file.remove("progress.log")
  
}

multiDFP <- function(RNAFinal, RNAPatientsFinal, datasetName, skipFactor = 3, zeta = 0.5, piVal = 0.5, overlapping = 1, filterGenes = TRUE, saveData = TRUE, core = 1) {
  totalCalls <- length(skipFactor) * length(zeta) * length(overlapping)
  i <- 0
  for (sf in skipFactor) {
    for (z in zeta) {
      for (o in overlapping) {
        i <- i + 1
        print(paste("Call ", i, " of ", totalCalls, "...", sep = ""))
        getDFP(RNAFinal = RNAFinal, RNAPatientsFinal = RNAPatientsFinal, datasetName = datasetName, skipFactor = sf, zeta = z, piVal = piVal, overlapping = o, filterGenes = filterGenes, saveData = saveData, core = core)
      }
    }
  }
}
