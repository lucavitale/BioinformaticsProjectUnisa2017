library(parallel)

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
