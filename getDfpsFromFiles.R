getDfpsFromFiles <- function(idx) {
  res <- list()
  for(i in idx) {
    fn <- paste("paramList", i, ".Rdata", sep = "")
    load(fn)
    res[[i]] <- list(
      skipFactor = paramList$skipFactor,
      zeta = paramList$zeta,
      piVal = paramList$piVal,
      overlapping = paramList$overlapping,
      dfps = paramList$dfps,
      numGenes = nrow(paramList$dfps)
    )
  }
  return(res)
}

getDfpResults <- function(idx){
  df <- data.frame()
  for(i in idx) {
    fn <- paste("paramList", i, ".Rdata", sep = "")
    load(fn)
    df = rbind(df,data.frame(paramList$skipFactor, paramList$zeta, paramList$piVal, paramList$overlapping, nrow(paramList$dfps),fn))
  }
  colnames(df) <- c("skipFactor","zeta","piVal","overlapping","numGenes","fileName")
  return(df)
}
