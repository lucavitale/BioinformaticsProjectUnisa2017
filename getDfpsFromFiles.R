getDfpsFromFiles <- function(idx) {
  res <- list()
  for(i in idx) {
    fn <- paste("paramList", i, ".Rdata", sep = "")
    load(fn)
    res[i] <- list(
      skipFactor = paramList$skipFactor,
      zeta = paramList$zeta,
      piVal = paramList$piVal,
      overlapping = paramList$overlapping,
      dfps = paramList$dfps,
      numGenes = nrow(paramList$dfps)
    )
    print(paste(i, "done."))
  }
  return(res)
}