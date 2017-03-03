genGraphGridSearch <- function(corTh, miRNATh, walktrapStep){
  dat = list()
  i <- 0
  pb <- txtProgressBar(0, length(corTh) * length(miRNATh) * length(walktrapStep), 0, style = 3)
  for(a in corTh) {
    for(co in miRNATh){
      for(step in walktrapStep){
        dat[length(dat)+1] <- list(mirnaTh = a, corTh = co, step = step, res = generateGraph(a, co, step, FALSE))
        i <- i + 1
        setTxtProgressBar(pb, i)
      }
      
    }
  }
  close(pb)
}

genGraphGridSearch(seq(0.1,0.3,0.1), seq(0.6,0.9,0.1), 2^(1:9))