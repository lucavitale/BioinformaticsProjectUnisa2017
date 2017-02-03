plotPatients <- function(ncol = 50, mixed = FALSE, maxy = 500000) {
  st1 <- which(RNAPatientsFinal$stageClass == "stage i")
  st2 <- which(RNAPatientsFinal$stageClass == "stage ii")
  st3 <- which(RNAPatientsFinal$stageClass == "stage iii")
  
  #ncol = 50
  if (mixed) {
    par(mfrow=c(1,1))
  } else {
    par(mfrow=c(3,1))
  }
  #maxy = 500000
  
  s1 <- sample(st1, ncol)
  s2 <- sample(st2, ncol)
  s3 <- sample(st3, ncol)
  
  if (mixed) {
    Sfull <- c(s1,s2,s3)
  }
  
  #apply(RNAReduced[s1], 2, boxplot, ylim = c(0, maxy), col = 2)
  #apply(RNAReduced[s2], 2, boxplot, ylim = c(0, maxy), col = 3)
  #apply(RNAReduced[s3], 2, boxplot, ylim = c(0, maxy), col = 4)
  
  if (mixed) {
    boxplot(RNAReduced[Sfull], ylim = c(0, maxy), col = c(rep(2, ncol), rep(3, ncol), rep(4, ncol)))
  } else {
    boxplot(RNAReduced[s1], ylim = c(0, maxy), col = 2)
    boxplot(RNAReduced[s2], ylim = c(0, maxy), col = 3)
    boxplot(RNAReduced[s3], ylim = c(0, maxy), col = 4)
  }
}