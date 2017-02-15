getTrainingPatients <- function(patientsData, testFraction) {
  idx1 <- which(patientsData$stageClass == "stage i")
  idx2 <- which(patientsData$stageClass == "stage ii")
  idx3 <- which(patientsData$stageClass == "stage iii")
  
  tr1 <- sample(idx1, length(idx1)*(1-testFraction))
  tr2 <- sample(idx2, length(idx2)*(1-testFraction))
  tr3 <- sample(idx3, length(idx3)*(1-testFraction))
  
  return(sort(c(tr1,tr2,tr3)))
}

trainingIdx <- getTrainingPatients(RNAPatientsFinal, 0.3)
save(trainingIdx, file = "trainingIdx.RData")