varThreshold <- 10

rnanormvar <- rnanorm[apply(rnanorm, 1, var) > 10,]
idx = which(rf$importance>0.6)
rf2 <- randomForest(x = t(rnanormvar[idx,]),y = RNAPatientsFinal$stageClass)
