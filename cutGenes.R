varThreshold <- 10

rnanormvar <- rnanorm[apply(rnanorm, 1, var) > 10,]
