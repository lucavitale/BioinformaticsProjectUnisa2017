library(e1071)
#nvar = 1000

#index <- 1:nrow(Rsvm)
#testindex <- sample(index, trunc(length(index)/3))
#testset <- Rsvm[testindex,c(1:nvar,ncol(Rsvm))]
#trainset <- Rsvm[-testindex,c(1:nvar,ncol(Rsvm))]

trainset <- as.data.frame(t(rnanorm[rownames(paramList$dfps),trainingIdx]))
trainset$classes <- RNAPatientsFinal[trainingIdx,"stageClass"]
testset <- as.data.frame(t(rnanorm[rownames(paramList$dfps),-trainingIdx]))
testset$classes <- RNAPatientsFinal[-trainingIdx,"stageClass"]

svm.model <- svm(classes ~ ., data = trainset, cost = 1e12, gamma = 1, kernel="radial")
svm.pred <- predict(svm.model, testset[,-(ncol(testset))])
table(pred = svm.pred, true = testset[,(ncol(testset))])