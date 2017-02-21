library(e1071)
#nvar = 1000

#index <- 1:nrow(Rsvm)
#testindex <- sample(index, trunc(length(index)/3))
#testset <- Rsvm[testindex,c(1:nvar,ncol(Rsvm))]
#trainset <- Rsvm[-testindex,c(1:nvar,ncol(Rsvm))]

cl13 <- which(RNAPatientsFinal$stageClass %in% c("stage i", "stage iii"))
cl2 <- which(RNAPatientsFinal$stageClass %in% c("stage ii"))

trainset <- as.data.frame(t(rnanorm[rownames(paramList$dfps),intersect(trainingIdx,cl13)]))
trainset$classes <- RNAPatientsFinal[intersect(trainingIdx,cl13),"stageClass"]
testset <- as.data.frame(t(rnanorm[rownames(paramList$dfps),-c(intersect(trainingIdx,cl13),cl2)]))
testset$classes <- RNAPatientsFinal[-c(intersect(trainingIdx,cl13),cl2),"stageClass"]

svm.model <- svm(classes ~ ., data = trainset, cost = 1)
svm.pred <- predict(svm.model, testset[,-(ncol(testset))])
table(pred = svm.pred, true = testset[,(ncol(testset))])
