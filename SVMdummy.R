nvar = 1000

index <- 1:nrow(Rsvm)
testindex <- sample(index, trunc(length(index)/3))
testset <- Rsvm[testindex,c(1:nvar,ncol(Rsvm))]
trainset <- Rsvm[-testindex,c(1:nvar,ncol(Rsvm))]
svm.model <- svm(classes ~ ., data = trainset, cost = Inf, gamma = 1)
svm.pred <- predict(svm.model, testset[,-(nvar+1)])
table(pred = svm.pred, true = testset[,nvar+1])
