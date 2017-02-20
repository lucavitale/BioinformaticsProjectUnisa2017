confusionMatrixFiles <- list.files(path = ".",pattern="*conf_mat.txt")
a <- list()
for (i in 1:length(confusionMatrixFiles)){
  a[[i]] <- read.table(confusionMatrixFiles[i],header=F,sep = " ")
}

b = list()
j=1
for (i in 1:length(confusionMatrixFiles)){
  if(sum(a[[1]] == a[[i]]) != 9){
    b[[j]] <- a[[i]]
    j=j+1
  }
}