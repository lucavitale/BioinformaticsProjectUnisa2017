####################################################################
### THIS SCRIPT FINDS THE GENES FOR EACH PATHWAY THAT ARE NEEDED ###
### TO REACH THE 80% OF THE CLASS SEPARABILITY CUMULATIVE SCORE  ###
####################################################################

library(readr)
path_class = "tcgaMiRNACommunities"
paths = list.files(paste("./",path_class,"/",sep=""))

for(path in paths){
  cat(path,"\n")
  path_mat = read.table(paste(path_class,"/",path,sep=""))
  path_name = unlist(strsplit(x = path,split = ".txt"))
  
  #ognuna di queste matrici è grande 6 x il numero di geni. 6 sono tutte le coppie di classi perchè la SVM multiclasse
  #in sklearn funzione SVC è implementata con uno schema uno-vs-uno
  #http://scikit-learn.org/stable/modules/generated/sklearn.svm.SVC.html
  fold0 = read.table(paste(path_class,"_res/",path_name,"_ranked_genes_fold0.txt",sep=""), quote="\"")
  fold1 = read.table(paste(path_class,"_res/",path_name,"_ranked_genes_fold1.txt",sep=""), quote="\"")
  fold2 = read.table(paste(path_class,"_res/",path_name,"_ranked_genes_fold2.txt",sep=""), quote="\"")
  # fold3 = read.table(paste(path_class,"_res/",path_name,"_ranked_genes_fold3.txt",sep=""), quote="\"")
  # fold4 = read.table(paste(path_class,"_res/",path_name,"_ranked_genes_fold4.txt",sep=""), quote="\"")
  colnames(fold0) = colnames(fold1) = colnames(fold2) = colnames(path_mat)#colnames(fold3) = colnames(fold4)=  colnames(path_mat)
  
  fold0 = abs(fold0)
  fold1 = abs(fold1)
  fold2 = abs(fold2)
  # fold3 = abs(fold3)
  # fold4 = abs(fold4)
  
  for(i in 1:nrow(fold0)){
    fold0[i,] = fold0[i,]/sum(fold0[i,])
  }
  
  for(i in 1:nrow(fold1)){
    fold1[i,] = fold1[i,]/sum(fold1[i,])
  }
  
  for(i in 1:nrow(fold2)){
    fold2[i,] = fold2[i,]/sum(fold2[i,])
  }
  
  # for(i in 1:nrow(fold3)){
  #   fold3[i,] = fold3[i,]/sum(fold3[i,])
  # }
  # 
  # for(i in 1:nrow(fold4)){
  #   fold4[i,] = fold4[i,]/sum(fold4[i,])
  # }
  
  fold0 = colMeans(fold0)
  fold1 = colMeans(fold1)
  fold2 = colMeans(fold2)
  # fold3 = colMeans(fold3)
  # fold4 = colMeans(fold4)
  # 
  fold14 = rbind(fold0,fold1,fold2)#,fold3,fold4)
  finalRank = colMeans(fold14)
  finalRank = finalRank[order(finalRank,decreasing = TRUE)]
  
  barplot(finalRank)
  plot(cumsum(finalRank))
  ng = sum(cumsum(finalRank) < 0.8)
  ggenes = names(finalRank)[1:ng]
  small_mat = path_mat[,ggenes]
  
  #### THE PATHWAY WITH THE RIGHT GENES WILL BE SAVED IN THE DIRECTORY allRet_path
  write.table(small_mat,paste("allRed_path/",path,sep=""))
}