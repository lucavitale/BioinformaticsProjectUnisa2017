##### FIND OPTIMAL GENES

library(TopKLists)
path_class = "all78"
paths = list.files(paste("./",path_class,"_path/",sep=""))

for(path in paths){
  cat(path,"\n")
  path_mat = read.table(paste(path_class,"_path/",path,sep=""))
  path_name = unlist(strsplit(x = path,split = ".txt"))
  
  #ognuna di queste matrici è grande 6 x il numero di geni. 6 sono tutte le coppie di classi perchè la SVM multiclasse
  #in sklearn funzione SVC è implementata con uno schema uno-vs-uno
  #http://scikit-learn.org/stable/modules/generated/sklearn.svm.SVC.html
  fold0 = read.table(paste(path_class,"_res/",path_name,"_ranked_genes_fold0.txt",sep=""), quote="\"")
  fold1 = read.table(paste(path_class,"_res/",path_name,"_ranked_genes_fold1.txt",sep=""), quote="\"")
  fold2 = read.table(paste(path_class,"_res/",path_name,"_ranked_genes_fold2.txt",sep=""), quote="\"")
  fold3 = read.table(paste(path_class,"_res/",path_name,"_ranked_genes_fold3.txt",sep=""), quote="\"")
  fold4 = read.table(paste(path_class,"_res/",path_name,"_ranked_genes_fold4.txt",sep=""), quote="\"")
  colnames(fold0) = colnames(fold1) = colnames(fold2) = colnames(fold3) = colnames(fold4)=  colnames(path_mat)
  
  LL = list()
  
  ListMat = matrix("",nrow=nrow(fold0),ncol=ncol(fold0))
  
  for(i in 1:6){
    L0 = colnames(fold0)[order(fold0[i,]^2,decreasing=TRUE)]
    L1 = colnames(fold1)[order(fold1[i,]^2,decreasing=TRUE)]
    L2 = colnames(fold2)[order(fold2[i,]^2,decreasing=TRUE)]
    L3 = colnames(fold3)[order(fold3[i,]^2,decreasing=TRUE)]
    L4 = colnames(fold4)[order(fold4[i,]^2,decreasing=TRUE)]

    L = list(L0,L1,L2,L3,L4)
    BRes = Borda(L)
    LL[[i]] = BRes$TopK$mean
    ListMat[i,] = BRes$TopK$mean
  }
  
  BRes2 = Borda(LL)
  ListMat = rbind(ListMat,BRes2$TopK$mean)
  
  write.table(BRes2$TopK$mean,file = paste(path_class,"_res/",path_name,"_genes_final_rank.txt",sep=""),sep="\t",row.names = FALSE,col.names=FALSE,quote=FALSE)
  write.table(ListMat,file = paste(path_class,"_res/",path_name,"_genes_all_rank.txt",sep=""),sep="\t",row.names = FALSE,col.names=FALSE,quote=FALSE)
}

#### FIND OPTIMAL PATHWAY

path_class = "c6"
paths = list.files(paste(path_class,"_path/",sep=""))
path_name = unlist(strsplit(x = paths,split = ".txt"))

n_classes = 4
n_pathways = length(paths)

fold0 = read.table(paste(path_class,"_res/ranked_pathways_fold0.txt",sep=""), quote="\"")
fold1 = read.table(paste(path_class,"_res/ranked_pathways_fold1.txt",sep=""), quote="\"")
fold2 = read.table(paste(path_class,"_res/ranked_pathways_fold2.txt",sep=""), quote="\"")

nomi_colonne = c()
for(path in path_name){
  nomi_colonne = c(nomi_colonne,paste(path,"_class_",1:4,sep=""))
}

colnames(fold0) = colnames(fold1) = colnames(fold2) = nomi_colonne

LL = list()

ListMat = matrix("",nrow=nrow(fold0),ncol=ncol(fold0))

for(i in 1:6){
  L0 = colnames(fold0)[order(fold0[i,]^2,decreasing=TRUE)]
  L1 = colnames(fold1)[order(fold1[i,]^2,decreasing=TRUE)]
  L2 = colnames(fold2)[order(fold2[i,]^2,decreasing=TRUE)]
  
  L = list(L0,L1,L2)
  BRes = Borda(L)
  LL[[i]] = BRes$TopK$mean
  ListMat[i,] = BRes$TopK$mean
}

BRes2 = Borda(LL)
ListMat = rbind(ListMat,BRes2$TopK$mean)

write.table(BRes2$TopK$mean,file = paste(path_class,"_res/pathways_classes_final_rank.txt",sep=""),sep="\t",row.names = FALSE,col.names=FALSE,quote=FALSE)
write.table(ListMat,file = paste(path_class,"_res/pathways_classes_all_rank.txt",sep=""),sep="\t",row.names = FALSE,col.names=FALSE,quote=FALSE)
