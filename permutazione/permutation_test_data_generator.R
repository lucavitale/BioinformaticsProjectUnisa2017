#### PERMUTATION TEST
dataset = read.table("tcga_breast/RNASeq_full.txt",sep="\t",header=T)
load("tcgaDRPResults/paramList96.Rdata")
DFP_genes = rownames(paramList$dfps)

n_perm_test = 1000

##### THE PERMUTATION TEST WILL BE RUN ONLY ON THE GENES THAT COMES FROM THE DFP ANALYSIS ##### 
X = dataset[DFP_genes,]
load("goodPathway.Rdata")
goodPathway=paste(goodPathway,".txt",sep="")

files <- list.files("tcgaPathways")
found = c()
for (f in 1:length(files)){
  if(files[f] %in% goodPathway){
    cat("moving file\n")
    found = c(found,f)
    file.rename(from= paste("tcgaPathways/",files[f],sep=""),to= paste("allRed_path/",files[f],sep="")) 
  }
}

pathways = list.files("allRed_path/")
path_class = "allRed"

dimensions = c() #NUMBER OF GENES IN EACH PATHWAYS
for(pat in pathways){
  path_mat = read.table(paste(path_class,"_path/",pat,sep=""))
  dimensions = c(dimensions,ncol(path_mat))
}

####################### CREATE DATASETS FOR PERMUTATION TEST #######################
pb = txtProgressBar(min = 1,max = n_perm_test,style=3)
for(i in 1:n_perm_test){
  dir.create(paste("perm_test/permutation_test_",i,"/",sep=""))
  for(j in 1:length(dimensions)){
    dir.create(paste("perm_test/permutation_test_",i,"/all_path",sep=""))
    dir.create(paste("perm_test/permutation_test_",i,"/all_res",sep=""))
    
    pj = t(X[rownames(X)[sample(x = 1:nrow(X),size = dimensions[j])],])
    write.table(pj,file=paste("perm_test/permutation_test_",i,"/all_path/path_",j,".txt",sep=""),quote=FALSE)
  }
  setTxtProgressBar(pb,i)
}
close(pb)

####################### RETRIEVE PERMUTATION TEST RESULTS #######################
all_78_avg_acc = read.table("allRed_res/test_avg_accuracy.txt",header=TRUE,sep=",")
accuracies = c()
nPerm = n_perm_test
pb = txtProgressBar(min = 1,max = n_perm_test,style=3)
for(i in 1:nPerm){
  avg_acc = read.table(paste("perm_test/permutation_test_",i,"/all_res/test_avg_accuracy.txt",sep=""),header=TRUE,sep=",")
  accuracies = cbind(accuracies,avg_acc$avg_acc)
  setTxtProgressBar(pb,i)
}
close(pb)
percs = c()
for(i in 1:nrow(avg_acc)){
  percs = c(percs,sum(all_78_avg_acc$avg_acc[i] > accuracies[i,]) / nPerm)
}
percs[all_78_avg_acc$avg_acc>0.78]
View(cbind(all_78_avg_acc[percs>0.95,],percs[percs>0.95]))
M = cbind(all_78_avg_acc,percs)
View(M)
#write.table(M,file="perm_test_res_red_pathways.txt",sep="\t",quote = FALSE)
MM = cbind(all_78_avg_acc[percs>0.95,],percs[percs>0.95])
View(MM)
#write.table(MM,file="significant_pathways_after_perm_test.txt",sep="\t",quote = FALSE)
