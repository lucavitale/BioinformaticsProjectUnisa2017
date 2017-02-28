##### FUNCTIONS DEFINITIONS

# this functions return the dimensions of the pathways that are in the folder "path_class"
get_pathway_dimensions = function(path_class){
  pathways = list.files(paste(path_class,"_path/",sep=""))

  dimensions = c() #NUMBER OF GENES IN EACH PATHWAYS
  for(pat in pathways){
    path_mat = read.table(paste(path_class,"_path/",pat,sep=""))
    dimensions = c(dimensions,ncol(path_mat))
  }
  return(dimensions)
}

# get_pathway_dimensions = function(path_class){
#   pathways = list.files(paste(path_class,"/",sep=""))
# 
#   dimensions = c() #NUMBER OF GENES IN EACH PATHWAYS
#   for(pat in pathways){
#     path_mat = read.table(paste(path_class,"/",pat,sep=""))
#     dimensions = c(dimensions,ncol(path_mat))
#   }
#   return(dimensions)
# }

#this functions create the dataset for the permutation test. In the perm_test_dir it creates n_perm_test folders.
# each folder contains a directory named all_path (that contains a file for each pathway)
# and a directory named all_res (that will contain the results of the classification task performed with the py functions)

create_dataset_for_permutation_test = function(X,perm_test_dir="perm_test",n_perm_test=100,dimensions){
  pb = txtProgressBar(min = 1,max = n_perm_test,style=3)
  for(i in 1:n_perm_test){
    dir.create(paste(perm_test_dir,"/permutation_test_",i,"/",sep=""))
    for(j in 1:length(dimensions)){
      dir.create(paste(perm_test_dir,"/permutation_test_",i,"/all_path",sep=""))
      dir.create(paste(perm_test_dir,"/permutation_test_",i,"/all_res",sep=""))
      
      pj = t(X[rownames(X)[sample(x = 1:nrow(X),size = dimensions[j])],])
      write.table(pj,file=paste(perm_test_dir,"/permutation_test_",i,"/all_path/path_",j,".txt",sep=""),quote=FALSE)
    }
    setTxtProgressBar(pb,i)
  }
  close(pb)
}

get_permutation_tests_accuracies = function(perm_test_dir="perm_test",n_perm_test=100){
  accuracies = c()
  nPerm = n_perm_test
  pb = txtProgressBar(min = 1,max = n_perm_test,style=3)
  for(i in 1:nPerm){
    avg_acc = read.table(paste(perm_test_dir,"/permutation_test_",i,"/all_res/test_avg_accuracy.txt",sep=""),header=TRUE,sep=",")
    avg_acc <- avg_acc[order(as.integer(substring(avg_acc$X,6))),]
    accuracies = cbind(accuracies,avg_acc$avg_acc)
    setTxtProgressBar(pb,i)
  }
  close(pb)
  return(accuracies)
}

evaluate_perm_test_results = function(real_accuracies, perm_test_accuracies,n_perm_test = n_perm_test){
  percs = c()
  for(i in 1:length(real_accuracies)){
    percs = c(percs,sum(real_accuracies[i] > perm_test_accuracies[i,]) / n_perm_test)
  }
  return(percs)
}


#### PERMUTATION TEST
n_perm_test = 500

##################################################################################
#   BEFORE THE PERMUTATION TEST THS SCRIPT stacking_level_1.py MUST BE EXECUTED  #
##################################################################################

#LOAD THE DATASET
load("../gbm_dataset.RData") 

##### THE PERMUTATION TEST WILL BE RUN ONLY ON THE GENES THAT COMES FROM THE DFP ANALYSIS ##### 
#LOAD THE DFP ANALYSIS RESULTS
DFP_genes = read.table("../discriminant_genes.txt",header = TRUE,sep=" ")
X = dataset[rownames(DFP_genes),]

####################### CREATE DATASETS FOR PERMUTATION TEST #######################

dimensions = get_pathway_dimensions("all")
create_dataset_for_permutation_test(X,perm_test_dir="perm_test",n_perm_test = n_perm_test,dimensions=dimensions)

################################################################################
#   HERE IN BETWEEN THE SCRIPT stacking_level_1_perm_test.py MUST BE EXECUTED  #
################################################################################

####################### RETRIEVE PERMUTATION TEST RESULTS #######################
# Average Accuracy of the pathway that contains only the genes that reach 80% of the cumulative score
# evaluated on the linear svm scores assigned to each gene
# all_78_avg_acc = read.table("allRed_res/test_avg_accuracy.txt",header=TRUE,sep=",")
# dimensions_red = get_pathway_dimensions("allRed")

# Average Accuracy of the original pathway 
load("goodPathway.Rdata")
test_avg_accuracy = read.table("tcgaPathways_res/test_avg_accuracy.txt",header=TRUE,sep=",")
all_avg_acc <- test_avg_accuracy[which(paste(test_avg_accuracy$X,".txt",sep="") %in% goodPathway),]
#all_avg_acc = read.table("all_res/test_avg_accuracy.txt",header=TRUE,sep=",")



dimensions_all = get_pathway_dimensions("allRed")
#dimensions_all <- dimensions_all[as.integer(row.names(all_avg_acc))]

#Reading accuracies obtained on the permuted pathways
#accuracies_red = get_permutation_tests_accuracies(perm_test_dir="perm_test_Red",n_perm_test=n_perm_test)
accuracies_all = get_permutation_tests_accuracies(perm_test_dir="perm_test",n_perm_test=n_perm_test)

#Evaluating the percentage of times that the original accuracy is higher that the one obtained on the permuted pathways
# perm_res_red = evaluate_perm_test_results(real_accuracies = all_78_avg_acc$avg_acc, perm_test_accuracies=accuracies_red,n_perm_test = n_perm_test)
perm_res_all = evaluate_perm_test_results(real_accuracies = all_avg_acc$avg_acc, perm_test_accuracies=accuracies_all,n_perm_test = n_perm_test)

M = cbind(all_avg_acc[,c("X","avg_acc")],perm_res_all,dimensions_all)#,
          # all_78_avg_acc[,c("avg_acc")],perm_res_red,dimensions_red)

colnames(M) = c("Pathway","Avg_Accuracy_All_Genes","Perm_Test_Result_All_Genes","Original_Pathway_Size")#,
                          # "Avg_Accuracy_Red_Genes","Perm_Test_Result_Red_Genes","Reduced_Pathway_Size")
View(M)
write.table(M,file="perm_test_res_all_pathways.txt",sep="\t",quote = FALSE)

MM = M[which(M$Perm_Test_Result_Red_Genes>0.95),]
View(MM)
write.table(MM,file="significant_pathways_after_perm_test.txt",sep="\t",quote = FALSE)
