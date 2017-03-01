a=0;
for (i in 1:10){
  n_perm_test = i*100
  
  all_avg_acc = read.table("tcgaPathways_res/test_avg_accuracy.txt",header=TRUE,sep=",")
  #all_avg_acc <- test_avg_accuracy[which(paste(test_avg_accuracy$X,".txt",sep="") %in% goodPathway),]
  #all_avg_acc = read.table("all_res/test_avg_accuracy.txt",header=TRUE,sep=",")
  
  dimensions_all = get_pathway_dimensions("allRed")
  #dimensions_all <- dimensions_all[as.integer(row.names(all_avg_acc))]
  
  #Reading accuracies obtained on the permuted pathways
  #accuracies_red = get_permutation_tests_accuracies(perm_test_dir="perm_test_Red",n_perm_test=n_perm_test)
  accuracies_all = get_permutation_tests_accuracies(perm_test_dir="perm_test",n_perm_test=n_perm_test)
  
  #Evaluating the percentage of times that the original accuracy is higher that the one obtained on the permuted pathways
  # perm_res_red = evaluate_perm_test_results(real_accuracies = all_78_avg_acc$avg_acc, perm_test_accuracies=accuracies_red,n_perm_test = n_perm_test)
  perm_res_all = evaluate_perm_test_results(real_accuracies = all_avg_acc$avg_acc, perm_test_accuracies=accuracies_all, n_perm_test = n_perm_test, perm_dimensions = dimensions, all_dimensions = dimensions_all)
  
  M = cbind(all_avg_acc[,c("X","avg_acc")],perm_res_all,dimensions_all)#,
  # all_78_avg_acc[,c("avg_acc")],perm_res_red,dimensions_red)
  
  colnames(M) = c("Pathway","Avg_Accuracy_All_Genes","Perm_Test_Result_All_Genes","Original_Pathway_Size")#,
  # "Avg_Accuracy_Red_Genes","Perm_Test_Result_Red_Genes","Reduced_Pathway_Size")
  a[i] = sum(M$Perm_Test_Result_All_Genes>0.95)

}
names(a) = seq(100,1000,100)
barplot(a)
