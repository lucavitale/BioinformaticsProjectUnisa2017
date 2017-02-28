library(readr)

#average accuracy for each pathway
tac <- read.table("tcgaPathways_res/test_avg_accuracy.txt",sep=",",header = T)
#results of the enrichment analysis
load("allPathways.RData") 

# 
# CTD_diseases_pathways <- read_delim("~/bioinfo_luca_gioele/BioinformaticsProjectUnisa2017/CTD_diseases_pathways.tsv", 
#                                     "\t", escape_double = FALSE, trim_ws = TRUE)
# 
# breast_neoplasm_pw <- CTD_diseases_pathways[which(CTD_diseases_pathways$DiseaseName == "Breast Neoplasms"),]
# breast_neoplasm_pw <- breast_neoplasm_pw$PathwayName
# breast_neoplasm_pw <- unique(breast_neoplasm_pw)
# 
# sum(tac$X %in% breast_neoplasm_pw)
allPathways <- allPathways[order(allPathways$ID),]
tac2 = cbind(tac,allPathways[tac$X,c("ID","Count")])
colnames(tac2) = c("Pathway","Accuracy","Params","PathwayCompleteName","Size")
idx = grep(tac$X, pattern = "breast",ignore.case = TRUE)

colors = rep("black",nrow(tac2))
colors[idx] = "red"

plot(tac2$Accuracy,tac2$Size,xlab = "Accuracy",ylab="Size",col=colors)
abline(v = 0.8)
abline(v= 0.7)

idx_breast = grep(pattern = "breast",x = CTD_diseases_pathways$DiseaseName,ignore.case = TRUE)
submat = CTD_diseases_pathways[idx_breast,]
idx_2 = grep(pattern = "REACT",x = submat$PathwayID,ignore.case = TRUE)
View(submat[idx_2,])

goodPathway <- tac2$PathwayCompleteName[tac2$Accuracy>0.75]
save(goodPathway,file="goodPathway.Rdata")
