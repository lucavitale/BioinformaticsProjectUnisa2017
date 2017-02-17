load("RNANORM.RData")
load("RNAPatientsFinal")

load("reactomeFinal.Rdata")
load("keggFinal.Rdata")

load("trainingIdx.RData")


pathways <- rbind(reactomeFinal,keggFinal)
pathways <- pathways[,c("ID","geneID")]

names = row.names(pathways)

for(i in 1:nrow(pathways)){
  cur <- pathways[i,]
  genes <- strsplit(cur$geneID,"/")[[1]]
  data = t(rnanorm[genes,trainingIdx])
  write.table(data,file=paste("pathways/",names[i],".txt",sep=""),quote=F)
}

class <- as.data.frame(as.integer(RNAPatientsFinal[trainingIdx,]$stageClass))
names(class) <- "x"
row.names(class) <- rownames(RNAPatientsFinal[trainingIdx,])
write.table(class,file="pathways/labels",quote=F)
