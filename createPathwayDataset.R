oxfordGenes <- read.table("oxford/genes.txt",sep="\t",header = T)
oxfordClass <- read.table("oxford/patient_classes.txt",sep="\t",header = T)

best = TRUE

if (best) {
  load("bestPathways.Rdata")
  pathways <- bestPathways
} else {
  load("reactomeFinal.Rdata")
  load("keggFinal.Rdata")
  #load("biocartaFinal.Rdata")
  pathways <- rbind(reactomeFinal,keggFinal)
}

load("trainingIdx.RData")

pathways <- pathways[,c("ID","geneID")]

names = row.names(pathways)

for(i in 1:nrow(pathways)){
  cur <- pathways[i,]
  genes <- strsplit(cur$geneID,"/")[[1]]
  data = t(oxfordGenes[genes,trainingIdx])
  write.table(data,file=paste("oxfordPathways/",names[i],".txt",sep=""),quote=F)
}

class <- as.data.frame(as.integer(oxfordClass[trainingIdx,]$x))
names(class) <- "x"
row.names(class) <- rownames(oxfordClass[trainingIdx,])
write.table(class,file="oxfordPathways/labels",quote=F)
