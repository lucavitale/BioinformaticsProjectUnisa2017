oxfordGenes <- read.table("tcga_breast/RNASeq_full.txt",sep="\t",header = T)
oxfordClass <- read.table("tcga_breast/patient_classes.txt",sep="\t",header = T)

best = FALSE
allPath = TRUE

if(allPath){
  load("allPathways.Rdata")
  pathways <- allPathways
}else if (best) {
  load("bestPathways.Rdata")
  pathways <- bestPathways
} else {
  load("reactomeFinal.Rdata")
  load("keggFinal.Rdata")
  #load("biocartaFinal.Rdata")
  pathways <- rbind(reactomeFinal,keggFinal)
}

load("tcga_trainingIdx.RData")

pathways <- pathways[,c("ID","geneID")]

names = row.names(pathways)

for(i in 1:nrow(pathways)){
  cur <- pathways[i,]
  genes <- strsplit(cur$geneID,"/")[[1]]
  data = t(oxfordGenes[genes,trainingIdx])
  write.table(data,file=paste("tcgaPathways/",names[i],".txt",sep=""),quote=F)
}

class <- as.data.frame(as.integer(oxfordClass[trainingIdx,]$x))
names(class) <- "x"
row.names(class) <- rownames(oxfordClass[trainingIdx,])
write.table(class,file="tcgaPathways/labels",quote=F)
