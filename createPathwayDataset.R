oxfordGenes <- read.table("oxford/genes.txt",sep="\t",header = T)
oxfordClass <- read.table("oxford/patient_classes.txt",sep="\t",header = T)

load("reactomeFinal.Rdata")
load("keggFinal.Rdata")
#load("biocartaFinal.Rdata")

load("trainingIdx.RData")


pathways <- rbind(reactomeFinal,keggFinal)
pathways <- pathways[,c("ID","geneID")]

names = row.names(pathways)

for(i in 1:nrow(pathways)){
  cur <- pathways[i,]
  genes <- strsplit(cur$geneID,"/")[[1]]
  data = t(oxfordGenes[genes,trainingIdx])
  write.table(data,file=paste("oxfordPathways/",names[i],".txt",sep=""),quote=F)
}

class <- as.data.frame(as.integer(oxfordClass[trainingIdx,]$stageClass))
names(class) <- "x"
row.names(class) <- rownames(oxfordClass[trainingIdx,])
write.table(class,file="pathways/labels",quote=F)
