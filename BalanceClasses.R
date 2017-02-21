pathways2 <- list.files(path = "pathways/")

stage1 <- which(RNAPatientsFinal$stageClass[trainingIdx] == "stage i")
stage2 <- which(RNAPatientsFinal$stageClass[trainingIdx] == "stage ii")
stage3 <- which(RNAPatientsFinal$stageClass[trainingIdx] == "stage iii")

numToSelect <- min(length(stage1),length(stage2),length(stage3))

stage1r <- stage1[1:numToSelect]
stage2r <- stage2[1:numToSelect]
stage3r <- stage3[1:numToSelect]

index <- c(stage1r,stage2r,stage3r)
index <- sample(index,length(index),replace=F)
save(index, file="train_permutation_363_indices.Rdata")

for(i in (1:length(pathways2))){
  p <- read.table(paste("pathways/",pathways2[i],sep=""),header = T,sep="",fill = F)
  p <- p[index,]
  write.table(p,file=paste("pathways2/",pathways2[i],sep=""),quote=F)
}

l <- read.table(paste("pathways/",pathways2[which(pathways2=="labels")],sep=""),header = T,sep="",fill = F)  
l2 <- as.data.frame(l[index,])
rownames(l2) <- rownames(l)[index]
names(l2) <- "x"
write.table(l2,file=paste("pathways2/",pathways2[which(pathways2=="labels")],sep=""),quote=F)
