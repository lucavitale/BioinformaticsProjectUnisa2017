pathways2 <- list.files(path = "pathways/")

stage1 <- which(RNAPatientsFinal$stageClass[trainingIdx] == "stage i")
stage2 <- which(RNAPatientsFinal$stageClass[trainingIdx] == "stage ii")
stage3 <- which(RNAPatientsFinal$stageClass[trainingIdx] == "stage iii")

numToSelect <- min(length(stage1),length(stage2),length(stage3))

stage1r <- stage1[1:numToSelect]
stage2r <- stage2[1:numToSelect]
stage3r <- stage3[1:numToSelect]

for(i in (1:length(pathways2))){
  p <- read.table(paste("pathways/",pathways2[i],sep=""),header = T,sep="",fill = F)
  p <- p[c(stage1r,stage2r,stage3r),]
  write.table(p,file=paste("pathways2/",pathways2[i],sep=""),quote=F)
}
