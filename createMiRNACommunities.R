load("miRNA_communities_walk.Rdata")
load("tcga_trainingIdx.RData")
load("good_communities_idx.Rdata")

miRNA <- read.table("tcga_breast/miRNASeq.txt", header = T, sep = "\t")
miRNATraining <- miRNA[,trainingIdx]

good_communities <- good_communities_idx

folder <- "tcgaMiRNACommunities"

for (i in 1:length(good_communities)) {
  data <- t(miRNATraining[membership(bw_walktrap) == good_communities[i],])
  write.table(data, file = paste(folder, "/COMM_", good_communities[i], ".txt", sep = ""), quote = F)
}

tot <- 0
for (i in 1:length(good_communities)) {
  tot <- tot + sum(membership(bw_walktrap) == good_communities[i])
}