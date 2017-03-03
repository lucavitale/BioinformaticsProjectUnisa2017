load("miRNA_communities_walk.Rdata")
load("tcga_trainingIdx.RData")

miRNA <- read.table("tcga_breast/miRNASeq.txt", header = T, sep = "\t")
miRNATraining <- miRNA[,trainingIdx]

good_communities <- which(sizes(bw_walktrap) != 1)

folder <- "tcgaMiRNACommunities"

for (i in 1:length(good_communities)) {
  data <- t(miRNATraining[membership(bw_walktrap) == good_communities[i],])
  write.table(data, file = paste(folder, "/community_", good_communities[i], ".txt", sep = ""), quote = F)
}
