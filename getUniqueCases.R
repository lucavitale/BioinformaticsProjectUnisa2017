t <- read.table("file_matching.txt", header = TRUE, sep = "\t")
tRNA <- t[which(t$experimental_strategy=="RNA-Seq"),]
tmiRNA <- t[which(t$experimental_strategy=="miRNA-Seq"),]
tmerged <- intersect(unique(tRNA$cases_0_submitter_id), unique(tmiRNA$cases_0_submitter_id))
tRNA2 <- tRNA[which(tRNA$cases_0_submitter_id %in% tmerged),]
tmiRNA2 <- tmiRNA[which(tRNA$cases_0_submitter_id %in% tmerged),]
write.table(tRNA2, file = "RNA_patients.txt", sep = "\t")
write.table(tmiRNA2, file = "miRNA_patients.txt", sep = "\t")
