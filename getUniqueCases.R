t <- read.table("file_matching.txt", header = TRUE, sep = "\t")
tRNA <- t[which(t$experimental_strategy=="RNA-Seq"),]
tmiRNA <- t[which(t$experimental_strategy=="miRNA-Seq"),]
tmerged <- intersect(unique(tRNA$cases_0_submitter_id), unique(tmiRNA$cases_0_submitter_id))
