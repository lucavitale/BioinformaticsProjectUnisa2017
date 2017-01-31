t <- read.table("file_matching.txt", header = TRUE, sep = "\t")
tRNA <- t[which(t$experimental_strategy=="RNA-Seq"),]
tmiRNA <- t[which(t$experimental_strategy=="miRNA-Seq"),]
v1 <- unique(tRNA$cases_0_case_id)
v2 <- unique(tmiRNA$cases_0_case_id)
tmerged <- intersect(unique(tRNA$cases_0_case_id), unique(tmiRNA$cases_0_case_id))
