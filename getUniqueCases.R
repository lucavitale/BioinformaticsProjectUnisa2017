t <- read.table("file_matching.txt", header = TRUE, sep = "\t")
tRNA <- t[which(t$experimental_strategy=="RNA-Seq"),]
tmiRNA <- t[which(t$experimental_strategy=="miRNA-Seq"),]
tmerged <- intersect(unique(tRNA$cases_0_submitter_id), unique(tmiRNA$cases_0_submitter_id))
tRNA2 <- tRNA[which(tRNA$cases_0_submitter_id %in% tmerged),]
tmiRNA2 <- tmiRNA[which(tmiRNA$cases_0_submitter_id %in% tmerged),]
write.table(tRNA2, file = "RNA_patients.txt", sep = "\t")
write.table(tmiRNA2, file = "miRNA_patients.txt", sep = "\t")

RNA2 <- RNAByAliquot[,which(names(RNAByAliquot) %in% tRNA2$cases_0_samples_0_portions_0_analytes_0_aliquots_0_submitter_id)]
miRNA2 <- miRNAByAliquot[,which(names(miRNAByAliquot) %in% tmiRNA2$cases_0_samples_0_portions_0_analytes_0_aliquots_0_submitter_id)]