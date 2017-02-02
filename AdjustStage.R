## miRNA
miRNAPatients <- read.table("miRNA_patients.txt",sep="\t",header=T)
stage <- miRNAPatients$cases_0_diagnoses_0_tumor_stage

stage[which(stage=="stage ia")] <- "stage i"
stage[which(stage=="stage ib")] <- "stage i"

stage[which(stage=="stage iia")] <- "stage ii"
stage[which(stage=="stage iib")] <- "stage ii"

stage[which(stage=="stage iiia")] <- "stage iii"
stage[which(stage=="stage iiib")] <- "stage iii"
stage[which(stage=="stage iiic")] <- "stage iii"

miRNAPatients$stage <- stage

save(miRNAPatients, file="miRNA_patients.Rdata")


## RNA
RNAPatients <- read.table("RNA_patients.txt",sep="\t",header=T)
stage <- RNAPatients$cases_0_diagnoses_0_tumor_stage
unique(stage)

stage[which(stage=="stage ia")] <- "stage i"
stage[which(stage=="stage ib")] <- "stage i"

stage[which(stage=="stage iia")] <- "stage ii"
stage[which(stage=="stage iib")] <- "stage ii"

stage[which(stage=="stage iiia")] <- "stage iii"
stage[which(stage=="stage iiib")] <- "stage iii"
stage[which(stage=="stage iiic")] <- "stage iii"

RNAPatients$stage <- stage

save(RNAPatients, file="RNA_patients.Rdata")
