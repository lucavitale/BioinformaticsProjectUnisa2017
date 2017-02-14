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

## BRCA2norm adjustment
patientsBRCA2norm <- read.table("PatientsBRCA2norm",sep="\t",header=T)
stage <- patientsBRCA2norm$diagnoses_0_tumor_stage
unique(stage)

stage[which(stage=="stage ia")] <- "stage i"
stage[which(stage=="stage ib")] <- "stage i"

stage[which(stage=="stage iia")] <- "stage ii"
stage[which(stage=="stage iib")] <- "stage ii"

stage[which(stage=="stage iiia")] <- "stage iii"
stage[which(stage=="stage iiib")] <- "stage iii"
stage[which(stage=="stage iiic")] <- "stage iii"

patientsBRCA2norm$stage <- stage

patientsBRCA2norm <- patientsBRCA2norm[which(patientsBRCA2norm$stage %in% c("stage i","stage ii","stage iii")),]
patientsBRCA2norm <- droplevels(patientsBRCA2norm)

save(patientsBRCA2norm, file="BRCA2_RNA_patients.Rdata")
