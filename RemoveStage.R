# Eliminate invalid stages and stages with too few patients
# Also remove patients with duplicated files

# RNA
RNAPatientsFil <- RNAPatients[which(RNAPatients$cases_0_samples_0_sample_type=="Primary Tumor"),]
RNAPatientsFil <- RNAPatientsFil[- which(RNAPatientsFil$stage %in% c("stage iv","not reported","stage x")),]
dupid <- RNAPatientsFil[which(duplicated(RNAPatientsFil$cases_0_case_id)),]$cases_0_case_id
asd <- RNAPatientsFil[which(RNAPatientsFil$cases_0_case_id %in% dupid),]
RNAPatientsFil <- RNAPatientsFil[-which(RNAPatientsFil$cases_0_samples_0_portions_0_analytes_0_aliquots_0_submitter_id %in% asd$cases_0_samples_0_portions_0_analytes_0_aliquots_0_submitter_id),]
save(RNAPatientsFil, file = "RNAPatientsFil.Rdata")

#miRNA
miRNAPatientsFil <- miRNAPatients[which(miRNAPatients$cases_0_samples_0_sample_type=="Primary Tumor"),]
miRNAPatientsFil <- miRNAPatientsFil[- which(miRNAPatientsFil$stage %in% c("stage iv","not reported","stage x")),]
dupid <- miRNAPatientsFil[which(duplicated(miRNAPatientsFil$cases_0_case_id)),]$cases_0_case_id
asd <- miRNAPatientsFil[which(miRNAPatientsFil$cases_0_case_id %in% dupid),]
miRNAPatientsFil <- miRNAPatientsFil[-which(miRNAPatientsFil$cases_0_samples_0_portions_0_analytes_0_aliquots_0_submitter_id %in% asd$cases_0_samples_0_portions_0_analytes_0_aliquots_0_submitter_id),]
save(miRNAPatientsFil, file = "miRNAPatientsFil.Rdata")