for (i in 1:length(names(miRNAFinal))) {
  names(miRNAFinal)[i] <- as.character(miRNAPatientsFinal[which(names(miRNAFinal)[i]==miRNAPatientsFinal$cases_0_samples_0_portions_0_analytes_0_aliquots_0_submitter_id),]$cases_0_submitter_id)
}