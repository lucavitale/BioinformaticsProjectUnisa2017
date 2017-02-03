load("RNAPatientsFinal.Rdata")
load("RNAFinal.Rdata")

load("miRNAPatientsFinal.Rdata")
load("miRNAFinal.Rdata")

## Order
RNAFinal <- RNAFinal[,order(names(RNAFinal))]
miRNAFinal <- miRNAFinal[,order(names(miRNAFinal))]

RNAPatientsFinal <- RNAPatientsFinal[order(RNAPatientsFinal$cases_0_submitter_id),]
miRNAPatientsFinal <- miRNAPatientsFinal[order(miRNAPatientsFinal$cases_0_submitter_id),]

## Change columns names
colPat <- c("caseID","sampleType","daysTolastFollowUp","fileName","siteOfResectionOrBiopsy","fileID","morphology","primaryDiagnosis","submitterID","vitalStatus","aliquotsSubmitterID","daysToDeath","daysToBirth","ageAtDiagnosis","tissueOrOrganOfOrigin","tumorStage","experimentalStrategy","stageClass")
names(miRNAPatientsFinal) <- colPat
names(RNAPatientsFinal) <- colPat

save(RNAPatientsFinal,file = "RNAPatientsFinal.Rdata")
save(RNAFinal, file = "RNAFinal.Rdata")

save(miRNAPatientsFinal,file = "miRNAPatientsFinal.Rdata")
save(miRNAFinal,file = "miRNAFinal.Rdata")
