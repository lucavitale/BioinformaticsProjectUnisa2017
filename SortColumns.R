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
colPat <- c("caseID","ageAtDiagnosis","daysToBirth","daysToDeath","daysTolastFollowUp","morphology","primaryDiagnosis","siteOfResectionOrBiopsy","tissueOrOrganOfOrigin","tumorStage","vitalStatus","aliquotsSubmitterID","sampleType","submitterID","experimentalStrategy","fileID","fileName","stageClass")
names(miRNAPatientsFinal) <- colPat
names(RNAPatientsFinal) <- colPat

## Change row names
rownames(RNAPatientsFinal) <- RNAPatientsFinal$submitterID
rownames(miRNAPatientsFinal) <- miRNAPatientsFinal$submitterID 

save(RNAPatientsFinal,file = "RNAPatientsFinal.Rdata")
save(RNAFinal, file = "RNAFinal.Rdata")

save(miRNAPatientsFinal,file = "miRNAPatientsFinal.Rdata")
save(miRNAFinal,file = "miRNAFinal.Rdata")
