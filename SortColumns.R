load("RNAPatientsFinal.Rdata")
load("RNAFinal.Rdata")

load("miRNAPatientsFinal.Rdata")
load("miRNAFinal.Rdata")

RNAFinal <- RNAFinal[,order(names(RNAFinal))]
miRNAFinal <- miRNAFinal[,order(names(miRNAFinal))]

RNAPatientsFinal <- RNAPatientsFinal[order(RNAPatientsFinal$cases_0_submitter_id),]
miRNAPatientsFinal <- miRNAPatientsFinal[order(miRNAPatientsFinal$cases_0_submitter_id),]

save(RNAPatientsFinal,file = "RNAPatientsFinal.Rdata")
save(RNAFinal, file = "RNAFinal.Rdata")

save(miRNAPatientsFinal,file = "miRNAPatientsFinal.Rdata")
save(miRNAFinal,file = "miRNAFinal.Rdata")
