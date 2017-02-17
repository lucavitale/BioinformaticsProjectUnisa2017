library(ReactomePA)

load("paramList98.Rdata")
gene <- rownames(paramList$dfps)
if (!exists("RNAFinal")) {
  load("RNAFinal.Rdata")
}
universe <- rownames(RNAFinal)
ReactomeResults <- enrichPathway(gene, organism = "human", universe)
