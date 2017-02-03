library(DFP)

load("RNAFinal.Rdata")
load("RNAPatientsFinal.Rdata")

phenoData <- new("AnnotatedDataFrame",data = RNAPatientsFinal)
esRNA <- ExpressionSet(as.matrix(RNAFinal),phenoData = phenoData)

skipFactor <- 3 # Factor to skip odd values
zeta <- 0.5 # Threshold used in the membership functions to label the float values with a discrete value
piVal <- 0.9 # Percentage of values of a class to determine the fuzzy patterns
overlapping <- 2 # Determines the number of discrete labels


mfs <- calculateMembershipFunctions(esRNA, skipFactor); mfs[[1]]
plotMembershipFunctions(esRNA, mfs, featureNames(esRNA)[1:2])


dvs <- discretizeExpressionValues(esRNA, mfs, zeta, overlapping); dvs[1:4,1:10]
showDiscreteValues(dvs, featureNames(esRNA)[1:10],c("healthy", "AML-inv"))

fps <- calculateFuzzyPatterns(esRNA, dvs, piVal, overlapping); fps[1:30,]
showFuzzyPatterns(fps, "healthy")[21:50]


dfps <- calculateDiscriminantFuzzyPattern(esRNA, fps); dfps[1:5,]
plotDiscriminantFuzzyPattern(dfps, overlapping)

