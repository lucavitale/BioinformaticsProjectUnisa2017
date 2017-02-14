rnasel <- RNAFinal[rowSums(RNAFinal)!=0,]
rnanorm =log2(rnasel+8)

range(rnanorm)
# stage i vs stage i
plot((rnanorm[,3]+rnanorm[,7])/2, (rnanorm[,3]-rnanorm[,7]),xlim=c(0,26))
v# stage iii vs stage i
plot((rnanorm[,4]+rnanorm[,7])/2, (rnanorm[,4]-rnanorm[,7]),xlim=c(0,26))
# stage iii vs stage i
plot((rnanorm[,4]+rnanorm[,3])/2, (rnanorm[,4]-rnanorm[,3]),xlim=c(0,26))