BRCA <-  read.table("BRCA.mRNAseq_RPKM_log2.txt",sep="\t",header=T) ## logaritmo

BRCA2  <- read.table("BRCA.mRNAseq_RPKM.txt",sep="\t",header=T) ## RPKM

rownames(BRCA2) <- BRCA2[,1]
BRCA2 <- BRCA2[,-1]

BRCA2sel <- BRCA2[rowSums(BRCA2)!=0,]
BRCA2norm <- log2(BRCA2sel+8)

range(BRCA2norm)
plot((BRCA2norm[,3]+BRCA2norm[,7])/2, (BRCA2norm[,3]-BRCA2norm[,7]),ylim=c(-5,5),xlim=c(0,20))
plot((BRCA2norm[,4]+BRCA2norm[,7])/2, (BRCA2norm[,4]-BRCA2norm[,7]),ylim=c(-5,5),xlim=c(0,20))
plot((BRCA2norm[,4]+BRCA2norm[,3])/2, (BRCA2norm[,4]-BRCA2norm[,3]),ylim=c(-5,5),xlim=c(0,20))


#### ROW COUNTS 

BRCA3  <- read.table("BRCA.mRNAseq_raw_counts.txt",sep="\t",header=T) ## counts per million
rownames(BRCA3) <- BRCA3[,1]
BRCA3 <- BRCA3[,-1]


BRCA3sel <- BRCA3[rowSums(BRCA3)!=0,]


libsizes <- colSums(BRCA3sel)
size.factor <- libsizes/exp(mean(log(libsizes)))
BRCA3sel <- t(t(BRCA3sel)/size.factor)

BRCA3norm <- log2(BRCA3sel+8)

range(BRCA3norm)
plot((BRCA3norm[,3]+BRCA3norm[,7])/2, (BRCA3norm[,3]-BRCA3norm[,7]),ylim=c(-5,5),xlim=c(0,20))
plot((BRCA3norm[,4]+BRCA3norm[,7])/2, (BRCA3norm[,4]-BRCA3norm[,7]),ylim=c(-5,5),xlim=c(0,20))
plot((BRCA3norm[,4]+BRCA3norm[,3])/2, (BRCA3norm[,4]-BRCA3norm[,3]),ylim=c(-5,5),xlim=c(0,20))
