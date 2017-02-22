library(GSA)
library(clusterProfiler)

load("oxfordDFPResults/paramList60.Rdata")
gene <- rownames(paramList$dfps)
breast_genes <- read_csv("breast_genes.txt", col_names = FALSE)

breast_genes <- breast_genes$X1[breast_genes$X1 %in% rownames(paramList$fps)]

gene <- union(gene,c("BRCA1","BRCA2"))
gene <- union(gene,breast_genes)


gmtfile <- system.file("extdata", "msigdb.v5.2.symbols.gmt",package="clusterProfiler")
c2kegg <- read.gmt(gmtfile)

#sum(gene %in% c2kegg$gene)

egmt <- enricher(gene, TERM2GENE=c2kegg,pvalueCutoff = 0.05)
#head(egmt)

#View(egmt@result[grep(pattern = "KEGG",x = rownames(egmt@result),ignore.case = F),])

reactome <- egmt@result[grep(pattern = "reactome",x = rownames(egmt@result),ignore.case = T),]
kegg <- egmt@result[grep(pattern = "kegg",x = rownames(egmt@result),ignore.case = T),]
biocarta <- egmt@result[grep(pattern = "biocarta",x = rownames(egmt@result),ignore.case = T),]

#gmtfile <- system.file("extdata", "c6.all.v5.2.entrez.gmt",package="clusterProfiler")
#c6 <- read.gmt(gmtfile)
#c6 <- unique(c6$ont)


reactomeFinal <- reactome[which(reactome$Count >=10),]
keggFinal <- kegg[which(kegg$Count >=10),]
biocartaFinal <- biocarta[which(biocarta$Count >=10),]

#c6Final <- egmt@result[egmt@result$ID %in% c6 & egmt@result$Count >= 10,]

save(reactomeFinal,file="reactomeFinal.Rdata")
save(keggFinal,file="keggFinal.Rdata")
#save(biocartaFinal,file="biocartaFinal.Rdata")
#save(c6Final,file="c6Final.Rdata")

bestPathways <- egmt@result[1:30,]
bestPathways <- bestPathways[which(bestPathways$Count >= 10),]
save(bestPathways,file="bestPathways.Rdata")
