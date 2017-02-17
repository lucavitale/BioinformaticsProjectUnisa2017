library(GSA)
gmtfile <- system.file("extdata", "msigdb.v5.2.entrez.gmt",package="clusterProfiler")
c2kegg <- read.gmt(gmtfile)

sum(gene %in% c2kegg$gene)

egmt <- enricher(gene, TERM2GENE=c2kegg,pvalueCutoff = 0.05)
#head(egmt)

#View(egmt@result[grep(pattern = "KEGG",x = rownames(egmt@result),ignore.case = F),])

reactome <- egmt@result[grep(pattern = "reactome",x = rownames(egmt@result),ignore.case = T),]
kegg <- egmt@result[grep(pattern = "kegg",x = rownames(egmt@result),ignore.case = T),]
biocarta <- egmt@result[grep(pattern = "biocarta",x = rownames(egmt@result),ignore.case = T),]