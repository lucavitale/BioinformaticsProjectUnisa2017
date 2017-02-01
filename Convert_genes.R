library(org.Hs.eg.db)

convertIDs <- function( ids, from, to, db, ifMultiple=c("putNA", "useFirst")) {
  stopifnot( inherits( db, "AnnotationDb" ) )
  ifMultiple <- match.arg( ifMultiple )
  suppressWarnings( selRes <- AnnotationDbi::select(
    db, keys=ids, keytype=from, columns=c(from,to) ) )
  
  if ( ifMultiple == "putNA" ) {
    duplicatedIds <- selRes[ duplicated( selRes[,1] ), 1 ]
    selRes <- selRes[ ! selRes[,1] %in% duplicatedIds, ]
  }
  
  return( selRes[ match( ids, selRes[,1] ), 2 ] )
}

load("RNATable.Rdata")
ensembles <- rownames(assay(RNATable,1))

entrez <- convertIDs(ensembles, "ENSEMBL", "ENTREZID", org.Hs.eg.db, ifMultiple = "useFirst")
#symbol <- convertIDs(as.character(gc$V1), "ENSEMBL", "SYMBOL", org.Hs.eg.db, ifMultiple = "useFirst")

## Number of genes
print(length(unique(entrez)))

##Only entrez genes with value
filtered <- assay(RNATable,1)[- which(is.na(entrez)),]

##Change rownames ensemble to entrezID
rownames(filtered) <- entrez[- which(is.na(entrez))]

