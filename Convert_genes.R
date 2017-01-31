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

gc <- read.table("gene_codes.txt",sep="\t")
gc$V1

entrez <- convertIDs(as.character(gc$V1), "ENSEMBL", "ENTREZID", org.Hs.eg.db, ifMultiple = "useFirst")
symbol <- convertIDs(as.character(gc$V1), "ENSEMBL", "SYMBOL", org.Hs.eg.db, ifMultiple = "useFirst")

notnullEntrez <- length(entrez) - sum(is.na(entrez))
notnullSymbol <- length(symbol) - sum(is.na(symbol))
