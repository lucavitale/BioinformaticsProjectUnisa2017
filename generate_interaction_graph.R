library(parallel)
library(igraph)
source("NetProj.R")

load("perm_test_results.Rdata")
MM <- MM[-which(MM$Pathway == "GO_REGULATION_OF_CELL_CYCLE_PHASE_TRANSITION"),]
pathwayPath <- paste("tcgaPathways/", MM$Pathway, ".txt", sep="")

pathwayData <- lapply(pathwayPath, FUN = function(path) {read.table(path, header = T, sep = " ", check.names = F)})

IS <- sapply(1:length(pathwayData), FUN = function(j) {
  sapply(1:length(pathwayData), FUN = function(i, j) {
    X <- as.numeric(as.matrix(pathwayData[[i]]))
    Y <- as.numeric(as.matrix(pathwayData[[j]]))
    abs(mean(X) - mean(Y)) / (sd(X) + sd(Y))
  }, j)
})

IS <- 1 - IS / max(IS)
q <- quantile(IS, probs = 0.8)
IS[IS <= q] <- 0
rownames(IS) <- MM$Pathway
colnames(IS) <- MM$Pathway
#g <- graph_from_adjacency_matrix(IS, mode = "undirected", weighted = T, diag = F)
#plot(g, vertex.size = 5, vertex.label.cex = 0.5, edge.width = E(g)$weight)

load("cor_miRNA_RNA.Rdata")

q <- c()
q[1] <- quantile(cor_miRNA_RNA[cor_miRNA_RNA < 0], probs = 0.1)
q[2] <- quantile(cor_miRNA_RNA[cor_miRNA_RNA >= 0], probs = 0.9)

par(mar = c(5, 4, 4, 2) + 0.1) # default value
plot(density(cor_miRNA_RNA))
abline(v = q)

cor_miRNA_RNA[cor_miRNA_RNA >= q[1] & cor_miRNA_RNA <= q[2]] <- 0
cor_miRNA_RNA <- t(cor_miRNA_RNA)

#g2 <- graph_from_adjacency_matrix(cor_miRNA_RNA, mode = "undirected", weighted = T)
#plot(g2, vertex.size = 5, vertex.label.cex = 0.5, edge.width = abs(E(g2)$weight), edge.col = ifelse(E(g2)$weight < 0, "red", "green"))

cl <- makeCluster(8)
clusterExport(cl, varlist = c("pathwayData","cor_miRNA_RNA"))
miRNAPathwayAdj <- parSapply(cl, 1:length(pathwayData), function(j) {
  genesInPathway <- colnames(pathwayData[[j]])
  sapply(1:nrow(cor_miRNA_RNA), function(i,genesInPathway){
    genesInMiRNA <- names(cor_miRNA_RNA[i,cor_miRNA_RNA[i,]>0])
    
    conting <- matrix(rep(0,4), nrow = 2)
    conting[1,1] <- length(intersect(genesInMiRNA, genesInPathway))
    conting[2,2] <- length(setdiff(colnames(cor_miRNA_RNA), union(genesInPathway, genesInMiRNA)))
    conting[1,2] <- length(setdiff(genesInMiRNA, genesInPathway))
    conting[2,1] <- length(setdiff(genesInPathway, genesInMiRNA))
    
    # if(sum(conting) != ncol(cor_miRNA_RNA)) {
    #   stop("There is a bug :(")
    # }
    
    res <- fisher.test(conting)
    res$p.value
  }, genesInPathway)
})
pvals <- matrix(p.adjust(miRNAPathwayAdj, method = "fdr"), byrow = F, nrow = nrow(cor_miRNA_RNA))
pvals_bool <- pvals <= 0.01

miRNAPathwayAdj <- parSapply(cl, 1:length(pathwayData), function(j) {
  genesInPathway <- colnames(pathwayData[[j]])
  sapply(1:nrow(cor_miRNA_RNA), function(i,genesInPathway){
    genesInMiRNA <- names(cor_miRNA_RNA[i,cor_miRNA_RNA[i,]>0])
    
    conting <- matrix(rep(0,4), nrow = 2)
    conting[1,1] <- length(intersect(genesInMiRNA, genesInPathway))
    conting[2,2] <- length(setdiff(colnames(cor_miRNA_RNA), union(genesInPathway, genesInMiRNA)))
    conting[1,2] <- length(setdiff(genesInMiRNA, genesInPathway))
    conting[2,1] <- length(setdiff(genesInPathway, genesInMiRNA))
    
    # if(sum(conting) != ncol(cor_miRNA_RNA)) {
    #   stop("There is a bug :(")
    # }
    
    res <- fisher.test(conting)
    res$estimate
  }, genesInPathway)
})

stopCluster(cl)

miRNAPathwayAdj <- miRNAPathwayAdj * pvals_bool
miRNAPathwayAdj[is.nan(miRNAPathwayAdj)] <- 0

miRNAPathwayAdj <- miRNAPathwayAdj / max(miRNAPathwayAdj)

colnames(miRNAPathwayAdj) <- MM$Pathway
rownames(miRNAPathwayAdj) <- rownames(cor_miRNA_RNA)

#miRNAPathwayAdjBinary <- miRNAPathwayAdj
#miRNAPathwayAdjBinary[miRNAPathwayAdjBinary > 0] <- 1
#miRNAmiRNAAdj <- counted.binary.net.proj(miRNAPathwayAdjBinary)
miRNAmiRNAAdj <- weighted.net.proj(miRNAPathwayAdj)
#miRNAmiRNAAdj <- miRNAmiRNAAdj/max(miRNAmiRNAAdj)
q2 <- quantile(miRNAmiRNAAdj[miRNAmiRNAAdj > 0], probs = 0.8)
miRNAmiRNAAdj[miRNAmiRNAAdj < q2] <- 0

g2 <- graph_from_adjacency_matrix(miRNAmiRNAAdj, mode = "undirected", weighted = T, diag = F)
bw <- cluster_edge_betweenness(g2, directed = F)
save(bw, file="miRNA_communities.Rdata")
bw_walktrap <- walktrap.community(g2)
save(bw_walktrap, file="miRNA_communities_walk.Rdata")


####################################################
#                                                  #
#   ADD miRNA GROUPS INSTEAD OF SINGLE miRNAs!!!   #
#                                                  #
####################################################
sup <- cbind(IS, t(miRNAPathwayAdj))
inf <- cbind(miRNAPathwayAdj, miRNAmiRNAAdj)
#inf <- cbind(miRNAPathwayAdj, matrix(0, nrow=nrow(miRNAmiRNAAdj), ncol = nrow(miRNAmiRNAAdj)))
adjFinal <- rbind(sup, inf)

#IndexPathwayToRemove <- which(rownames(adjFinal) =="GO_REGULATION_OF_CELL_CYCLE_PHASE_TRANSITION")
#adjFinal<- adjFinal[-IndexPathwayToRemove,-IndexPathwayToRemove]

g <- graph_from_adjacency_matrix(adjFinal, mode = "undirected", weighted = T, diag = F)
colors <- c(rep("red",nrow(MM)), rep("green",nrow(cor_miRNA_RNA)))

# removing disconnected vertices
# toremove <- which(degree(g) == max(degree(g)))
# colors <- colors[-toremove]
# g <- delete.vertices(g, toremove)

toremove <- which(degree(g)<1)
colors <- colors[-toremove]
g <- delete.vertices(g, toremove)

par(mar = c(0,0,0,0))
plot(g, vertex.size = 5, vertex.label.cex = 0.5, edge.width = abs(E(g)$weight), vertex.color = colors,layout=layout.davidson.harel)
