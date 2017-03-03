# Library for bipartite network projections
# January 2013

# Function to implement a binary network projection from a bipartite heterogeneous network to a homogeneous network.
# Input:
# bip: adjacency matrix of a bipartite graph. Rows represent the set of nodes S, columns a set of different nodes T.
#      The matrix must be a binary matrix. 
# Output:
# the adjacency matrix of the nodes S obtained by binary network projection. The matrix has a number or rows and columns
# equal to nrow(bip) and may have rows/columns of zeros if the corresponding node have no common neighbhours with another 
# node in the original bipartite network

binary.net.proj <- function (bip) {
  d <- dim(bip);
  n.S <- d[1];
  n.T <- d[2];
  
  res<-apply(bip,2,function(x) return(which(x==1)));

  sets <- lapply(res,function(x) {
       if (length(x)>1) {
    	 y <- y2 <- combn(x,2);
		 y2[1,] <- y[2,];
		 y2[2,] <- y[1,];
		 z <- cbind(y,y2);
		 return(t(z));
	   }
	   return(NULL);
	   });

  if (is.null(rownames(bip)))
	rownames(bip)=1:n.S;

  proj<-matrix(integer(n.S*n.S), nrow=n.S);
  rownames(proj) <- colnames(proj) <- rownames(bip);

  for (i in 1:n.T) 
	if (!is.null(sets[[i]])) 
      proj[sets[[i]]] <- 1; 
  
  return(proj);
}




# Function to implement a counted binary network projection from a bipartite heterogeneous network to a homogeneous network.
# If two nodes shares n nodes in the bipartite network, the weight between them is n
# Input:
# bip: adjacency matrix of a bipartite graph. Rows represent the set of nodes S, columns a set of different nodes T.
#      The matrix must be a binary matrix. 
# Output:
# the adjacency matrix of the nodes S obtained by binary network projection. The matrix has a number or rows and columns
# equal to nrow(bip) and may have rows/columns of zeros if the corresponding node have no common neighbhours with another 
# node in the original bipartite network

counted.binary.net.proj <- function (bip) {
  d <- dim(bip);
  n.S <- d[1];
  n.T <- d[2];
  
  res<-apply(bip,2,function(x) return(which(x==1)));

  sets <- lapply(res,function(x) {
       if (length(x)>1) {
    	 y <- y2 <- combn(x,2);
		 y2[1,] <- y[2,];
		 y2[2,] <- y[1,];
		 z <- cbind(y,y2);
		 return(t(z));
	   }
	   return(NULL);
	   });

  if (is.null(rownames(bip)))
	rownames(bip)=1:n.S;

  proj<-matrix(integer(n.S*n.S), nrow=n.S);
  rownames(proj) <- colnames(proj) <- rownames(bip);

  for (i in 1:n.T) 
	if (!is.null(sets[[i]])) 
      proj[sets[[i]]] <- proj[sets[[i]]] + 1; 
	  
  return(proj);
}



# Function to implement a weighted network projection from a bipartite heterogeneous network to a homogeneous network.
# If two nodes shares a nodes in the bipartite network, the weight between them is the projected network is the
# average, max, min or the result of whatever supplied function between the weights of the edges in the bipartite network.
# Input:
# bip: adjacency matrix of a bipartite graph. Rows represent the set of nodes S, columns a set of different nodes T.
#      Bip is a real valued matrix 
# thresh : threshold below which the edges in the bipartite network are not considered (def: 0) 
# fun : function to compute the resulting weight in the projected network (def: mean). Other available functions are e.g. max, min
#       or any binary functions acting on a vector with 2 elements (i.e. the weights of the edges of the bipartite network).
# Output:
# the adjacency matrix of the nodes S obtained by weighted network projection. The matrix has a number or rows and columns
# equal to nrow(bip) and may have rows/columns of zeros if the corresponding node have no common neighbhours with another 
# node in the original bipartite network

weighted.net.proj <- function (bip, thresh=0, fun=mean) {
  d <- dim(bip);
  n.S <- d[1];
  n.T <- d[2];
  
  res<-apply(bip,2,function(x) return(which(x>thresh)));
  
  cnt <<- 0; # counter of the columns of bip
  
  sets <- lapply(res,function(x,fun) {
       cnt <<- cnt + 1;
       if (length(x)>1) {
    	 y <- y2 <- combn(x,2);
		 y2[1,] <- y[2,];
		 y2[2,] <- y[1,];
		 z <- t(cbind(y,y2));
		 w <- apply(z,1, function(xx,fun) {return(fun(c(bip[xx[1],cnt], bip[xx[2],cnt])))}, fun);
	     z <- cbind(z,w);
		 return(z);
	   }
	   return(NULL);
	   }, fun);

  if (is.null(rownames(bip)))
	rownames(bip)=1:n.S;

  proj<-matrix(numeric(n.S*n.S), nrow=n.S);
  rownames(proj) <- colnames(proj) <- rownames(bip);
  for (i in 1:n.T) 
	if (!is.null(sets[[i]])) {
	  part <- sets[[i]][,1:2];
      proj[part] <- pmax(proj[part], sets[[i]][,3]);	  
	}
	  
  return(proj);
}

