## Query miRNA-target interactions by combining expression data and putative miRNA-target
## interactions
querymiRTargetbinding <- function(ExpData, miRTarget, type = c("all", "miRNA", "target")) {
  
  colnames(ExpData) <- gsub("\\.", "-", colnames(ExpData))
  ExpDataNames <- colnames(ExpData)
  miRTarget <- as.matrix(miRTarget)
  if (type == "all") {
    miRTargetCandidate <- miRTarget[intersect(which(miRTarget[, 1] %in% ExpDataNames),
                                              which(miRTarget[, 2] %in% ExpDataNames)), ]
  } else if (type == "miRNA") {
    miRTargetCandidate <- miRTarget[which(miRTarget[, 1] %in% ExpDataNames), ]
  } else if (type == "target") {
    miRTargetCandidate <- miRTarget[which(miRTarget[, 2] %in% ExpDataNames), ]
  }
  
  return(miRTargetCandidate)
}


## Internal functions (calCMI, combpvalue, predCor) of hermes method are from the website:
## http://califano.c2b2.columbia.edu/hermes Copyright Columbia University in the City of New
## York. You may not use this file except in compliance with the License. You may obtain a copy
## of the License at http://califano.c2b2.columbia.edu/hermes-license

## calCMI: Calculate conditional mutual information
calCMI <- function(dmat) {

    # rank data
    N <- nrow(dmat)
    cdim <- ncol(dmat)
    idat <- apply(dmat, 2, order)

    ydat <- idat
    sapply(seq_len(cdim), function(i) ydat[idat[, i], i] <- seq_len(N))

    ddim <- 2^cdim
    dim2 <- 2 * cdim

    # init
    rm(dmat, idat)
    poc <- vector()
    kon <- vector()
    cmi <- 0
    npar <- 1
    poc[1] <- 1
    kon[1] <- N
    poradi <- seq_len(N)
    NN <- matrix(0, nrow = 1, ncol = ddim)
    marg <- matrix(0, nrow = 8 * ddim, ncol = dim2)
    marg[1, ] <- c(rep(1, cdim), rep(N, cdim))
    Imm <- matrix(c(0, 1), ncol = 1)  # Binary matrix for all combinations
    for (d in seq(2, cdim)) {
        Imm <- rbind(cbind(matrix(0, nrow = nrow(Imm), ncol = 1), Imm),
            cbind(matrix(1, nrow = nrow(Imm), ncol = 1), Imm))
    }

    chi2 <- c(0, 7.81, 13.9, 25, 42)
    run <- 0

    # partition
    while (npar > 0) {
        run <- run + 1
        apoc <- poc[npar]
        akon <- kon[npar]
        apor <- poradi[seq(apoc, akon)]
        Nex <- length(apor)
        margave <- floor((marg[npar, seq_len(cdim)] + marg[npar, seq(cdim + 1, dim2)])/2)
        J <- (ydat[apor, ] <= (matrix(1, nrow = Nex, ncol = 1) %*% margave)) * 1

        cI <- matrix(0, nrow = Nex, ncol = ddim)
        amarg <- matrix(1, nrow = ddim, ncol = 1) %*%
            subset(marg, subset = seq_len(nrow(marg)) %in% npar)

        for (d in seq_len(ddim)) {
            cI[, d] <- matrix(1, nrow = Nex, ncol = 1)
            for (k in seq_len(cdim)) {
                if (Imm[d, k] == 1) {
                  cI[, d] <- 1 * (cI[, d] & (!J[, k]))
                  amarg[d, k] <- margave[k] + 1
                } else {
                  cI[, d] <- 1 * (cI[, d] & J[, k])
                  amarg[d, k + cdim] <- margave[k]
                }
            }
        }

        NN <- colSums(cI)
        tst <- ddim * sum((NN - Nex/ddim * matrix(1, nrow = 1, ncol = ddim))^2)/Nex

        if ((tst > chi2[cdim]) | (run == 1)) {
            # decide partition or not
            npar <- npar - 1
            for (ind in seq_len(ddim)) {
                if (NN[ind] > ddim) {
                  npar <- npar + 1
                  akon <- apoc + NN[ind] - 1
                  poc[npar] <- apoc
                  kon[npar] <- akon
                  marg[npar, ] <- amarg[ind, ]
                  poradi[seq(apoc, akon)] <- apor[which(cI[, ind] != 0, arr.ind = TRUE)]
                  apoc <- akon + 1
                } else {
                  if (NN[ind] > 0) {
                    Nxx <- apply(amarg[ind, seq(cdim + 1, dim2)] - amarg[ind, seq_len(cdim)] + matrix(1,
                      nrow = 1, ncol = cdim), 2, prod)

                    Nz <- amarg[ind, 6] - amarg[ind, 3] + 1
                    Jx <- 1 * ((ydat[, 1] >= amarg[ind, 1]) & (ydat[, 1] <= amarg[ind, 4]))
                    Jy <- 1 * ((ydat[, 2] >= amarg[ind, 2]) & (ydat[, 2] <= amarg[ind, 5]))
                    Jz <- 1 * ((ydat[, 3] >= amarg[ind, 3]) & (ydat[, 3] <= amarg[ind, 6]))

                    Nxz <- sum(1 * (Jx & Jz))
                    Nyz <- sum(1 * (Jy & Jz))
                    cond <- (NN[ind] * Nz)/(Nxz * Nyz)
                    if (is.infinite(cond))
                      cond <- 1
                    if (cond == 0)
                      cond <- 1
                    cmi <- cmi + NN[ind] * log(cond)
                  }
                }
            }
        } else {

            Nxx <- apply(marg[npar, seq(cdim + 1, dim2)] - marg[npar, seq_len(cdim)] + matrix(1,
                nrow = 1, ncol = cdim), 2, prod)

            Nz <- marg[npar, 6] - marg[npar, 3] + 1
            Jx <- 1 * ((ydat[, 1] >= marg[npar, 1]) & (ydat[, 1] <= marg[npar, 4]))
            Jy <- 1 * ((ydat[, 2] >= marg[npar, 2]) & (ydat[, 2] <= marg[npar, 5]))
            Jz <- 1 * ((ydat[, 3] >= marg[npar, 3]) & (ydat[, 3] <= marg[npar, 6]))

            Nxz <- sum(1 * (Jx & Jz))
            Nyz <- sum(1 * (Jy & Jz))
            cond <- (Nex * Nz)/(Nxz * Nyz)
            if (is.infinite(cond))
                cond <- 1
            if (cond == 0)
                cond <- 1
            cmi <- cmi + Nex * log(cond)
            npar <- npar - 1
        }
    }
    # normalize
    cmi <- cmi/N
    return(cmi)
}

## combpvalue: Combine p-values using Fisher's method
combpvalue <- function(p_values) {

    # calculate chi-square statistic and combined p-value
    Q <- -2 * sum(log(p_values))
    Degree <- length(p_values)
    Combp <- pchisq(Q, 2 * Degree, lower.tail = FALSE)
    return(Combp)
}

## predCor: Predict competing endogenous RNAs via evidence for competition for miRNA regulation
## based on conditional mutual information (CMI) or partial pearson correlation (PPC)
predCor <- function(expr, num_perm, method = c("CMI", "PPC")) {

    exprRownames <- rownames(expr)
    expr <- as.matrix(expr)

    # generate random sequence
    num_sample <- ncol(expr)
    perm_seq <- t(sapply(seq_len(num_perm), function(i) sample(num_sample)))

    # identify mediator candidate
    num_cand <- nrow(expr) - 2

    # evaluate significance of triplet
    tri_id <- matrix(0, nrow = 2 * num_cand, ncol = 3)
    tri_cor <- matrix(0, nrow = 2 * num_cand, ncol = 1)
    tri_pval <- rep(1, 2 * num_cand)

    for (i in seq_len(2)) {
        for (j in seq_len(num_cand)) {

            idx_tar <- i  # target index
            idx_miR <- j + 2  # miRNA index
            idx_mod <- 2/i  # modulator index
            idx_tri <- (i - 1) * num_cand + j  # triplet index

            tri_id[idx_tri, ] <- c(exprRownames[idx_tar], exprRownames[idx_miR], exprRownames[idx_mod])

            # calculate Cor( target ; miRNA | modulator )
            data <- t(expr[c(idx_tar, idx_miR, idx_mod), ])

            # construct null distribution
            rand_exp <- lapply(seq_len(num_perm), function(i) expr[idx_mod, ][perm_seq[i, ]])
            nulldata <- lapply(seq_len(num_perm), function(i) cbind(t(expr[c(idx_tar, idx_miR),
                ]), rand_exp[[i]]))

            if (method == "CMI") {
                tri_cor[idx_tri, ] <- calCMI(data)
                null <- unlist(lapply(seq_len(num_perm), function(i) calCMI(nulldata[[i]])))
            } else if (method == "PPC") {
                tri_cor[idx_tri, ] <- corpcor::pcor.shrink(data, verbose = FALSE)[1, 2]
                null <- unlist(lapply(seq_len(num_perm), function(i) corpcor::pcor.shrink(nulldata[[i]],
                  verbose = FALSE)[1, 2]))
            }

            # calculate p-value
            tri_pval[idx_tri] <- max(1, sum(tri_cor[idx_tri, ] <= null))/num_perm
        }
    }

    # evaluate significance of interaction
    tri_idx <- order(tri_pval)
    tri_pval <- tri_pval[tri_idx]

    tri_id <- tri_id[tri_idx, ]
    tri_cor <- tri_cor[tri_idx]
    pcomb <- as.vector(sapply(seq_along(tri_pval), function(i) combpvalue(tri_pval[seq_len(i)])[1]))

    # identify final mediators
    min_pval <- min(pcomb)
    return(min_pval)
}

## Internal functions (parMM, graphWeights, recommendation, dtHybrid) of cernia method are from
## the website: https://github.com/dsardina/cernia Copyright 2016 Rosalba Giugno Licensed under
## the Apache License, Version 2.0 (the 'License'); you may not use this file except in
## compliance with the License. You may obtain a copy of the License at
## http://www.apache.org/licenses/LICENSE-2.0

## Unless required by applicable law or agreed to in writing, software distributed under the
## License is distributed on an 'AS IS' BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND,
## either express or implied. See the License for the specific language governing permissions and
## limitations under the License.

## Matrix multiplication in parallel
parMM <- function(cl, A, B) {

    if (!all(is.na(cl)) && is.object(cl)) {

        nA <- nrow(A)
        ncl <- length(cl)
        # Split an indice equally
        i <- seq_len(nA)
        if (ncl == 1) {
            splitIndices <- i
        } else {
            fuzz <- min((nA - 1)/1000, 0.4 * nA/ncl)
            breaks <- seq(1 - fuzz, nA + fuzz, length = ncl + 1)
            splitIndices <- structure(split(i, cut(i, breaks)), names = NULL)
        }

        # Split a matrix equally according to the rows of the matrix
        splitRows <- lapply(splitIndices, function(i) A[i, , drop = FALSE])

        fun <- get(as.character("rbind"), envir = .GlobalEnv, mode = "function")
        R <- do.call("fun", lapply(clusterApply(cl = cl, x = splitRows, get("%*%"), B), enquote))
    } else {
        R <- A %*% B
    }
    return(R)
}

## The first step of DT-Hybrid recommendation algorithm: generating the weights for each pair of
## target nodes
graphWeights <- function(n, m, A, lambda = 0.5, alpha = 0.5, S = NA, S1 = NA, cl = NA) {

    if (nrow(A) != n || ncol(A) != m) {
        stop("The matrix A should be an n by m matrix.")
    }

    has.similarity <- (!all(is.na(S)) && is.matrix(S) && !all(is.na(S1)) && is.matrix(S1))

    if (has.similarity) {
        if (nrow(S1) != m || ncol(S1) != m) {
            stop("The matrix S1 should be an m by m matrix.")
        }
        if (nrow(S) != n || ncol(S) != n) {
            stop("The matrix S should be an n by n matrix.")
        }
    }

    Ky <- diag(1/colSums(A))
    Ky[is.infinite(Ky) | is.na(Ky)] <- 0  #BugFix: 1/0=Infinite replaced with 0

    kx <- rowSums(A)
    Nx <- 1/(matrix(kx, nrow = n, ncol = n, byrow = TRUE)^(lambda) * matrix(kx, nrow = n, ncol = n,
        byrow = FALSE)^(1 - lambda))

    Nx[is.infinite(Nx) | is.na(Nx)] <- 0  #BugFix: 1/0=Infinite replaced with 0
    kx[is.infinite(kx) | is.na(kx)] <- 0  #BugFix: 1/0=Infinite replaced with 0

    W <- t(parMM(cl, A, Ky))
    W <- parMM(cl, A, W)
    W <- Nx * W
    rownames(W) <- rownames(A)
    colnames(W) <- rownames(A)

    if (has.similarity) {
        X5 <- parMM(cl, A, S1)
        X6 <- parMM(cl, X5, t(A))
        X7 <- parMM(cl, A, matrix(1, nrow = m, ncol = m))
        X8 <- parMM(cl, X7, t(A))
        S2 <- X6/X8
        W <- W * (1 + (alpha * S) + ((1 - alpha) * S2))
    }

    W[is.nan(W)] <- 0  #This should never happen
    return(W)
}

## The second step of DT-Hybrid recommendation algorithm: generating ecommendation scores of each
## RNA-RNA pair
recommendation <- function(A, lambda = 0.5, alpha = 0.5, S = NA, S1 = NA, cl = NA) {

    n <- nrow(A)
    m <- ncol(A)
    W <- graphWeights(n = n, m = m, A = A, lambda = lambda, alpha = alpha, S = S, S1 = S1, cl = cl)

    R <- parMM(cl, W, A)
    return(R)
}

## Make projection from bipartite network using DT-hybrid sources
dtHybrid <- function(miRTarget) {

    # Extract miRs and their targets
    mir <- unique(miRTarget[, 1])
    tar <- unique(miRTarget[, 2])

    # Create the matrix of the miRTarget
    A <- matrix(nrow = length(tar), ncol = length(mir), data = 0)
    colnames(A) <- mir
    rownames(A) <- tar

    for (i in seq_len(nrow(miRTarget))) {
        A[which(tar %in% as.character(miRTarget[i, 2])),
            which(mir %in% as.character(miRTarget[i, 1]))] <- 1
    }

    # Make projection from bipartite network using DT-hybrid sources
    cl <- makeCluster(detectCores() - 2)
    M <- recommendation(A, cl = cl)
    W <- graphWeights(nrow(M), ncol(M), M, cl = cl)

    stopCluster(cl)

    return(W)
}

## Internal function cluster from ProNet package 
## (https://github.com/cran/ProNet) with GPL-2 license.
cluster <- function(graph, method="MCL", expansion = 2, inflation = 2, 
                  hcmethod = "average", directed = FALSE, outfile = NULL, ...)
{
  #method<-match.arg(method)
  if(method=="FN"){
	  graph <- simplify(graph)
	  fc <- fastgreedy.community(graph, merges = TRUE, modularity = TRUE)
	  membership <- membership(fc)
	  if(!is.null(V(graph)$name)){
              names(membership) <- V(graph)$name
	  }	  
	  if(!is.null(outfile)){
	      cluster.save(cbind(names(membership),membership),outfile=outfile)		
	  }else{
	      return(membership)
	  }
  }else if(method=="LINKCOMM"){
	  edgelist <- get.edgelist(graph)
	  if(!is.null(E(graph)$weight)){
              edgelist <- cbind(edgelist,E(graph)$weight)
	  }
	  lc <- getLinkCommunities(edgelist,plot=FALSE,directed=directed,hcmethod=hcmethod)
	  if(!is.null(outfile)){
		  cluster.save(lc$nodeclusters,outfile=outfile)		  
	  }else{		  
		  return(lc$nodeclusters)
	  }
  }else if(method=="MCL"){
        adj <- matrix(rep(0,length(V(graph))^2),nrow=length(V(graph)),ncol=length(V(graph)))
        for(i in seq_along(V(graph))){
            neighbors <- neighbors(graph,v=V(graph)$name[i],mode="all")
            j <- match(neighbors$name,V(graph)$name,nomatch=0)
            adj[i,j] = 1
        }
        lc <- mcl(adj,addLoops=TRUE,expansion=expansion,inflation=inflation,allow1=TRUE,max.iter=100,ESM=FALSE)
        lc$name <- V(graph)$name
        lc$Cluster <- lc$Cluster
        
        if(!is.null(outfile)){
            cluster.save(cbind(lc$name,lc$Cluster),outfile=outfile)		
        }else{
        result <- lc$Cluster
        names(result) <- V(graph)$name
        return(result)
        }
  }else if(method=="MCODE"){
	  compx <- mcode(graph,vwp=0.9,haircut=T,fluff=T,fdt=0.1)
	  index <- which(!is.na(compx$score))
	  membership <- rep(0,vcount(graph))
	  for(i in seq_along(index)){
	      membership[compx$COMPLEX[[index[i]]]]<-i
	  }
	      if(!is.null(V(graph)$name)) names(membership)<-V(graph)$name	  
	  if(!is.null(outfile)){
		  cluster.save(cbind(names(membership),membership),outfile=outfile)
		  invisible(NULL)
	  }else{
		  return(membership)
	  }
  }
}

## Internal function cluster.save from ProNet package
## (https://github.com/cran/ProNet) with GPL-2 license.
cluster.save <- function(membership, outfile){
	wd <- dirname(outfile)
	wd <- ifelse(wd==".",paste(wd,"/",sep=""),wd)
	filename <- basename(outfile)
	if((filename=="")||(grepl(":",filename))){
		filename <- "membership.txt"
	}else if(grepl("\\.",filename)){
		filename <- sub("\\.(?:.*)",".txt", filename)
	}
	write.table(membership,file=paste(wd,filename,sep="/"),
              row.names=FALSE,col.names=c("node","cluster"),quote =FALSE)
}

## Internal function mcode.vertex.weighting from ProNet package
## (https://github.com/cran/ProNet) with GPL-2 license.
mcode.vertex.weighting<-function(graph, neighbors){	
	stopifnot(is.igraph(graph))  
	weight <- lapply(seq_len(vcount(graph)),
                 function(i){
		              subg<-induced.subgraph(graph,neighbors[[i]])
		              core<-graph.coreness(subg)
		              k<-max(core)
				          ### k-coreness
				          kcore<-induced.subgraph(subg,which(core==k))
				          if(vcount(kcore)>1){
					          if(any(is.loop(kcore))){
						          k*ecount(kcore)/choose(vcount(kcore)+1,2)						
					          }else{
						          k*ecount(kcore)/choose(vcount(kcore),2)
					          }
				          }else{
                                             0
				          }
				         }
			 )
  
	return(unlist(weight))
}

## Internal function mcode.find.complex from ProNet package
## (https://github.com/cran/ProNet) with GPL-2 license.
mcode.find.complex <- function(neighbors, neighbors.indx, vertex.weight,
                             vwp, seed.vertex, seen)
{
    
    res<-.C("complex",as.integer(neighbors),as.integer(neighbors.indx),
          as.single(vertex.weight),as.single(vwp),as.integer(seed.vertex),
          seen=as.integer(seen),COMPLEX=as.integer(rep(0,length(seen))), PACKAGE = "miRspongeR"
          )
	
	  return(list(seen=res$seen,COMPLEX=which(res$COMPLEX!=0)))
}

## Internal function mcode.find.complexex from ProNet package
## (https://github.com/cran/ProNet) with GPL-2 license.
mcode.find.complexex <- function(graph, neighbors, vertex.weight, vwp)
{
	seen<-rep(0,vcount(graph))

	neighbors<-lapply(neighbors,function(item){item[-1]})
	neighbors.indx<-cumsum(unlist(lapply(neighbors,length)))
	
	neighbors.indx<-c(0,neighbors.indx)
	neighbors<-unlist(neighbors)-1
	
	COMPLEX<-list()
	n<-1
        w.order<-order(vertex.weight,decreasing=TRUE)
	for(i in w.order){
		if(!(seen[i])){
			res<-mcode.find.complex(neighbors,neighbors.indx,vertex.weight,vwp,i-1,seen)
			if(length(res$COMPLEX)>1){
				COMPLEX[[n]]<-res$COMPLEX
				seen<-res$seen
				n<-n+1
			}
		}
	}	
	rm(neighbors)	
	return(list(COMPLEX=COMPLEX,seen=seen))
}

## Internal function mcode.fluff.complex from ProNet package
## (https://github.com/cran/ProNet) with GPL-2 license.
mcode.fluff.complex <- function(graph, vertex.weight, fdt=0.8, complex.g, seen)
{
	seq_complex.g<-seq_along(complex.g)
	for(i in seq_complex.g){
	    node.neighbor<-unlist(neighborhood(graph,1,complex.g[i]))
	    if(length(node.neighbor)>1){
                subg<-induced.subgraph(graph,node.neighbor)
            if(graph.density(subg, loops=FALSE)>fdt){
                complex.g<-c(complex.g,node.neighbor)
       }
		 }
	}
  
	return(unique(complex.g))
}

## Internal function mcode.post.process from ProNet package
## (https://github.com/cran/ProNet) with GPL-2 license.
mcode.post.process<-function(graph, vertex.weight, haircut, fluff, fdt=0.8,
                             set.complex.g, seen)
{
	indx<-unlist(lapply(set.complex.g,
                      function(complex.g){
		          if(length(complex.g)<=2)
			      0
			  else
		              1
			  }
		    ))
	set.complex.g<-set.complex.g[indx!=0]
	set.complex.g<-lapply(set.complex.g,
                        function(complex.g){
			    coreness<-graph.coreness(induced.subgraph(graph,complex.g))						
			    if(fluff){
				complex.g<-mcode.fluff.complex(graph,vertex.weight,fdt,complex.g,seen)
			    if(haircut){
				## coreness needs to be recalculated
				coreness<-graph.coreness(induced.subgraph(graph,complex.g))
				complex.g<-complex.g[coreness>1]
				}
				}else if(haircut){
				complex.g<-complex.g[coreness>1]
				}
				return(complex.g)
				})
	set.complex.g<-set.complex.g[lapply(set.complex.g,length)>2]
	return(set.complex.g)
}

## Internal function mcode from ProNet package
## (https://github.com/cran/ProNet) with GPL-2 license.
mcode <- function(graph, vwp=0.5, haircut=FALSE, fluff=FALSE, fdt=0.8, loops=TRUE)
{
	stopifnot(is.igraph(graph))
	if(vwp>1 | vwp <0){
            stop("vwp must be between 0 and 1")
	}
	if(!loops){
            graph<-simplify(graph,remove.multiple=FALSE,remove.loops=TRUE)
	}
	neighbors<-neighborhood(graph,1)
	W<-mcode.vertex.weighting(graph,neighbors)
	res<-mcode.find.complexex(graph,neighbors=neighbors,vertex.weight=W,vwp=vwp)
	COMPLEX<-mcode.post.process(graph,vertex.weight=W,haircut=haircut,fluff=fluff,
                              fdt=fdt,res$COMPLEX,res$seen)		
	score<-unlist(lapply(COMPLEX,
                       function(complex.g){
		           complex.g<-induced.subgraph(graph,complex.g)
			   if(any(is.loop(complex.g)))
			   score<-ecount(complex.g)/choose(vcount(complex.g)+1,2)*vcount(complex.g)
			   else
			   score<-ecount(complex.g)/choose(vcount(complex.g),2)*vcount(complex.g)
			   return(score)
			}
		    ))
	order_score<-order(score,decreasing=TRUE)
	return(list(COMPLEX=COMPLEX[order_score],score=score[order_score]))
}

## Utility methods for identifying miRNA sponge interactions For input expression data, the
## columns are genes and the rows are samples.  For input miRTarget, the miRNA-target
## interactions could be miRNA-mRNA, miRNA-lncRNA, miRNA-circRNA, miRNA-pseudogene, etc.  For
## input mres, each row contains five elements: Mirna, Target, energy, gap_l, gap_r.

## 1. miRHomology
# Original version
miRHomology <- function(miRTarget, minSharedmiR = 3, padjustvaluecutoff = 0.01, padjustmethod = "BH") {
  
  miRTarget <- as.matrix(miRTarget)
  
  m1 <- nrow(miRTarget)
  n1 <- ncol(miRTarget)
  
  miR <- miRTarget[, 1]
  tar <- miRTarget[, 2]
  
  miRSym <- unique(miR)
  targetSym <- unique(tar)
  
  m2 <- length(targetSym)
  
  # Initialize variables
  ceRInt <- matrix(NA, m2 * (m2 - 1)/2, 2)
  C <- matrix(NA, m2 * (m2 - 1)/2, 2)
  
  for (i in seq_len(m2 - 1)) {
    for (j in seq(i + 1, m2)) {
      
      Interin1 <- miRTarget[which(miRTarget[, 2] %in% targetSym[i]), 1]
      Interin2 <- miRTarget[which(miRTarget[, 2] %in% targetSym[j]), 1]
      
      M1 <- length(Interin1)
      M2 <- length(Interin2)
      M3 <- length(intersect(Interin1, Interin2))
      M4 <- length(miRSym)
      M5 <- 1 - phyper(M3 - 1, M2, M4 - M2, M1)
      
      if (M3 >= minSharedmiR & M5 < padjustvaluecutoff) {
        
        ceRInt[(i - 1) * m2 + j - sum(seq_len(i)), 1] <- targetSym[i]
        ceRInt[(i - 1) * m2 + j - sum(seq_len(i)), 2] <- targetSym[j]
        
        C[(i - 1) * m2 + j - sum(seq_len(i)), 1] <- M3
        C[(i - 1) * m2 + j - sum(seq_len(i)), 2] <- M5
        
      }
      
    }
  }
  
  # Extract RNA-RNA pair based on the homology of sharing miRNAs
  C[, 2] <- p.adjust(C[, 2], method = padjustmethod)
  index <- which(C[, 2] < padjustvaluecutoff)
  
  ceRInt <- ceRInt[index, ]
  C <- C[index, ]
  
  if (is.vector(C)) {
    miRHomologyceRInt <- c(ceRInt, C)
    names(miRHomologyceRInt) <- c("sponge_1", "sponge_2", "#shared miRNAs", "p.adjusted_value of shared miRNAs")
  } else {
    miRHomologyceRInt <- cbind(ceRInt, C)
    colnames(miRHomologyceRInt) <- c("sponge_1", "sponge_2", "#shared miRNAs", "p.adjusted_value of shared miRNAs")
  }
  
  return(miRHomologyceRInt)
}

# Parallel version
miRHomology_parallel <- function(miRTarget, minSharedmiR = 3, padjustvaluecutoff = 0.01, padjustmethod = "BH") {
  
  miRTarget <- as.matrix(miRTarget)
  
  m1 <- nrow(miRTarget)
  n1 <- ncol(miRTarget)
  
  miR <- miRTarget[, 1]
  tar <- miRTarget[, 2]
  
  miRSym <- unique(miR)
  targetSym <- unique(tar)
  
  m2 <- length(targetSym)
  index <- t(utils::combn(m2, 2))  
  
  Res <- foreach(i = seq_len(nrow(index))) %dopar% {    
    
    Interin1 <- miRTarget[which(miRTarget[, 2] %in% targetSym[index[i, 1]]), 1]
    Interin2 <- miRTarget[which(miRTarget[, 2] %in% targetSym[index[i, 2]]), 1]
    
    M1 <- length(Interin1)
    M2 <- length(Interin2)
    M3 <- length(intersect(Interin1, Interin2))
    M4 <- length(miRSym)
    M5 <- 1 - phyper(M3 - 1, M2, M4 - M2, M1)
    
    if (M3 >= minSharedmiR & M5 < padjustvaluecutoff) {
      
      C1 <- targetSym[index[i, 1]]
      C2 <- targetSym[index[i, 2]]
      C3 <- M3
      C4 <- M5
      
    } else {
      
      C1 <- NA; C2 <- NA; C3 <- NA; C4 <- NA 
      
    }
    
    tmp <- c(C1, C2, C3, C4)    
    return(tmp)
  }
  
  Res <- do.call(rbind, Res)
  
  # Extract RNA-RNA pair based on the homology of sharing miRNAs
  Res[, 4] <- p.adjust(as.numeric(Res[, 4]), method = padjustmethod)
  
  miRHomologyceRInt <- Res[which(as.numeric(Res[, 4]) < padjustvaluecutoff), ]         
  
  if (is.vector(miRHomologyceRInt)) {      
    names(miRHomologyceRInt) <- c("sponge_1", "sponge_2", "#shared miRNAs", "p.adjusted_value of shared miRNAs")
  } else {      
    colnames(miRHomologyceRInt) <- c("sponge_1", "sponge_2", "#shared miRNAs", "p.adjusted_value of shared miRNAs")
  }
  
  return(miRHomologyceRInt)
}


## 2. Positive Correlation (PC) method
# Original version
pc <- function(miRTarget, ExpData, consider.miR.expr = "TRUE", minSharedmiR = 3, poscorcutoff = 0, padjustvaluecutoff = 0.01,
               padjustmethod = "BH") {
  
  if (consider.miR.expr == "TRUE") {
    miRTargetCandidate <- querymiRTargetbinding(ExpData, miRTarget, type = "all")
    miRTargetCandidate <- as.matrix(miRTargetCandidate)
  } else {
    miRTargetCandidate <- querymiRTargetbinding(ExpData, miRTarget, type = "target")
    miRTargetCandidate <- as.matrix(miRTargetCandidate)
  }
  m1 <- nrow(miRTargetCandidate)
  n1 <- ncol(miRTargetCandidate)
  
  colnames(ExpData) <- gsub("\\.", "-", colnames(ExpData))
  ExpDataNames <- colnames(ExpData)
  
  miR <- miRTargetCandidate[, 1]
  tar <- miRTargetCandidate[, 2]
  
  miRSym <- unique(miR)
  targetSym <- unique(tar)
  
  m2 <- length(targetSym)
  
  # Initialize variables
  ceRInt <- matrix(NA, m2 * (m2 - 1)/2, 2)
  C <- matrix(NA, m2 * (m2 - 1)/2, 4)
  
  for (i in seq_len(m2 - 1)) {
    for (j in seq(i + 1, m2)) {
      
      Interin1 <- miRTargetCandidate[which(miRTargetCandidate[, 2] %in% targetSym[i]), 1]
      Interin2 <- miRTargetCandidate[which(miRTargetCandidate[, 2] %in% targetSym[j]), 1]
      
      M1 <- length(Interin1)
      M2 <- length(Interin2)
      M3 <- length(intersect(Interin1, Interin2))
      M4 <- length(miRSym)
      M5 <- 1 - phyper(M3 - 1, M2, M4 - M2, M1)
      
      if (M3 >= minSharedmiR & M5 < padjustvaluecutoff) {
        
        ceRInt[(i - 1) * m2 + j - sum(seq_len(i)), 1] <- targetSym[i]
        ceRInt[(i - 1) * m2 + j - sum(seq_len(i)), 2] <- targetSym[j]
        
        tarExpIdx1 <- which(ExpDataNames %in% targetSym[i])
        tarExpIdx2 <- which(ExpDataNames %in% targetSym[j])
        
        # Calculate Pearson correlation of each RNA-RNA pair
        M6 <- cor.test(ExpData[, tarExpIdx1], ExpData[, tarExpIdx2])$estimate
        M7 <- cor.test(ExpData[, tarExpIdx1], ExpData[, tarExpIdx2])$p.value
        
        C[(i - 1) * m2 + j - sum(seq_len(i)), 1] <- M3
        C[(i - 1) * m2 + j - sum(seq_len(i)), 2] <- M5
        C[(i - 1) * m2 + j - sum(seq_len(i)), 3] <- M6
        C[(i - 1) * m2 + j - sum(seq_len(i)), 4] <- M7
      }
    }
  }
  
  # Extract positive correlated RNA-RNA pairs.
  C[, 2] <- p.adjust(C[, 2], method = padjustmethod)
  C[, 4] <- p.adjust(C[, 4], method = padjustmethod)
  index <- which(C[, 2] < padjustvaluecutoff & C[, 3] > poscorcutoff & C[, 4] < padjustvaluecutoff)
  
  ceRInt <- ceRInt[index, ]
  C <- C[index, ]
  
  if (is.vector(C)) {
    PCceRInt <- c(ceRInt, C)
    names(PCceRInt) <- c("sponge_1", "sponge_2", "#shared miRNAs", "p.adjusted_value of shared miRNAs",
                         "correlation", "p.adjusted_value of positive correlation")
  } else {
    PCceRInt <- cbind(ceRInt, C)
    colnames(PCceRInt) <- c("sponge_1", "sponge_2", "#shared miRNAs", "p.adjusted_value of shared miRNAs",
                            "correlation", "p.adjusted_value of positive correlation")
  }
  
  return(PCceRInt)
}

# Parallel version
pc_parallel <- function(miRTarget, ExpData, consider.miR.expr = "TRUE", minSharedmiR = 3, poscorcutoff = 0, padjustvaluecutoff = 0.01,
                        padjustmethod = "BH") {
  
  if (consider.miR.expr == "TRUE") {
    miRTargetCandidate <- querymiRTargetbinding(ExpData, miRTarget, type = "all")
    miRTargetCandidate <- as.matrix(miRTargetCandidate)
  } else {
    miRTargetCandidate <- querymiRTargetbinding(ExpData, miRTarget, type = "target")
    miRTargetCandidate <- as.matrix(miRTargetCandidate)
  }
  m1 <- nrow(miRTargetCandidate)
  n1 <- ncol(miRTargetCandidate)
  
  colnames(ExpData) <- gsub("\\.", "-", colnames(ExpData))
  ExpDataNames <- colnames(ExpData)
  
  miR <- miRTargetCandidate[, 1]
  tar <- miRTargetCandidate[, 2]
  
  miRSym <- unique(miR)
  targetSym <- unique(tar)
  
  m2 <- length(targetSym)
  index <- t(utils::combn(m2, 2))    
  
  Res <- foreach(i = seq_len(nrow(index))) %dopar% {	    
    
    Interin1 <- miRTargetCandidate[which(miRTargetCandidate[, 2] %in% targetSym[index[i, 1]]), 1]
    Interin2 <- miRTargetCandidate[which(miRTargetCandidate[, 2] %in% targetSym[index[i, 2]]), 1]
    
    M1 <- length(Interin1)
    M2 <- length(Interin2)
    M3 <- length(intersect(Interin1, Interin2))
    M4 <- length(miRSym)
    M5 <- 1 - phyper(M3 - 1, M2, M4 - M2, M1)
    
    if (M3 >= minSharedmiR & M5 < padjustvaluecutoff) {
      
      C1 <- targetSym[index[i, 1]]
      C2 <- targetSym[index[i, 2]]
      
      tarExpIdx1 <- which(ExpDataNames %in% targetSym[index[i, 1]])
      tarExpIdx2 <- which(ExpDataNames %in% targetSym[index[i, 2]])
      
      # Calculate Pearson correlation of each RNA-RNA pair
      M6 <- cor.test(ExpData[, tarExpIdx1], ExpData[, tarExpIdx2])$estimate
      M7 <- cor.test(ExpData[, tarExpIdx1], ExpData[, tarExpIdx2])$p.value
      
      C3 <- M3
      C4 <- M5
      C5 <- M6
      C6 <- M7
    } else {
      C1 <- NA; C2 <- NA; C3 <- NA; C4 <- NA; C5 <- NA; C6 <- NA 	       
    }
    
    tmp <- c(C1, C2, C3, C4, C5, C6)    
    return(tmp)
  }
  
  Res <- do.call(rbind, Res)
  # Extract positive correlated RNA-RNA pairs.
  Res[, 4] <- p.adjust(as.numeric(Res[, 4]), method = padjustmethod)
  Res[, 6] <- p.adjust(as.numeric(Res[, 6]), method = padjustmethod)
  
  PCceRInt <- Res[which(as.numeric(Res[, 4]) < padjustvaluecutoff & as.numeric(Res[, 5]) > poscorcutoff & as.numeric(Res[, 6]) < padjustvaluecutoff), ]         
  
  if (is.vector(PCceRInt)) {
    names(PCceRInt) <- c("sponge_1", "sponge_2", "#shared miRNAs", "p.adjusted_value of shared miRNAs",
                         "correlation", "p.adjusted_value of positive correlation")
  } else {
    colnames(PCceRInt) <- c("sponge_1", "sponge_2", "#shared miRNAs", "p.adjusted_value of shared miRNAs",
                            "correlation", "p.adjusted_value of positive correlation")
  }
  
  return(PCceRInt)
}

## 3. Sensitivity Partial Pearson Correlation (SPPC) method
# Original version
sppc <- function(miRTarget, ExpData, minSharedmiR = 3, poscorcutoff = 0, padjustvaluecutoff = 0.01,
                 padjustmethod = "BH", senscorcutoff = 0.3) {
  
  miRTargetCandidate <- querymiRTargetbinding(ExpData, miRTarget, type = "all")
  miRTargetCandidate <- as.matrix(miRTargetCandidate)
  m1 <- nrow(miRTargetCandidate)
  n1 <- ncol(miRTargetCandidate)
  
  colnames(ExpData) <- gsub("\\.", "-", colnames(ExpData))
  ExpDataNames <- colnames(ExpData)
  
  miR <- miRTargetCandidate[, 1]
  tar <- miRTargetCandidate[, 2]
  
  miRSym <- unique(miR)
  targetSym <- unique(tar)
  
  m2 <- length(targetSym)
  
  # Initialize variables
  ceRInt <- matrix(NA, m2 * (m2 - 1)/2, 2)
  C <- matrix(NA, m2 * (m2 - 1)/2, 5)
  
  for (i in seq_len(m2 - 1)) {
    for (j in seq(i + 1, m2)) {
      
      Interin1 <- miRTargetCandidate[which(miRTargetCandidate[, 2] %in% targetSym[i]), 1]
      Interin2 <- miRTargetCandidate[which(miRTargetCandidate[, 2] %in% targetSym[j]), 1]
      
      M1 <- length(Interin1)
      M2 <- length(Interin2)
      M3 <- length(intersect(Interin1, Interin2))
      M4 <- length(miRSym)
      M5 <- 1 - phyper(M3 - 1, M2, M4 - M2, M1)
      
      if (M3 >= minSharedmiR & M5 < padjustvaluecutoff) {
        
        ceRInt[(i - 1) * m2 + j - sum(seq_len(i)), 1] <- targetSym[i]
        ceRInt[(i - 1) * m2 + j - sum(seq_len(i)), 2] <- targetSym[j]
        
        tarExpIdx1 <- which(ExpDataNames %in% targetSym[i])
        tarExpIdx2 <- which(ExpDataNames %in% targetSym[j])
        miRExpIdx <- which(ExpDataNames %in% intersect(Interin1, Interin2))
        
        # Calculate sensitivity correlation of each RNA-RNA pair
        M6 <- cor.test(ExpData[, tarExpIdx1], ExpData[, tarExpIdx2])$estimate
        M7 <- cor.test(ExpData[, tarExpIdx1], ExpData[, tarExpIdx2])$p.value
        M8 <- M6 - corpcor::pcor.shrink(cbind(ExpData[, tarExpIdx1], ExpData[, tarExpIdx2],
                                              ExpData[, miRExpIdx]), verbose = FALSE)[1, 2]
        
        C[(i - 1) * m2 + j - sum(seq_len(i)), 1] <- M3
        C[(i - 1) * m2 + j - sum(seq_len(i)), 2] <- M5
        C[(i - 1) * m2 + j - sum(seq_len(i)), 3] <- M6
        C[(i - 1) * m2 + j - sum(seq_len(i)), 4] <- M7
        C[(i - 1) * m2 + j - sum(seq_len(i)), 5] <- M8
      }
    }
  }
  
  # Extract RNA-RNA pairs with sensitivity correlation more than senscorcutoff.
  C[, 2] <- p.adjust(C[, 2], method = padjustmethod)
  C[, 4] <- p.adjust(C[, 4], method = padjustmethod)
  index <- which(C[, 2] < padjustvaluecutoff & C[, 3] > poscorcutoff & C[, 4] < padjustvaluecutoff &
                   C[, 5] > senscorcutoff)
  
  ceRInt <- ceRInt[index, ]
  C <- C[index, ]
  
  if (is.vector(C)) {
    SPPCceRInt <- c(ceRInt, C)
    names(SPPCceRInt) <- c("sponge_1", "sponge_2", "#shared miRNAs", "p.adjusted_value of shared miRNAs",
                           "correlation", "p.adjusted_value of positive correlation", "sensitivity correlation")
  } else {
    SPPCceRInt <- cbind(ceRInt, C)
    colnames(SPPCceRInt) <- c("sponge_1", "sponge_2", "#shared miRNAs", "p.adjusted_value of shared miRNAs",
                              "correlation", "p.adjusted_value of positive correlation", "sensitivity correlation")
  }
  
  return(SPPCceRInt)
}

# Parallel version
sppc_parallel <- function(miRTarget, ExpData, minSharedmiR = 3, poscorcutoff = 0, padjustvaluecutoff = 0.01,
                          padjustmethod = "BH", senscorcutoff = 0.3) {
  
  miRTargetCandidate <- querymiRTargetbinding(ExpData, miRTarget, type = "all")
  miRTargetCandidate <- as.matrix(miRTargetCandidate)
  m1 <- nrow(miRTargetCandidate)
  n1 <- ncol(miRTargetCandidate)
  
  colnames(ExpData) <- gsub("\\.", "-", colnames(ExpData))
  ExpDataNames <- colnames(ExpData)
  
  miR <- miRTargetCandidate[, 1]
  tar <- miRTargetCandidate[, 2]
  
  miRSym <- unique(miR)
  targetSym <- unique(tar)
  
  m2 <- length(targetSym)
  index <- t(utils::combn(m2, 2))
  
  Res <- foreach(i = seq_len(nrow(index))) %dopar% {	    
    
    Interin1 <- miRTargetCandidate[which(miRTargetCandidate[, 2] %in% targetSym[index[i, 1]]), 1]
    Interin2 <- miRTargetCandidate[which(miRTargetCandidate[, 2] %in% targetSym[index[i, 2]]), 1]
    
    M1 <- length(Interin1)
    M2 <- length(Interin2)
    M3 <- length(intersect(Interin1, Interin2))
    M4 <- length(miRSym)
    M5 <- 1 - phyper(M3 - 1, M2, M4 - M2, M1)
    
    if (M3 >= minSharedmiR & M5 < padjustvaluecutoff) {
      
      C1 <- targetSym[index[i, 1]]
      C2 <- targetSym[index[i, 2]]
      
      tarExpIdx1 <- which(ExpDataNames %in% targetSym[index[i, 1]])
      tarExpIdx2 <- which(ExpDataNames %in% targetSym[index[i, 2]])          
      miRExpIdx <- which(ExpDataNames %in% intersect(Interin1, Interin2))
      
      # Calculate sensitivity correlation of each RNA-RNA pair
      M6 <- cor.test(ExpData[, tarExpIdx1], ExpData[, tarExpIdx2])$estimate
      M7 <- cor.test(ExpData[, tarExpIdx1], ExpData[, tarExpIdx2])$p.value
      M8 <- M6 - corpcor::pcor.shrink(cbind(ExpData[, tarExpIdx1], ExpData[, tarExpIdx2],
                                            ExpData[, miRExpIdx]), verbose = FALSE)[1, 2]
      
      C3 <- M3
      C4 <- M5
      C5 <- M6
      C6 <- M7
      C7 <- M8
      
    } else {
      C1 <- NA; C2 <- NA; C3 <- NA; C4 <- NA; C5 <- NA; C6 <- NA; C7 <- NA 	       
    }
    
    tmp <- c(C1, C2, C3, C4, C5, C6, C7)    
    return(tmp)
  }
  
  Res <- do.call(rbind, Res)
  
  # Extract RNA-RNA pairs with sensitivity correlation more than senscorcutoff.
  Res[, 4] <- p.adjust(as.numeric(Res[, 4]), method = padjustmethod)
  Res[, 6] <- p.adjust(as.numeric(Res[, 6]), method = padjustmethod)
  
  SPPCceRInt <- Res[which(as.numeric(Res[, 4]) < padjustvaluecutoff & as.numeric(Res[, 5]) > poscorcutoff & 
                            as.numeric(Res[, 6]) < padjustvaluecutoff & as.numeric(Res[, 7]) > senscorcutoff), ]         
  
  if (is.vector(SPPCceRInt)) {
    names(SPPCceRInt) <- c("sponge_1", "sponge_2", "#shared miRNAs", "p.adjusted_value of shared miRNAs",
                           "correlation", "p.adjusted_value of positive correlation", "sensitivity correlation")
  } else {
    colnames(SPPCceRInt) <- c("sponge_1", "sponge_2", "#shared miRNAs", "p.adjusted_value of shared miRNAs",
                              "correlation", "p.adjusted_value of positive correlation", "sensitivity correlation")
  }
  
  return(SPPCceRInt)
}

## 4. Partial Pearson Correlation (PPC) method
# Original version
ppc <- function(miRTarget, ExpData, minSharedmiR = 3, num_perm = 100, padjustvaluecutoff = 0.01,
                padjustmethod = "BH") {
  
  miRTargetCandidate <- querymiRTargetbinding(ExpData, miRTarget, type = "all")
  miRTargetCandidate <- as.matrix(miRTargetCandidate)
  m1 <- nrow(miRTargetCandidate)
  n1 <- ncol(miRTargetCandidate)
  
  colnames(ExpData) <- gsub("\\.", "-", colnames(ExpData))
  ExpDataNames <- colnames(ExpData)
  
  miR <- miRTargetCandidate[, 1]
  tar <- miRTargetCandidate[, 2]
  
  miRSym <- unique(miR)
  targetSym <- unique(tar)
  
  m2 <- length(targetSym)
  
  # Initialize variables
  ceRInt <- matrix(NA, m2 * (m2 - 1)/2, 2)
  C <- matrix(NA, m2 * (m2 - 1)/2, 3)
  
  for (i in seq_len(m2 - 1)) {
    for (j in seq(i + 1, m2)) {
      
      Interin1 <- miRTargetCandidate[which(miRTargetCandidate[, 2] %in% targetSym[i]), 1]
      Interin2 <- miRTargetCandidate[which(miRTargetCandidate[, 2] %in% targetSym[j]), 1]
      
      M1 <- length(Interin1)
      M2 <- length(Interin2)
      M3 <- length(intersect(Interin1, Interin2))
      M4 <- length(miRSym)
      M5 <- 1 - phyper(M3 - 1, M2, M4 - M2, M1)
      
      if (M3 >= minSharedmiR & M5 < padjustvaluecutoff) {
        
        ceRInt[(i - 1) * m2 + j - sum(seq_len(i)), 1] <- targetSym[i]
        ceRInt[(i - 1) * m2 + j - sum(seq_len(i)), 2] <- targetSym[j]
        
        tarExpIdx1 <- which(ExpDataNames %in% targetSym[i])
        tarExpIdx2 <- which(ExpDataNames %in% targetSym[j])
        miRExpIdx <- which(ExpDataNames %in% intersect(Interin1, Interin2))      
        
        inputdata <- t(cbind(ExpData[, tarExpIdx1], ExpData[, tarExpIdx2], ExpData[, miRExpIdx]))
        
        C[(i - 1) * m2 + j - sum(seq_len(i)), 1] <- M3
        C[(i - 1) * m2 + j - sum(seq_len(i)), 2] <- M5
        C[(i - 1) * m2 + j - sum(seq_len(i)), 3] <- predCor(inputdata, num_perm, method = "PPC")
        
      }
    }
  }
  
  # Extract significant RNA-RNA pairs.
  C[, 2] <- p.adjust(C[, 2], method = padjustmethod)
  C[, 3] <- p.adjust(C[, 3], method = padjustmethod)
  index <- which(C[, 2] < padjustvaluecutoff & C[, 3] < padjustvaluecutoff)
  
  ceRInt <- ceRInt[index, ]
  C <- C[index, ]
  
  if (is.vector(C)) {
    PPCceRInt <- c(ceRInt, C)
    names(PPCceRInt) <- c("sponge_1", "sponge_2", "#shared miRNAs", "p.adjusted_value of shared miRNAs",
                          "p.adjusted_value of RNA competition")    
  } else {
    PPCceRInt <- cbind(ceRInt, C)
    colnames(PPCceRInt) <- c("sponge_1", "sponge_2", "#shared miRNAs", "p.adjusted_value of shared miRNAs",
                             "p.adjusted_value of RNA competition")
  }
  
  return(PPCceRInt)
}

# Parallel version
ppc_parallel <- function(miRTarget, ExpData, minSharedmiR = 3, num_perm = 100, padjustvaluecutoff = 0.01,
                         padjustmethod = "BH") {
  
  miRTargetCandidate <- querymiRTargetbinding(ExpData, miRTarget, type = "all")
  miRTargetCandidate <- as.matrix(miRTargetCandidate)
  m1 <- nrow(miRTargetCandidate)
  n1 <- ncol(miRTargetCandidate)
  
  colnames(ExpData) <- gsub("\\.", "-", colnames(ExpData))
  ExpDataNames <- colnames(ExpData)
  
  miR <- miRTargetCandidate[, 1]
  tar <- miRTargetCandidate[, 2]
  
  miRSym <- unique(miR)
  targetSym <- unique(tar)
  
  m2 <- length(targetSym)
  index <- t(utils::combn(m2, 2))
  
  Res <- foreach(i = seq_len(nrow(index)), .export = c("predCor", "combpvalue")) %dopar% {	    
    
    Interin1 <- miRTargetCandidate[which(miRTargetCandidate[, 2] %in% targetSym[index[i, 1]]), 1]
    Interin2 <- miRTargetCandidate[which(miRTargetCandidate[, 2] %in% targetSym[index[i, 2]]), 1]      
    
    M1 <- length(Interin1)
    M2 <- length(Interin2)
    M3 <- length(intersect(Interin1, Interin2))
    M4 <- length(miRSym)
    M5 <- 1 - phyper(M3 - 1, M2, M4 - M2, M1)
    
    if (M3 >= minSharedmiR & M5 < padjustvaluecutoff) {
      
      C1 <- targetSym[index[i, 1]]
      C2 <- targetSym[index[i, 2]]
      
      tarExpIdx1 <- which(ExpDataNames %in% targetSym[index[i, 1]])
      tarExpIdx2 <- which(ExpDataNames %in% targetSym[index[i, 2]])
      miRExpIdx <- which(ExpDataNames %in% intersect(Interin1, Interin2))
      
      inputdata <- t(cbind(ExpData[, tarExpIdx1], ExpData[, tarExpIdx2], ExpData[, miRExpIdx]))
      
      C3 <- M3
      C4 <- M5
      C5 <- predCor(inputdata, num_perm, method = "PPC")
    } else {
      C1 <- NA; C2 <- NA; C3 <- NA; C4 <- NA; C5 <- NA 	       
    }
    
    tmp <- c(C1, C2, C3, C4, C5)    
    return(tmp)
  }
  
  Res <- do.call(rbind, Res)
  
  # Extract significant RNA-RNA pairs.
  Res[, 4] <- p.adjust(as.numeric(Res[, 4]), method = padjustmethod)
  Res[, 5] <- p.adjust(as.numeric(Res[, 5]), method = padjustmethod)
  
  PPCceRInt <- Res[which(as.numeric(Res[, 4]) < padjustvaluecutoff & as.numeric(Res[, 5]) < padjustvaluecutoff), ]
  
  if (is.vector(PPCceRInt)) {
    names(PPCceRInt) <- c("sponge_1", "sponge_2", "#shared miRNAs", "p.adjusted_value of shared miRNAs",
                          "p.adjusted_value of RNA competition")
  } else {
    colnames(PPCceRInt) <- c("sponge_1", "sponge_2", "#shared miRNAs", "p.adjusted_value of shared miRNAs",
                             "p.adjusted_value of RNA competition")
  }
  
  return(PPCceRInt)
}

## 5. Hermes method
# Original version
hermes <- function(miRTarget, ExpData, minSharedmiR = 3, num_perm = 100, padjustvaluecutoff = 0.01,
                   padjustmethod = "BH") {
  
  miRTargetCandidate <- querymiRTargetbinding(ExpData, miRTarget, type = "all")
  miRTargetCandidate <- as.matrix(miRTargetCandidate)
  m1 <- nrow(miRTargetCandidate)
  n1 <- ncol(miRTargetCandidate)
  
  colnames(ExpData) <- gsub("\\.", "-", colnames(ExpData))
  ExpDataNames <- colnames(ExpData)
  
  miR <- miRTargetCandidate[, 1]
  tar <- miRTargetCandidate[, 2]
  
  miRSym <- unique(miR)
  targetSym <- unique(tar)
  
  m2 <- length(targetSym)
  
  # Initialize variables
  ceRInt <- matrix(NA, m2 * (m2 - 1)/2, 2)
  C <- matrix(NA, m2 * (m2 - 1)/2, 3)
  
  for (i in seq_len(m2 - 1)) {
    for (j in seq(i + 1, m2)) {
      
      Interin1 <- miRTargetCandidate[which(miRTargetCandidate[, 2] %in% targetSym[i]), 1]
      Interin2 <- miRTargetCandidate[which(miRTargetCandidate[, 2] %in% targetSym[j]), 1]
      
      M1 <- length(Interin1)
      M2 <- length(Interin2)
      M3 <- length(intersect(Interin1, Interin2))
      M4 <- length(miRSym)
      M5 <- 1 - phyper(M3 - 1, M2, M4 - M2, M1)
      
      if (M3 >= minSharedmiR & M5 < padjustvaluecutoff) {
        
        ceRInt[(i - 1) * m2 + j - sum(seq_len(i)), 1] <- targetSym[i]
        ceRInt[(i - 1) * m2 + j - sum(seq_len(i)), 2] <- targetSym[j]
        
        tarExpIdx1 <- which(ExpDataNames %in% targetSym[i])
        tarExpIdx2 <- which(ExpDataNames %in% targetSym[j])
        miRExpIdx <- which(ExpDataNames %in% intersect(Interin1, Interin2))
        
        inputdata <- t(cbind(ExpData[, tarExpIdx1], ExpData[, tarExpIdx2], ExpData[, miRExpIdx]))
        
        C[(i - 1) * m2 + j - sum(seq_len(i)), 1] <- M3
        C[(i - 1) * m2 + j - sum(seq_len(i)), 2] <- M5
        C[(i - 1) * m2 + j - sum(seq_len(i)), 3] <- predCor(inputdata, num_perm, method = "CMI")
      }
    }
  }
  
  # Extract significant RNA-RNA pairs.
  C[, 2] <- p.adjust(C[, 2], method = padjustmethod)
  C[, 3] <- p.adjust(C[, 3], method = padjustmethod)
  index <- which(C[, 2] < padjustvaluecutoff & C[, 3] < padjustvaluecutoff)
  
  ceRInt <- ceRInt[index, ]
  C <- C[index, ]
  
  if (is.vector(C)) {
    HermesceRInt <- c(ceRInt, C)
    names(HermesceRInt) <- c("sponge_1", "sponge_2", "#shared miRNAs", "p.adjusted_value of shared miRNAs",
                             "p.adjusted_value of RNA competition")
  } else {
    HermesceRInt <- cbind(ceRInt, C)
    colnames(HermesceRInt) <- c("sponge_1", "sponge_2", "#shared miRNAs", "p.adjusted_value of shared miRNAs",
                                "p.adjusted_value of RNA competition")
  }
  
  return(HermesceRInt)
}

# Parallel version
hermes_parallel <- function(miRTarget, ExpData, minSharedmiR = 3, num_perm = 100, padjustvaluecutoff = 0.01,
                            padjustmethod = "BH") {
  
  miRTargetCandidate <- querymiRTargetbinding(ExpData, miRTarget, type = "all")
  miRTargetCandidate <- as.matrix(miRTargetCandidate)
  m1 <- nrow(miRTargetCandidate)
  n1 <- ncol(miRTargetCandidate)
  
  colnames(ExpData) <- gsub("\\.", "-", colnames(ExpData))
  ExpDataNames <- colnames(ExpData)
  
  miR <- miRTargetCandidate[, 1]
  tar <- miRTargetCandidate[, 2]
  
  miRSym <- unique(miR)
  targetSym <- unique(tar)
  
  m2 <- length(targetSym)
  index <- t(utils::combn(m2, 2))
  
  Res <- foreach(i = seq_len(nrow(index)), .export = c("predCor", "combpvalue", "calCMI")) %dopar% {	    
    
    Interin1 <- miRTargetCandidate[which(miRTargetCandidate[, 2] %in% targetSym[index[i, 1]]), 1]
    Interin2 <- miRTargetCandidate[which(miRTargetCandidate[, 2] %in% targetSym[index[i, 2]]), 1]      
    
    M1 <- length(Interin1)
    M2 <- length(Interin2)
    M3 <- length(intersect(Interin1, Interin2))
    M4 <- length(miRSym)
    M5 <- 1 - phyper(M3 - 1, M2, M4 - M2, M1)
    
    if (M3 >= minSharedmiR & M5 < padjustvaluecutoff) {
      
      C1 <- targetSym[index[i, 1]]
      C2 <- targetSym[index[i, 2]]
      
      tarExpIdx1 <- which(ExpDataNames %in% targetSym[index[i, 1]])
      tarExpIdx2 <- which(ExpDataNames %in% targetSym[index[i, 2]])
      miRExpIdx <- which(ExpDataNames %in% intersect(Interin1, Interin2))
      
      inputdata <- t(cbind(ExpData[, tarExpIdx1], ExpData[, tarExpIdx2], ExpData[, miRExpIdx]))
      
      C3 <- M3
      C4 <- M5
      C5 <- predCor(inputdata, num_perm, method = "CMI")
    } else {
      C1 <- NA; C2 <- NA; C3 <- NA; C4 <- NA; C5 <- NA 	       
    }
    
    tmp <- c(C1, C2, C3, C4, C5)    
    return(tmp)
  }
  
  Res <- do.call(rbind, Res)
  
  # Extract significant RNA-RNA pairs.
  Res[, 4] <- p.adjust(as.numeric(Res[, 4]), method = padjustmethod)
  Res[, 5] <- p.adjust(as.numeric(Res[, 5]), method = padjustmethod)
  
  HermesceRInt <- Res[which(as.numeric(Res[, 4]) < padjustvaluecutoff & as.numeric(Res[, 5]) < padjustvaluecutoff), ]
  
  if (is.vector(HermesceRInt)) {
    names(HermesceRInt) <- c("sponge_1", "sponge_2", "#shared miRNAs", "p.adjusted_value of shared miRNAs",
                             "p.adjusted_value of RNA competition")
  } else {
    colnames(HermesceRInt) <- c("sponge_1", "sponge_2", "#shared miRNAs", "p.adjusted_value of shared miRNAs",
                                "p.adjusted_value of RNA competition")
  }
  
  return(HermesceRInt)
}

## 6. MuTaME method
# Original version
muTaME <- function(miRTarget, mres, minSharedmiR = 3, padjustvaluecutoff = 0.01, padjustmethod = "BH",
                   scorecutoff = 0.5) {
  
  miRTarget <- as.matrix(miRTarget)
  m1 <- nrow(miRTarget)
  n1 <- ncol(miRTarget)
  
  miR <- miRTarget[, 1]
  tar <- miRTarget[, 2]
  
  miRSym <- unique(miR)
  targetSym <- unique(tar)
  
  m2 <- length(targetSym)
  
  # Initialize variables
  ceRInt <- matrix(NA, m2 * (m2 - 1)/2, 2)
  C <- matrix(NA, m2 * (m2 - 1)/2, 8)
  
  for (i in seq_len(m2 - 1)) {
    for (j in seq(i + 1, m2)) {
      
      Interin1 <- miRTarget[which(miRTarget[, 2] %in% targetSym[i]), 1]
      Interin2 <- miRTarget[which(miRTarget[, 2] %in% targetSym[j]), 1]
      cm <- intersect(Interin1, Interin2)
      SharedMREs <- mres[mres[, 2] %in% c(targetSym[i], targetSym[j]) & mres[, 1] %in% cm, ]
      
      M1 <- length(Interin1)
      M2 <- length(Interin2)
      M3 <- length(intersect(Interin1, Interin2))
      M4 <- length(miRSym)
      M5 <- 1 - phyper(M3 - 1, M2, M4 - M2, M1)
      
      if (M3 >= minSharedmiR & M5 < padjustvaluecutoff & nrow(SharedMREs) > 0) {
        
        ceRInt[(i - 1) * m2 + j - sum(seq_len(i)), 1] <- targetSym[i]
        ceRInt[(i - 1) * m2 + j - sum(seq_len(i)), 2] <- targetSym[j]
        
        C[(i - 1) * m2 + j - sum(seq_len(i)), 1] <- M3
        C[(i - 1) * m2 + j - sum(seq_len(i)), 2] <- M5
        
        # Score 1 for the fraction of coomon miRNAs
        C[(i - 1) * m2 + j - sum(seq_len(i)), 3] <- log(M3/min(M1, M2))
        
        # Score 2 for the density of the MREs for all shared miRNAs
        C[(i - 1) * m2 + j - sum(seq_len(i)), 4] <- sum(sapply(cm, function(miR) {
          MREs <- SharedMREs[SharedMREs[, 1] == miR, ]
          if (nrow(MREs) <= 0) return(1)
          MREslr <- MREs[, c("gap_l", "gap_r")]
          D <- abs(max(MREslr[, 1]) - min(MREslr[, 2]))
          return(log(nrow(MREs)/D))
        }))
        
        # Score 3 for the distribution of MREs of the putative RNA-RNA pairs
        C[(i - 1) * m2 + j - sum(seq_len(i)), 5] <- sum(sapply(cm, function(miR) {
          positions <- SharedMREs[SharedMREs[, 1] == miR, c("gap_l", "gap_r")]
          if (nrow(positions) <= 0) return(1)
          return(log(abs(max(positions[, 1]) - min(positions[, 2]))^2/sum((positions[, 2] -
                                                                             positions[, 1])^2)))
        }))
        
        # Score 4 for the relation between the overall number of MREs for a putative miRNA sponge,
        # compared with the number of miRNAs that yield these MREs
        B <- nrow(SharedMREs)
        if (B == length(unique(SharedMREs[, 1]))) {
          C[(i - 1) * m2 + j - sum(seq_len(i)), 6] <- log(1/B)
        } else {
          C[(i - 1) * m2 + j - sum(seq_len(i)), 6] <- log((B - length(unique(SharedMREs[,
                                                                                        1])) - 1)/B)
        }
        
        C[(i - 1) * m2 + j - sum(seq_len(i)), 7] <- C[(i - 1) * m2 + j - sum(seq_len(i)),
                                                      3] + C[(i - 1) * m2 + j - sum(seq_len(i)), 4] + C[(i - 1) * m2 + j - sum(seq_len(i)),
                                                                                                        5] + C[(i - 1) * m2 + j - sum(seq_len(i)), 6]
      }
    }
  }
  
  # Extract RNA-RNA pair based on four scores.
  ceRInt <- ceRInt[apply(ceRInt, 1, function(x) !all(is.na(x))), ]
  C <- C[apply(C, 1, function(x) !all(is.na(x))), ]
  C[, 8] <- (C[, 7] - min(C[, 7]))/(max(C[, 7]) - min(C[, 7]))
  
  C[, 2] <- p.adjust(C[, 2], method = padjustmethod)
  index <- which(C[, 2] < padjustvaluecutoff & C[, 8] > scorecutoff)
  
  ceRInt <- ceRInt[index, ]
  C <- C[index, ]
  
  if (is.vector(C)) {
    MuTaMEceRInt <- c(ceRInt, C)
    names(MuTaMEceRInt) <- c("sponge_1", "sponge_2", "#shared miRNAs", "p.adjusted_value of shared miRNAs",
                             "Score 1", "Score 2", "Score 3", "Score 4", "Combined score", "Normalized score")
  } else {
    MuTaMEceRInt <- cbind(ceRInt, C)
    colnames(MuTaMEceRInt) <- c("sponge_1", "sponge_2", "#shared miRNAs", "p.adjusted_value of shared miRNAs",
                                "Score 1", "Score 2", "Score 3", "Score 4", "Combined score", "Normalized score")
  }
  return(MuTaMEceRInt)
}

# Parallel version
muTaME_parallel <- function(miRTarget, mres, minSharedmiR = 3, padjustvaluecutoff = 0.01, padjustmethod = "BH",
                            scorecutoff = 0.5) {
  
  miRTarget <- as.matrix(miRTarget)
  m1 <- nrow(miRTarget)
  n1 <- ncol(miRTarget)
  
  miR <- miRTarget[, 1]
  tar <- miRTarget[, 2]
  
  miRSym <- unique(miR)
  targetSym <- unique(tar)
  
  m2 <- length(targetSym)
  index <- t(utils::combn(m2, 2))
  
  Res <- foreach(i = seq_len(nrow(index))) %dopar% {	    
    
    Interin1 <- miRTarget[which(miRTarget[, 2] %in% targetSym[index[i, 1]]), 1]
    Interin2 <- miRTarget[which(miRTarget[, 2] %in% targetSym[index[i, 2]]), 1] 
    
    cm <- intersect(Interin1, Interin2)
    SharedMREs <- mres[mres[, 2] %in% c(targetSym[index[i, 1]], targetSym[index[i, 2]]) & mres[, 1] %in% cm, ]
    
    M1 <- length(Interin1)
    M2 <- length(Interin2)
    M3 <- length(intersect(Interin1, Interin2))
    M4 <- length(miRSym)
    M5 <- 1 - phyper(M3 - 1, M2, M4 - M2, M1)
    
    if (M3 >= minSharedmiR & M5 < padjustvaluecutoff & nrow(SharedMREs) > 0) {
      
      C1 <- targetSym[index[i, 1]]
      C2 <- targetSym[index[i, 2]]
      
      C3 <- M3
      C4 <- M5
      
      # Score 1 for the fraction of coomon miRNAs
      C5 <- log(M3/min(M1, M2))
      
      # Score 2 for the density of the MREs for all shared miRNAs
      C6 <- sum(sapply(cm, function(miR) {
        MREs <- SharedMREs[SharedMREs[, 1] == miR, ]
        if (nrow(MREs) <= 0) return(1)
        MREslr <- MREs[, c("gap_l", "gap_r")]
        D <- abs(max(MREslr[, 1]) - min(MREslr[, 2]))
        return(log(nrow(MREs)/D))
      }))
      
      # Score 3 for the distribution of MREs of the putative RNA-RNA pairs
      C7 <- sum(sapply(cm, function(miR) {
        positions <- SharedMREs[SharedMREs[, 1] == miR, c("gap_l", "gap_r")]
        if (nrow(positions) <= 0) return(1)
        return(log(abs(max(positions[, 1]) - min(positions[, 2]))^2/sum((positions[, 2] -
                                                                           positions[, 1])^2)))
      }))
      
      # Score 4 for the relation between the overall number of MREs for a putative miRNA sponge,
      # compared with the number of miRNAs that yield these MREs
      B <- nrow(SharedMREs)
      if (B == length(unique(SharedMREs[, 1]))) {
        C8 <- log(1/B)
      } else {
        C8 <- log((B - length(unique(SharedMREs[, 1])) - 1)/B)
      }
      
      C9 <- C5 + C6 + C7 + C8
    } else {
      C1 <- NA; C2 <- NA; C3 <- NA; C4 <- NA; C5 <- NA; C6 <- NA; C7 <- NA; C8 <- NA; C9 <- NA 	       
    }
    
    tmp <- c(C1, C2, C3, C4, C5, C6, C7, C8, C9)    
    return(tmp)
  }
  
  Res <- do.call(rbind, Res)
  
  # Extract RNA-RNA pair based on four scores.
  Res <- Res[apply(Res, 1, function(x) !all(is.na(x))), ]
  C10 <- (as.numeric(Res[, 9]) - min(as.numeric(Res[, 9])))/(max(as.numeric(Res[, 9])) - min(as.numeric(Res[, 9])))
  Res <- cbind(Res, C10)
  Res[, 4] <- p.adjust(as.numeric(Res[, 4]), method = padjustmethod)
  
  MuTaMEceRInt <- Res[which(as.numeric(Res[, 4]) < padjustvaluecutoff & as.numeric(Res[, 10]) > scorecutoff), ]    
  
  if (is.vector(MuTaMEceRInt)) {
    names(MuTaMEceRInt) <- c("sponge_1", "sponge_2", "#shared miRNAs", "p.adjusted_value of shared miRNAs",
                             "Score 1", "Score 2", "Score 3", "Score 4", "Combined score", "Normalized score")
  } else {
    colnames(MuTaMEceRInt) <- c("sponge_1", "sponge_2", "#shared miRNAs", "p.adjusted_value of shared miRNAs",
                                "Score 1", "Score 2", "Score 3", "Score 4", "Combined score", "Normalized score")
  }
  return(MuTaMEceRInt)
}

## 7. CERNIA method
# Original version
cernia <- function(miRTarget, ExpData, mres, consider.miR.expr = "TRUE", minSharedmiR = 3, poscorcutoff = 0, 
                   padjustvaluecutoff = 0.01, padjustmethod = "BH", scorecutoff = 0.5) {
  
  if (consider.miR.expr == "TRUE") {
    miRTargetCandidate <- querymiRTargetbinding(ExpData, miRTarget, type = "all")
    miRTargetCandidate <- as.matrix(miRTargetCandidate)
  } else {
    miRTargetCandidate <- querymiRTargetbinding(ExpData, miRTarget, type = "target")
    miRTargetCandidate <- as.matrix(miRTargetCandidate)
  }
  m1 <- nrow(miRTargetCandidate)
  n1 <- ncol(miRTargetCandidate)
  
  colnames(ExpData) <- gsub("\\.", "-", colnames(ExpData))
  ExpDataNames <- colnames(ExpData)
  
  miR <- miRTargetCandidate[, 1]
  tar <- miRTargetCandidate[, 2]
  
  miRSym <- unique(miR)
  targetSym <- unique(tar)
  
  m2 <- length(targetSym)
  
  # Initialize variables
  ceRInt <- matrix(NA, m2 * (m2 - 1)/2, 2)
  C <- matrix(NA, m2 * (m2 - 1)/2, 11)
  
  for (i in seq_len(m2 - 1)) {
    for (j in seq(i + 1, m2)) {
      
      Interin1 <- miRTargetCandidate[which(miRTargetCandidate[, 2] %in% targetSym[i]), 1]
      Interin2 <- miRTargetCandidate[which(miRTargetCandidate[, 2] %in% targetSym[j]), 1]
      cm <- intersect(Interin1, Interin2)
      SharedMREs <- mres[mres[, 2] %in% c(targetSym[i], targetSym[j]) & mres[, 1] %in% cm, ]
      
      M1 <- length(Interin1)
      M2 <- length(Interin2)
      M3 <- length(intersect(Interin1, Interin2))
      M4 <- length(miRSym)
      M5 <- 1 - phyper(M3 - 1, M2, M4 - M2, M1)
      
      tarExpIdx1 <- which(ExpDataNames %in% targetSym[i])
      tarExpIdx2 <- which(ExpDataNames %in% targetSym[j])
      
      M6 <- cor.test(ExpData[, tarExpIdx1], ExpData[, tarExpIdx2])$estimate
      
      if (M3 >= minSharedmiR & M5 < padjustvaluecutoff & M6 > poscorcutoff &
          nrow(SharedMREs) > 0) {
        
        ceRInt[(i - 1) * m2 + j - sum(seq_len(i)), 1] <- targetSym[i]
        ceRInt[(i - 1) * m2 + j - sum(seq_len(i)), 2] <- targetSym[j]
        
        C[(i - 1) * m2 + j - sum(seq_len(i)), 1] <- M3
        C[(i - 1) * m2 + j - sum(seq_len(i)), 2] <- M5
        
        # Score 1 for the fraction of coomon miRNAs
        C[(i - 1) * m2 + j - sum(seq_len(i)), 3] <- log(M3/min(M1, M2))
        
        # Score 2 for the density of the MREs for all shared miRNAs
        C[(i - 1) * m2 + j - sum(seq_len(i)), 4] <- sum(sapply(cm, function(miR) {
          MREs <- SharedMREs[SharedMREs[, 1] == miR, ]
          if (nrow(MREs) <= 0) return(1)
          MREslr <- MREs[, c("gap_l", "gap_r")]
          D <- abs(max(MREslr[, 1]) - min(MREslr[, 2]))
          return(log(nrow(MREs)/D))
        }))
        
        # Score 3 for the distribution of MREs of the putative RNA-RNA pairs
        C[(i - 1) * m2 + j - sum(seq_len(i)), 5] = sum(sapply(cm, function(miR) {
          positions <- SharedMREs[SharedMREs[, 1] == miR, c("gap_l", "gap_r")]
          if (nrow(positions) <= 0) return(1)
          return(log(abs(max(positions[, 1]) - min(positions[, 2]))^2/sum((positions[, 2] -
                                                                             positions[, 1])^2)))
        }))
        
        # Score 4 for the relation between the overall number of MREs for a putative miRNA sponge,
        # compared with the number of miRNAs that yield these MREs
        B <- nrow(SharedMREs)
        if (B == length(unique(SharedMREs[, 1]))) {
          C[(i - 1) * m2 + j - sum(seq_len(i)), 6] <- log(1/B)
        } else {
          C[(i - 1) * m2 + j - sum(seq_len(i)), 6] <- log((B - length(unique(SharedMREs[,
                                                                                        1])) - 1)/B)
        }
        
        # Score 5 for the density of the hybridization energies related to MREs for all shared miRNAs
        SharedMREs <- mres[mres[, 2] %in% c(targetSym[i], targetSym[j]) & mres[, 1] %in%
                             cm, ]
        C[(i - 1) * m2 + j - sum(seq_len(i)), 7] <- sum(sapply(cm, function(miR) {
          MREs <- SharedMREs[SharedMREs[, 1] == miR, ]
          if (nrow(MREs) <= 0) return(1)
          MREslr <- MREs[, c("gap_l", "gap_r")]
          D <- abs(max(MREslr[, 1]) - min(MREslr[, 2]))
          return(log(sum(abs(MREs[, 3]))/D))
        }))
        
        # Score 6 for the DT-Hybrid recommendation scores
        cerna_recommendations <- dtHybrid(miRTargetCandidate)
        C[(i - 1) * m2 + j - sum(seq_len(i)), 8] <- cerna_recommendations[targetSym[i],
                                                                          targetSym[j]]
        
        # Score 7 for the pairwise Peason correlation between putative RNA-RNA pair expression data
        C[(i - 1) * m2 + j - sum(seq_len(i)), 9] <- log(M6)
        
        C[(i - 1) * m2 + j - sum(seq_len(i)), 10] <- C[(i - 1) * m2 + j - sum(seq_len(i)),
                                                       3] + C[(i - 1) * m2 + j - sum(seq_len(i)), 4] + C[(i - 1) * m2 + j - sum(seq_len(i)),
                                                                                                         5] + C[(i - 1) * m2 + j - sum(seq_len(i)), 6] + C[(i - 1) * m2 + j - sum(seq_len(i)),
                                                                                                                                                           7] + C[(i - 1) * m2 + j - sum(seq_len(i)), 8] + C[(i - 1) * m2 + j - sum(seq_len(i)),
                                                                                                                                                                                                             9]
      }
    }
  }
  
  # Extract RNA-RNA pair based on seven scores
  ceRInt <- ceRInt[apply(ceRInt, 1, function(x) !all(is.na(x))), ]
  C <- C[apply(C, 1, function(x) !all(is.na(x))), ]
  C[, 11] <- (C[, 10] - min(C[, 10]))/(max(C[, 10]) - min(C[, 10]))
  
  C[, 2] <- p.adjust(C[, 2], method = padjustmethod)
  index <- which(C[, 2] < padjustvaluecutoff & C[, 11] > scorecutoff)
  
  ceRInt <- ceRInt[index, ]
  C <- C[index, ]
  
  if (is.vector(C)) {
    CERNIAceRInt <- c(ceRInt, C)
    names(CERNIAceRInt) <- c("sponge_1", "sponge_2", "#shared miRNAs", "p.adjusted_value of shared miRNAs",
                             "Score 1", "Score 2", "Score 3", "Score 4", "Score 5", "Score 6", "Score 7", "Combined score",
                             "Normalized score")
  } else {
    CERNIAceRInt <- cbind(ceRInt, C)
    colnames(CERNIAceRInt) <- c("sponge_1", "sponge_2", "#shared miRNAs", "p.adjusted_value of shared miRNAs",
                                "Score 1", "Score 2", "Score 3", "Score 4", "Score 5", "Score 6", "Score 7", "Combined score",
                                "Normalized score")
  }
  
  return(CERNIAceRInt)
  
}

# Parallel version
cernia_parallel <- function(miRTarget, ExpData, mres, consider.miR.expr = "TRUE", minSharedmiR = 3, poscorcutoff = 0, 
                            padjustvaluecutoff = 0.01, padjustmethod = "BH", scorecutoff = 0.5) {
  
  if (consider.miR.expr == "TRUE") {
    miRTargetCandidate <- querymiRTargetbinding(ExpData, miRTarget, type = "all")
    miRTargetCandidate <- as.matrix(miRTargetCandidate)
  } else {
    miRTargetCandidate <- querymiRTargetbinding(ExpData, miRTarget, type = "target")
    miRTargetCandidate <- as.matrix(miRTargetCandidate)
  }
  m1 <- nrow(miRTargetCandidate)
  n1 <- ncol(miRTargetCandidate)
  
  colnames(ExpData) <- gsub("\\.", "-", colnames(ExpData))
  ExpDataNames <- colnames(ExpData)
  
  miR <- miRTargetCandidate[, 1]
  tar <- miRTargetCandidate[, 2]
  
  miRSym <- unique(miR)
  targetSym <- unique(tar)
  
  m2 <- length(targetSym)
  index <- t(utils::combn(m2, 2))
  
  Res <- foreach(i = seq_len(nrow(index)), .export = c("dtHybrid", "makeCluster", "detectCores", "recommendation", 
                                                       "graphWeights", "parMM", "stopCluster", "clusterApply")) %dopar% {	    
                                                         
                                                         Interin1 <- miRTargetCandidate[which(miRTargetCandidate[, 2] %in% targetSym[index[i, 1]]), 1]
                                                         Interin2 <- miRTargetCandidate[which(miRTargetCandidate[, 2] %in% targetSym[index[i, 2]]), 1]
                                                         
                                                         cm <- intersect(Interin1, Interin2)
                                                         SharedMREs <- mres[mres[, 2] %in% c(targetSym[index[i, 1]], targetSym[index[i, 2]]) & mres[, 1] %in% cm, ]
                                                         
                                                         M1 <- length(Interin1)
                                                         M2 <- length(Interin2)
                                                         M3 <- length(intersect(Interin1, Interin2))
                                                         M4 <- length(miRSym)
                                                         M5 <- 1 - phyper(M3 - 1, M2, M4 - M2, M1)
                                                         
                                                         tarExpIdx1 <- which(ExpDataNames %in% targetSym[index[i, 1]])
                                                         tarExpIdx2 <- which(ExpDataNames %in% targetSym[index[i, 2]])
                                                         
                                                         M6 <- cor.test(ExpData[, tarExpIdx1], ExpData[, tarExpIdx2])$estimate
                                                         
                                                         if (M3 >= minSharedmiR & M5 < padjustvaluecutoff & M6 > poscorcutoff &
                                                             nrow(SharedMREs) > 0) {
                                                           
                                                           C1 <- targetSym[index[i, 1]]
                                                           C2 <- targetSym[index[i, 2]]
                                                           
                                                           C3 <- M3
                                                           C4 <- M5
                                                           
                                                           # Score 1 for the fraction of coomon miRNAs
                                                           C5 <- log(M3/min(M1, M2))
                                                           
                                                           # Score 2 for the density of the MREs for all shared miRNAs
                                                           C6 <- sum(sapply(cm, function(miR) {
                                                             MREs <- SharedMREs[SharedMREs[, 1] == miR, ]
                                                             if (nrow(MREs) <= 0) return(1)
                                                             MREslr <- MREs[, c("gap_l", "gap_r")]
                                                             D <- abs(max(MREslr[, 1]) - min(MREslr[, 2]))
                                                             return(log(nrow(MREs)/D))
                                                           }))
                                                           
                                                           # Score 3 for the distribution of MREs of the putative RNA-RNA pairs
                                                           C7 = sum(sapply(cm, function(miR) {
                                                             positions <- SharedMREs[SharedMREs[, 1] == miR, c("gap_l", "gap_r")]
                                                             if (nrow(positions) <= 0) return(1)
                                                             return(log(abs(max(positions[, 1]) - min(positions[, 2]))^2/sum((positions[, 2] -
                                                                                                                                positions[, 1])^2)))
                                                           }))
                                                           
                                                           # Score 4 for the relation between the overall number of MREs for a putative miRNA sponge,
                                                           # compared with the number of miRNAs that yield these MREs
                                                           B <- nrow(SharedMREs)
                                                           if (B == length(unique(SharedMREs[, 1]))) {
                                                             C8 <- log(1/B)
                                                           } else {
                                                             C8 <- log((B - length(unique(SharedMREs[, 1])) - 1)/B)
                                                           }
                                                           
                                                           # Score 5 for the density of the hybridization energies related to MREs for all shared miRNAs
                                                           SharedMREs <- mres[mres[, 2] %in% c(targetSym[index[i, 1]], targetSym[index[i, 2]]) & mres[, 1] %in%
                                                                                cm, ]
                                                           C9 <- sum(sapply(cm, function(miR) {
                                                             MREs <- SharedMREs[SharedMREs[, 1] == miR, ]
                                                             if (nrow(MREs) <= 0) return(1)
                                                             MREslr <- MREs[, c("gap_l", "gap_r")]
                                                             D <- abs(max(MREslr[, 1]) - min(MREslr[, 2]))
                                                             return(log(sum(abs(MREs[, 3]))/D))
                                                           }))
                                                           
                                                           # Score 6 for the DT-Hybrid recommendation scores
                                                           cerna_recommendations <- dtHybrid(miRTargetCandidate)
                                                           C10 <- cerna_recommendations[targetSym[index[i, 1]], targetSym[index[i, 2]]]
                                                           
                                                           # Score 7 for the pairwise Pearson correlation between putative RNA-RNA pair expression data
                                                           C11 <- log(M6)
                                                           
                                                           C12 <- C5 + C6 + C7 + C8 + C9 + C10 + C11
                                                         } else {
                                                           C1 <- NA; C2 <- NA; C3 <- NA; C4 <- NA; C5 <- NA; C6 <- NA; C7 <- NA; C8 <- NA; C9 <- NA; C10 <- NA; C11 <- NA; C12 <- NA 	       
                                                         }
                                                         
                                                         tmp <- c(C1, C2, C3, C4, C5, C6, C7, C8, C9, C10, C11, C12)    
                                                         return(tmp)
                                                       }
  
  Res <- do.call(rbind, Res)
  
  # Extract RNA-RNA pair based on seven scores
  Res <- Res[apply(Res, 1, function(x) !all(is.na(x))), ]
  C13 <- (as.numeric(Res[, 12]) - min(as.numeric(Res[, 12])))/(max(as.numeric(Res[, 12])) - min(as.numeric(Res[, 12])))
  Res <- cbind(Res, C13)
  Res[, 4] <- p.adjust(as.numeric(Res[, 4]), method = padjustmethod)
  
  CERNIAceRInt <- Res[which(as.numeric(Res[, 4]) < padjustvaluecutoff & as.numeric(Res[, 13]) > scorecutoff), ]   
  
  if (is.vector(CERNIAceRInt)) {
    names(CERNIAceRInt) <- c("sponge_1", "sponge_2", "#shared miRNAs", "p.adjusted_value of shared miRNAs",
                             "Score 1", "Score 2", "Score 3", "Score 4", "Score 5", "Score 6", "Score 7", "Combined score",
                             "Normalized score")
  } else {
    colnames(CERNIAceRInt) <- c("sponge_1", "sponge_2", "#shared miRNAs", "p.adjusted_value of shared miRNAs",
                                "Score 1", "Score 2", "Score 3", "Score 4", "Score 5", "Score 6", "Score 7", "Combined score",
                                "Normalized score")
  }
  
  return(CERNIAceRInt)
  
}

## 8. Sparse Partial correlation ON Gene Expression (SPONGE) method
sponge_parallel <- function(miRTarget, ExpData, minSharedmiR = 3, poscorcutoff = 0, padjustvaluecutoff = 0.01,
                            padjustmethod = "BH", senscorcutoff = 0, null_model) {
  
  miRTargetCandidate <- querymiRTargetbinding(ExpData, miRTarget, type = "all")
  miRTargetCandidate <- as.matrix(miRTargetCandidate)
  m1 <- nrow(miRTargetCandidate)
  n1 <- ncol(miRTargetCandidate)
  
  colnames(ExpData) <- gsub("\\.", "-", colnames(ExpData))
  ExpDataNames <- colnames(ExpData)
  
  miR <- miRTargetCandidate[, 1]
  tar <- miRTargetCandidate[, 2]
  
  miRSym <- unique(miR)
  targetSym <- unique(tar)
  
  m2 <- length(targetSym)
  index <- t(utils::combn(m2, 2))    
   
  Res <- foreach(i = seq_len(nrow(index)), .packages = c("SPONGE")) %dopar% {	    
    
    Interin1 <- miRTargetCandidate[which(miRTargetCandidate[, 2] %in% targetSym[index[i, 1]]), 1]
    Interin2 <- miRTargetCandidate[which(miRTargetCandidate[, 2] %in% targetSym[index[i, 2]]), 1]
    
    M1 <- length(Interin1)
    M2 <- length(Interin2)
    M3 <- length(intersect(Interin1, Interin2))
    M4 <- length(miRSym)
    M5 <- 1 - phyper(M3 - 1, M2, M4 - M2, M1)
    
    if (M3 >= minSharedmiR & M5 < padjustvaluecutoff) {
      
      C1 <- targetSym[index[i, 1]]
      C2 <- targetSym[index[i, 2]]
      
      tarExpIdx1 <- which(ExpDataNames %in% targetSym[index[i, 1]])
      tarExpIdx2 <- which(ExpDataNames %in% targetSym[index[i, 2]])          
      miRExpIdx <- which(ExpDataNames %in% intersect(Interin1, Interin2))
      
      # Calculate sensitivity correlation of each RNA-RNA pair
      M6 <- cor.test(ExpData[, tarExpIdx1], ExpData[, tarExpIdx2])$estimate
      M7 <- cor.test(ExpData[, tarExpIdx1], ExpData[, tarExpIdx2])$p.value
      M8 <- corpcor::pcor.shrink(cbind(ExpData[, tarExpIdx1], ExpData[, tarExpIdx2],
                                       ExpData[, miRExpIdx]), verbose = FALSE)[1, 2]
      M9 <- M6 - M8
      
      C3 <- M3
      C4 <- M5
      C5 <- M6
      C6 <- M7
      C7 <- M8
      C8 <- M9
    } else {
      C1 <- NA; C2 <- NA; C3 <- NA; C4 <- NA; C5 <- NA; C6 <- NA; C7 <- NA; C8 <- NA 	       
    }
    
    tmp <- c(C1, C2, C3, C4, C5, C6, C7, C8)    
    return(tmp)
  }
  
  Res <- do.call(rbind, Res)
  
  # Extract RNA-RNA pairs with sensitivity correlation more than senscorcutoff.
  Res[, 4] <- p.adjust(as.numeric(Res[, 4]), method = padjustmethod)
  Res[, 6] <- p.adjust(as.numeric(Res[, 6]), method = padjustmethod)
  
  Res <- Res[which(as.numeric(Res[, 4]) < padjustvaluecutoff & as.numeric(Res[, 5]) > poscorcutoff & 
                     as.numeric(Res[, 6]) < padjustvaluecutoff & as.numeric(Res[, 8]) > senscorcutoff), ] 
  
  sponge_result <- data.frame(geneA = Res[, 1], 
                              geneB = Res[, 2], 
                              df = as.numeric(Res[, 3]), 
                              cor = as.numeric(Res[, 5]), 
                              pcor = as.numeric(Res[, 7]), 
                              mscor = as.numeric(Res[, 8]))
  
  sponge_result_null_model <- sponge_compute_p_values(sponge_result = sponge_result, null_model = null_model)
  sponge_result_fdr <- sponge_result_null_model[which(sponge_result_null_model$p.adj < padjustvaluecutoff), ]
  Res_fdr <- Res[which(sponge_result_null_model$p.adj < padjustvaluecutoff), ]
  
  SPONGEceRInt <- cbind(Res_fdr, sponge_result_fdr[, 8])
  
  if (is.vector(SPONGEceRInt)) {
    names(SPONGEceRInt) <- c("sponge_1", "sponge_2", "#shared miRNAs", "p.adjusted_value of shared miRNAs",
                             "correlation", "p.adjusted_value of positive correlation", "partial correlation", 
                             "sensitivity correlation", "p.adjusted_value of sensitivity correlation")
  } else {
    colnames(SPONGEceRInt) <- c("sponge_1", "sponge_2", "#shared miRNAs", "p.adjusted_value of shared miRNAs",
                                "correlation", "p.adjusted_value of positive correlation", "partial correlation", 
                                "sensitivity correlation", "p.adjusted_value of sensitivity correlation")
  }
  
  return(SPONGEceRInt)
}

## 9. Integrate method for miRNA sponge interactions by integrating different methods.
integrateMethod <- function(Interlist, Intersect_num) {
  
  if (length(Interlist) >= 2 & length(Interlist) >= Intersect_num) {
    Combcase <- t(utils::combn(length(Interlist), Intersect_num))
    Combnum <- dim(Combcase)[1]
    Integrate_Inter <- list()
    
    for (i in seq_len(Combnum)) {
      Interin <- do.call(rbind, lapply(Combcase[i, ], function(i) Interlist[[i]]))
      Interin_paste <- paste(Interin[, 1], Interin[, 2], sep = "-")
      Interin_table <- table(Interin_paste)
      Interin_names <- names(Interin_table)[which(Interin_table == Intersect_num)]
      Integrate_Inter[[i]] <- Interin[which(Interin_paste %in% Interin_names), ]
    }
    
    Integrate_res <- unique(do.call(rbind, Integrate_Inter))
    return(Integrate_res)
  } else {
    stop("Please check your input!\n")
  }
  
}

## Consolidating 15 functions: miRHomology, miRHomology_parallel, pc, pc_parallel, sppc, sppc_parallel,  
## ppc, ppc_parallel, hermes, hermes_parallel, muTaME, muTaME_parallel, cernia, cernia_parallel and sponge_parallel.
spongeMethod <- function(miRTarget, 
                         ExpData = NULL, 
                         mres = NULL, 
                         consider.miR.expr = "TRUE", 
                         minSharedmiR = 3, 
                         poscorcutoff = 0,
                         num_perm = 100, 
                         padjustvaluecutoff = 0.01, 
                         padjustmethod = "BH", 
                         senscorcutoff = 0.3, 
                         scorecutoff = 0.5, 
                         null_model,
                         method = c("miRHomology", "miRHomology_parallel", 
                                    "pc", "pc_parallel", "sppc", "sppc_parallel", 
                                    "ppc", "ppc_parallel", 
                                    "hermes", "hermes_parallel", 
                                    "muTaME", "muTaME_parallel", 
                                    "cernia", "cernia_parallel", 
                                    "sponge_parallel"),
                         num.cores = 2) {
  
  if (method == "miRHomology_parallel" | method == "pc_parallel" | method == "sppc_parallel" |
      method == "ppc_parallel" | method == "hermes_parallel" | method == "muTaME_parallel" |
      method == "cernia_parallel" | method == "sponge_parallel") {
    # get number of cores to run
    cores <- makeCluster(num.cores)
    registerDoParallel(cores)
  }
  
  if (method == "miRHomology") {
    ceRInt <- miRHomology(miRTarget, minSharedmiR, padjustvaluecutoff, padjustmethod)
    
  } else if (method == "miRHomology_parallel"){
    ceRInt <- miRHomology_parallel(miRTarget, minSharedmiR, padjustvaluecutoff, padjustmethod)
    
  } else if (method == "pc") {
    ceRInt <- pc(miRTarget, ExpData, consider.miR.expr, minSharedmiR, poscorcutoff,
                 padjustvaluecutoff, padjustmethod)
    
  } else if (method == "pc_parallel") {
    ceRInt <- pc_parallel(miRTarget, ExpData, consider.miR.expr, minSharedmiR, poscorcutoff,
                          padjustvaluecutoff, padjustmethod)
    
  } else if (method == "sppc") {
    ceRInt <- sppc(miRTarget, ExpData, minSharedmiR, poscorcutoff,
                   padjustvaluecutoff, padjustmethod, senscorcutoff)
    
  } else if (method == "sppc_parallel") {
    ceRInt <- sppc_parallel(miRTarget, ExpData, minSharedmiR, poscorcutoff,
                            padjustvaluecutoff, padjustmethod, senscorcutoff)
    
  } else if (method == "ppc") {
    ceRInt <- ppc(miRTarget, ExpData, minSharedmiR, num_perm, padjustvaluecutoff,
                  padjustmethod)
    
  } else if (method == "ppc_parallel") {
    ceRInt <- ppc_parallel(miRTarget, ExpData, minSharedmiR, num_perm, padjustvaluecutoff,
                           padjustmethod) 
    
  } else if (method == "hermes") {
    ceRInt <- hermes(miRTarget, ExpData, minSharedmiR, num_perm, padjustvaluecutoff,
                     padjustmethod)
    
  } else if (method == "hermes_parallel") {
    ceRInt <- hermes_parallel(miRTarget, ExpData, minSharedmiR, num_perm, padjustvaluecutoff,
                              padjustmethod)  
    
  } else if (method == "muTaME") {
    ceRInt <- muTaME(miRTarget, mres, minSharedmiR, padjustvaluecutoff, 
                     padjustmethod, scorecutoff)
    
  } else if (method == "muTaME_parallel") {
    ceRInt <- muTaME_parallel(miRTarget, mres, minSharedmiR, padjustvaluecutoff, 
                              padjustmethod, scorecutoff)
    
  } else if (method == "cernia") {
    ceRInt <- cernia(miRTarget, ExpData, mres, consider.miR.expr, minSharedmiR, poscorcutoff,
                     padjustvaluecutoff, padjustmethod, scorecutoff)
    
  } else if (method == "cernia_parallel") {
    ceRInt <- cernia_parallel(miRTarget, ExpData, mres, consider.miR.expr, minSharedmiR, poscorcutoff,
                              padjustvaluecutoff, padjustmethod, scorecutoff)
    
  } else if (method == "sponge_parallel") {        
    ceRInt <- sponge_parallel(miRTarget, ExpData, minSharedmiR, poscorcutoff,
                              padjustvaluecutoff, padjustmethod, senscorcutoff, null_model)	
    
  }
  
  if (method == "miRHomology_parallel" | method == "pc_parallel" | method == "sppc_parallel" |
      method == "ppc_parallel" | method == "hermes_parallel" | method == "muTaME_parallel" |
      method == "cernia_parallel" | method == "sponge_parallel") {
    # shut down the workers
    stopCluster(cores)
    stopImplicitCluster()
  }
  
  return(ceRInt)
}

## Inferring sample-specific miRNA sponge interaction network using sample 
## control variable strategy. 
sponge_sample_specific <- function(miRTarget, 
                                   ExpData = NULL, 
                                   mres = NULL, 
                                   consider.miR.expr = "TRUE", 
                                   minSharedmiR = 3, 
                                   poscorcutoff = 0,
                                   num_perm = 100, 
                                   padjustvaluecutoff = 0.01, 
                                   padjustmethod = "BH", 
                                   senscorcutoff = 0.3, 
                                   scorecutoff = 0.5, 
                                   null_model,
                                   method = c("pc", "pc_parallel", "sppc", 
                                              "sppc_parallel", "ppc", "ppc_parallel", 
                                              "hermes", "hermes_parallel", "cernia", 
                                              "cernia_parallel", "sponge_parallel"),
                                   num.cores = 2) {    
  
  # get number of cores to run
  cores <- makeCluster(num.cores)
  registerDoParallel(cores)    
  
  if (method == "pc") {
    ceRInt.all <- pc(miRTarget, ExpData, consider.miR.expr, minSharedmiR, poscorcutoff,
                     padjustvaluecutoff, padjustmethod)	
    
    ceRInt.single <- foreach(i = seq(nrow(ExpData)), .export = c("pc", "querymiRTargetbinding")) %dopar% {            
      pc(miRTarget, ExpData[-i, ], consider.miR.expr, minSharedmiR, poscorcutoff,  padjustvaluecutoff, padjustmethod)
    }	
    
  } else if (method == "pc_parallel") {
    ceRInt.all <- pc_parallel(miRTarget, ExpData, consider.miR.expr, minSharedmiR, poscorcutoff,
                              padjustvaluecutoff, padjustmethod)
    
    ceRInt.single <- foreach(i = seq(nrow(ExpData)), .export = c("pc_parallel", "querymiRTargetbinding")) %dopar% {            
      pc_parallel(miRTarget, ExpData[-i, ], consider.miR.expr, minSharedmiR, poscorcutoff,
                  padjustvaluecutoff, padjustmethod)
    }
    
  } else if (method == "sppc") {
    ceRInt.all <- sppc(miRTarget, ExpData, minSharedmiR, poscorcutoff,
                       padjustvaluecutoff, padjustmethod, senscorcutoff)
    
    ceRInt.single <- foreach(i = seq(nrow(ExpData)), .export = c("sppc", "querymiRTargetbinding")) %dopar% {            
      sppc(miRTarget, ExpData[-i, ], minSharedmiR, poscorcutoff,
           padjustvaluecutoff, padjustmethod, senscorcutoff)
    }
    
  } else if (method == "sppc_parallel") {
    ceRInt.all <- sppc_parallel(miRTarget, ExpData, minSharedmiR, poscorcutoff,
                                padjustvaluecutoff, padjustmethod, senscorcutoff)
    
    ceRInt.single <- foreach(i = seq(nrow(ExpData)), .export = c("sppc_parallel", "querymiRTargetbinding")) %dopar% {            
      sppc_parallel(miRTarget, ExpData[-i, ], minSharedmiR, poscorcutoff,
                    padjustvaluecutoff, padjustmethod, senscorcutoff)
    }
    
  } else if (method == "ppc") {
    ceRInt.all <- ppc(miRTarget, ExpData, minSharedmiR, num_perm, padjustvaluecutoff,
                      padjustmethod)
    
    ceRInt.single <- foreach(i = seq(nrow(ExpData)), .export = c("ppc", "querymiRTargetbinding", "predCor", "combpvalue")) %dopar% {            
      ppc(miRTarget, ExpData[-i, ], minSharedmiR, num_perm, padjustvaluecutoff,
          padjustmethod)
    }
    
  } else if (method == "ppc_parallel") {
    ceRInt.all <- ppc_parallel(miRTarget, ExpData, minSharedmiR, num_perm, padjustvaluecutoff,
                               padjustmethod)
    
    ceRInt.single <- foreach(i = seq(nrow(ExpData)), .export = c("ppc_parallel", "querymiRTargetbinding", "predCor", "combpvalue")) %dopar% {            
      ppc_parallel(miRTarget, ExpData[-i, ], minSharedmiR, num_perm, padjustvaluecutoff,
                   padjustmethod)
    }
    
  } else if (method == "hermes") {
    ceRInt.all <- hermes(miRTarget, ExpData, minSharedmiR, num_perm, padjustvaluecutoff,
                         padjustmethod)
    
    ceRInt.single <- foreach(i = seq(nrow(ExpData)), .export = c("hermes", "querymiRTargetbinding", "predCor", "combpvalue", "calCMI")) %dopar% {            
      hermes(miRTarget, ExpData[-i, ], minSharedmiR, num_perm, padjustvaluecutoff,
             padjustmethod)
    }
    
  } else if (method == "hermes_parallel") {
    ceRInt.all <- hermes_parallel(miRTarget, ExpData, minSharedmiR, num_perm, padjustvaluecutoff,
                                  padjustmethod)
    
    ceRInt.single <- foreach(i = seq(nrow(ExpData)), .export = c("hermes_parallel", "querymiRTargetbinding", "predCor", "combpvalue", "calCMI")) %dopar% {            
      hermes_parallel(miRTarget, ExpData[-i, ], minSharedmiR, num_perm, padjustvaluecutoff,
                      padjustmethod)
    }
    
  } else if (method == "cernia") {
    ceRInt.all <- cernia(miRTarget, ExpData, mres, consider.miR.expr, minSharedmiR, poscorcutoff,
                         padjustvaluecutoff, padjustmethod, scorecutoff)
    
    ceRInt.single <- foreach(i = seq(nrow(ExpData)), .export = c("cernia", "querymiRTargetbinding", "dtHybrid", "makeCluster", 
                                                                 "detectCores", "recommendation", "graphWeights", "parMM", "stopCluster", "clusterApply")) %dopar% {            
                                                                   cernia(miRTarget, ExpData[-i, ], mres, consider.miR.expr, minSharedmiR, poscorcutoff,
                                                                          padjustvaluecutoff, padjustmethod, scorecutoff)
                                                                 }
    
  } else if (method == "cernia_parallel") {
    ceRInt.all <- cernia_parallel(miRTarget, ExpData, mres, consider.miR.expr, minSharedmiR, poscorcutoff,
                                  padjustvaluecutoff, padjustmethod, scorecutoff)
    
    ceRInt.single <- foreach(i = seq(nrow(ExpData)), .export = c("cernia_parallel", "querymiRTargetbinding", "dtHybrid", "makeCluster", 
                                                                                              "detectCores", "recommendation", "graphWeights", "parMM", "stopCluster", "clusterApply")) %dopar% {            
                                                                                                cernia_parallel(miRTarget, ExpData[-i, ], mres, consider.miR.expr, minSharedmiR, poscorcutoff,
                                                                                                                padjustvaluecutoff, padjustmethod, scorecutoff)
                                                                                              }
    
  } else if (method == "sponge_parallel") {
    ceRInt.all <- sponge_parallel(miRTarget, ExpData, minSharedmiR, poscorcutoff,
                                  padjustvaluecutoff, padjustmethod, senscorcutoff, null_model)
    
    ceRInt.single <- foreach(i = seq(nrow(ExpData)), .packages = c("SPONGE"), .export = c("sponge_parallel", "querymiRTargetbinding")) %dopar% {            
      sponge_parallel(miRTarget, ExpData[-i, ], minSharedmiR, poscorcutoff,
                      padjustvaluecutoff, padjustmethod, senscorcutoff, null_model)
    }
    
  }
  
  # shut down the workers
  stopCluster(cores)
  stopImplicitCluster()    
  
  ceRInt.all.graph <- make_graph(c(t(ceRInt.all[, 1:2])), directed = FALSE)
  ceRInt.single.graph <- lapply(seq(ceRInt.single), function(i) make_graph(c(t(ceRInt.single[[i]][, 1:2])), directed = FALSE))
  ceRInt <- lapply(seq(ceRInt.single), function(i) as_data_frame((ceRInt.all.graph %m% ceRInt.single.graph[[i]]) %u% 
                                                                   (ceRInt.single.graph[[i]] %m% ceRInt.all.graph)))
  names(ceRInt) <- rownames(ExpData)
  
  return(ceRInt)
}

## Identifying sample-sample correlation network in terms of 
## sample-specific miRNA sponge networks.
sample_cor_network <- function(ceRNet, 
                               genes_num, 
                               method = "Simpson",
                               simcutoff = 0.5,
                               padjustvaluecutoff = 0.01,
                               padjustmethod = "BH",
                               num.cores = 2){
  
  # get number of cores to run
  cores <- makeCluster(num.cores)
  registerDoParallel(cores)
  
  ceRNet.graph <- lapply(seq(ceRNet), function(i) graph_from_data_frame(ceRNet[[i]], directed = FALSE))
  m <- length(ceRNet.graph)
  index <- t(utils::combn(m, 2))
  background <- genes_num * (genes_num - 1)/2
  
  Res <- foreach(i = seq_len(nrow(index)), .packages = c("igraph")) %dopar% {
    net1 <- ceRNet.graph[[index[i, 1]]]
    net2 <- ceRNet.graph[[index[i, 2]]]
    overlap <- length(E(net1 %s% net2))
    if (method == "Simpson") {
      sim <- overlap/min(length(E(net1)), length(E(net2)))
    } else if (method == "Jaccard") {	    
      sim <- overlap/length(E(net1 %u% net2))
    } else if (method == "Lin") {
      sim <- 2 * overlap/(length(E(net1)) + length(E(net2)))
    }
    sim.pvalue <- 1 - phyper(overlap - 1, length(E(net2)), background - length(E(net2)), length(E(net1)))
    sample_1 <- names(ceRNet)[index[i, 1]]
    sample_2 <- names(ceRNet)[index[i, 2]]
    tmp <- c(sample_1, sample_2, sim, sim.pvalue)    
    return(tmp)
  }
  
  # shut down the workers
  stopCluster(cores)
  stopImplicitCluster() 
  
  Res <- do.call(rbind, Res)
  
  Res[, 4] <- p.adjust(as.numeric(Res[, 4]), method = padjustmethod)    
  Res <- Res[which(as.numeric(Res[, 3]) > simcutoff & as.numeric(Res[, 4]) < padjustvaluecutoff), ]    
  
  if (is.vector(Res)) {
    names(Res) <- c("sample_1", "sample_2", "similarity", "p.adjusted_value of similarity")
  } else {
    colnames(Res) <- c("sample_1", "sample_2", "similarity", "p.adjusted_value of similarity")
  }
  
  return(Res)
}

## Validation of computationally predicted miRNA sponge interactions.
spongeValidate <- function(spongenetwork, directed = FALSE, Groundtruth) {
  
  spongenetwork_graph <- graph_from_data_frame(spongenetwork, directed = directed)
  Groundtruth_graph <- graph_from_data_frame(Groundtruth, directed = directed)
  Validated_interactions <- as_data_frame(spongenetwork_graph %s% Groundtruth_graph)
  colnames(Validated_interactions) <- c("sponge_1", "sponge_2")
  
  return(Validated_interactions)
}

## netModule function for identifying miRNA sponge modules from network. Possible methods
## include FN, MCL, LINKCOMM, MCODE, betweenness, infomap, prop, eigen, louvain, walktrap.
netModule <- function(spongenetwork, method = "MCL", directed = FALSE, modulesize = 3, save = FALSE) {
  
  if (method == "FN" | method == "MCL" | method == "MCODE") {
    spongenetwork_Cluster <- cluster(graph_from_data_frame(spongenetwork, directed = directed),
                                     method = method, directed = directed)
  } else if (method == "LINKCOMM") {
    edgelist <- get.edgelist(graph_from_data_frame(spongenetwork, directed = directed))
    spongenetwork_Cluster <- getLinkCommunities(edgelist, directed = directed)$nodeclusters
  }
  
  if (method == "FN" | method == "MCL") {
    spongenetwork_Cluster_result <- lapply(seq_len(max(spongenetwork_Cluster)),
                                           function(i) rownames(as.matrix(spongenetwork_Cluster))[which(spongenetwork_Cluster == i)])
    size <- unlist(lapply(seq_len(max(spongenetwork_Cluster)), function(i) length(spongenetwork_Cluster_result[[i]])))
    spongenetwork_Cluster_result <- lapply(which(size >= modulesize), function(i) spongenetwork_Cluster_result[[i]])
  } else if (method == "LINKCOMM") {
    spongenetwork_Cluster_result <- lapply(seq_len(max(c(spongenetwork_Cluster$cluster))), 
                                           function(i) as.character(spongenetwork_Cluster$node[which(c(spongenetwork_Cluster$cluster) == i)]))
    size <- unlist(lapply(seq_len(max(c(spongenetwork_Cluster$cluster))), function(i) length
                          (spongenetwork_Cluster_result[[i]])))
    spongenetwork_Cluster_result <- lapply(which(size >= modulesize), function(i) spongenetwork_Cluster_result[[i]])
  } else if (method == "MCODE") {
    spongenetwork_Cluster <- spongenetwork_Cluster + 1
    spongenetwork_Cluster_result <- lapply(seq_len(max(spongenetwork_Cluster)),
                                           function(i) rownames(as.matrix(spongenetwork_Cluster))[which(spongenetwork_Cluster == i)])
    size <- unlist(lapply(seq_len(max(spongenetwork_Cluster)), function(i) length(spongenetwork_Cluster_result[[i]])))
    spongenetwork_Cluster_result <- lapply(which(size >= modulesize), function(i) spongenetwork_Cluster_result[[i]])
  } else if (method == "betweenness") {
    spongenetwork_Cluster_result <- cluster_edge_betweenness(graph_from_data_frame(spongenetwork, directed = directed))
    size <- unlist(lapply(seq_len(length(spongenetwork_Cluster_result)), function(i) length(spongenetwork_Cluster_result[[i]])))
    spongenetwork_Cluster_result <- lapply(which(size >= modulesize), function(i) spongenetwork_Cluster_result[[i]])
  } else if (method == "infomap") {
    spongenetwork_Cluster_result <- cluster_infomap(graph_from_data_frame(spongenetwork, directed = directed))
    size <- unlist(lapply(seq_len(length(spongenetwork_Cluster_result)), function(i) length(spongenetwork_Cluster_result[[i]])))
    spongenetwork_Cluster_result <- lapply(which(size >= modulesize), function(i) spongenetwork_Cluster_result[[i]])
  } else if (method == "prop") {
    spongenetwork_Cluster_result <- cluster_label_prop(graph_from_data_frame(spongenetwork, directed = directed))
    size <- unlist(lapply(seq_len(length(spongenetwork_Cluster_result)), function(i) length(spongenetwork_Cluster_result[[i]])))
    spongenetwork_Cluster_result <- lapply(which(size >= modulesize), function(i) spongenetwork_Cluster_result[[i]])
  } else if (method == "eigen") {
    spongenetwork_Cluster_result <- cluster_leading_eigen(graph_from_data_frame(spongenetwork, directed = directed))
    size <- unlist(lapply(seq_len(length(spongenetwork_Cluster_result)), function(i) length(spongenetwork_Cluster_result[[i]])))
    spongenetwork_Cluster_result <- lapply(which(size >= modulesize), function(i) spongenetwork_Cluster_result[[i]])
  } else if (method == "louvain") {
    spongenetwork_Cluster_result <- cluster_louvain(graph_from_data_frame(spongenetwork, directed = directed))
    size <- unlist(lapply(seq_len(length(spongenetwork_Cluster_result)), function(i) length(spongenetwork_Cluster_result[[i]])))
    spongenetwork_Cluster_result <- lapply(which(size >= modulesize), function(i) spongenetwork_Cluster_result[[i]])
  } else if (method == "walktrap") {
    spongenetwork_Cluster_result <- cluster_walktrap(graph_from_data_frame(spongenetwork, directed = directed))
    size <- unlist(lapply(seq_len(length(spongenetwork_Cluster_result)), function(i) length(spongenetwork_Cluster_result[[i]])))
    spongenetwork_Cluster_result <- lapply(which(size >= modulesize), function(i) spongenetwork_Cluster_result[[i]])
  }
  
  if (save) {
    fileName <- paste("spongenetwork_Cluster_", method, ".txt", sep = "")
    for (i in seq_len(length(spongenetwork_Cluster_result))) {
        cat(c(i, "\t", length(spongenetwork_Cluster_result[[i]])), file = fileName, sep = "", append = TRUE)
        for (j in seq_len(length(spongenetwork_Cluster_result[[i]]))) {
          cat(c("\t", spongenetwork_Cluster_result[[i]][j]), file = fileName, sep = "", append = TRUE)
        }
        if (i != length(spongenetwork_Cluster_result)) {
          cat("\n", file = fileName, sep = "", append = TRUE)
        }
    }
  }
  
  return(spongenetwork_Cluster_result)
  
}

## Disease enrichment analysis of modules
moduleDEA <- function(Modulelist, OrgDb = "org.Hs.eg.db", ont = "DO", padjustvaluecutoff = 0.05,
    padjustedmethod = "BH") {

    entrezIDs <- lapply(seq_along(Modulelist), function(i) bitr(Modulelist[[i]], fromType = "SYMBOL",
        toType = "ENTREZID", OrgDb = OrgDb)$ENTREZID)

    entrezIDs <- lapply(seq_along(Modulelist), function(i) as.character(entrezIDs[[i]]))

    enrichDOs <- lapply(seq_along(Modulelist), function(i) enrichDO(entrezIDs[[i]], ont = ont, pvalueCutoff = padjustvaluecutoff,
        pAdjustMethod = padjustedmethod))

    enrichDGNs <- lapply(seq_along(Modulelist), function(i) enrichDGN(entrezIDs[[i]], pvalueCutoff = padjustvaluecutoff,
        pAdjustMethod = padjustedmethod))

    enrichNCGs <- lapply(seq_along(Modulelist), function(i) enrichNCG(entrezIDs[[i]], pvalueCutoff = padjustvaluecutoff,
        pAdjustMethod = padjustedmethod))

    return(list(enrichDOs, enrichDGNs, enrichNCGs))
}

## Functional GO, KEGG and Reactome enrichment analysis of modules
moduleFEA <- function(Modulelist, ont = "BP", KEGGorganism = "hsa", Reactomeorganism = "human",
    OrgDb = "org.Hs.eg.db", padjustvaluecutoff = 0.05, padjustedmethod = "BH") {

    entrezIDs <- lapply(seq_along(Modulelist), function(i) bitr(Modulelist[[i]], fromType = "SYMBOL",
        toType = "ENTREZID", OrgDb = OrgDb)$ENTREZID)

    entrezIDs <- lapply(seq_along(Modulelist), function(i) as.character(entrezIDs[[i]]))

    enrichGOs <- lapply(seq_along(Modulelist), function(i) enrichGO(entrezIDs[[i]], OrgDb = OrgDb,
        ont = ont, pvalueCutoff = padjustvaluecutoff, pAdjustMethod = padjustedmethod))

    enrichKEGGs <- lapply(seq_along(Modulelist), function(i) enrichKEGG(entrezIDs[[i]], organism = KEGGorganism,
        pvalueCutoff = padjustvaluecutoff, pAdjustMethod = padjustedmethod))

    enrichReactomes <- lapply(seq_along(Modulelist), function(i) enrichPathway(entrezIDs[[i]], organism = Reactomeorganism,
        pvalueCutoff = padjustvaluecutoff, pAdjustMethod = padjustedmethod))

    return(list(enrichGOs, enrichKEGGs, enrichReactomes))

}

## Survival analysis of modules
moduleSurvival <- function(Modulelist, ExpData, SurvData, devidePercentage = 0.5, plot = FALSE) {
  
  colnames(ExpData) <- gsub("\\.", "-", colnames(ExpData))
  ExpDataNames <- colnames(ExpData)
  
  myfit <- list()
  LogRank <- list()
  
  for (i in seq_along(Modulelist)) {
    Interin_Data <- cbind(SurvData[, seq(2, 3)], ExpData[, which(ExpDataNames %in% Modulelist[[i]])])
    Interin_Data <- na.omit(Interin_Data)
    
    try_mm <- try(coxph(survival::Surv(time, status) ~ ., data = data.frame(Interin_Data)),
                  silent = TRUE)
    if ("try-error" %in% class(try_mm))
      next
    
    mm <- coxph(survival::Surv(time, status) ~ ., data = data.frame(Interin_Data))
    
    Risk_score <- predict(mm, newdata = data.frame(Interin_Data), type = "risk")
    
    group <- rep("NA", dim(Interin_Data)[1])
    group[Risk_score > quantile(Risk_score, probs = devidePercentage)] <- "High"
    group[Risk_score <= quantile(Risk_score, probs = devidePercentage)] <- "Low"
    
    Data <- cbind(Interin_Data[, seq_len(2)], group)
    myfit[[i]] <- survfit(survival::Surv(time, status) ~ group, data = Data)
    
    sdf <- survdiff(survival::Surv(time, status) ~ group, data = Data)
    sdf.p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
    HR <- (sdf$obs[1]/sdf$exp[1])/(sdf$obs[2]/sdf$exp[2])
    HRlow95 <- exp(log(HR) - qnorm(0.975) * sqrt(1/sdf$exp[1] + 1/sdf$exp[2]))
    HRup95 <- exp(log(HR) + qnorm(0.975) * sqrt(1/sdf$exp[1] + 1/sdf$exp[2]))
    
    LogRank[[i]] <- c(sdf$chisq, sdf.p.val, HR, HRlow95, HRup95)
  }
  
  if (plot) {
    for (i in seq_along(myfit)) {
      if (!is.null(LogRank[[i]])) {
        dev.new()
        plot(myfit[[i]], lty = 1, col = c("red", "green"), main = paste("Module", i), xlab = "Time (Months)",
             ylab = "Probability of survival")
        
        legend("topright", legend = c("High risk group", "Low risk group"), lty = seq_len(2),
               col = c("red", "green"))
      }
    }
  }
  
  LogRank_res <- do.call(rbind, LogRank)
  
  if (length(myfit) >= 1) {
    colnames(LogRank_res) <- c("Chi-square", "p-value", "HR", "HRlow95", "HRup95")
    names(LogRank) <- seq_along(Modulelist)
    LogRank[sapply(LogRank, is.null)] <- NULL
    rownames(LogRank_res) <- paste("Module", names(LogRank))
  }
  
  return(LogRank_res)
}

