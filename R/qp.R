# qp package - this R code implements part of the R qp package described
# in R. Castelo and A. Roverato. A robust procedure for Gaussian graphical
# model search from microarray data with p larger than n, Journal of Machine
# Learning Research, accepted for publication.
#
# Copyright (C) 2006 R. Castelo and A. Roverato
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, you can obtain one via WWW at
# http://www.gnu.org/copyleft/gpl.html, or by writing to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.





# function: qp.search
# purpose: carry out the search for the q-order partial correlation
#          graph
# parameters: S - sample covariance matrix
#             N - sample size
#             q - partial-correlation order
#             T - number of tests per adjacency
#             significance - significance level of each test
#             binary - flag set to TRUE when using the C implementation
# return: a list with two members, a matrix A with the acceptance test
#         counts for each adjacency and the number T of tests per adjacency 

qp.search <- function(S, N, q=0, T=500, significance=0.05, binary=TRUE) {

  if (binary == TRUE) {
    return(list(A=qp.fast.search(S, N, q, T, significance),T=T))
  }

  n.var <- nrow(S)
  if (q < 0 || q > n.var-2)
    stop(paste("q=",q," > n.var-2=",n.var-2))

  A <- matrix(0,n.var,n.var)
  ppct <- -1
  k <- 0
  t <- n.var*(n.var-1)/2;
  for (i in 1:(n.var-1)) {
    for(j in (i+1):n.var) {
      A[i,j] <- qp.edge.prob(S, N, i, j, q, T, significance,FALSE)
      k <- k + 1
      pct <- floor((k*100)/t)
      if (pct != ppct) {
        if (pct %% 10 == 0) {
          cat(pct)
        } else {
          cat(".")
        }
        ppct = pct
      }
    }
  }

  return(list(A=A,T=T))
}



# function: qp.edge.prob
# purpose: calculate the probability of the edge as the number of tests
#          that accept the null hypothesis of independence given the
#          q-order conditionals
# parameters: S - sample covariance matrix
#             N - sample size
#             i - vertex
#             j - vertex
#             q - partial-correlation order
#             T - number of tests to perform
#             significance - significance level of each test
#             binary - flag set to TRUE when using the C implementation

qp.edge.prob <- function(S, N, i=1, j=2, q=0, T=500, significance=0.05, binary=TRUE) {

  if (binary == TRUE) {
    return(qp.fast.edge.prob(S, N, i, j, q, T, significance));
  }

  n.var  <- nrow(S)
  pop    <- (1:n.var)[-c(i, j)]

  thr     <- qt(p=1-(significance/2),df=N-q-2,lower.tail=TRUE,log.p=FALSE)
  lambda  <- c()
  for (k in 1:T) {
    sp <- sample(pop, q, rep=F)
    cit <- qp.ci.test(S, N, i, j, sp, binary=FALSE)
    lambda  <- c(lambda,abs(cit$t.value))
  }

  return(sum(lambda<thr))
}



# function: qp.ci.test
# purpose: perform a test for conditional independence between variables
#          indexed by i and j given the conditioning set Q
# parameters: S - sample covariance matrix
#             N - sample size
#             i - vertex
#             j - vertex
#             Q - vector of vertices forming the conditioning set
#             binary - flag set to TRUE when using the C implementation
# return: a list with two members, the t-statistic value and the p-value
#         on rejecting the null hypothesis of independence

qp.ci.test <- function(S, N, i=1, j=2, Q=c(), binary=TRUE) {

  if (binary == TRUE) {
    return(qp.fast.ci.test(S, N, i, j, Q));
  }

  q       <- length(Q)
  Mmar    <- S[c(i, j, Q), c(i, j, Q)]
  S11     <- Mmar[1,1]
  S12     <- Mmar[1,-1]
  S21     <- Mmar[-1,1]
  S22     <- Mmar[-1,-1]
  S22inv  <- solve(S22)
  betahat <- S12 %*% S22inv[,1]
  sigma   <- sqrt((S11 - S12 %*% S22inv %*% S21) * (N - 1) / (N - q - 2))
  se      <- sigma * sqrt(S22inv[1,1] / (N - 1))
  t.value <- betahat / se
  p.value <- 2*(1-pt(abs(t.value), N - q - 2))

  return(list(t.value=t.value,p.value=p.value))
}



# function: qp.analyse
# purpose: perform a first exploratory analysis of the results of qp.search
# parameters: qp.output - output of qp.search
#             threshold - threshold for edge removal
#             largest.clique - calculate the dimension of the largest clique
#             plot.image - flag that when set it plots the incidence matrix
#                          resulting of applying the threshold on the
#                          non-rejection rates
#             exact.calculation - flag that when set to TRUE the exact maximum clique
#                                 size is calculated and when set to FALSE a lower
#                                 bound is calculated instead
#             approximation.iterations - number of iterations performed to calculate
#                                        the lower bound on the clique number of
#                                        each graph (exact.calculation is FALSE)
# return: a matrix with the nr. of selected edges, the nr. of edges of the
#         complete graph and the percentage of selected edges. When
#         largest.clique=TRUE it gives also the size of the largest clique
#         and when plot.image=TRUE it plots the incidence matrix

qp.analyse <- function(qp.output, threshold, largest.clique=TRUE, plot.image=TRUE,
                       exact.calculation=FALSE, approximation.iterations=100) {
  M <- qp.output$A
  grayscale <- gray(1:99 / 99)

  M <- M + t(M) # make M symmetric

  if (is.null(qp.output$T)) { # the resulting incidence matrix has no loops
    diag(M) <- max(M)
  } else {
    diag(M) <- qp.output$T
  }

  Mflip <- M[nrow(M):1,] # the image function belows works bottom-up rows

  if (is.null(qp.output$T)) {
    Mflip <- (Mflip*100) / max(Mflip)
    M <- (M*100) / max(M)
  } else {
    Mflip <- (Mflip*100) / qp.output$T
    M <- (M*100) / qp.output$T
  }

  threshold <- threshold*100
  n.par <- (sum(Mflip <= threshold)-nrow(Mflip))/2
  Mflip[Mflip > threshold] <- 100
  Mflip[Mflip <= threshold] <- 0
  n.par.sat <- nrow(Mflip)*(nrow(Mflip)-1)/2
  outp <- matrix(c(n.par, n.par.sat, round(n.par*100/n.par.sat, 0)))
  rownames(outp) <- c("numb. of selected edges",
                      "edges of complete graph",
                      "percentage of sel. edges")

  if (plot.image) { 
    image(list(x=c(1:nrow(Mflip)), y=c(1:nrow(Mflip)), z=Mflip), col=grayscale, axes=F)
    nr <- nrow(Mflip)
    nr7<-round(nr/7, 0)
    atx=seq(nr7, nr-nr7+1, by=nr7)
    axis(3, at=atx)
    axis(2, at=nr-atx+1, lab=atx)
    box()
  }

  if (largest.clique == TRUE) {
    label <- "exact maximum clique size"
    if (exact.calculation == FALSE) {
      label <- "lower bound on the maximum clique size"
    }
    I <- xor(M <= threshold,diag(1,nrow(M))==1) 
    dimnames(I) <- list(1:length(I[,1]),1:length(I[1,]))
    maxsize <- qp.clique.number(I,exact.calculation,approximation.iterations)
    outp <- rbind(outp,maxsize)
    rownames(outp) <- c(rownames(outp)[1:(length(outp[,1])-1)],label)
  }

  return(outp)
}



# function: qp.clique
# purpose: plot the relationship between non-rejection rate and maximum clique size
# parameters: qp.output - output from qp.search
#             N - sample size
#             threshold.lim - range of threshold values
#             breaks - either a number of threshold bins or a vector of threshold
#                      breakpoints
#             plot.image - flag that when set it makes a plot of the result
#             exact.calculation - flag that when set to TRUE the exact maximum clique
#                                 size is calculated and when set to FALSE a lower
#                                 bound is calculated instead
#             approximation.iterations - number of iterations performed to calculate
#                                        the lower bound on the clique number of
#                                        each graph (exact.calculation is FALSE)
# return: a list with two members, the threshold on the non-rejection rate that provides
#         the maximum clique size that is strictly smaller than the sample size N, and
#         the resulting maximum clique size

qp.clique <- function(qp.output, N, threshold.lim=c(0,1), breaks=5, plot.image=TRUE,
                      exact.calculation=FALSE,approximation.iterations=100) {
  if (length(breaks) == 1) {
    len <- threshold.lim[2] - threshold.lim[1]
    br <- seq(threshold.lim[1],threshold.lim[2],by=len/breaks)
  } else {
    br <- breaks
    threshold.lim = range(br)
  }

  M <- qp.output$A
  M=M+t(M) # make M symmetric

  if (is.null(qp.output$T)) {
    M <- M / max(M)
  } else {
    M <- M / qp.output$T
  }

  maxclqszeunderN <- 0
  thrmaxclqszeunderN <- 0
  mpctedclqsze <- matrix(rep(0,length(br)*2),nrow=length(br),ncol=2)
  n.var <- nrow(M)
  n.par.sat <- n.var*(n.var-1)/2

  for (i in 1:length(br)) {
    threshold <- br[i]
    n.par <- (sum(M <= threshold) - n.var)/2
    I <- xor(M <= threshold,diag(1,n.var)==1) 
    dimnames(I) <- list(1:length(I[,1]),1:length(I[1,]))
    maxsize <- qp.clique.number(I,exact.calculation,approximation.iterations)
    mpctedclqsze[i,] <- c(round(n.par*100/n.par.sat,digits=0),maxsize)
    if (maxsize > maxclqszeunderN && maxsize < N) {
      maxclqszeunderN <- maxsize
      thrmaxclqszeunderN <- threshold
    }
  }

  linetype <- 1
  label <- "exact maximum clique size"
  if (exact.calculation == FALSE) {
    linetype <- 2
    label <- "lower bound on the maximum clique size"
  }

  if (plot.image == TRUE) {
    plot(br,mpctedclqsze[,2],type="o",xlim=threshold.lim,
         ylim=range(0,N,mpctedclqsze[,2]),lwd=2,lty=linetype,
         xlab="threshold",ylab="maximum clique size",
         main="maximum clique size as function of threshold")
    abline(h=N,col="red",lwd=2)
    text(br,mpctedclqsze[,2],lab=paste(mpctedclqsze[,1],"%",sep=""),pos=1,cex=.7)
    legend(min(threshold.lim),max(N,mpctedclqsze[,2]),label,lty=linetype,lwd=2,pch=1)
  }

  return(list(threshold=thrmaxclqszeunderN,size=maxclqszeunderN))
}



# function: qp.clique.number
# purpose: calculate the size of the largest maximal clique in the given graph
# parameters: I - incidence matrix of the graph
#             exact.calculation - flag that when set to TRUE the exact maximum clique
#                                 size is calculated and when set to FALSE a lower
#                                 bound is calculated instead
#             approximation.iterations - number of iterations performed to calculate
#                                        the lower bound on the clique number of
#                                        each graph (exact.calculation is FALSE)
# return: the size of the largest maximal clique in the given graph, also known as
#         its clique number

qp.clique.number <- function(I,exact.calculation=FALSE,approximation.iterations=100) {
  n.var <- nrow(I)

  if (exact.calculation == TRUE) {

    I <- I == 1 # make sure we get a boolean matrix
    clqlst <- qp.get.cliques(I,binary=TRUE)
    clique.number <- sort(as.numeric(as.matrix(lapply(clqlst,length))),decreasing=TRUE)[1]

  } else {

    clique.number <- 0
    I <- I + 0 # make sure we get a 0-1 matrix
    deg <- sort(rowsum(I,rep(1,n.var)),index.return=TRUE,decreasing=TRUE) # order by degree

    for (i in 1:approximation.iterations) {

      pdeg <- deg$ix
      if (i %% n.var + 1 > 1) {
        sset <- sample(1:n.var,i %% n.var + 1,replace=FALSE) # we alter the order of the ranking
        ssetelem <- pdeg[sset]                               # by degree with increasing levels
        ssetelem <- sample(ssetelem)                         # of randomness cyclically
        pdeg[sset] <- ssetelem
      }
      clq <- c(pdeg[1])
      j <- 2
      for (j in 2:n.var) {
        v <- pdeg[j]
        clq2 <- c(clq,v)
        if (sum(I[clq2,clq2]) == length(clq2)*length(clq2)-length(clq2)) {
          clq <- clq2
        }
      }
                                                                                            
      if (length(clq) > clique.number) {
        clique.number <- length(clq)
      }

    }

  }

  return(clique.number)
}



# function: qp.hist
# purpose: output a histogram of the non-rejection rates
# parameters: qp.output - output from qp.search
#             IC - inverse correlation matrix from the generative graph
#             prob - flag that when set to true the histograms show
#                    densities, otherwise they show absolute frequencies (counts)
# return: none

qp.hist <- function(qp.output, IC=NULL, prob=FALSE) {
  G <- qp.output$A
  N.loop <- qp.output$T
  # all
  x=G[col(G)>row(G)]
  y=x/N.loop
  xl=range(y)
  if(is.null(IC)){
    hist(y,  col="yellow", main="all estimated\nnon-rejection rates", xlab="non-rejection rate", prob=prob)
  } else {
    # only beta
    T=G
    T[IC==0]=-1
    xbeta=T[col(T)>row(T)]
    xbeta=xbeta[xbeta!=-1]
    ybeta=xbeta/N.loop
    # not beta
    T=G
    T[IC!=0]=-1
    xnotbeta=T[col(T)>row(T)]
    xnotbeta=xnotbeta[xnotbeta!=-1]
    ynotbeta=xnotbeta/N.loop
    # plots
    split.screen(c(2, 2))
    screen(1)
    H=hist(y,  prob=prob, col="yellow", xlim=xl, main="all estimated\nnon-rejection rates", plot=F)
    if(prob){
       yl=NULL
    }else{
       yl=c(0, max(H$counts))
    }
    hist(y,  prob=prob, col="yellow", xlim=xl, ylim=yl, main="all estimated\nnon-rejection rates", xlab="non-rejection rates")
    screen(2)
    boxplot(ybeta,ynotbeta, col=c("red", "cyan"),ylab="non-rejection rate",naxt="n")
    axis(1,at=1,"present edges")
    axis(1,at=2,"missing edges")
    screen(3)
    hist(ybeta, prob=prob,  col="yellow", xlim=xl, ylim=yl, main="present edges\nnon-rejection rates", xlab="non-rejection rates",breaks=length(H$breaks))
    screen(4)
    hist(ynotbeta,  prob=prob, col="yellow", xlim=xl, ylim=yl, main="missing edges\nnon-rejection rates", xlab="non-rejection rates",breaks=length(H$breaks))
    close.screen(all=T)
  }
}



# function: qp.graph
# purpose: output the incidence matrix resulting of thresholding
#          the output of qp.search
# parameters: qp.output - output of qp.search
#             threshold - threshold for edge removal
# return: incidence matrix

qp.graph <- function(qp.output,threshold) {

  M <- qp.output$A
  M=M+t(M) # make M symmetric
  n.var <- nrow(M)

  if (is.null(qp.output$T)) {
    M <- M / max(M)
  } else {
    M <- M / qp.output$T
  }

  I <- xor(M <= threshold,diag(1,n.var)==1) 

  return(I)
}



# function: qp.matrix.image
# purpose: make an image plot of the absolute value of an inverse
#          correlation matrix and reports the number of edges of
#          the corresponding independence graph
# parameters: M - the matrix
#             col - flag that when set to NULL a grey scale is used
#             plot - flag that when set to FALSE doesn't plot the matrix
# return: number of present edges, number of edges in the complete
#         graph, percentage of present edges

qp.matrix.image <- function(M, col=NULL, plot=TRUE) {

  if (is.null(col)) col <- gray(99:1 / 99)
  M<-M[nrow(M):1,]
  if (plot){
    image(list(x=c(1:nrow(M)), y=c(1:nrow(M)), z=abs(M)), col=col, axes=F)
    nr <- nrow(M)
    nr7<-round(nr/7, 0)
    atx=seq(nr7, nr-nr7+1, by=nr7)
    axis(3, at=atx)
    axis(2, at=nr-atx+1, lab=atx)
    box()
  }
  n.par <- (sum(M!=0)-nrow(M))/2
  n.par.sat <- nrow(M)*(nrow(M)-1)/2
  outp <- matrix(c(n.par, n.par.sat, round(n.par*100/n.par.sat, 0)))
  rownames(outp) <- c("numb. of present edges", "edges of complete graph", "percentage of pres. edges")

  return(outp)
}



# function: qp.get.cliques (exported) and extend (internal)
# purpose: find the set of (maximal) cliques of an undirected graph.
#          this function implements the algorithm in
#
#          Bron, C. and Kerbosch, J.
#          Finding all cliques of an undirected graph,
#          Communications of the ACM, 16:575-577, 1973.
# parameters: I - incidence matrix of the graph
#             binary - flag set to TRUE when using the C implementation.
#                      note that because the algorithm is recursive (through
#                      the 'extend' function) the R code implementation may
#                      easily reach the maximum number of recursive calls
# return: a list of maximal cliques

qp.get.cliques <- function(I,binary=TRUE) {

  if (binary == TRUE) {
    return(qp.fast.getcliques(I))
  }

  N <- dim(I)[1]
  I[diag(N) == 1] <- TRUE
  ALL <- 1:N
  compsub <- rep(NA,N)
  c <- 0
  pos <- NA
 
  r <- extend(I,compsub,c,pos,ALL,0,N,list())
  r$clqlst
}

extend <- function(connected, compsub, c, pos, old, ne, ce, clqlst) {
                                                                                                
  new <- 1:ce
  minnod <- ce
  nod <- 0
  i <- 1
                                                                                                
  while (i <= ce && minnod != 0) {
    p <- old[i]
    count <- 0
    j <- ne
    j <- j + 1
    while (j <= ce && count < minnod) {
      if (!connected[p,old[j]]) {
        count <- count + 1
        pos <- j
      }
      j <- j + 1
    }
    if (count < minnod) {
      fixp <- p
      minnod <- count
      if (i <= ne) {
        s <- pos
      } else {
        s <- i
        nod <- 1
      }
    }
    i <- i + 1
  }

  nod <- minnod + nod
  while (nod >= 1) { # BEGIN BACKTRACKCYCLE
    p <- old[s]
    old[s] <- old[ne+1]
    sel <- old[ne+1] <- p
    newne <- 0
    i <- 1
    while (i <= ne) {
      if (connected[sel,old[i]]) {
        newne <- newne + 1
        new[newne] <- old[i]
      }
      i <- i + 1
    }
    newce <- newne
    i <- ne + 1
    i <- i + 1
    while (i <= ce) {
      if (connected[sel,old[i]]) {
        newce <- newce + 1
        new[newce] <- old[i]
      }
      i <- i + 1
    }
    c <- c + 1
    compsub[c] <- sel
    if (newce == 0) {
      clq <- c()
      for (loc in 1:c) {
        clq <- c(clq,compsub[loc])
      }
      clqlst[[length(clqlst)+1]] <- clq
    } else {
      if (newne < newce) {
        r <- extend(connected,compsub,c,pos,new,newne,newce,clqlst)
        compsub <- r$compsub
        c <- r$c
        pos <- r$pos
        clqlst <- r$clqlst
      }
    }
    c <- c - 1 # REMOVE FROM COMPSUB
    ne <- ne + 1 # ADD TO NOT
    if (nod > 1) { # SELECT A CANDIDATE DISCONNECTED TO THE FIXED POINT
      s <- ne
      repeat {
        s <- s + 1
        if (!connected[fixp,old[s]]) break
      }
    } # END SELECTION
    nod <- nod - 1
  } # END BACKTRACKCYCLE
                                                                                                
  list(compsub=compsub,c=c,pos=pos,clqlst=clqlst)
}


#############################################
# ENTRY POINTS OF THE C CODE OF THE PACKAGE #
#############################################

qp.fast.search <- function(S, N, q=0, T=500, significance=0.05){
  return(.Call("qp_fast_search",S,as.integer(N),as.integer(q),as.integer(T),as.double(significance)))
}

qp.fast.edge.prob <- function(S, N, i, j, q, T, significance=0.05){
  return(.Call("qp_fast_edge_prob",S,as.integer(N),as.integer(i),as.integer(j),as.integer(q),as.integer(T),as.double(significance)))
}

qp.fast.ci.test <- function(S, N, i, j, C=c()){
  return(.Call("qp_fast_ci_test",S,as.integer(N),as.integer(i),as.integer(j),C))
}

qp.fast.getcliques <- function(I) {
  return(.Call("qp_fast_get_cliques",I))
}

######################
# INTERNAL FUNCTIONS #
######################

# function: qp.d.scale
# purpose: scales a PD matrix
# parameters: V - the matrix
# return: the scaled matrix

qp.d.scale <- function(V){
  d <- 1/sqrt(diag(V))
                                                                                                  
  return(V*outer(d, d))
}
