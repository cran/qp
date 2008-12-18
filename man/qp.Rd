\name{qp}
\alias{qp}
\title{
The package 'qp': summary information
}
\description{
This package has been replaced by 'qpgraph' available at the
Bioconductor project website (http://www.bioconductor.org)
}
\section{Data}{
This synthetic data forms still part of this package.
\itemize{
  \item \code{\link{jmlr06data}} synthetic data used in the referenced article
}
}

\section{Functions}{
\itemize{
  \item \code{qp.search} calculates the estimates of the non-rejection
        rates for every pair of variables. Now \code{qpNrr} in \code{qpgraph}.
  \item \code{qp.edge.prob} calculates the estimate of the non-rejection
        rate for a particular pair of variables, this function is also called
        by \code{qp.search}. Now \code{qpEdgeNrr} in \code{qpgraph}.
  \item \code{qp.ci.test} performs a test for conditional independence. Now \code{qpCItest} in \code{qpgraph}.
  \item \code{qp.analyse} provides some exploratory analyses on the
        output of \code{qp.search}. This function does not exist anymore.
  \item \code{qp.clique} calculates the maximum clique size as a function
        of the minimum threshold on the non-rejection rate for removing an edge. Now \code{qpClique} in \code{qpgraph}.
  \item \code{qp.hist} shows a histogram of the estimated non-rejection
        rates obtained through \code{qp.search}. Now \code{qpHist} in \code{qpgraph}.
  \item \code{qp.graph} returns the qp-graph, in the form of an incidence
        matrix, resulting of thresholding the non-rejection rates in the output
        of \code{qp.search}. Now \code{qpGraph} in \code{qpgraph}.
  \item \code{qp.matrix.image} makes an image plot of the absolute value
        of an inverse correlation matrix. This function does not exist anymore.
  \item \code{qp.get.cliques} finds the set of cliques of an undirected
        graph. Now \code{qpGetCliques} in \code{qpgraph}.
}
The package provides an implementation of the procedures described by
Castelo and Roverato (2006) and is a contribution to the gR-project
described by Lauritzen (2002).
}
\section{Authors}{
Robert Castelo, Departament de Ci\`encies Experimentals i de la Salut,
Universitat Pompeu Fabra, Barcelona, Spain.

Alberto Roverato, Dipartimento di Scienze Statistiche,
Universit\`a di Bologna, Italy.
}
\references{
Lauritzen, S. L. (2002). gRaphical Models in R.
\emph{R News}, 3(2)39.

Castelo, R. and Roverato, A. (2006). A robust procedure for Gaussian graphical model
search from microarray data with p larger than n, \emph{J. Mach. Learn. Res.}, 7:2621-2650
}
\keyword{multivariate}
\keyword{models}
\keyword{graphs}
