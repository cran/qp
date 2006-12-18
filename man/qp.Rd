\name{qp}
\alias{qp}
\title{
The package 'qp': summary information
}
\description{
This package provides functions for implementing the q-order partial-correlation
graph search algorithm, q-partial, or qp, algorithm for short.
The qp algorithm is a robust procedure for structure learning of undirected
Gaussian graphical Markov models (UGGMMs) from "small n, large p" data, that
is, multivariate normal data coming from a number of random variables p larger
than the number of multidimensional data points n as in the case of, e.g.,
microarray data.
}
\section{Data}{
\itemize{
  \item \code{\link{jmlr06data}} synthetic data used in the referenced article
}
}

\section{Functions}{
\itemize{
  \item \code{\link{qp.search}} calculates the estimates of the non-rejection
        rates for every pair of variables
  \item \code{\link{qp.edge.prob}} calculates the estimate of the non-rejection
        rate for a particular pair of variables, this function is also called
        by \code{\link{qp.search}}
  \item \code{\link{qp.ci.test}} performs a test for conditional independence
  \item \code{\link{qp.analyse}} provides some exploratory analyses on the
        output of \code{\link{qp.search}}
  \item \code{\link{qp.clique}} calculates the maximum clique size as a function
        of the minimum threshold on the non-rejection rate for removing an edge
  \item \code{\link{qp.hist}} shows a histogram of the estimated non-rejection
        rates obtained through \code{\link{qp.search}}
  \item \code{\link{qp.graph}} returns the qp-graph, in the form of an incidence
        matrix, resulting of thresholding the non-rejection rates in the output
        of \code{\link{qp.search}}
  \item \code{\link{qp.matrix.image}} makes an image plot of the absolute value
        of an inverse correlation matrix
  \item \code{\link{qp.get.cliques}} finds the set of cliques of an undirected
        graph
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
