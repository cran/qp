\name{qp.graph}
\alias{qp.graph}
\title{
Incidence matrix of the qp-graph
}
\description{
Using the output of \code{\link{qp.search}} this function returns
the qp-graph, in the form of an incidence matrix, resulting of thresholding
the non-rejection rates in the output of \code{\link{qp.search}}
}
\usage{
qp.graph(qp.output, threshold)

}
\arguments{
  \item{qp.output}{output of \code{\link{qp.search}}}
  \item{threshold}{threshold on the non-rejection rate}
}
\references{Castelo, R. and Roverato, A. (2006). A robust procedure for
Gaussian graphical model search from microarray data with p larger than n,
\emph{J. Mach. Learn. Res.}, 7:2621-2650}

\author{Robert Castelo and Alberto Roverato}
\seealso{\code{\link{qp.search}}}
\examples{
data(jmlr06data)
I <- qp.graph(qp.out.bd5.N20.q10,threshold=0.9)
sprintf("the graph has \%.0f edges\n",sum(I)/2)
}
\keyword{models}
\keyword{multivariate}
\keyword{graphs}
