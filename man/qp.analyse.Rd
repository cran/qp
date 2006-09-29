\name{qp.analyse}
\alias{qp.analyse}
\title{
Performs some exploratory analyses on the q-partial graph
}
\description{
Using the output of \code{\link{qp.search}} this function provides some
exploratory analyses on the resulting q-partial graph.
}
\usage{
qp.analyse(qp.output, threshold, largest.clique=TRUE, plot.image=TRUE)

}
\arguments{
  \item{qp.output}{output of qp.search}
  \item{threshold}{threshold on the minimum non-rejection rate required
                   for edge removal}
  \item{largest.clique}{when this flag is set to TRUE it calculates the
                        dimension of the largest clique}
  \item{plot.image}{when this flag is set it plots the incidence matrix
                    resulting of thresholding the non-rejection rates}
}
\details{
Returns an object of the class matrix showing the number of selected edges,
the number of edges of the complete graph and the percentage of selected
edges. When largest.clique=TRUE it gives also the size of the largest clique
and when plot.image=TRUE it plots the incidence matrix resulting of
thresholding the non-rejection rates.

Beware that the combination or arguments largest.clique=TRUE and a
threshold > 0.95 may result in a very long time of computation due to the
NP-hardness of the problem of calculating the size of the largest clique.
}
\references{Castelo, R. and Roverato, A. (2006). A robust procedure for Gaussian
graphical model search from microarray data with p larger than n, \emph{J. Mach.
Learn. Res., accepted}}

\author{Robert Castelo and Alberto Roverato}
\seealso{\code{\link{qp.search}}}
\examples{
data(jmlr06data)
qp.analyse(qp.out.bd5.N20.q10,threshold=0.9,largest.clique=TRUE)
}
\keyword{models}
\keyword{multivariate}
\keyword{graphs}
