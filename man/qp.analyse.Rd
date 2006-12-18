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
qp.analyse(qp.output, threshold, largest.clique=TRUE,
           plot.image=TRUE, exact.calculation=FALSE,
           approximation.iterations=100)

}
\arguments{
  \item{qp.output}{output of qp.search}
  \item{threshold}{threshold on the minimum non-rejection rate required
                   for edge removal}
  \item{largest.clique}{when this flag is set to TRUE it calculates the
                        size of the largest clique}
  \item{plot.image}{when this flag is set it plots the incidence matrix
                    resulting of thresholding the non-rejection rates}
  \item{exact.calculation}{when this flag is set to \code{TRUE}, the exact
                           maximum clique size is calculated and when set
                           to \code{FALSE} a lower bound is calculated
                           instead. It applies only when
                           \code{largest.clique=TRUE}}
  \item{approximation.iterations}{number of iterations performed to calculate
                                  the lower bound on the clique number of
                                  each graph. It applies only when
                                  \code{largest.clique=TRUE} and \\
                                  \code{exact.calculation=FALSE}}
}
\details{
Returns an object of the class matrix showing the number of selected edges,
the number of edges of the complete graph and the percentage of selected
edges. When largest.clique=TRUE it gives also the size of the largest clique
and when plot.image=TRUE it plots the incidence matrix resulting of
thresholding the non-rejection rates.

Beware that setting \code{largest.clique=TRUE} and \code{exact.calculation=TRUE}
when giving breakpoints between 0.95 and 1.0 (which may result into very dense
graphs) can lead to a very long time of computation due to the NP-completeness
of the problem of calculating the size of the largest clique which is therefore
bounded by an exponential growth of the running time as function of the
graph density (cf.~Pardalos and Xue, 1994).

The lower bound on the maximum clique size is calculated by ranking the
vertices by their connectivity degree, put the first vertex in a set and
go through the rest of the ranking adding those vertices to the set that
form a clique with the vertices currently within the set. Once the entire
ranking has been examined a large clique should have been built and hopefully
the largest one. This process is repeated a number of times
(\code{approximation.iterations}) each of which the ranking is altered with
increasing levels of randomness acyclically (altering 1 to $p$ vertices and
again). Larger values of \code{approximation.iterations} should provide
tighter lower bounds and eventually the exact maximum clique size (the clique
number).
}
\references{
Castelo, R. and Roverato, A. (2006). A robust procedure for Gaussian
graphical model search from microarray data with p larger than n,
\emph{J. Mach. Learn. Res.}, 7:2621-2650

Pardalos, P.M. and Xue, J. (1994). The maximum clique problem,
\emph{J. Global Optim.}, 4:301-328
}

\author{Robert Castelo and Alberto Roverato}
\seealso{\code{\link{qp.search}}, \code{\link{qp.clique}}}
\examples{
data(jmlr06data)
qp.analyse(qp.out.bd5.N20.q10,threshold=0.9,largest.clique=TRUE)
}
\keyword{models}
\keyword{multivariate}
\keyword{graphs}
