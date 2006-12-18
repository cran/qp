\name{qp.search}
\alias{qp.search}
\title{
Matrix of non-rejection rates
}
\description{
Calculates the estimates of the non-rejection rates for every pair
of variables.
}
\usage{
qp.search(S, N, q=0, T=500, significance=0.05, binary=TRUE)
}
\arguments{
  \item{S}{sample variance-covariance matrix}
  \item{N}{sample size}
  \item{T}{number of tests per adjacency}
  \item{q}{partial-correlation order}
  \item{significance}{significance level of each test}
  \item{binary}{flag to switch to the compiled C code}
}
\details{
By default \code{binary=TRUE} and the compiled and faster C code
corresponding to this function will be executed. If \code{binary=FALSE}
is set, then the R code will be executed.
}
\value{
  \item{A}{matrix with the acceptance test counts}
  \item{T}{number of tests per adjacency (copied from the input parameter)}
}
\references{Castelo, R. and Roverato, A. (2006). A robust procedure for
Gaussian graphical model search from microarray data with p larger than n,
\emph{J. Mach. Learn. Res.}, 7:2621-2650}
\author{Robert Castelo and Alberto Roverato}
\seealso{\code{\link{qp.edge.prob}}, \code{\link{qp.analyse}}, \code{\link{qp.hist}}}
\examples{
data(jmlr06data)
S <- S.bd5.N20
N <- 20
q <- 6
T <- 100
qp.out <- qp.search(S,N,q,T)
}
\keyword{models}
\keyword{multivariate}
\keyword{graphs}
