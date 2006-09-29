\name{qp.edge.prob}
\alias{qp.edge.prob}
\title{
Estimate of the non-rejection rate
}
\description{
Calculates the estimate of the non-rejection rate for a pair of variables,
that is, the proportion of conditional independence tests
that accept the null hypothesis of zero partial correlation given the
q-order conditionals.
}
\usage{
qp.edge.prob(S, N, i=1, j=2, q=0, T=500, significance=0.05, binary=TRUE)
}
\arguments{
  \item{S}{sample variance-covariance matrix}
  \item{N}{sample size}
  \item{i}{index of one variable}
  \item{j}{index of another variable}
  \item{q}{partial-correlation order}
  \item{T}{number of tests per adjacency}
  \item{significance}{significance level of each test}
  \item{binary}{flag to switch to the compiled C code}
}
\details{
By default \code{binary=TRUE} and the compiled and faster C code
corresponding to this function will be executed. If \code{binary=FALSE}
is set, then the R code will be executed.
}
\references{Castelo, R. and Roverato, A. (2006). A robust procedure for
Gaussian graphical model search from microarray data with p larger than
n, \emph{J. Mach. Learn. Res, accepted}}
\author{Robert Castelo and Alberto Roverato}
\seealso{\code{\link{qp.search}},\code{\link{qp.ci.test}}}
\examples{
data(jmlr06data)
S <- S.bd5.N20
N <- 20
q <- 6
T <- 100

qp.edge.prob(S,N,i=3,j=4,q,T)
}
\keyword{models}
\keyword{multivariate}
\keyword{graphs}
