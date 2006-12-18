\name{qp.ci.test}
\alias{qp.ci.test}
\title{
Conditional independence test
}
\description{
Performs a test for conditional independence between variables
indexed by i and j given the conditioning set Q
}
\usage{
qp.ci.test(S, N, i=1, j=2, Q=c(), binary=TRUE)
}
\arguments{
  \item{S}{sample variance-covariance matrix}
  \item{N}{sample size}
  \item{i}{index of one variable}
  \item{j}{index of another variable}
  \item{Q}{conditioning set}
  \item{binary}{flag to switch to the compiled C code}
}
\details{
By default \code{binary=TRUE} and the compiled and faster C code
corresponding to this function will be executed. If \code{binary=FALSE}
is set, then the R code will be executed.
}
\value{
  \item{t.value}{the t-statistic value}
  \item{p.value}{the p-value on rejecting the null hypothesis of
                 conditional independence}
}
\examples{
data(jmlr06data)
S <- S.bd5.N20
N <- 20
qp.ci.test(S,N,i=3,j=4,Q=c(5,6,7))
}
\references{Castelo, R. and Roverato, A. (2006). A robust procedure for
Gaussian graphical model search from microarray data with p larger than
n, \emph{J. Mach. Learn. Res.}, 7:2621-2650}
\author{Robert Castelo and Alberto Roverato}
\seealso{\code{\link{qp.edge.prob}}}
\keyword{models}
\keyword{multivariate}
\keyword{graphs}
