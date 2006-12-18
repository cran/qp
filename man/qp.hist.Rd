\name{qp.hist}
\alias{qp.hist}
\title{
Histogram of the non-rejection rates
}
\description{
Using the output of \code{\link{qp.search}} this function plots
the histogram of the estimated non-rejection rates. When the inverse
correlation matrix from the generative graph is given, it provides
additional plots of information.
}
\usage{
qp.hist(qp.output, IC=NULL, prob=FALSE)

}
\arguments{
  \item{qp.output}{output of qp.search}
  \item{IC}{inverse correlation matrix from the generative graph}
  \item{prob}{when this flag is set to TRUE the histograms show densities,
              otherwise they show absolute frequencies.}
}
\references{Castelo, R. and Roverato, A. (2006). A robust procedure for
Gaussian graphical model search from microarray data with p larger than n,
\emph{J. Mach. Learn. Res.}, 7:2621-2650}

\author{Robert Castelo and Alberto Roverato}
\seealso{\code{\link{qp.search}}}
\examples{
data(jmlr06data)

# if are working with synthetic data and have the IC matrix

qp.hist(qp.out.bd5.N20.q10,IC.bd5,prob=TRUE)


# otherwise, we just look at the non-rejection rate distribution

qp.hist(qp.out.bd5.N20.q10,NULL,prob=TRUE)

}
\keyword{models}
\keyword{multivariate}
\keyword{graphs}
