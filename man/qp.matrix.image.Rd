\name{qp.matrix.image}
\alias{qp.matrix.image}
\title{
Image of an inverse correlation matrix
}
\description{
Makes an image plot of the absolute value of an inverse correlation
matrix and reports the number of edges of the corresponding
independence graph, the total number of adjacencies of the graph and
the percentage of edges respect to this total number of adjacencies
}
\usage{
qp.matrix.image(M, col=NULL, plot=TRUE)
}
\arguments{
  \item{M}{the matrix to make the image plot}
  \item{col}{flag that when set to \code{NULL} the gray scale is used in the plot}
  \item{plot}{when this flag is set to \code{TRUE} it plots the function}
}
\details{
Returns an object of the class matrix containing the the number of edges of
the corresponding independence graph, the total number of adjacencies of the
graph and the percentage of edges respect to this total number of adjacencies.

When \code{plot=TRUE} it plots the partial correlation coefficients as a
matrix.
}
\references{Castelo, R. and Roverato, A. (2006). A robust procedure for
Gaussian graphical model search from microarray data with p larger than n,
\emph{J. Mach. Learn. Res.}, 7:2621-2650}

\author{Robert Castelo and Alberto Roverato}
\examples{
data(jmlr06data)
qp.matrix.image(IC.bd5)
}
\keyword{models}
\keyword{multivariate}
\keyword{graphs}
