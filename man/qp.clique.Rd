\name{qp.clique}
\alias{qp.clique}
\title{
Relationship between non-rejection rate and maximum clique size
}
\description{
Using the output of \code{\link{qp.search}} this function calculates
the maximum clique size as a function of the minimum threshold on the
non-rejection rate for removing an edge
}
\usage{
qp.clique(qp.output, N, threshold.lim=c(0,1), breaks=5, plot.image=TRUE)

}
\arguments{
  \item{qp.output}{output of \code{\link{qp.search}}}
  \item{N}{sample size}
  \item{threshold.lim}{range of the non-rejection rate threshold on which
                       calculate the funcion}
  \item{breaks}{one of:
    \itemize{
      \item a vector giving the breakpoints along the range defined by
            \code{threshold.lim},
      \item a single number giving the number of equidistant breakpoints
            that divide the range defined by \code{threshold.lim}.
    }
  }
  \item{plot.image}{when this flag is set to \code{TRUE}, the
                    \code{qp.clique} plot is produced}
}
\details{
The \code{qp.clique} plot provides information on the graphs
potentially selected by specifying different values of the threshold.
Every circle in the plot corresponds to a graph and has three values
associated with it: the threshold value used to construct the graph
(horizontal axis); the number of vertices of the largest clique of the
graph (vertical axis); the percentage of present edges in the graph
(number inside the plot, beside the circle). Furthermore, adjacent
circles are joined by a line and the dotted horizontal line corresponds
to the sample size N.

Beware that setting breakpoints between 0.95 and 1.0 may result in a very
long time of computation due to the NP-hardness of the problem of
calculating the size of the largest clique.
}
\value{
  \item{threshold}{threshold on the non-rejection rate that provides the
                   maximum clique size that is strictly smaller than the
                   sample size \code{N}}
  \item{size}{maximum clique size strictly smaller than the sample size
              \code{N}}
}
\references{Castelo, R. and Roverato, A. (2006). A robust procedure for Gaussian
graphical model search from microarray data with p larger than n, \emph{J. Mach.
Learn. Res., accepted}}

\author{Robert Castelo and Alberto Roverato}
\seealso{\code{\link{qp.search}}}
\examples{
data(jmlr06data)
qp.clique(qp.out.bd5.N20.q10,20)
}
\keyword{models}
\keyword{multivariate}
\keyword{graphs}
