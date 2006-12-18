\name{qp.get.cliques}
\alias{qp.get.cliques}
\title{
Cliques of an undirected graph
}
\description{
It finds the set of cliques, i.e. maximal complete subsets of vertices,
of an undirected graph taken as an incidence matrix.
}
\usage{
qp.get.cliques(I, binary=TRUE)

}
\arguments{
  \item{I}{incidence matrix}
  \item{binary}{flag to switch to the compiled C code}
}
\details{
It uses the algorithm described in Bron and Kerbosch (1973) and returns
a list where each member is a vector of vertices forming a clique in the
given graph. Beware that the problem of finding the set of cliques is
NP-complete and the time of computation of this algorithm grows
exponentially in the graph density (number of actual edges over the total
number of adjacencies).
}
\references{Castelo, R. and Roverato, A. (2006). A robust procedure for
Gaussian graphical model search from microarray data with p larger than n,
\emph{J. Mach. Learn. Res.}, 7:2621-2650

Bron, C. and Kerbosch, J (1973). Finding all cliques of an undirected graph,
\emph{Commun. ACM}, 16:575--577
}

\author{Robert Castelo and Alberto Roverato}
\seealso{\code{\link{qp.graph}}, \code{\link{qp.clique}}}
\examples{
data(jmlr06data)
I <- qp.graph(qp.out.bd5.N20.q10,threshold=0.9)
cliquelist <- qp.get.cliques(I)
sprintf("the graph has \%d cliques\n",length(cliquelist))
}
\keyword{models}
\keyword{multivariate}
\keyword{graphs}
