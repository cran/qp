\name{jmlr06data}
\docType{data}
\alias{jmlr06data}
\alias{IC.bd5}
\alias{IC.bd20}
\alias{S.bd5.N20}
\alias{S.bd5.N50}
\alias{S.bd5.N150}
\alias{S.bd20.N20}
\alias{S.bd20.N50}
\alias{S.bd20.N150}
\alias{qp.out.bd5.N20.q10}
\alias{qp.out.bd20.N20.q10}
\title{Synthetic data from the article by Castelo and Roverato (2006)}
\usage{data(jmlr06data)}
\description{Synthetic data generated from two graphs with 150 vertices,
$G_1$ and $G_2$. In $G_1$ the boundary of every vertex is at most 5, while
in $G_2$ the boundary of every vertext is at most 20}
\format{
  \tabular{ll}{
    IC.bd5:\tab inverse correlation matrix encoding the independence
                structure of $G_1$\cr
    IC.bd20:\tab inverse correlation matrix encoding the independence
                structure of $G_2$\cr
    S.bd5.N20:\tab sample covariance matrix from a sample of size 20
                   drawn from a normal\cr\tab distribution with mean 0
                   and inverse correlation matrix IC.bd5\cr
    S.bd5.N50:\tab sample covariance matrix from a sample of size 50
                   drawn from a normal\cr\tab distribution with mean 0
                   and inverse correlation matrix IC.bd5\cr
    S.bd5.N150:\tab sample covariance matrix from a sample of size 150
                   drawn from a normal\cr\tab distribution with mean 0
                   and inverse correlation matrix IC.bd5\cr
    S.bd20.N20:\tab sample covariance matrix from a sample of size 20
                   drawn from a normal\cr\tab distribution with mean 0
                   and inverse correlation matrix IC.bd20\cr
    S.bd20.N50:\tab sample covariance matrix from a sample of size 50
                   drawn from a normal\cr\tab distribution with mean 0
                   and inverse correlation matrix IC.bd20\cr
    S.bd20.N150:\tab sample covariance matrix from a sample of size 150
                   drawn from a normal\cr\tab distribution with mean 0
                   and inverse correlation matrix IC.bd20\cr
    qp.out.bd5.N20.q10:\tab output from \code{qp.search} applied
                            to S.bd5.N20 with q=10 and T=500\cr
    qp.out.bd20.N20.q10:\tab output from \code{qp.search} applied
                            to S.bd20.N20 with q=10 and T=500\cr
  }
}
\references{Castelo, R. and Roverato, A. (2006). A robust procedure for
Gaussian graphical model search from microarray data with p larger than n,
\emph{J. Mach. Learn. Res.}, 7:2621-2650}
\keyword{datasets}
