\name{theta.k.landscape}
\alias{theta.k.landscape}
\title{Calculate theta using the number of alleles}

\description{
  Calculate theta using number of alleles from a landscape.}
\usage{
 theta.k.mat <- theta.k.landscape(rland)
}
\arguments{
  \item{rland}{the Rmetasim landscape object}
}
\details{Uses routines in the package 'ape'}
\value{
A matrix with num loci columns and num populations rows.  Each element
reflects the estimated theta for that population x locus combination  }
\seealso{
theta.h.landscape, theta.s.landscape
}
\examples{
  exampleland <- new.example.landscape()
  exampleland <- simulate.landscape(exampleland, 4)

  library(ape)
  theta.k.mat <- theta.k.landscape(exampleland)
  theta.k.mat
}
\keyword{misc}
