\name{theta.h.landscape}
\alias{theta.h.landscape}
\title{Calculate theta using heterozygosity}

\description{Calculate theta from a landscape based upon
  heterozygosity.}

\usage{
 theta.h.mat <- theta.h.landscape(rland)
}
\arguments{
  \item{rland}{the Rmetasim landscape object}
}
\details{Uses routines in the package 'ape'}
  
\value{
A matrix with num loci columns and num populations rows.  Each element
reflects the estimated theta for that population x locus combination  }
\seealso{
theta.k.landscape, theta.s.landscape
}
\examples{

  exampleland <- new.example.landscape()
  exampleland <- sim.landscape(exampleland, 4)

  library(ape)
  theta.h.mat <- theta.h.landscape(exampleland)
  theta.h.mat
}
\keyword{misc}
