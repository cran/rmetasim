\name{theta.s.landscape}
\alias{theta.s.landscape}
\title{Calculate theta using segregating sites}

\description{Calculate theta from a landscape based upon
  the number of seqregating sites.}

\usage{
 theta.s.mat <- theta.s.landscape(rland)
}
\arguments{
  \item{rland}{the Rmetasim landscape object}
}
\details{Uses routines in the package 'ape'}
  
\value{
A matrix with num loci columns and num populations rows.  Each element
reflects the estimated theta for that population x locus combination  }
\seealso{
theta.k.landscape, theta.h.landscape
}
\examples{
  exampleland <- new.example.landscape()
  exampleland <- simulate.landscape(exampleland, 4)

  library(ape)
  theta.s.mat <- theta.s.landscape(exampleland)
  theta.s.mat
}
\keyword{misc}
