\name{new.local.demo}
\alias{new.local.demo}
\title{Create a Local Demography}

\description{
  Create a local demography for an Rmetasim Landscape object}
\usage{
  ## must be called AFTER integer, switch, and float params have been created
  ## S, R, and M matricies must be square matricies of size rland$intparam$stages by rland$intparam$stages
  rland <- new.local.demo(rland,S=Smatrix,R=Rmatrix,M=Mmatrix)}
\arguments{
  \item{rland}{partially created landscape object, required}
  \item{S}{Survivablity matrix for demograpy, required}
  \item{R}{female Reproduction matrix for demography, required}
  \item{M}{Male reporduction matrix for demography, required}
}

\examples{
  
  exampleS <- matrix(c(0.1, 0, 0.5, 0.3), nrow = 2)
  exampleR <- matrix(c(0, 1.1, 0, 0), nrow = 2)
  exampleM <- matrix(c(0, 0, 0, 1), nrow = 2)
  
  exampleland <- new.landscape.empty()
  exampleland <- new.intparam.land(exampleland, s=2)
  exampleland <- new.floatparam.land(exampleland)
  exampleland <- new.switchparam.land(exampleland)
  exampleland <- new.local.demo(exampleland,exampleS,exampleR,exampleM)

  exampleland$demography$localdem

  rm(exampleS)
  rm(exampleR)
  rm(exampleM)
  rm(exampleland)
}
\keyword{misc}
