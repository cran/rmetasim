\name{simulate.landscape}
\alias{simulate.landscape}
\title{Simulate a landscape}

\description{
  Simulate a Rmetasim landscape for a number of generations.}
\usage{
  rland <- simulate.landscape(rland=l,numit=100,seed=-1,compress=FALSE)
}
\arguments{
  \item{rland}{the Rmetasim landscape object}
  \item{numit}{the number of generations/iterations to simulate, note that
    landscapes will not run past the rland$intparam$totalgens value}
  \item{seed}{The default value of seed uses the seed set in the
    calling environment.  Any other value for seed uses 'set.seed()' to
    reset the random number generator.  simulate.landscape uses the RNG
    selected by the
    calling environment.}
  \item{compress}{If true, simulate.landscape executes a survival and carrying
    capacity step before returning.  In demographies with high
    reproductive potential, this can significantly reduce the size of R
    objects returned}
}

\examples{
  exampleland <- new.example.landscape()
  exampleland <- simulate.landscape(exampleland, 4)
  exampleland
  rm(exampleland)
}
\keyword{misc}
