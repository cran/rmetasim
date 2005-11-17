\name{sim.landscape}
\alias{sim.landscape}
\title{Run a simulation for a single landscape through time}

\description{
  Simulate a Rmetasim landscape for a number of generations.}
\usage{
  rland <- sim.landscape(rland=l,numit=100,seed=-1,compress=FALSE,adj.lambda=0)
}
\arguments{
  \item{rland}{the Rmetasim landscape object}
  \item{numit}{the number of generations/iterations to simulate, note that
    landscapes will not run past the rland$intparam$totalgens value}
  \item{seed}{The default value of seed uses the seed set in the
    calling environment.  Any other value for seed uses 'set.seed()' to
    reset the random number generator.  sim.landscape uses the RNG
    selected by the
    calling environment.}
  \item{compress}{If true, sim.landscape executes a survival and carrying
    capacity step before returning.  In demographies with high
    reproductive potential, this can significantly reduce the size of R
    objects returned}
  \item{adj.lambda}{Tries to apply a correction to population growth
  that makes the observed growth rate more closely approximate that
  predicted from standard analysis of the Lefkovitch matrices}
}

\examples{
  exampleland <- new.example.landscape()
  exampleland <- sim.landscape(exampleland, 4)
  exampleland
  rm(exampleland)
}
\keyword{misc}
