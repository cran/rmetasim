\name{new.floatparam.land}
\alias{new.floatparam.land}
\title{Create a set of floating point parameters}

\description{
  Create a set of floating point parameters for a Rmetasim landscape.}
\usage{
  ## must be called AFTER new.landscape.empty()
  rland <- new.floatparam.land(rland,s=0)}
\arguments{
  \item{rland}{skeletion of landscape object, required}
  \item{s}{selfing (default=0), the selfing rate of the species}
}
\examples{
  ## Defaults
  exampleland <- new.landscape.empty()
  exampleland <- new.floatparam.land(exampleland)
  exampleland$floatparam

  ## .5 selfing rate
  exampleland <- new.landscape.empty()
  exampleland <- new.floatparam.land(exampleland,s=0.5)
  exampleland$floatparam

  rm(exampleland)
}
\keyword{misc}
