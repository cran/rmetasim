\name{Components of Simulations}
\alias{advance.landscape}
\alias{carry.landscape}
\alias{extinct.landscape}
\alias{reproduce.landscape}
\alias{survive.landscape}
\title{Code components to simulate a landscape}

\description{
  These functions can be used to construct custom simulations of
  landscapes.  Each conducts only a single generations worth of change}
\usage{
  rland <- advance.landscape(rland=l)
  rland <- carry.landscape(rland=l)
  rland <- extinct.landscape(rland=l)
  rland <- reproduce.landscape(rland=l)
  rland <- survive.landscape(rland=l)
}
\arguments{
  \item{rland}{the Rmetasim landscape object}
}
\details{
  \code{advance.landscape()} merely advances the generation counter and selects
  the new generations demographic conditions if such conditions can
  vary.  The other functions implement carrying capacity, local
  extinction, reproduction, and survival/growth, respectively.  The
  function \code{sim.landscape()} bundles the functionality of these
  components into a single function.
}

\seealso{
  \code{sim.landscape}
  }
\examples{
  exampleland <- new.example.landscape()
  exampleland <- sim.landscape(exampleland, 4)
  exampleland
  rm(exampleland)
}
\keyword{misc}
