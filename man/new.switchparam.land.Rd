\name{new.switchparam.land}
\alias{new.switchparam.land}
\title{Create a set of boolean parameters}

\description{
  Create a set of boolean (1 or 0) parameters for a Rmetasim landscape.}
\usage{
  ## must be called AFTER new.landscape.empty()
  rland <- new.switchparam.land(rland,re=1,rd=1,mp=1)}
\arguments{
  \item{rland}{skeletion of landscape object, required}
  \item{re}{randepoch (default=0), 1=randomly pick a new epoch (from the
    epochs listed in the landscape) after an epoch completes, 0=epochs
    are chosen in order}
  \item{rd}{randdemo (default=0), 1=randomly choose a demography (from the demographies listed in the landscape) for each subpopulation, 0=demographies are assigned in order}
  \item{mp}{multp (defualt=1), 1=multiple paternity,0=entire families from a single mating}
}

\examples{
  ## Defaults
  exampleland <- new.landscape.empty()
  exampleland <- new.switchparam.land(exampleland)
  exampleland$switchparam

  ## Random epochs, random demographies, and no multiple paternity
  exampleland <- new.landscape.empty()
  exampleland <- new.switchparam.land(exampleland,re=1,rd=1,mp=0)
  exampleland$switchparam

  rm(exampleland)
}
\keyword{misc}
