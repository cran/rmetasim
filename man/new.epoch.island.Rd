\name{new.epoch.island}
\alias{new.epoch.island}
\title{Create an Island Model Epoch}

\description{
  Create an epoch for a Rmetasim landscape object that follows the
  conventions of the island model}
\usage{
  ## must be called AFTER integer, switch, and float params have
  ## been created and after the demography has been created

  rland <- new.epoch.island(rland,s,sfrom,sto,m,mfrom,mto,f,ffrom,fto,
                            epochprob=1,startgen=0,extinct=NULL,carry=NULL,localprob=NULL)
}
\arguments{
  \item{rland}{partially created landscape object, required}
  \item{s}{Rate of migration by individuals}
  \item{sfrom}{vector of lifecycle stages in which the individuals
    migrate from their current subpopulation, must be
    rland$intparam$stages in length}
  \item{sto}{vector of lifecycle stages in which the individuals arrive
    in a new subpopulation, must be rland$intparam$stages in length}
  \item{m}{Rate of migration by newborns}
  \item{mfrom}{vector of lifecycle stages in which the newborns leave
    their subpopulation, must be rland$intparam$stages in length}
  \item{mto}{vector of lifecycle stages in which the newborns arrive
    at their new subpopulation, must be rland$intparam$stages in length}
  \item{f}{Rate of sperm/pollen movement}
  \item{ffrom}{vector of lifecycle stages in which the males produce
    sperm/pollen, must be rland$intparam$stages in length}
  \item{fto}{???}
  \item{epochprob}{(default=1) probability of choosing this epoch randomly if randepoch==1}
  \item{startgen}{(default=0) generation in which this epoch starts}
  \item{extinct}{(default=NULL) vector of extinction probabilities per
    generation for each subpopulation, must be rland$intparam$habitats
    in length, passing NULL gives a 0\% probability of extinction to each subpopulation}
    \item{carry}{(default=NULL) vector of carrying capacities for each subpopulation, must be rland$intparam$habitats in length, passing NULL gives a 1000 individual carrying capacity to each subpopulation}
    \item{localprob}{(default=NULL) vector of probabilites for choosing local demographies, must be length(rland$demography$localdem) in length, passing NULL gives each demography an equal probability}
}

  \examples{
  
  exampleS <- matrix(c(0.1, 0, 0.5, 0.3), nrow = 2)
  exampleR <- matrix(c(0, 1.1, 0, 0), nrow = 2)
  exampleM <- matrix(c(0, 0, 0, 1), nrow = 2)
  
  exampleland <- new.landscape.empty()
  exampleland <- new.intparam.land(exampleland, s=2, h=2)
  exampleland <- new.floatparam.land(exampleland)
  exampleland <- new.switchparam.land(exampleland)
  exampleland <- new.local.demo(exampleland,exampleS,exampleR,exampleM)

  exampleland <- new.epoch.island(exampleland,.1,c(0,1),c(0,1),.3,c(0,1),c(1,0),.5,c(1,0),c(0,1))

  exampleland$demography$epochs[[1]]

  rm(exampleland)
}
\keyword{misc}
