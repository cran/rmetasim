\name{read.landscape}
\alias{read.landscape}
\title{Load a landscape from a file}

\description{
  Load a Rmetasim landscape from a file into an R object.}
\usage{
  rland <- read.landscape(fn="mylandscape.lnd")}
\arguments{
  \item{fn}{the path and name of the file containing the landscape}
}

\examples{
  ## Needs write access to the current directory, files created!!
  write.landscape(new.example.landscape(), "exampleland.lnd")
  exampleland <- read.landscape("exampleland.lnd")
  exampleland
  rm(exampleland)
}
\keyword{misc}
