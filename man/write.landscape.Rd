\name{write.landscape}
\alias{write.landscape}
\title{Save a landscape to a file}

\description{
  Save a Rmetasim landscape object to a file}
\usage{
  rland <- write.landscape(rland=l,fn="mylandscape.lnd")
}
\arguments{
  \item{rland}{the Rmetasim landscape object}
  \item{fn}{the path and name of the file to save the landscape to}
}

\examples{
  ## Needs write access to the current directory, files created!!
  exampleland <- new.example.landscape()
  write.landscape(exampleland, "exampleland.lnd")

  rm(exampleland)
}
\keyword{misc}
