\name{compress.landscape}
\alias{compress.landscape}
\title{Function to resolve inconsistencies within a landscape, deprecated}

\description{
  Deprecated, 'clean.landscape()' does the same.}
  \usage{
  rland <- compress.landscape(rland)
}
\arguments{
  \item{rland}{the Rmetasim landscape object}
}

\examples{
  exampleland <- new.example.landscape()
  exampleland <- simulate.landscape(exampleland, 4)
  exampleland.clean <- compress.landscape(exampleland)
  rm(exampleland)
}
\keyword{misc}
