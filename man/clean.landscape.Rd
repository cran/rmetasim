\name{clean.landscape}
\alias{clean.landscape}
\title{Function to resolve inconsistencies within a landscape}

\description{Converts a landscape to internal format and back.  This can
resolve inconsistencies in a 'hand-built' landscape}
  
  \usage{
  rland <- clean.landscape(rland)
}
\arguments{
  \item{rland}{the Rmetasim landscape object}
}

\examples{
  exampleland <- new.example.landscape()
  exampleland <- simulate.landscape(exampleland, 4)
  exampleland.clean <- clean.landscape(exampleland)
  rm(exampleland)
}
\keyword{misc}
