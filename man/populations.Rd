\name{populations}
\alias{populations}
\title{return a vector of population IDs from a landscape}

\description{return a vector of population IDs from a landscape
}
\usage{
  populations(Rland)
}
\arguments{
  \item{Rland}{the Rmetasim landscape object}
}


\value{
  a vector
}
\details{  Returns a vector of length \code{dim(rland$individuals)[1]} where rland is a
  landscape object.  The vector classifies individuals into populations
  (or habitats)  }

\seealso{
  locus
  }
}
\examples{
  exampleland <- new.example.landscape()
  exampleland <- simulate.landscape(exampleland, 4)
  hist(populations(exampleland),main="Distribution of population size in landscape")
  rm(exampleland)
}
\keyword{misc}
