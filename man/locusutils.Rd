\name{Genetic locus utilities}
\alias{locusvec}
\alias{ploidy}
\title{Provides convenient summarys of locus characteristics}

\description{return a vector of information about loci
}
\usage{
  locusvec(Rland)
  ploidy(Rland)
}
\arguments{
  \item{Rland}{the Rmetasim landscape object}
}


\value{
  vector
}

\details{\code{locusvec()} returns a vector that gives the numeric locus
ids for the columns in the \code{rland$individuals} matrix.
\code{ploidy()} returns a vector that gives the ploidy of each locus.}
}
\examples{
  exampleland <- new.example.landscape()
  exampleland <- sim.landscape(exampleland, 4)
  locusvec(exampleland)
  ploidy(exampleland)
}
\keyword{misc}
