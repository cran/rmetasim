\name{landscape.democol}
\alias{landscape.democol}
\title{return a vector of population IDs from a landscape}

\description{return a vector of population IDs from a landscape
}
\usage{
  landscape.democol()
}


\value{ a scalar integer representing the largest column of
  demographic information in a landscape's individuals object }

\details{Useful to write functions that will be insensitive to some changes in the individuals object (mainly addition of non-genetic information)  }

\seealso{
  landscape.locus
  }
}
\examples{
	landscape.democol()
}
\keyword{misc}
