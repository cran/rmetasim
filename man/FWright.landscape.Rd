\name{FWright.landscape}
\alias{FWright.landscape}
\title{Calculates Wright's Fixation index at each locus in each population}

\description{
  Calculate F
  }
\usage{
  Fmat <- FWright.landscape(rland)
}
\arguments{
  \item{rland}{the Rmetasim landscape object}
}
\value{
  A matrix of population x locus measures of Wright's F
}
\seealso{
obs.het.landscape, exp.het.landscape, allelefreq.landscape, Fst.landscape
}
\examples{
  exampleland <- new.example.landscape()
  FWright.landscape(exampleland)
  exampleland <- simulate.landscape(exampleland, 4)
  FWright.landscape(exampleland)
  rm(exampleland)
}
\keyword{misc}
