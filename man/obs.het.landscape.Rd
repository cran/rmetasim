\name{obs.het.landscape}
\alias{obs.het.landscape}
\title{Calculate observed heterozygosity}

\description{
  Calculate observed heterozygosity from a landscape}
\usage{
  hetmat <- obs.het.landscape(rland)
}
\arguments{
  \item{rland}{the Rmetasim landscape object}
}
\value{
A matrix with num loci columns and num populations rows.  Each element
reflects the observed heterozygosity for that population x locus combination  }

\seealso{
exp.het.landscape, Fst.landscape
}

\examples{
  exampleland <- new.example.landscape()
  exampleland <- sim.landscape(exampleland, 4)
  obshet <- obs.het.landscape(exampleland)
  rm(exampleland)
}
\keyword{misc}
