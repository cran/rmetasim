\name{exp.het.landscape}
\alias{exp.het.landscape}
\title{Calculate expected heterozygosity}

\description{
  Calculate expected heterozygosity from a landscape}
\usage{
  hetmat <- exp.het.landscape(rland)
}
\arguments{
  \item{rland}{the Rmetasim landscape object}
}
\value{
A matrix with num loci columns and num populations rows.  Each element
reflects the expected heterozygosity for that population x locus combination  }
\details{
  Calculates the expected heterozygosity in each population:
  \deqn{1-\Sigma_{i_k} p_i^2}{1 - sum(p^2)}
  where \eqn{p} is a vector of allele frequencies for a locus in a population.
}
\seealso{
obs.het.landscape, Fst.landscape
}
\examples{
  exampleland <- new.example.landscape()
  exampleland <- simulate.landscape(exampleland, 4)
  exphet <- exp.het.landscape(exampleland)
  rm(exampleland)
}
\keyword{misc}
