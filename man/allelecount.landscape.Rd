\name{allelecount.landscape}
\alias{allelecount.landscape}
\title{Calculate allele numbers (frequency in the statistical sense) at each locus in each population}

\description{
  Calculate allele counts}
\usage{
  hetmat <- exp.het.landscape(rland,tbl.out=F)
}
\arguments{
  \item{rland}{the Rmetasim landscape object}
  \item{tbl.out}{return as a (three-dimensional) table if TRUE.  If
  FALSE, return as a dataframe with categorical variables denoting the
  locus, population and allele.}
}
\value{
  Depends on the value of tbl.out.  See above.
}
\seealso{
allelefreq.landscape, obs.het.landscape, exp.het.landscape, Fwright.landscape, Fst.landscape
}
\examples{
  exampleland <- new.example.landscape()
  exampleland <- simulate.landscape(exampleland, 4)
  allelefreq.landscape(exampleland,tbl.out=TRUE)
  allelefreq.landscape(exampleland,tbl.out=FALSE)
  rm(exampleland)
}
\keyword{misc}
