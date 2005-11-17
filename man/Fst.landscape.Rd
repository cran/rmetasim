\name{Fst.landscape}
\alias{Fst.landscape}
\title{Calculates population structure statistic for the entire landscape}

\description{
  Calculate Fst for each allele at each locus in the landscape.  If
  verb is set to TRUE, the function prints average Fst for loci and
  overall.
}
\usage{
  Fstmat <- Fst.landscape(rland,verb=F)
}
\arguments{
  \item{rland}{the Rmetasim landscape object}
  \item{verb}{determines whether there is verbose output}
}
\value{
A matrix with num alleles columns and num loci rows.  Each element
reflects the value of Fst for that allelexlocus combination.  NA is
assigned to alleles that are not present at a locus (either no longer
or ever)}
\details{
  Calculates Fst based upon the ratio of variance in allele frequency
  across subpopulations to the total variance in that allele's
  frequency.  Does not calculate Wright's other statistics.
}
\seealso{
obs.het.landscape, exp.het.landscape, FWright.landscape
}
\examples{
  exampleland <- new.example.landscape()
  exampleland <- sim.landscape(exampleland, 4)
  Fst <- Fst.landscape(exampleland,verb=TRUE)
  Fst
  rm(exampleland,Fst)
}

\keyword{misc}
