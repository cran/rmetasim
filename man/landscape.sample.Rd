\name{landscape.sample}
\alias{landscape.sample}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{simulates sampling for genetics on the landscape}
\description{

Randomly pulls a max of \code{ns} individuals from a max of \code{np}
populations and returns a landscape object that could be used for
further simulation, but is usually used for analyses and summary
statistics calculatiuons

}
\usage{
landscape.sample(rland, np = NULL, ns = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rland}{ landscape object }
  \item{np}{ number populations }
  \item{ns}{ number samples per population }
}
\value{
landscape object
}

\examples{
	l <- landscape.new.example()
	l <- landscape.simulate(l,1)
	l.samp <- landscape.sample(l,np=3,ns=24)
	landscape.amova.pairwise(l.samp)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

\keyword{ misc }% __ONLY ONE__ keyword per line
