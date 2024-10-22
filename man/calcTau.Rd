% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ploidyUtils.R
\name{calcTau}
\alias{calcTau}
\title{find tumor fraction from copy number and normalized read depth}
\usage{
calcTau(cn, NRD)
}
\arguments{
\item{cn}{copy number}

\item{NRD}{normalized read depth}
}
\description{
find tumor fraction from copy number and normalized read depth
}
\examples{
calcTau(NRD=2.9,cn=3)
}
\seealso{
Other copyNumberCalcs: 
\code{\link{calcCopyNumber}()},
\code{\link{calcNrd}()},
\code{\link{calcTumorRatio}()}
}
\concept{copyNumberCalcs}
\keyword{internal}
