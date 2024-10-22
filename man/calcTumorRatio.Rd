% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ploidyUtils.R
\name{calcTumorRatio}
\alias{calcTumorRatio}
\title{find tumor fraction from the ratio between two copy number peaks and their corresponding read depth}
\usage{
calcTumorRatio(rd1, rd2, cn1, cn2)
}
\arguments{
\item{rd1, rd2}{read depth of two peaks}

\item{cn1, cn2}{copy number  of two peaks}
}
\description{
find tumor fraction from the ratio between two copy number peaks and their corresponding read depth
}
\examples{
calcTumorRatio(rd1=1.2,rd2=2,cn1=2,cn2=4)
calcTumorRatio(rd1=2, rd2=3.5, cn1=1,cn2=2) # 58163 = .155
}
\seealso{
Other copyNumberCalcs: 
\code{\link{calcCopyNumber}()},
\code{\link{calcNrd}()},
\code{\link{calcTau}()}
}
\concept{copyNumberCalcs}
\keyword{internal}
