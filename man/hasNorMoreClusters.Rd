% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ploidyUtils.R
\name{hasNorMoreClusters}
\alias{hasNorMoreClusters}
\title{determine if the peak has a minimum number of clusters}
\usage{
hasNorMoreClusters(
  hetScoreDensityResult,
  N,
  heterozygosityScoreThreshold,
  minObservations = 20
)
}
\arguments{
\item{hetScoreDensityResult}{hetScore values returned from density function}

\item{N}{minimum number of clusters}

\item{heterozygosityScoreThreshold}{peaks with a hetScore mode above this value are considered heterozygous, typically 0.98,
but may vary depending on NGS library quality and preparation}

\item{minObservations}{required number of input values to determine number of clusters present}
}
\description{
determine if the peak has a minimum number of clusters
}
\keyword{internal}
