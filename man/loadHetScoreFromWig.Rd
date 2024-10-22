% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/heterozygosityScore.R
\name{loadHetScoreFromWig}
\alias{loadHetScoreFromWig}
\title{Load hetScore as data.frame from the binned wig file}
\usage{
loadHetScoreFromWig(wigFile)
}
\arguments{
\item{wigFile}{full path to Heterozygosity Score per bin wig file as created by \code{calculateHetScore}}
}
\value{
data.frame with seqnames, start, end, strand, score columns. Note that the strand column is not set (all values will be '*')
}
\description{
Load hetScore as data.frame from the binned wig file
}
