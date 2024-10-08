% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bmdUtils.R
\name{linearCoordinates}
\alias{linearCoordinates}
\title{Coordinate system for multiple chromosomes placed linearly one after the other}
\usage{
linearCoordinates(chroms, chromStart, chromEnd)

linearCoordinates(chroms, chromStart, chromEnd)
}
\arguments{
\item{chroms}{A vector of reference sequence numbers for the chromosomes being used}

\item{chromStart}{A vector of all chromosome start positions in linear coordinates.}

\item{chromEnd}{A vector of all chromosome positions in linear coordinates}
}
\value{
\link{linearCoordinates} class
}
\description{
The coordinate system takes a list of chromosomes and allocates vectors of their start and end positions
in the linear coordinate space. All the coordinate vectors are designed to
be indexable directly with chromosome number, so you can do e.g. \code{coords@chromStart[7]}
to get 7-th chromosome start.
}
\details{
The bin coordinates are calculated using \link{binnedPosStart} and \link{binnedPosEnd}.
See those functions for more information.

Because the bins cannot cleanly cut chromosomes whose lengths are not divisible by 10000, there are following issues: the start
bin can contain some values from the end of previous chromosome. The end bin will be pure, but values
at the end of the chromosomes are chopped off and merged with the next chromosome. Up to 9999 values at the very end
of the array are chopped off by the 10K array.

This object can be returned directly from the internal rgdObject by \link{getLinearCoordinates}.

For each start/end pair it is guaranteed that start<=end
}
\section{Slots}{

\describe{
\item{\code{chroms}}{A vector of reference sequence numbers for the chromosomes being used}

\item{\code{maxcn}}{Maximum chromosome number. Often used when allocating arrays.}

\item{\code{totalLength}}{The last coordinate utilized by this system (max chromEnd)}

\item{\code{totalLength1K}}{The last coordinate utilized by this system (max chromEnd) for 1K binned arrays}

\item{\code{totalLength10K}}{The last coordinate utilized by this system (max chromEnd) for 10K binned arrays}

\item{\code{chromStart}}{A vector of chromosome start positions in linear coordinates}

\item{\code{chromEnd}}{A vector of chromosome end positions in linear coordinates (inclusive)}

\item{\code{chromStart1K}}{A vector of chromosome start positions for 1K binned arrays.
The 1K bin that contains \code{chromStart}.}

\item{\code{chromEnd1K}}{A vector of chromosome end positions for 1K binned arrays. The 1K bin
just before the 1K bin containing \code{chromEnd}. This ensures that 1K array segments
do not overlap.}

\item{\code{chromStart10K}}{A vector of chromosome start positions for 10K binned arrays.
The 10K bin that contains \code{chromStart}.}

\item{\code{chromEnd10K}}{A vector of chromosome end positions for 10K binned arrays.
The 10K bin just before the 10K bin containing \code{chromEnd}. This ensures that 10K array segments
do not overlap.}
}}

\note{
We do not support non-consecutive chromosome indices and fail if an attempt is made to create
such coordinate system.
}
\seealso{
Other coordinates: 
\code{\link{bimaToLinear}()},
\code{\link{binnedPosEnd}()},
\code{\link{binnedPosStart}()},
\code{\link{getLinearCoordinates}()},
\code{\link{linearToBima}()}
}
\concept{coordinates}
