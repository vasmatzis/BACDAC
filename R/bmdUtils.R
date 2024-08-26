#' Get object describing linear coordinate space for main chromosomes
#'
#' @param chromosomes List of main chromosomes (technically .sva files) to form the coordinate system from
#' @return \link{linearCoordinates}
#'
#' @family coordinates
#'
#' @export
getLinearCoordinates <- function(chromosomes=1:24) {
  uniqueChroms <- sort(unique(chromosomes))

  # Test that we have consecutive chromosome numbers. Our algorithms depend on this now.
  if(max(diff(uniqueChroms))!=1) {
    stop("The sequence of chromosomes contains holes. This is not supported yet.", uniqueChroms)
  }

  chstrt   <- rgdObject$linearGenomeCoords[, 'start']
  chend    <- rgdObject$linearGenomeCoords[, 'end']

  if(sum(chstrt>chend,na.rm = TRUE)!=0) {
    stop("There is a chromosome that starts after it ends")
  }

  return(linearCoordinates(chroms=uniqueChroms, chromStart=chstrt, chromEnd=chend))
}


#' Coordinate system for multiple chromosomes placed linearly one after the other
#'
#' The coordinate system takes a list of chromosomes and allocates vectors of their start and end positions
#' in the linear coordinate space. All the coordinate vectors are designed to
#' be indexable directly with chromosome number, so you can do e.g. \code{coords@@chromStart[7]}
#' to get 7-th chromosome start.
#'
#' The bin coordinates are calculated using \link{binnedPosStart} and \link{binnedPosEnd}.
#' See those functions for more information.
#'
#' Because the bins cannot cleanly cut chromosomes whose lengths are not divisible by 10000, there are following issues: the start
#' bin can contain some values from the end of previous chromosome. The end bin will be pure, but values
#' at the end of the chromosomes are chopped off and merged with the next chromosome. Up to 9999 values at the very end
#' of the array are chopped off by the 10K array.
#'
#' This object can be returned directly from the internal rgdObject by \link{getLinearCoordinates}.
#'
#' For each start/end pair it is guaranteed that start<=end
#'
#' @note We do not support non-consecutive chromosome indices and fail if an attempt is made to create
#'  such coordinate system.
#'
#' @slot chroms A vector of reference sequence numbers for the chromosomes being used
#' @slot maxcn Maximum chromosome number. Often used when allocating arrays.
#'
#' @slot totalLength The last coordinate utilized by this system (max chromEnd)
#' @slot totalLength1K The last coordinate utilized by this system (max chromEnd) for 1K binned arrays
#' @slot totalLength10K The last coordinate utilized by this system (max chromEnd) for 10K binned arrays
#'
#' @slot chromStart A vector of chromosome start positions in linear coordinates
#' @slot chromEnd A vector of chromosome end positions in linear coordinates (inclusive)
#'
#' @slot chromStart1K A vector of chromosome start positions for 1K binned arrays.
#' The 1K bin that contains \code{chromStart}.
#' @slot chromEnd1K A vector of chromosome end positions for 1K binned arrays. The 1K bin
#' just before the 1K bin containing \code{chromEnd}. This ensures that 1K array segments
#' do not overlap.
#'
#' @slot chromStart10K A vector of chromosome start positions for 10K binned arrays.
#' The 10K bin that contains \code{chromStart}.
#' @slot chromEnd10K A vector of chromosome end positions for 10K binned arrays.
#' The 10K bin just before the 10K bin containing \code{chromEnd}. This ensures that 10K array segments
#' do not overlap.
#'
#' @export
linearCoordinates <- setClass("linearCoordinates",
                              slots=c(
                                "chroms"="numeric",
                                "maxcn"="numeric",
                                "totalLength"="numeric",
                                "totalLength1K"="numeric",
                                "totalLength10K"="numeric",
                                "chromStart"="numeric",
                                "chromStart1K"="numeric",
                                "chromStart10K"="numeric",
                                "chromEnd"="numeric",
                                "chromEnd1K"="numeric",
                                "chromEnd10K"="numeric"
                              )
)

#' Linear coordinates class constructor
#'
#' @param chroms A vector of reference sequence numbers for the chromosomes being used
#' @param chromStart A vector of all chromosome start positions in linear coordinates.
#' @param chromEnd A vector of all chromosome positions in linear coordinates
#'
#' @return \link{linearCoordinates} class
#'
#' @family coordinates
#'
#' @export
linearCoordinates <- function(chroms, chromStart, chromEnd) {
  lengths <- chromEnd-chromStart+1
  minIndex <- chroms[which.min(lengths[chroms])]
  minLength <- lengths[minIndex]
  if(minLength<10000) {
    stop(sprintf("We cannot represent chromosome %d which is shorter than 10000 base pairs (%d).", chroms[minIndex], minLength))
  }

  chromStart1K   <- binnedPosStart(chromStart, 1000)
  chromStart10K  <- binnedPosStart(chromStart, 10000)
  chromEnd1K     <- binnedPosEnd(chromEnd, 1000)
  chromEnd10K    <- binnedPosEnd(chromEnd, 10000)
  maxcn          <- max(chroms)
  totalLength    <- unname(chromEnd[which.max(chroms)])
  totalLength1K  <- max(chromEnd1K,na.rm=TRUE)
  totalLength10K <- max(chromEnd10K,na.rm=TRUE)
  return(new("linearCoordinates",
             chroms=chroms,
             maxcn=maxcn,
             totalLength=totalLength,
             totalLength1K=totalLength1K,
             totalLength10K=totalLength10K,
             chromStart=chromStart,
             chromStart1K=chromStart1K,
             chromStart10K=chromStart10K,
             chromEnd=chromEnd,
             chromEnd1K=chromEnd1K,
             chromEnd10K=chromEnd10K
  ))
}

#' For given position and bin size, return the binned position
#'
#' This is done for consistency. There are many off by one errors
#' one can make when doing this.
#'
#' This function assumes the following binning for
#' 1K and 10K bins and will return the \strong{bin starts}:
#'
#' \tabular{ccc}{
#' bin# \tab 1K \tab 10K \cr
#' 1 \tab    \strong{1}..1000 \tab     \strong{1}..10000 \cr
#' 2 \tab \strong{1001}..2000 \tab \strong{10001}..20000 \cr
#' 3 \tab \strong{2001}..3000 \tab \strong{20001}..30000 \cr
#' \dots \tab \tab \cr
#' }
#'
#' @param pos Vector of positions
#' @param binSize Size of the bin
#'
#' @return The bin which contains the starting position.
#'
#' @family coordinates
#'
#' @export
binnedPosStart <- function(pos, binSize=1000) {
  return(trunc((pos-1)/binSize)+1)
}


#' Translate from BIMA coordinates to linear genome.
#'
#' Translate a list of BIMA chromosomes (svaNumbers) and BIMA positions to the
#' linear genome coordinates. The linear genome starts at position 1 and
#' consists of all sequences in linear order.
#'
#' @param rgd Reference Genome Descriptor
#' @param svaNumber BIMA chromosome number (as present in the SVA file)
#' @param svaPos Position within the sequence (1-based)
#' @param binSize How to round the coordinates. The result will be divided by this number.
#' Note - recreates a bug in George's code at the moment,
#' so the result can be up to binSize away from what you'd expect.
#'
#' @return Linear genome position (1-based). For negative \code{svaPos}
#' return negative number.
#'
#' @note Fix the bug with binSize rounding in two passes
#'
#' @note \code{bimaToLinear(rgd, svaNumber, abs(svaPos), binSize=10000)} replaces
#' the original \code{conver10Kserial(svaPos, svaNumber)}
#'
#' @family coordinates
#'
#' @export
bimaToLinear <- function(rgd, svaNumber, svaPos, binSize=1)
{
  # checkSvaNumber(svaNumber)

  # store signs of the positions and take absolute values
  posSign <- as.numeric(svaPos>=0)*2-1 # Convert >=0 to 1, <0 to -1
  svaPos <- pmax(0, abs(svaPos)-1) # Clamp because inconsistency of BIMA

  # Convert BIMA coordinates to global ones
  globalPos <- trunc((rgd$globalCoordinates$svaOffset[svaNumber+1]-1)/binSize)+trunc(svaPos/binSize)+1

  return(globalPos * posSign)
}


#' Translate linear genome coordinates back to reference genome coordinates
#'
#' This is an inverse function to \link{bimaToLinear}.
#'
#' @param rgd Reference Genome Descriptor
#' @param globalPos linear genome position (1-based)
#' @param binSize Size of bin that was originally used
#'
#' #' @return A list with two elements:
#' \tabular{ll}{
#' \code{svaNumber} \tab reference sequence number \cr
#' \code{svaPos} \tab position within the reference sequence (1-based) \cr
#' }
#'
#' @family coordinates
#'
#' @export
linearToBima <- function(rgd, globalPos, binSize=1)
{
  # store signs of the positions and take absolute values
  posSign <- as.numeric(globalPos>=0)*2-1 # Convert >=0 to 1, <0 to -1

  gp <- abs(globalPos)-1

  # The offset part of the calculation
  svaOffsets <- trunc((rgd$globalCoordinates$svaOffset-1)/binSize)

  # Deal with NAs when broken sequence of svas
  names(svaOffsets) <- seqFwd(from=0, to=length(svaOffsets)-1)
  svaOffsets <- svaOffsets[!is.na(svaOffsets)]

  # Find where our positions fall
  intervals <- findInterval(gp, svaOffsets)

  # Determine relative positions within those reference sequences
  svaPos <- unname((gp - svaOffsets[intervals])*binSize) + 1

  svaNumber <- as.integer(names(svaOffsets)[intervals])

  # Re-apply signs
  svaPos <- svaPos * posSign

  result <- list(
    svaNumber = svaNumber,
    svaPos = svaPos
  )
  return(result)
}

#' For given position and bin size, return the binned position
#'
#' This is done for consistency. There are many off by one errors
#' one can make when doing this.
#'
#' This function assumes the following bining for
#' 1K and 10K bins and will return the \strong{bin ends}:
#'
#' \tabular{ccc}{
#' bin# \tab 1K \tab 10K \cr
#' 1 \tab    1..\strong{1000} \tab     1..\strong{10000} \cr
#' 2 \tab 1001..\strong{2000} \tab 10001..\strong{20000} \cr
#' 3 \tab 2001..\strong{3000} \tab 20001..\strong{30000} \cr
#' \dots \tab \tab \cr
#' }
#'
#' @param pos Vector of positions
#' @param binSize Size of the bin
#'
#' @return The bin just before the bin containing the ending position (last fully populated bin)
#'
#' @family coordinates
#'
#' @export
binnedPosEnd <- function(pos, binSize=1000) {
  return(trunc(pos/binSize))
}


#' Forward only \code{:} replacement
#'
#' @param from first sequence item
#' @param to last sequence item
#'
#' @return Just like \code{from:to} except if \code{to < from}, the returned
#' sequence is empty.
#'
#' @export
seqFwd <- function(from=1, to=1) {
  if(from<=to) {
    return(from:to)
  } else {
    return(vector(mode = "integer", length = 0))
  }
}

#' Convert chromosome integer to a character string
#'
#' allosome 23 will be reported as X, 24 as Y, 25 as M.
#' Option to include "chr" prefix, if FALSE and if the 'chr' prefix is present, 'chr' will be removed.
#'
#' @param ichr number vector, for example:  2 or c(2,5,23)
#' @param withChrPrefix if TRUE, will add "chr" to the front of every ichr value
#'
#' @return character vector of ichr
#'
convertChromToCharacter <- function(ichr,withChrPrefix=FALSE){
  result <- vector(mode = "numeric", length = length(ichr))

  # remove the 'chr' first to simplify the tests, and to make sure all are consistant
  prefixPresent <- (substr(ichr,1,3) =='chr')
  ichr[prefixPresent] <- substr(ichr[prefixPresent],4,nchar(ichr[prefixPresent]))
  autosomes <- 1:22
  allosomes <- 23:24
  mitochondria <- 25


  nkey <- c(which(ichr %in% autosomes))
  xkey <- c(which(ichr==allosomes[1]), which(ichr=='X'))
  ykey <- c(which(ichr==allosomes[2]), which(ichr=='Y'))
  mkey <- c(which(ichr==mitochondria), which(ichr=='M'))
  #vkey=which(ichr==virus)

  # do as a separate step from above, to prevent premature conversion to character
  result[nkey] <- as.character(ichr[nkey])
  result[xkey] <- 'X'
  result[ykey] <- 'Y'
  result[mkey] <- 'M'

  if(withChrPrefix){
    result <- paste0('chr',result)
  }
  return(result)
}

#' Get centromere positions, will return NULL if no centromeres exist (i.e. in mouse, rat)
#'
#' @return Centromere position array. Two-dimensions - first is chromosome number, the other separates start and end positions.
#'
#' @family ideogram
#'
#'
#' @export
getCentromerePositions <- function() {
  cytobands <- ideogram[['cytoBands']]

  xind <- 23
  yind <- 24
  maxChrInd <- 24 # Index of the last chromosome, not last sva

  centro <- cytobands[which(cytobands$gieStain == "acen"),]
  if(nrow(centro)>0){
    centroArray <- array(data=0, dim=c(maxChrInd,2), dimnames=list(chromosome=seq_len(maxChrInd), c("start", "end")))
    #Loop finds and saves the start and end of the centromeres
    for(ict in seq(1,nrow(centro),by=2)){
      chrN  <-  substring(centro$chrom[ict],4)
      if(chrN=="X") { chrN <- xind }
      if(chrN=="Y") { chrN <- yind }
      chrN  <-  as.integer(chrN)
      centroArray[chrN, "start"] <- centro[(ict),2]
      centroArray[chrN, "end"] <- centro[(ict+1),3]
    }
  }else{
    centroArray <- NULL
  }
  return(centroArray)
}
