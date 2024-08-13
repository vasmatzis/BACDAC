#' Get object describing linear coordinate space for main chromosomes
#'
#' @param rgdObject sysdata object GRCh38 reference genome descriptor used for translating some chromosome names to numbers
#' @param chromosomes List of main chromosomes (technically .sva files) to form the coordinate system from
#' @return \link{linearCoordinates}
#'
#' @family coordinates
#'
#' @export
getLinearCoordinates <- function(rgdObject, chromosomes=1:24) {
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
#' @param rgdObject Reference genome descriptor, if NULL assumes autosomes= chr1-22 and allosomes= chr23,24
#' @param withChrPrefix if TRUE, will add "chr" to the front of every ichr value
#'
#' @return character vector of ichr
#'
convertChromToCharacter <- function(ichr,rgdObject=NULL, withChrPrefix=FALSE){
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
#' @param ideogram sysdata object with cytoband descriptions for GRCh38
#'
#' @return Centromere position array. Two-dimensions - first is chromosome number, the other separates start and end positions.
#'
#' @family ideogram
#'
#'
#' @export
getCentromerePositions <- function(ideogram) {
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
