#' Draw a linear genome plot
#'
#' Displays the read-depth array in a linear fashion. Can report CNV as all black or use color calls.
#' If sampleId is provided, title annotations will be added. But should be specified in order to retrieve other data (segmentation, allelicSegData, normalPeakMethod),
#' unless those are passed in too.
#' Does not display structural variants.
#'
#' @param wsz Window size to use for binning when drawing the genome plot
#' @param yLimQuantile To prevent outliers from dominating the plot, we set the y axis to go up to
#'                     this quantile (number 0-1, 1=100\%, 0.5=median)
#' @param segmentation provides the color coding for each segment (dot), if it is not provided, try to load it. If no segmentation then color the segments gray
#' @param nrd plot y axis with normalize read depth to 2 for the median frequency
#' @param ... Parameters passed onto the actual plot command
#'
#' @inheritParams commonParameters
#'
#' @example inst/examples/linearGenomePlotExample.R
#'
#' @export
linearGenomePlot <- function( postProcessingDir, rgd, wsz=30000, yLimQuantile=0.99, sampleId=NULL, noDelAmpDetection = FALSE, segmentation=NULL,
                              nrd=FALSE, allelicSegData=NULL, normalPeakMethod=NULL, ...) {
  # wsz=30000; yLimQuantile=0.99; noDelAmpDetection = FALSE; segmentation=NULL; nrd=FALSE; allelicSegData=NULL;normalPeakMethod=NULL; sampleId=getSampleId(19015)
  if(R.version$major>=4){
    # with R 4.0 the default color palette changed, making reds and blues in the genome plot not true red and blue anymore.
    # specifying the palette as the original color scheme to avoid these muted reds and blues in the genome plot, on exit return palette to the default
    palette(value = 'R3')
    on.exit(palette(value='default'),add = TRUE)
  }

  rgdObject <- loadRgd(rgd)
  sampleAlias <- rgdObject$folderId

  coords <- getLinearCoordinates(rgd=rgdObject) # Linear coordinate system over the genome
  postProcessingDir <- rcfPrefix(postProcessingDir)

  cnvDetectFrequency <- loadCnvBinned(postProcessingDir)


  if(cnvDetectFrequency[['version']]>=2) {
    frq30 <- cnvBinnedToLegacy(shrinkCnvBinned(cnvDetectFrequency, wsz=wsz))
  } else {
    frq30 <- cnvDetectFrequency
    wsz <- 30000 # The window size is fixed
  }

  y <- frq30[['frq']]
  x <- frq30[['wdnsMSK']]

  if(nrd){
    expectedNormalBin <- getExpectedNormalBin(cnvDetectFrequency)
    medfrq2N <- expectedNormalBin * wsz # medfrq of all the 2N level which may or may not be the dominate level
    # medfrq <- median(y) #  medfrq of all the data, will be above 2N level if the dominate level is 3N or higher

    y <- y/(2*medfrq2N)
    ylabel <- "Normalized Read Depth (NRD)"
  }else{
    ylabel <- "Counts per window"
  }

  # TODO: this was using DELAMPar10saveFile which is basically replaced by segmentation... need to convert to that and stop making DELAMPar10saveFile in cnvDetect
  # segmentation provides the color coding for each segment (dot), if it is not provided, try to load it. If no segmentation then color the segments gray
  segmentationFile <-getTypedFile("segmentation",dir = postProcessingDir,values = list(sampleId = sampleId),legacy = TRUE )
  if(is.null(segmentation)){
    if(file.exists(segmentationFile@path)) {
      segmentation <- bmd.read.csv(segmentationFile)
    }
  }
  if(is.null(normalPeakMethod)){
    cnvMetadata <- readMetadata(segmentationFile)
    if(!is.null(cnvMetadata[['normalPostProcessingDir']])) {
      normalPeakMethod <- cnvMetadata[['normalPeakMethod']]
    }
  }

  if(!is.null(segmentation)){
    delAmp   <- makeLegacyDelAmp(segmentation,   rgdObject = rgdObject, coords = coords)
  }else{
    delAmp  <-  NULL
  }

  colorVector <- makeCNVcolorVector(postProcessingDir,x=x, noDelAmpDetection, delAmp)


  if(length(x)>0) {

    chromosomesToDisplay <- union(svaAutosomes(rgdObject), svaAllosomes(rgdObject))
    chrCharacters <- convertChromToCharacter(coords@chroms,rgdObject=rgdObject)

    xlim <- c(0, binnedPosEnd(max(coords@chromEnd[chromosomesToDisplay]), wsz))
    #ylim=c(0,max(y)/2)
    extraYAxis <- 1.3 # Fudge factor - add a little bit of extra space on top
    ylim <- c(0, quantile(x=y, probs=yLimQuantile)/yLimQuantile*extraYAxis) # Stretch the y axis, assuming linear distribution of data


    plot(x, y, pch=".", xaxs="i", xlim=xlim, xaxt="n", ylim=ylim,
         ylab=ylabel,xlab="",cex.axis=1,col=colorVector) # Plot the summed array
    title(xlab='chromosome', line=0)

    # add purple lines for LOH segments
    # load allelicSegData if possible
    if(is.null(allelicSegData)){
      if(is.null(normalPeakMethod)){
        segmentationFile <-getTypedFile("segmentation",dir = postProcessingDir,values = list(sampleId = sampleId),legacy = TRUE )
        cnvMetadata <- readMetadata(segmentationFile)
        if(!is.null(cnvMetadata[['normalPostProcessingDir']])) {
          normalPeakMethod <- cnvMetadata[['normalPeakMethod']]
        }
      }


      if(!is.null(normalPeakMethod)){
        if(normalPeakMethod=='ploidyBased'){
          ploidySegments <- getPloidySegments(postProcessingDir, sampleId, stopIfNoPloidySegments=stopIfNoPloidySegments) # load this separately to do all the checks and stuff
          if(all(!is.na(ploidySegments))){
            hetScores_StarCloudDataInfo <- getTypedFile("hetScores_StarCloudData", dir = postProcessingDir, values = list(sampleId = sampleId),legacy = TRUE)
            starCloudData <- loadRdata(file=hetScores_StarCloudDataInfo@path, fileLabel = 'output from the star-cloud plot' )
            allelicSegData <- allelicCNV(starLookUp = starCloudData$starLookUp, segmentDataIn = ploidySegments)
          }else{
            loginfo('%s: ploidySegments is.na cannot add allelicSegData to linear plot.',sampleId)
          }
        }
      }

    }

    if(!is.null(allelicSegData)){
      # we don't want to include the 1N segments
      minorZeroSegments <- allelicSegData[which(allelicSegData$minor==0 & allelicSegData$major>=2),]
      if(nrow(minorZeroSegments)>0){
        # convert to linear coordinates
        linPosStart <- abs(bimaToLinear(rgd=rgdObject,  svaNumber=minorZeroSegments$chr, svaPos=minorZeroSegments$start) )
        linPosEnd   <- abs(bimaToLinear(rgd=rgdObject,  svaNumber=minorZeroSegments$chr, svaPos=minorZeroSegments$end) )
        # convert to binned - linear coordinates
        linBinStart <- binnedPosStart(linPosStart, wsz)
        linBinEnd <- binnedPosEnd(linPosEnd, wsz)
        # add allele=0 segments to plot and data.frame
        segments(x0 = linBinStart, y0 = 0, x1 = linBinEnd, y1 = 0, col = 'purple', lwd=3)
      }
    }



    # Display the chromosome start as a line between the consecutive windows
    chromosomeStarts <- coords@chromStart[coords@chroms]/wsz
    chromosomeEnds   <- coords@chromEnd[coords@chroms]/wsz
    axis(side=3, at=chromosomeEnds, labels=NA, lwd=0, lwd.ticks = 1, tck = 1, col='darkgray') # Draw ticks
    axis(side=3, at=chromosomeStarts+(chromosomeEnds-chromosomeStarts)/2, line = -2, labels = chrCharacters, cex.axis=0.85, lwd=0, padj=0) # Draw labels
  } else {
    plot(x=c(), y=c(), xlim=c(0,1), ylim=c(0,1), type="n", xlab="", ylab="", ...)

    ## Draw a cross signaling we had no data
    lines(x=c(0, 1), y=c(1, 0), type="l")
    lines(x=c(0, 1), y=c(0, 1), type="l")
  }
  # ## annotate, if sampleId is provided
  # if(!is.null(sampleId)) {
  #   title(main= sampleId,    adj= 0, line=1)
  #   title(main= sampleAlias, adj= 1, line=1)
  # }


}


#' Provide the CNV color coding for display on the genome plot (U-shaped and linear)
#'
#' @param x x-axis data from the CNV binned data, aka frq30[['wdnsMSK']]
#' @param delAmp cnvIntervals converted to the legacy delAmp format using \code{makeLegacyDelAmp}
#'
#' @inheritParams commonParameters
#'
#' @return a vector of color codes to reflect the CNV calls: 2=red/loss, 4=blue/gains, 8=gray/normal; or 1=black/no color coding
#'
makeCNVcolorVector <- function(postProcessingDir, x, noDelAmpDetection, delAmp=NULL){
  if(!noDelAmpDetection && !is.null(delAmp)) {

    lda <- length(delAmp)
    delAmpCol <- array(0,length(delAmp))
    delAmpCol[subset(seqFwd(1,lda),delAmp < 0)] <- 1
    delAmpCol[subset(seqFwd(1,lda),delAmp > 1)] <- 3

    colorVector <- 1+pmax(delAmpCol[x*3-2],
                          delAmpCol[x*3-1],
                          delAmpCol[x*3+0],
                          delAmpCol[x*3+1],
                          delAmpCol[x*3+2])
    colorVector[which(colorVector==1)] <- 8
  } else {
    # We do not have color data, just use gray
    colorVector <- rep(8, length(x))
  }
  return(colorVector)
}


#' for compatibility with genomePlot. Make George's old array from new structure
#' @param delGainInfo cnvIntervals loaded via \link{typedFile-class}
#' @param rgdObject reference genome descriptor object loaded via \code{loadRgd()}
#' @param coords linear genome coordinates, loaded via \code{getLinearCoordinates(rgd=rgdObject)}
#'
makeLegacyDelAmp <- function(delGainInfo, rgdObject, coords) {
  wsz <- 10000
  # Resulting array
  DELAMPar10save <- array(1, coords@totalLength10K)
  attr(DELAMPar10save, 'wsz') <- wsz

  for (i in seq_len(nrow(delGainInfo))) {
    chrom_number <- delGainInfo[i, "chr"]
    edge1 <- delGainInfo[i, "start"]
    edge2 <- delGainInfo[i, "end"]
    interpretation <- delGainInfo[i, "cnvState"]

    binStart <-
      binnedPosStart(bimaToLinear(rgdObject, svaNumber = chrom_number, svaPos = edge1),
                     wsz)
    binEnd <-
      binnedPosEnd(bimaToLinear(rgdObject, svaNumber = chrom_number, svaPos = edge2),
                   wsz)

    if (interpretation == 2) {
      value <- 1 # normal
    } else if (interpretation < 2) {
      value <- -1 # loss
    } else if (interpretation > 2) {
      value <- 2 # gain
    }

    DELAMPar10save[binStart:binEnd] <- value
  }

  return(DELAMPar10save)
}
