#' Draw a linear genome plot
#'
#' Displays the read-depth array in a linear fashion. Can report CNV as all gray or use color calls.
#' If sampleId is provided, title annotations will be added.
#'
#' @param readDepthBinSize size of each window (bin) default is 30 kb
#' @param noDelAmpDetection do not color code deletions and gains in genome plot
#' @param gainColor color to use for gains in linear genome plot, default is blue
#' @param lossColor color to use for losses in linear genome plot, default is red
#' @param zebraStrips option to have alternating gray/white background for chromosome delineation
#' @param yAxisLimits option to set y axis to user specified limits, e.g. to align with the constellation plot
#' @param ... Parameters passed onto the actual plot command
#' @inheritParams commonParameters
#' @examples
#' sampleId='TCGA-14-1402-02A_ds'
#' inputDir <- system.file('extdata', package = "BACDAC")
#' outputDir=tempdir()
#' segmentationFile <- file.path(inputDir, paste0(sampleId, '_segmentation.csv'))
#' segmentation <- read.csv(segmentationFile, comment.char = '#') # chr, start, end, rd per
#' segmentation <- checkSegmentation(segmentation)
#' thirtyKbFile=file.path(inputDir, paste0(sampleId,'_','readDepthPer30kbBin.Rds'))
#' readDepthPer30kbBin = readRDS(file=thirtyKbFile )
#'
#' op <- par(mfrow=c(3,1),mai=c(.25,0.5, 0.3,0.25), mgp=c(2, .5, 0))
#' # default copy number color coding
#' linearGenomePlot(readDepthPer30kbBin=readDepthPer30kbBin,sampleId=sampleId,segmentation=segmentation)
#' # example with no copy number color coding
#' linearGenomePlot(readDepthPer30kbBin=readDepthPer30kbBin,sampleId=sampleId,noDelAmpDetection=TRUE,segmentation=segmentation)
#' # alternate copy number color coding, no sample id label
#' linearGenomePlot(readDepthPer30kbBin=readDepthPer30kbBin,sampleId=NULL,segmentation=segmentation,
#' gainColor = 'red', lossColor= 'blue' )
#'
#' @export
linearGenomePlot <- function( readDepthPer30kbBin, readDepthBinSize=30000, sampleId=NULL, alternateId=NULL, segmentation=NULL,
                              allelicSegments=NULL, noDelAmpDetection = FALSE, gainColor = 'blue', lossColor= 'red', zebraStrips=FALSE,
                              yAxisLimits=NULL, ...) {
  # readDepthBinSize=30000; yLimQuantile=0.99; noDelAmpDetection = FALSE;  gainColor = 'blue'; lossColor= 'red';zebraStrips=FALSE; alternateId=NULL

  if(R.version$major>=4){
    # with R 4.0 the default color palette changed, making reds and blues in the genome plot not true red and blue anymore.
    # specifying the palette as the original color scheme to avoid these muted reds and blues in the genome plot, on exit return palette to the default
    palette(value = 'R3')
    on.exit(palette(value='default'),add = TRUE)
  }

  mainChroms <- 1:24
  # We skip Y chromosome because hetScore does not make much sense there
  mainChromsNoY <- 1:23
  chromsToPlot = mainChroms

  coords <- getLinearCoordinates(mainChroms)

  x <- readDepthPer30kbBin$goodWindowArray
  y = readDepthPer30kbBin$readDepthArray
  ylabel <- paste("Counts per", readDepthBinSize/1000, 'kb bin')


  if(!is.null(segmentation)){
    delAmp   <- makeLegacyDelAmp(segmentation, coords = coords)
  }else{
    delAmp  <-  NULL
  }

  if(noDelAmpDetection){
    # We do not have color data, just use gray
    colorVector <- grayColorVector(x)
  }else{
    colorVector <- makeCNVcolorVector(x=x, delAmp,gainColor = gainColor, lossColor= lossColor)
  }

  if(length(x)>0) {
    chrCharacters <- convertChromToCharacter(coords@chroms)

    xlimit <- c(0, binnedPosEnd(max(coords@chromEnd[chromsToPlot]), readDepthBinSize))
    #ylimit=c(0,max(y)/2)
    if(is.null(yAxisLimits)){
      yLimQuantile=.99 # To prevent an outlier from dominating the plot, we set the max y axis to this quantile (number 0-1, 1=100%, 0.5=median)
      extraYAxis <- 1.3 # Fudge factor - add a little bit of extra space on top
      ylimit <- c(0, quantile(x=y, probs=yLimQuantile)/yLimQuantile*extraYAxis) # Stretch the y axis, assuming linear distribution of data
    }else{
      ylimit=yAxisLimits
    }

    plot(x, y, pch=".", xaxs="i", xlim=xlimit, xaxt="n", ylim=ylimit,
         ylab=ylabel,xlab="",cex.axis=1.3,col=colorVector,cex.lab=1.3) # Plot the summed array
    title(xlab='chromosome',cex.lab=1.4)

    ## option for Zebra bars for the chromosomes. Draw first so the dots can go over
    if(zebraStrips){
      for(i in seq_len(length(chromsToPlot))) {
        if(i%%2==1){
          rect(xleft = coords@chromStart[i],
               ybottom= par('usr')[3],
               xright = coords@chromEnd[i],
               ytop   = par('usr')[4],
               col=rgb(0.5,0.5,0.5,alpha=0.15), border=NA)
        }
      }
    }

    # Display the chromosome start as a line between the consecutive windows
    chromosomeStarts <- coords@chromStart[coords@chroms]/readDepthBinSize
    chromosomeEnds   <- coords@chromEnd[coords@chroms]/readDepthBinSize
    axis(side=3, at=chromosomeEnds, labels=NA, lwd=0, lwd.ticks = 1, tck = 1, col='darkgray') # Draw ticks
    axis(side=3, at=chromosomeStarts+(chromosomeEnds-chromosomeStarts)/2, line = -1.75, labels = chrCharacters, cex.axis=1.2, lwd=0, padj=0) # Draw labels


    # add purple lines for LOH segments
    # load segmentation TODO:  check for "minor" and "major" in segmentation file?
    if(!is.null(allelicSegments)){
      # we don't want to include the 1N segments
      minorZeroSegments <- allelicSegments[which(allelicSegments$minor_copy_number==0 & allelicSegments$major_copy_number>=2),]
      if(nrow(minorZeroSegments)>0){
        # convert to linear coordinates
        linPosStart <- abs(bimaToLinear(svaNumber=minorZeroSegments$chr, svaPos=minorZeroSegments$start) )
        linPosEnd   <- abs(bimaToLinear(svaNumber=minorZeroSegments$chr, svaPos=minorZeroSegments$end) )
        # convert to binned - linear coordinates
        linBinStart <- binnedPosStart(linPosStart, readDepthBinSize)
        linBinEnd <- binnedPosEnd(linPosEnd, readDepthBinSize)
        # add allele=0 segments to plot and data.frame
        segments(x0 = linBinStart, y0 = ylimit[1], x1 = linBinEnd, y1 = ylimit[1], col = 'purple', lwd=3)

        legend('topright', inset=c(0, .04),
               legend = c('2N+LOH'),col = c('purple'), lty=c(1), pch=c(NA), cex=1.3)
        #legend = c('2N+LOH', 'gain', 'loss'), col = c('purple',  gainColor, lossColor),lty=c(1,NA,NA), pch=c(NA, ".", "."), cex=.95, pt.cex=2)
      }
    }



  } else {
    plotEmptyLinearGenomePlot(chromsToPlot=chromsToPlot,coords=coords)
  }

  # ## annotate, if sampleId is provided
  if(!is.null(sampleId)) {
    title(main= sampleId,    adj= 0, line=1)
  }
  if(!is.null(alternateId)) {
    title(main= alternateId, adj= 1, line=1)
  }

}


#' create an empty linear genome plot to fill a space when data is not available
#'
#' @param chromsToPlot vector of chromosome (as integers) to plot
#' @param coords object describing linear coordinate space for the chromosomes
#' @keywords internal
plotEmptyLinearGenomePlot <- function(chromsToPlot,coords){
  # Get a plot started
  maxX <- max(coords@chromEnd[coords@maxcn])
  plot(x=0, y=0, type="n",
       xaxs="i",
       xlim=c(1,maxX), ylim=c(0, 1),
       xaxt="n", yaxt="n",
       xlab='', ylab="")
  ## annotations
  chrCharacters <- convertChromToCharacter(chromsToPlot)
  axis(side=3, at=coords@chromEnd[chromsToPlot], labels=NA, lwd=0, lwd.ticks = 1, tck = 1, col='lightgray') # Draw ticks
  axis(side=3, at=(coords@chromEnd[chromsToPlot]+coords@chromStart[chromsToPlot])/2, line = -2, labels = chrCharacters, cex.axis=0.85, lwd=0, padj=0) # Draw labels
  title(xlab='chromosome', line=0)

}



#' Provide the color coding for display on the genome plot
#'
#' @param x linear genome x-axis, windows, from the read-depth binned data
#' @param delAmp segmentation converted to the legacy delAmp format using \code{makeLegacyDelAmp}
#' @param gainColor color to use for gains in linear genome plot, default is blue
#' @param lossColor color to use for losses in linear geneome plot, default is red
#'

#' @return a vector of color codes to reflect the CNV calls: 2=red/loss, 4=blue/gains, 8=gray/normal; or 1=black/no color coding
#' @keywords internal
makeCNVcolorVector <- function(x, delAmp=NULL, gainColor = 'blue', lossColor= 'red'){
  if(!is.null(delAmp)) {

    lda <- length(delAmp)
    delAmpCol <- array(0,length(delAmp))
    delAmpCol[subset(seqFwd(1,lda),delAmp < 0)] <- 1
    delAmpCol[subset(seqFwd(1,lda),delAmp > 1)] <- 3

    # method to make the gains = 4, losses= 2, and normal = 1
    #  these values corresponds to the vector of colors in the default colors palette(), "black" "red" "green3" "blue" "cyan" "magenta" "yellow" "gray" ish
    colorVector <- 1+pmax(delAmpCol[x*3-2],
                          delAmpCol[x*3-1],
                          delAmpCol[x*3+0],
                          delAmpCol[x*3+1],
                          delAmpCol[x*3+2])
    colorVector[which(colorVector==1)] <- 8
  } else {
    # We do not have color data, just use gray which is 8
    colorVector <- grayColorVector(x)
  }

  # change color scheme if user has specified a non-default color for gains and losses
  # TODO: test different options
  if(gainColor!='blue'){
    colorVector[which(colorVector==4)] = gainColor
  }
  if(lossColor!='red'){
    colorVector[which(colorVector==2)] = lossColor
  }
  return(colorVector)
}



#' Provide gray color coding only for display on the genome plot
#'
#' @param x linear genome x-axis, windows, from the read-depth binned data
#' @keywords internal
grayColorVector=function(x){
  # We do not have color data, just use gray
  colorVector <- rep(8, length(x))
  return(colorVector)
}




#' for compatibility with genomePlot. Make old array from new structure: 1=normal, -1=loss, 2=gain
#' @param segmentation read depth data.frame with required columns: chr, start, end, rd; optional: cnvState for color coded linear genome plot
#' @param coords linear genome coordinates, loaded via \code{getLinearCoordinates()}
#'
makeLegacyDelAmp <- function(segmentation, coords) {
  wsz <- 10000
  # Initialize array data = 1, length = genomeLength/10k
  DELAMPar10save <- array(1, coords@chromEnd10K[24]) # was coords@totalLength10K but that included contigs etc. so shorten it to the mainChroms
  attr(DELAMPar10save, 'wsz') <- wsz

  for (i in seq_len(nrow(segmentation))) {
    # segmentation[i,]
    chrom_number <- segmentation[i, "chr"]
    edge1 <- segmentation[i, "start"]
    edge2 <- segmentation[i, "end"]
    if('cnvState' %in% names(segmentation)){
      interpretation <- segmentation[i, "cnvState"]
    }else{
      # do know know state so assign it to 2 aka 'normal'
      interpretation <- 2
    }

    binStart <- binnedPosStart(bimaToLinear(svaNumber = chrom_number, svaPos = edge1),
                               wsz)
    binEnd <- binnedPosEnd(bimaToLinear(svaNumber = chrom_number, svaPos = edge2),
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
