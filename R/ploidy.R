#' Use density on 100 kb binned read depth to determine read depth distribution (peaks)
#'
#' peaks are returned in rank order of read depth with highest frequency read depth peak listed first
#'
#' @param sampleId sample Identifier
#' @param alternateId secondary sample identifier
#' @param readDepthPer100kbBin  read depth for 100kb bins, with bin index in linear coordinates
#' @param segmentation identified regions of the genome with constant read depth
#' @param segmentationBinSize bin size used for the read depth in the segmentation data
#' @param wszPeaks window size for peaks by density
#' @param grabDataPercentManual portion of main peak data to grab, other peaks will be scaled based on read depth (x location) set to -1 to use main Peak width
#' @param origMaxPercentCutoffManual peaks smaller than this portion of the max peak are not considered; set to -1 to use default value
#' @param qualityPostNorm used to set max allowed value for grabDataPercent
#'
#' @inheritParams commonParameters
#'
peaksByDensity <-function(sampleId,readDepthPer100kbBin, segmentation, segmentationBinSize=30000, wszPeaks = 100000, grabDataPercentManual= -1, origMaxPercentCutoffManual=-1,
                          addAreaLinesToPlot=FALSE,qualityPostNorm=NULL,pause=FALSE, omitAnnotations=FALSE,alternateId=NULL){
  # peaksByDensity() provided by Jamie, tweaked by Roman, and then Sarah, plotting added by Sarah
  # peaksByDensity(sampleId,  cnvBinnedData, segmentation=segmentation, wszPeaks = 100000, grabDataPercentManual= 0.08,pause=F);
  #  grabDataPercentManual= -1; segmentation=segmentation; segmentationBinSize=30000; wszPeaks = 100000; addAreaLinesToPlot=F; pause=F;origMaxPercentCutoffManual=-1;  qualityPostNorm=NULL;omitAnnotations=FALSE

  coords <- getLinearCoordinates(chromosomes = 1:numChroms)

  ### get frequency array ---------------
  frqToUse <- readDepthPer100kbBin$readDepthArray
  wdnsToUse <- readDepthPer100kbBin$goodWindowArray


  ### do density
  bandwidth <- bw.nrd(frqToUse)
  meanFrqFactor <- 4  # originally this was 8, but down below it was 4, assigning a variable to be consistent.
  denTempOrig <- density(frqToUse,bw=bandwidth,from=1,to=meanFrqFactor*mean(frqToUse),n=4096) # choose n to be a factor of 2;
  # TODO: use densityTo(mean(frqToUse))??
  origMaxPeakY <- denTempOrig$y[which.max(denTempOrig$y)]
  origMaxPeakX <- denTempOrig$x[which.max(denTempOrig$y)]


  ###  one peak check ----------
  # see if 'origMaxPercentCutoff' can be lowered in order to get more than one peak
  # sometimes this works to pull in a really low level peak, not necessary for ploidy but necessary to get a tumor percent, see: 19039, 19006
  # don't want to pull in too many peaks and end up with subclones though, hence
  #  1) we only do this check when there is only 1 peak to begin with and 2) we don't split the peak here like we do below
  # keep in mind this new peak also won't have many data points for the het Score density test.

  if(origMaxPercentCutoffManual > 0){  #
    origMaxPercentCutoff <- origMaxPercentCutoffManual
  }else{
    origMaxPercentCutoffInitial <- 0.025 # default value

    # Find the peaks
    tpTemp <- pastecs::turnpoints(denTempOrig$y)
    peaksXTemp  <- denTempOrig$x[tpTemp$pos[which(tpTemp$peaks & tpTemp$points > origMaxPercentCutoffInitial*origMaxPeakY)]]

    # lower the cutoff if applicable
    if(length(peaksXTemp)==1){ # only one peak
      origMaxPercentCutoffTemp <- origMaxPercentCutoffInitial
      while(length(peaksXTemp)==1 & origMaxPercentCutoffTemp >= 0.016){ # choose 0.016 so that after subtracting 0.006 you are not lower than 0.01
        origMaxPercentCutoffTemp <- origMaxPercentCutoffTemp-0.006
        peaksXTemp  <- denTempOrig$x[tpTemp$pos[which(tpTemp$peaks & tpTemp$points > origMaxPercentCutoffTemp*origMaxPeakY)]]
      }

      # we successfully found another peak
      if(length(peaksXTemp)>1){
        logwarn('origMaxPercentCutoffTemp lowered to %s to find another peak, now have %i-ish peaks',origMaxPercentCutoffTemp,length(peaksXTemp))
        origMaxPercentCutoff <- origMaxPercentCutoffTemp
      }else{
        origMaxPercentCutoff <- origMaxPercentCutoffInitial
      }
    }else{
      origMaxPercentCutoff <- origMaxPercentCutoffInitial
    }
  }


  ##############'

  # convert rd to the proper window size, match window size for peak determination
  if(!'rd' %in% names(segmentation) ){
    stop('rd column missing from segmentation file')  # could be converted from 'nrd' but shouldn't have to
  }else{
    segmentation[,'rd_wszPeaks']  <- segmentation[,'rd'] * wszPeaks / segmentationBinSize  # rd was stored using segmentationBinSize=30000, so convert for use with wszPeaks
  }


  ### start of plot -------------------
  # Assume a peak that is greater than origMaxPercentCutoff the size of the main peak has value and we want to consider it for the digital peak algorithm.
  logdebug("origMaxPercentCutoff: %s",origMaxPercentCutoff) # Jamie's original heuristic was 2%

  # if(!skipExtras){
  op <- par(mfrow=c(2,1),mar=c(2.75, 3.5, 2, 1.5),mgp=c(1.5, 0.5,0))

  # plot stuff
  xMaxPlot <- min(2.75*mean(frqToUse),quantile(frqToUse,probs = .99))*1.07
  xMinPlot <- min(20000,              quantile(frqToUse,probs = .0051))*0.9

  plot(denTempOrig$x, denTempOrig$y, type='n',  xlim=c(xMinPlot, xMaxPlot), xlab = '', ylab = '');
  title(xlab=paste0('read depth per window (',wszPeaks/1000,' kb)')) #(passing mask)
  mtext(3, text=c(sampleId, alternateId),adj=c(0,1),line=0.5)
  mtext(3, text='Peak Rank',adj=0, line=0, cex=.7)


  polygon(denTempOrig$x, denTempOrig$y, col='gray92',border='gray92')
  abline(h=origMaxPercentCutoff*origMaxPeakY, col='orange')
  mtext(text = origMaxPercentCutoff, side=2, line=-2, at = origMaxPercentCutoff*origMaxPeakY, cex=.75, col='orange', las=2,padj = -1 )
  # }

  # find grabDataPercent
  # by using the main peak width if there is more than one peak (dataBasedGrabDataPercent by getting the pit distance from the main peak)
  # or
  # using a value that was manually provided
  if(grabDataPercentManual < 0){  #

    if(!is.null(qualityPostNorm)){
      dataBasedGrabDataPercentMax <- ifelse(qualityPostNorm < 2.5, 0.15, 0.17)
    }else{
      logdebug('need qualityPostNorm to be able to determined dataBasedGrabDataPercentMax, will set this max to 0.15')
      dataBasedGrabDataPercentMax <- 0.15
    }
    #Find the peaks for the unclassified data
    tpOrig <- pastecs::turnpoints(denTempOrig$y)
    peaksX  <- denTempOrig$x[tpOrig$pos[which(tpOrig$peaks & tpOrig$points > origMaxPercentCutoff*origMaxPeakY)]]
    peaksY  <- denTempOrig$y[tpOrig$pos[which(tpOrig$peaks & tpOrig$points > origMaxPercentCutoff*origMaxPeakY)]]

    if(length(peaksX)>=1){ # at least one peak

      # (pitScoresX  <- denTempOrig$x[pastecs::extract(tpOrig, no.tp = FALSE, peak = FALSE, pit = TRUE)] )
      # (pitScoresY  <- denTempOrig$y[pastecs::extract(tpOrig, no.tp = FALSE, peak = FALSE, pit = TRUE)])
      # allPeaksX <- dx[tp$pos[which (tp$peaks)]]
      # allPeaksY <- dy[tp$pos[which (tp$peaks)]]
      mainPeakReadDepth    <- peaksX[which.max(peaksY)]

      allPitsX <- denTempOrig$x[tpOrig$pos[which (tpOrig$pits)]] # aka pitScoresX
      #allPitsY <- denTempOrig$y[tpOrig$pos[which (tpOrig$pits)]] # aka pitScoresY

      pitL <- max(allPitsX[(allPitsX < mainPeakReadDepth* .99)],0) # 0 if there are no left peaks, *.99 just in case there is an adjacent peak, don't pick it i.e. MD66263
      pitR <- min(allPitsX[(allPitsX > mainPeakReadDepth*1.01)])
      # abline(v=c(pitL,pitR),col=c(2,3))
      pitR-pitL

      ratioRight <- -1 + pitR/mainPeakReadDepth
      ratioLeft  <-  1 - pitL/mainPeakReadDepth
      # make sure the grabDataPercent is not bigger than a given max, or smaller than 0.04)
      whoaDoNotGoSmallerThanThis <- 0.04
      dataBasedGrabDataPercent <- max(whoaDoNotGoSmallerThanThis,
                                      min(ratioRight*.95,ratioLeft*.95, dataBasedGrabDataPercentMax) )

      grabDataPercent <- round(dataBasedGrabDataPercent,3)
      loginfo('grabDataPercent: %s based on mainPeak width',grabDataPercent)
    }else{
      grabDataPercent <- 0.08 # default
      loginfo('grabDataPercent: %s default for one peak',grabDataPercent)
    }
  }else{
    grabDataPercent <- grabDataPercentManual
    loginfo('grabDataPercent: %s from manual input',grabDataPercent)
  }

  ### initialize necessary lists and parameters for the peak determination loop.
  i <- 0
  colNew <- rep(0,length(frqToUse))
  peakReadDepthList <- NULL
  peakHeightList <- NULL
  peakInt <- 0
  minimX <- 0

  cnColors <- getCNcolors()

  scaledGrabDataPercentPerPeak <- data.frame()
  while(length(minimX)>0) {
    #Take the points that have been deemed aberrant and have not yet been sorted into a peak (colNew value = 0)
    denTemp <- density(frqToUse[which(colNew==0)],bw=bandwidth,from=1,to=meanFrqFactor*mean(frqToUse),n=4096)
    dx <- denTemp$x
    dy <- denTemp$y * ((length(which(colNew==0)))/(length(colNew)))
    dx[which.max(dy)]

    #Find the peaks for the unclassified data
    tp <- pastecs::turnpoints(dy)

    #Get all of them that seem interesting (i.e. big enough)
    #Final heuristic. Assume a peak that is greater than origMaxPercentCutoff the size of the original peak has value. Smaller and not sure.
    minimX <- dx[tp$pos[which(tp$peaks & tp$points > (origMaxPercentCutoff-0.001)*origMaxPeakY)]] # subtract 0.001 to use a slightly smaller value to account for previous extraction # 66257
    minimY <- dy[tp$pos[which(tp$peaks & tp$points > (origMaxPercentCutoff-0.001)*origMaxPeakY)]]

    # added plot stuff
    i <- i+1
    # if(!skipExtras){
    points(dx, dy, col=cnColors[i], type='l')

    # abline(v=minimX[which.max(minimY)], col=cnColors[i], lwd=3, lty='dashed') # max peak
    if(addAreaLinesToPlot)abline(v=minimX, col=cnColors[i])                           # all peaks
    polygon(dx, dy, col=cnColors[i],border=cnColors[i])
    # }

    #If we seem to have found some peaks find the biggest one and add to the lists
    if (length(minimX)>0) {
      newpeakReadDepth    <- minimX[which.max(minimY)]
      newpeakHeight <- minimY[which.max(minimY)]
      peakReadDepthList <- c(peakReadDepthList,newpeakReadDepth)
      peakHeightList <- c(peakHeightList,newpeakHeight)
      #The new peak will be classified in the DELAMPar10Test array
      peakInt <- peakInt+1

      # grabDataPercent is a hardcoded value, but a single value will be too big for higher coverage (right) peaks and too small for lower coverage (left) peaks
      # scale the value so it decreases slightly for higher coverage peaks, thus the range for the data selected is still expanding, but at a slower rate.
      mainPeakReadDepth <- peakReadDepthList[1] # the mainPeak is the one detected and recorded first
      scaledGrabDataPercent <- grabDataPercent/sqrt(newpeakReadDepth/mainPeakReadDepth)

      scaledGrabDataPercentPerPeak <- rbind(scaledGrabDataPercentPerPeak, cbind(peakRank=peakInt,scaledGrabDataPercent))

      # if(!skipExtras){
      mtext(peakInt, side=3, line=0, at = newpeakReadDepth)
      points(newpeakReadDepth,newpeakHeight,pch=21,bg='black',col='white')

      # extra lines which tend to make the plot cluttered but can be useful for debugging
      if(addAreaLinesToPlot){
        mtext(paste(peakInt,round(scaledGrabDataPercent,3), sep=': '), side=3, line= -peakInt, adj=1)

        # lines to illustrate what data will be grabbed by "whichChoose"
        abline(v=newpeakReadDepth*(1+scaledGrabDataPercent), col=cnColors[i], lty='dashed', lwd=0.8)
        abline(v=newpeakReadDepth*(1-scaledGrabDataPercent), col=cnColors[i], lty='dashed', lwd=0.8)
      }
      # }

      # TODO: if X and Y are part of frqToUse,  then 'whichChoose' or maybe cnvSeq$rd needs to be adjusted because
      # currently there is a bug such that X and Y aren't extracted from the peaks correctly.
      whichChoose <- unlist(sapply(which(abs(segmentation$rd_wszPeaks-newpeakReadDepth)/newpeakReadDepth < scaledGrabDataPercent),function(x)
        which(wdnsToUse >=  binnedPosStart(coords@chromStart[segmentation[x,'chr']], binSize = wszPeaks) + segmentation[x,'start']/wszPeaks &
                wdnsToUse <=  binnedPosStart(coords@chromStart[segmentation[x,'chr']], binSize = wszPeaks) + segmentation[x,'end']  /wszPeaks ) ))

      colNewKeys <- which(colNew[whichChoose]==0)
      if(length(colNewKeys) > 100){
        colNew[whichChoose][colNewKeys] <- peakInt  #  only if colNew was still zero, change it to peakInt
        table(colNew)
      }else{
        # remove that last peak, it is too small, no good, removed it and exit
        peakReadDepthList <- peakReadDepthList[-peakInt]
        peakHeightList    <- peakHeightList[-peakInt]
        scaledGrabDataPercentPerPeak <- scaledGrabDataPercentPerPeak[-peakInt, ]
        # if(!skipExtras){
        mtext(text ='-------------', side=3, line= -peakInt, adj=1)
        # }


        break # jump out of while loop, not enough new data in peak to continue searching for peaks
      }
    }

    if(peakInt >= 15){
      # prevent a non-stop while loop, if too many peaks are found, exit out of while loop
      # digital grid while loop will also stop after 'numPeaks' loops
      break
    }

    if(interactive() & pause){BBmisc::pause()}
  }

  # annotate plot with the approx number of fragments
  # if(!skipExtras){
  if(!omitAnnotations){
    mtext(1, text='~fragments',adj=1, line=-3)
    mtext(1, text='chr1-22 ',adj=1, line=-2)
    mtext(1, text=paste(round(sum(frqToUse)/1000000,1), 'mil') ,adj=1, line=-1)
  }
  # }


  # peakReadDepthList_per1bp = X positions of peaks (normalized as if window size was 1), this is so the output is consistent no matter what wsz was used, don't need to know what wsz was used later to use the values
  # normalize the peakHeightList so biggest peak has a height of one
  return(list(peakReadDepthList_per1bp = peakReadDepthList/wszPeaks,
              peakHeightList_maxIsOne  = peakHeightList/max(peakHeightList),
              scaledGrabDataPercentPerPeak    = scaledGrabDataPercentPerPeak))
}

#' each peak gets a group of grid heights based on the mode of the peak
#'
#' @inheritParams commonParameters
#'
assignGridHeights <- function(peakInfo, n00, minGridHeight){
  gridModes <- round(peakInfo$peakReadDepth_normX*n00)
  gridHeights <- array(0,c(10*n00))
  nPeaks <- nrow(peakInfo)
  # NOTE: making the grid with simple 0/1 windows can be very sensitive to noise
  # consider using e.g. gaussian kernel of larger size (approximated below) instead of 0/1 windows

  # assign values to gridPositions based on
  #   1) the height of the peak (bonus)
  #   2) distance from the peak: directly on the peak = 1, decreasing as you move away from the peak, (.95, .75.5)
  #   3) min value is no less than minGridHeight, the center of the peak is slightly higher

  # Using an artificially growing digital binary grid
  fillInGrid <- function(gHeight, gridModes, bonus, minGridHeight){
    gHeight[gridModes - 4] <- pmax(.50 * bonus, minGridHeight)
    gHeight[gridModes - 3] <- pmax(.75 * bonus, minGridHeight)
    gHeight[gridModes - 2] <- pmax(.95 * bonus, minGridHeight)
    gHeight[gridModes - 1] <- pmax(.99 * bonus, minGridHeight)
    gHeight[gridModes]     <- pmax( 1  * bonus, minGridHeight+0.025)
    gHeight[gridModes + 1] <- pmax(.99 * bonus, minGridHeight)
    gHeight[gridModes + 2] <- pmax(.95 * bonus, minGridHeight)
    gHeight[gridModes + 3] <- pmax(.75 * bonus, minGridHeight)
    gHeight[gridModes + 4] <- pmax(.50 * bonus, minGridHeight)
    return(gHeight)
  }

  ## the two methods for filling in the gridPositions
  # makeTestPlots: will plot two graphs for comparing the differences in methods
  if(TRUE){
    makeTestPlots <- FALSE
    if(makeTestPlots){
      ## fill in gridModes at once, which will fill in left to right, so right edge of peak will overwrite left edge of next peak if they overlap
      gridHeights <- fillInGrid(gridHeights, gridModes, peakInfo[, 'bonus'], minGridHeight)
      plot(which(gridHeights > 0), gridHeights[which(gridHeights > 0)])
      ## remove first and last grid mode for the first peak, to force a tighter fit to the first peak; it should have a smaller distribution in read depth
      gridHeights[96] <- 0
      gridHeights[104] <- 0
      points(which(gridHeights > 0), gridHeights[which(gridHeights > 0)], pch=19)
      grid(ny=NA)
    }
    #
    # ## fill in lower peaks first which will fill in bottom to top, giving preference to bigger peaks, so bigger peak will overwrite smaller peak if they overlap
    if(makeTestPlots) plot(which(gridHeights > 0), gridHeights[which(gridHeights > 0)], type = 'n')
    for(i in nPeaks:1){
      iKey <- which(peakInfo$rankByHeight==i)
      iGridMode <- gridModes[iKey]
      gridHeights <- fillInGrid(gridHeights, iGridMode,peakInfo[iKey, 'bonus'], minGridHeight)
      if(makeTestPlots){ points(which(gridHeights > 0), gridHeights[which(gridHeights > 0)],col=iKey)}

      ## remove first and last grid mode for the first peak, to force a tighter fit to the first peak; it should have a smaller distribution in read depth
      if(iKey==1){
        gridHeights[96] <- 0
        gridHeights[104] <- 0
      }
    }
    if(makeTestPlots){
      points(which(gridHeights > 0), gridHeights[which(gridHeights > 0)], pch=19)
      grid(ny=NA)
    }
  }

  if(FALSE){
    gridHeights <- assignGridHeights(peakInfo, n00, minGridHeight )
    nonZeroGridCoords <- which(gridHeights > 0)
    gridHeights[which(gridHeights > 0)]
  }

  return(gridHeights)
}



#' fit a group of peaks to a digital grid
#'
#' @param peakInfo summary table of info for each peak found in \code{peaksByDensity}
#' @param gridHeights group of 8-10 grid heights set for each peak
#' @param gridIteration iteration count for digital grid
#' @param numOfGridCoordsToTest number of grid coordinates to test for the first digital peak
#' @param previousPeriod period used in the previous \code{gridIteration}
#' @param bonusFactor factor for assigning bonuses to the peaks
#' @param iterationStats stats collected during each \code{gridIteration}
#' @param omitAnnotations should extra annotation be included in digital grid plot
#'
#' @inheritParams commonParameters
#'
#' @return dPeaks grid coordinates of the digital peaks, where NA is a peak not on the digital grid, and nCopyPeaks_dig-first digital peak set to 1
digitalGrid <- function(peakInfo, gridHeights,
                        dPeaksCutoff, penaltyCoefForAddingGrids,
                        n00, bonusFactor=8, sampleId=NULL, alternateId=NULL,
                        gridIteration, minPeriodManual=minPeriodManual,maxPeriodManual=maxPeriodManual,
                        numOfGridCoordsToTest=NULL,
                        previousPeriod=NULL,
                        minGridHeight,iterationStats,omitAnnotations=FALSE){
  # dPeaksCutoff=0.01; minGridHeight=0.2; penaltyCoefForAddingGrids=0.049; n00=100;bonusFactor=8; previousPeriod=NULL;omitAnnotations=FALSE

  gridModes       <- round(peakInfo$peakReadDepth_normX*n00)
  peakRD_normX    <- peakInfo$peakReadDepth_normX
  peakRD_1bp      <- peakInfo$peakReadDepth_1bp
  peakHeights     <- peakInfo$peakHeight
  peakHeightRanks <- peakInfo$rankByHeight
  bonus           <- peakInfo$bonus
  numPeaks        <- nrow(peakInfo)
  mainPeakIndex   <- which(peakInfo$rankByHeight==1)


  nonZeroGridCoords <- which(gridHeights > 0)
  maxGridCoord  <- max(which(gridHeights > 0))
  minGridCoord  <- min(which(gridHeights > 0))
  gridHeights[which(gridHeights > 0)]
  gridHeights[gridModes]




  ### find if gridModes align with the digital grid
  # fit digital grid
  fitGridPeriod <- function(gridHeights,numOfGridCoordsToTest, minPeriod, maxPeriodToTry=NULL,plotDigDebug=FALSE){
    # plotDigDebug=TRUE

    maxGridCoord <- max(which(gridHeights > 0))
    minGridCoord <- min(which(gridHeights > 0))
    gShiftsToTry <- which(gridHeights>0)[1:numOfGridCoordsToTest]
    minPeriodForLoop <- minPeriod

    # maximize the number of peaks hit while minimizing the number of periods
    # maximize the number of peaks hit while maximizing the size of the period. (TODO: would this be a reason to have periodsToTry decreasing rather than increasing?)

    # initialize
    iCol <- 0  # for plotDebug
    maxSum <- -Inf
    maxShift <- -1
    maxPeriod <- -1
    logdebug('gShift - gPeriod - sumD | scoreD - numPeaksHit - numPeriods - penaltyMaxRo | penaltyMax | penaltyRZ | penaltyToUse')


    # gShift = starting point for the grid
    # gPeriod = the distance between each grid
    for(gShift in gShiftsToTry ) {

      # periodsToTry: define the sequence of periods to try
      #TODO: periodsToTry -- originally it was decreasing, what if it is increasing?
      # if(minPeriodForLoop < gShift){
      # periodsToTry=gShift:minPeriodForLoop
      # }else{
      #   periodsToTry= seqFwd(from=minPeriodForLoop, to=(maxGridCoord-minGridCoord))
      # }

      if( (maxGridCoord-minGridCoord) > minPeriodForLoop ){
        toPeriod <- (maxGridCoord-minGridCoord)
        if(!is.null(maxPeriodToTry)) toPeriod <- min(maxPeriodToTry, toPeriod)
        periodsToTry <- seqFwd(from=minPeriodForLoop, to=toPeriod)
      }else{
        periodsToTry <- minPeriodForLoop
      }

      for(gPeriod in periodsToTry) {
        bins0    <- gShift + seq(from = 0, to = maxGridCoord - gShift, by = gPeriod) # bins0: proposed grid lines

        digtGrid <- gridHeights*0 # re-initialize
        digtGrid[bins0] <- 1  # digtGrid[bins0]*gridHeights[bins0]

        # Score is positive on hits, but we penalize too fine-grained grids
        # sumD = score: how many gridlines overlap a red dot minus the penalty
        maxNumPeriodsRo  <- round( (maxGridCoord  -minGridCoord)/gPeriod,1) # does not depend on gShift (start point), not as sensitive to the gPeriod
        maxNumPeriods  <- ((maxGridCoord  -minGridCoord)/ gPeriod) # does not depend on gShift (start point)
        posNumPeriods <- ((maxGridCoord - gShift) / gPeriod )     # depends on gShift, will be smaller if gShift is not the first dot, especially smaller if the first peak is skipped
        penaltyMaxRo <- maxNumPeriodsRo  * penaltyCoefForAddingGrids  #
        penaltyMax <-  maxNumPeriods * penaltyCoefForAddingGrids  #
        penaltyRZ <-   posNumPeriods * penaltyCoefForAddingGrids  # favors skipping the first peak because it will be smaller

        penaltyToUse <- penaltyMax

        numPeaksHit <-  sum((digtGrid*gridHeights)>0)
        numPeriods  <-  length(bins0)

        sumD <- sum(digtGrid*gridHeights) - penaltyToUse
        scoreD <- round(sumD/numPeaksHit,4)
        if(sumD > maxSum) {
          maxScore <- scoreD
          maxSum <- sumD
          maxSumShift <- gShift
          maxSumPeriod <- gPeriod
          # TODO: should the minPeriodForLoop always be minPeriod or reset to gPeriod??
          # minPeriodForLoop <-gPeriod # reset the min period so it can not be smaller than the previously found period
          minPeriodForLoop <-minPeriod # reset the min period
          # output the values each time a new max score is found

          # print(paste('maxNumPeriods:', round(maxNumPeriods,3),'possNumPeriods:', round(possNumPeriods,3) ))
          logdebug(paste(gShift,'----',gPeriod,'----',round(maxSum,4), '|', round(maxScore,4), '----', numPeaksHit, '----', numPeriods, '----',round(penaltyMaxRo,4),'|',round(penaltyMax,4),'|',round(penaltyRZ,4), '|',round(penaltyToUse,4))   )

          if(plotDigDebug){
            # make a new plot since the previous one is probably getting cluttered
            if(iCol %% 8 ==0){
              yMin <- -bonusFactor/10
              plot(x=gridModes,
                   y=peakHeights,
                   main=paste('Digital grid iteration:',gridIteration),
                   xlim=c(min(gridModes)-10, max(gridModes)+10),
                   ylim=c(yMin, 1), cex=0.65, xlab='grid coordinate',
                   pch=19)
              points(x=which(gridHeights>0), y=yMin + gridHeights[gridHeights>0] * 0.1, cex=0.3, pch=19, col='red')
            }

            # add the lines
            iCol <- iCol+1
            abline(v=bins0, col=iCol, lwd=iCol*.5)
            points(x=which(gridHeights>0), y=yMin + gridHeights[gridHeights>0] * 0.1, cex=0.3, pch=19, col='red')

            BBmisc::pause()
          }

        }
      }
    }

    return(list(maxShift=maxSumShift,maxGridCoord=maxGridCoord, maxPeriod=maxSumPeriod,maxSum=maxSum))
  }

  # Translate the best digital grid result (maxSum, maxShift, maxPeriod) to a list of labeled peaks
  translateBinsToPeaks <- function(){
    if(TRUE){
      bins0 <- maxShift+seq(0, maxGridCoord - maxShift, maxPeriod)
      digtGrid <- gridHeights*0
      digtGrid[bins0] <- 1
      sumD <- sum(digtGrid*gridHeights)
      print(table(digtGrid*gridHeights))
      dPeaks <- which(digtGrid*gridHeights >= dPeaksCutoff)  # TODO: dPeaksCutoff-what is the correct value? minGridHeight? will depend on how gridHeights is calculated
      wd1 <- which(digtGrid==1)
      digtGrid[wd1] <- digtGrid[wd1]+seq(1,length(wd1))-1
      nCopyPeaks_dig <- digtGrid[dPeaks]

      maxDigPeakOnGrid  <-  which.max(gridHeights[dPeaks]); # Note: this isn't guaranteed to be the "main" peak but it will be the biggest digital peak, depends on on the grid falls on the red dots see 28040
      nCopyMaxDigPeakOnGrid  <-  nCopyPeaks_dig[maxDigPeakOnGrid]
    }
    return(list(dPeaks=dPeaks,nCopyPeaks_dig=nCopyPeaks_dig,nCopyMaxDigPeakOnGrid=nCopyMaxDigPeakOnGrid ))
  }

  # Plot digital grid results
  plotDigitalGrid <- function(omitAnnotations){
    if(TRUE){
      yMin <- -bonusFactor/10
      plot(x=gridModes,
           y=peakHeights,
           main="",
           xlim=c(min(gridModes)-30, max(gridModes)+30),
           ylim=c(yMin, 1), cex=0.65, xlab='grid coordinate',
           yaxt="n",
           pch=19)
      if(!omitAnnotations){
        title(main=paste('Digital grid iteration:',gridIteration) )
        if(!is.null(sampleId)){mtext(3, text=sampleId,adj=0)}
        if(!is.null(alternateId)){mtext(3, text=alternateId,adj=1)}
      }
      axis(2, at=c(0,.5, 1), labels= c(0,.5, 1))
      abline(h=0.01,col="black")

      # annotate the grid dots range of values
      maxGridMode <- max(gridHeights[gridHeights>0])
      abline(h= (yMin +maxGridMode*.1),col="red")
      abline(h=yMin,col="red")
      mtext(side=4, at=(yMin + c(0,maxGridMode*.1)),text = c(0,maxGridMode), col='red', las=1, line=0.5)

      points(x=which(gridHeights>0), y=yMin + gridHeights[gridHeights>0] * 0.1, cex=0.3, pch=19, col='red')
      for(igrd in 0:100) {
        lines(rep(maxShift+igrd*maxPeriod,2),  c(yMin,1),col="grey")       # final gridlines # abline(v=maxShift+c(0:100)*maxPeriod,  col="grey")
      }
      for(igrd in -100:100) {
        lines(rep(maxShift+igrd*maxPeriod/2,2),c(yMin,(yMin+.02)),col="grey")   # halfway between final gridlines
      }

      if(!omitAnnotations){
        text(dPeaks,rep((yMin-.02),length(nCopyPeaks_dig)),labels=nCopyPeaks_dig)
        mtext(3, text=paste('penaltyCoefForAddingGrids:',penaltyCoefForAddingGrids),adj=0, line=-1)
        mtext(3, text=paste0('numOfGridCoordsToTest: ',numOfGridCoordsToTest),adj=0, line=-2)
        mtext(3, text=paste('minPeriod:',minPeriod),adj=0, line=-3)
        mtext(3, text=paste('period:',maxPeriod),adj=0, line=-4)

        legend('topright', legend=c('gridHeights>0', 'Peaks'), col=c('red', 'black'), pch = 19, pt.cex=c(0.3, 0.5))
      }
    }
  }


  # check for multiple peaks, needed to be able to do these calculates
  if(numPeaks > 1 ){
    # NOTE: peak 1 has a width of 7 while the other peaks have a width of 9
    firstToSecPeakDis   <-  (gridModes[2]   ) - (gridModes[1]   )-1  # CENTER of secondPeak to CENTER of firstPeak
    firstToSecPeakSpan  <-  (gridModes[2] +4) - (gridModes[1] -3)+1  # end of secondPeak to start of firstPeak  aka the span between the two peaks
    firstToSecPeakGap   <-  (gridModes[2] -4) - (gridModes[1] +3)-1  # start of secondPeak to end of firstPeak  aka the gap between the two peaks

    SecToThirdPeakGap   <-  (gridModes[3] -4) - (gridModes[2] +3)-1  # start of thirdPeak to end of secondPeak  aka the gap between the two peaks

    firstToThirdPeakGap <-  (gridModes[3] -4) - (gridModes[1] +3)  # start of thirdPeak to end of firstPeak  aka the gap between the two peaks
    firstToForthPeakGap <-  (gridModes[4] -4) - (gridModes[1] +3)  # start of forthPeak to end of firstPeak  aka the gap between the two peaks
    firstToFifthPeakGap <-  (gridModes[5] -4) - (gridModes[1] +3)  # start of fifthPeak to end of firstPeak  aka the gap between the two peaks
    twoBiggestPeaksGap  <- abs( (gridModes[which(peakHeightRanks==1)]) - (gridModes[which(peakHeightRanks==2)]) )-8

    if(mainPeakIndex!=1 & mainPeakIndex!=2){
      firstToMainPeakGap <- (gridModes[mainPeakIndex] -4) - (gridModes[1] +3)  # end of firstPeak  to start of mainPeak
      firstToMainPeakSpan <- (gridModes[mainPeakIndex] +4) - (gridModes[1] -3) # start of firstPeak  to end of mainPeak
    }else{
      firstToMainPeakGap <- NA
    }

    # firstToSecPeakDis: CENTER of secondPeak to CENTER of firstPeak
    # minPeriod: 10 is the minimum allowed because there are 9 possible gridHeights per peak
    #            set to firstToSecPeakGap, but is this big enough
    #            set to firstToSecPeakDis, but this is not small enough
    #            when firstToSecPeakDis is less than 5, use different criteria for minPeriod
    #                                              5 is a bit of a heuristic, see 28034,
    # if firstToSecPeakDis is super small, 1) high tumor and 2nd peak is a subclone that we should skip (58092, 58147, 58150)
    #                                      2) low tumor and it is hard to separate the peaks, going to struggle no matter what (40005)
    #                                      3) are there any other reasons?

  }else{
    firstToSecPeakDis <- 20 # doesn't really matter
    firstToSecPeakGap <- 1  # doesn't really matter
  }

  if(firstToSecPeakGap < 0){               # 1 means they are next to each other but not overlapping
    numGridCoordsFirstpeak <- 7 + firstToSecPeakGap
  }else{
    numGridCoordsFirstpeak <- 7
  }


  absMinPeriod <- 9
  # absMinPeriod=numGridCoordsFirstpeak+1 # this is not good becasue then you can get two grids per peak in the later peaks

  #' decision tree for maxPeriodToTry, do not allow period to get bigger than this, set to firstToMainPeakGap when: 5 <= mainPeakIndex < 7
  setMaxPeriodToTry <- function(mainPeakIndex, maxGridCoord,  minGridCoord, maxPeriodManual){
    if(maxPeriodManual>0){
      maxPeriodToTry <- maxPeriodManual
    }else{
      if(mainPeakIndex >= 4){
        maxPeriodToTry <- firstToMainPeakSpan
      }else{
        maxPeriodToTry <- maxGridCoord  - minGridCoord
      }
      loginfo('do digital grid with maxPeriodToTry:%s',maxPeriodToTry)
    }


    return(maxPeriodToTry)
  }


  checkThePeakGaps <- function(previousPeriod,firstToSecPeakGap,firstToThirdPeakGap,firstToForthPeakGap,firstToFifthPeakGap){
    if(previousPeriod < firstToSecPeakGap){ # this may be true if main peak is peak 4-6, which has a special minPeriod setting
      minPeriod  <-  firstToSecPeakGap
      extendToPeak <- 2
    }else if(previousPeriod < firstToThirdPeakGap){
      minPeriod  <-  firstToThirdPeakGap
      extendToPeak <- 3
    }else if(!is.na(firstToForthPeakGap)){
      if(previousPeriod < firstToForthPeakGap){
        minPeriod  <-  firstToForthPeakGap
        extendToPeak <- 4
      }else if(!is.na(firstToFifthPeakGap)){
        if(previousPeriod < firstToFifthPeakGap){
          minPeriod  <-  firstToFifthPeakGap
          extendToPeak <- 5
        }else{
          logwarn('need to decide what to do now,  previousPeriod > %s (firstToFifthPeakGap)', firstToFifthPeakGap)
          return(NA)
        }

      }else{
        logwarn('need to decide what to do now, firstToFifthPeakGap is %s', firstToFifthPeakGap)
        return(NA)
      }

    }else{
      logwarn('need to decide what to do now, firstToForthPeakGap is %s', firstToForthPeakGap)
      return(NA)
    }
    loginfo('redo digital grid with (bigger) minPeriod:%i, extending to peak %i',minPeriod,extendToPeak)
    return(minPeriod)
  }

  #' decision tree for minPeriod
  #' @param absMinPeriod the smallest minPeriod allowed
  #'
  #' @inheritParams commonParameters
  setMinPeriod <- function(previousPeriod,minPeriodManual,absMinPeriod,firstToSecPeakGap,firstToThirdPeakGap,firstToForthPeakGap,firstToFifthPeakGap){
    if(minPeriodManual < 0){

      if(gridIteration==1){

        #default value
        minPeriod <- max(absMinPeriod, firstToSecPeakGap)

        # special situations to consider
        if(firstToSecPeakDis <= 10 ){
          # use the default minPeriod, and skip the other checks.


          # force bigger period, skip second peak even if first red dot of first peak is used

          # TODO consider what happens if the first peak is skipped, then trying to skip the second peak is not relevant,
          #      but how do you know if the first peak will be skipped or not? ie PT58251

          # consider using firstToSecPeakSpan only if first peak is use, but smaller if other peaks are used
          # or only open up minPeriod if the numOfGridCoordsToTest is limited to the first peak
          # or ...
          # proposedNumOfGridCoordsToTest= setNumOfGridCoordsToTest(previousPeriod) # calling it "proposed" because this is determined in another step
          # if(proposedNumOfGridCoordsToTest>=numGridCoordsFirstpeak){
          #   minPeriod <- firstToSecPeakSpan
          # }
        }else if( mainPeakIndex>=4 &
                  mainPeakIndex< 7){ # aka 4, 5 or 6;

          # there are lots of early peaks, allow the first peak to be a sub-clone
          # i.e. PT58158 first peak is 1N but main peak is the 6th peak so don't open it too far
          minPeriod <- max(absMinPeriod*2, min(minPeriod,SecToThirdPeakGap))
          # make sure minPeriod isn't smaller than absMinPeriod plus a little bit
          # absMinPeriod does not work for 64648 because the third peak is a subClone and the distance (11) is too small
          # this is required for 43027 because the first peak(1N) is shifted from all the others for some reason and need to use the 2N to 3N peak distance
          # would only want to do this if the first peak is a subClone and too far away from the second peak ie:

        }else if( mainPeakIndex>=7){
          # use the default minPeriod
        }
        loginfo('do digital grid with minPeriod:%i',minPeriod)
      }else if(gridIteration==2){
        if(numOfGridCoordsToTest == numGridCoordsFirstpeak){
          # be conservative, will increase just the numOfGridCoordsToTest first, leave minPeriod the same
          minPeriod  <-  previousPeriod                            # keep the minPeriod the same as previous iteration
        }else{
          minPeriod <- checkThePeakGaps(previousPeriod,firstToSecPeakGap,firstToThirdPeakGap,firstToForthPeakGap,firstToFifthPeakGap)
        }
      }else{  # this will be a x.1+ or x.01+ iteration implemented when nCopyMaxDigPeakOnGrid >= 5 or !isMainPeakDigital
        if(numOfGridCoordsToTest <= numGridCoordsFirstpeak){
          # be conservative, will increase just the numOfGridCoordsToTest first, leave minPeriod the same
          minPeriod  <-  previousPeriod                            # keep the minPeriod the same as previous iteration
        }else{
          minPeriod <- checkThePeakGaps(previousPeriod,firstToSecPeakGap,firstToThirdPeakGap,firstToForthPeakGap,firstToFifthPeakGap) # <---- this is where it is failing because it doesn't see
        }
      }
    }else{
      minPeriod <- minPeriodManual
    }

    return(minPeriod)
  }


  setNumOfGridCoordsToTest <- function(previousPeriod){
    if(gridIteration==1){
      # default value
      numOfGridCoordsToTest <-   numGridCoordsFirstpeak  # default, test only grids in the first peak, which has 7 (or fewer if there is overlap) positions, other peaks have 9
      # special situations to consider
      if(firstToSecPeakDis <= 10 ){
        numOfGridCoordsToTest <- length(nonZeroGridCoords)    # do not force first peak to be the first digital peak, allow other peaks to be the first digital peak
      }else if( mainPeakIndex>=4 &
                mainPeakIndex< 7){ # aka 4, 5 or 6

        # there are lots of early peaks, allow the first peak to be a sub-clone
        # PT58158 first peak is 1N but main peak is the 6th peak so don't open it too far
        numOfGridCoordsToTest <- numGridCoordsFirstpeak+9    # do not force first peak to be the first digital peak, allow second peak to be the first digital peak


      }else if( mainPeakIndex>=7){
        # there are lots of early peaks, allow the first peak to be a sub-clone
        # PT58158 first peak is 1N but main peak is the 6th peak so don't open it too far
        numOfGridCoordsToTest <-   numGridCoordsFirstpeak+9+9    # do not force first peak to be the first digital peak, allow other peaks to be the first digital peak
      }
      loginfo(' do digital grid with numOfGridCoordsToTest: %i', numOfGridCoordsToTest)

    }else if(gridIteration==2){
      if(numOfGridCoordsToTest == numGridCoordsFirstpeak){
        numOfGridCoordsToTest <- length(nonZeroGridCoords)    # do not force first peak to be the first digital peak, allow other peaks to be the first digital peak
        loginfo(' redo digital grid with (more) numOfGridCoordsToTest: %i', numOfGridCoordsToTest)
      }
    }else{  # this will be a x.1+ or x.01+ iteration implemented when nCopyMaxDigPeakOnGrid >= 5 or !isMainPeakDigital
      # be conservative, try increasing just the numOfGridCoordsToTest first, leave minPeriod the same
      if(numOfGridCoordsToTest == numGridCoordsFirstpeak){
        numOfGridCoordsToTest <- numGridCoordsFirstpeak+9    # do not force first peak to be the first digital peak, allow second peak to be the first digital peak
      }else{
        # max out numOfGridCoordsToTest and bump up minPeriod
        numOfGridCoordsToTest <- length(nonZeroGridCoords)    # do not force first peak to be the first digital peak, allow other peaks to be the first digital peak
      }
      loginfo(' redo digital grid with (more) numOfGridCoordsToTest: %i', numOfGridCoordsToTest)

    }
    return(numOfGridCoordsToTest)
  }

  # not dependent on iteration
  maxPeriodToTry <- setMaxPeriodToTry(mainPeakIndex,maxGridCoord, minGridCoord,maxPeriodManual=maxPeriodManual)


  # first fit
  if(gridIteration==1){
    #1) setMinPeriod uses numOfGridCoordsToTest for some of its logic so find this first
    minPeriod <- setMinPeriod(previousPeriod=NULL, minPeriodManual=minPeriodManual,absMinPeriod=absMinPeriod, firstToSecPeakGap = firstToSecPeakGap) # find this before
    #2)
    numOfGridCoordsToTest <- setNumOfGridCoordsToTest(previousPeriod=NULL)

    out <-  fitGridPeriod(gridHeights,numOfGridCoordsToTest, minPeriod, maxPeriodToTry=maxPeriodToTry,plotDigDebug = FALSE)
    maxShift <- out$maxShift; maxGridCoord <- out$maxGridCoord; maxPeriod <- out$maxPeriod; maxSum <- out$maxSum

    outPeaks  <-  translateBinsToPeaks()
    dPeaks <- outPeaks$dPeaks;   nCopyPeaks_dig <- outPeaks$nCopyPeaks_dig;  nCopyMaxDigPeakOnGrid <- outPeaks$nCopyMaxDigPeakOnGrid

    iterationStats <- plyr::rbind.fill(iterationStats,
                                       data.frame(iteration=gridIteration, numOfGridCoordsToTest=numOfGridCoordsToTest,
                                                  minPeriod=minPeriod,period=maxPeriod, numPeaks=numPeaks,numDigPeaks=length(nCopyPeaks_dig),
                                                  gridScore=round(maxSum/length(nCopyPeaks_dig),3)))
    # if(!skipExtras){
    plotDigitalGrid(omitAnnotations=omitAnnotations)
    # }
  }


  # next fit
  # first bump up numOfGridCoordsToTest, next time around bump up minPeriod
  if(gridIteration >= 2){
    #1) setMinPeriod uses numOfGridCoordsToTest for some of its logic so find this first
    minPeriod <- setMinPeriod(previousPeriod=previousPeriod,minPeriodManual=minPeriodManual,absMinPeriod=absMinPeriod,firstToSecPeakGap,firstToThirdPeakGap,firstToForthPeakGap,firstToFifthPeakGap)

    # no new solution is possible, exit and use previous iteration
    if(is.na(minPeriod)){
      return(list(nCopyPeaks_dig=NA)) # so the test on the other side does spew a warning message
    }

    #2)
    numOfGridCoordsToTest <- setNumOfGridCoordsToTest(previousPeriod=previousPeriod)

    out <-  fitGridPeriod(gridHeights,numOfGridCoordsToTest, minPeriod)
    maxShift <- out$maxShift; maxGridCoord <- out$maxGridCoord; maxPeriod <- out$maxPeriod; maxSum <- out$maxSum

    outPeaks  <-  translateBinsToPeaks()
    dPeaks <- outPeaks$dPeaks;   nCopyPeaks_dig <- outPeaks$nCopyPeaks_dig;  nCopyMaxDigPeakOnGrid <- outPeaks$nCopyMaxDigPeakOnGrid

    iterationStats <- plyr::rbind.fill(iterationStats,
                                       data.frame(iteration=gridIteration, numOfGridCoordsToTest=numOfGridCoordsToTest,
                                                  minPeriod=minPeriod,period=maxPeriod, numPeaks=numPeaks, numDigPeaks=length(nCopyPeaks_dig),
                                                  gridScore=round(maxSum/length(nCopyPeaks_dig),3)))
    plotDigitalGrid(omitAnnotations=omitAnnotations)
  }




  ### test: if max digital grid height has cn >= 5, redo fit-----
  ### usually this will be the main peak, but not necessarily, depends on how the grid falls on the peaks: if the grid overlaps a dot on the 2nd largest peak that is higher than a dot on the main peak i.e. 28040
  if(nCopyMaxDigPeakOnGrid >= 5 ){
    gridIteration <- gridIteration + 0.1

    if(TRUE){
      loginfo('nCopyMaxDigPeakOnGrid=%i redo digital grid',nCopyMaxDigPeakOnGrid)
      minPeriod <- setMinPeriod(previousPeriod=maxPeriod,minPeriodManual=minPeriodManual,absMinPeriod=absMinPeriod,firstToSecPeakGap,firstToThirdPeakGap,firstToForthPeakGap,firstToFifthPeakGap)
      numOfGridCoordsToTest <- setNumOfGridCoordsToTest(previousPeriod=maxPeriod)
    }else{
      # be conservative, try increasing just the numOfGridCoordsToTest first, leave minPeriod the same
      if(numOfGridCoordsToTest == numGridCoordsFirstpeak){
        numOfGridCoordsToTest <- numGridCoordsFirstpeak+9    # do not force first peak to be the first digital peak, allow other peaks to be the first digital peak
        loginfo('nCopyMaxDigPeakOnGrid=%i redo digital grid with (more) numOfGridCoordsToTest: %i',nCopyMaxDigPeakOnGrid,numOfGridCoordsToTest)
      }else{
        # max out numOfGridCoordsToTest and bump up minPeriod
        numOfGridCoordsToTest <- length(nonZeroGridCoords)    # do not force first peak to be the first digital peak, allow other peaks to be the first digital peak

        if(maxPeriod < firstToThirdPeakGap){
          minPeriod  <-  firstToThirdPeakGap
        }else if(maxPeriod < firstToForthPeakGap){
          minPeriod  <-  firstToForthPeakGap
        }else{
          stop(sprintf('need to decide what to do now, previous period (maxPeriod) > %s firstToForthPeakGap ',firstToForthPeakGap))
        }
        loginfo('nCopyMaxDigPeakOnGrid=%i redo digital grid with (bigger) minPeriod: %i and (more) numOfGridCoordsToTest: %i',nCopyMaxDigPeakOnGrid,minPeriod,numOfGridCoordsToTest)
      }
    }

    # another fit
    out <-  fitGridPeriod(gridHeights,numOfGridCoordsToTest, minPeriod)
    maxShift <- out$maxShift; maxGridCoord <- out$maxGridCoord; maxPeriod <- out$maxPeriod; maxSum <- out$maxSum

    outPeaks  <-  translateBinsToPeaks()
    dPeaks <- outPeaks$dPeaks;   nCopyPeaks_dig <- outPeaks$nCopyPeaks_dig;  nCopyMaxDigPeakOnGrid <- outPeaks$nCopyMaxDigPeakOnGrid

    iterationStats <- plyr::rbind.fill(iterationStats,
                                       data.frame(iteration=gridIteration, numOfGridCoordsToTest=numOfGridCoordsToTest,
                                                  minPeriod=minPeriod,period=maxPeriod, numPeaks=numPeaks, numDigPeaks=length(nCopyPeaks_dig),
                                                  gridScore=round(maxSum/length(nCopyPeaks_dig),3)))
    plotDigitalGrid(omitAnnotations=omitAnnotations)
  }


  ### test: if main peak is not on digital grid, redo fit -----
  mainPeakGridMarkers <- (which.max(gridHeights)-4):(which.max(gridHeights)+4)
  if(any(dPeaks %in% mainPeakGridMarkers)){
    isMainPeakDigital <- TRUE
  }else{
    isMainPeakDigital <- FALSE
    logwarn('main peak is not a digital peak')
  }

  if(!isMainPeakDigital){
    gridIteration <- gridIteration + 0.01

    if(TRUE){
      # TODO: should minPeriod be adjusted?
      #minPeriod = setMinPeriod(previousPeriod=maxPeriod)
      numOfGridCoordsToTest <- setNumOfGridCoordsToTest(previousPeriod=maxPeriod)
    }else{
      numOfGridCoordsToTest <- length(nonZeroGridCoords)    # do not force first peak to be the first digital peak, allow other peaks to be the first digital peak
    }

    loginfo('isMainPeakDigital=%s redo digital grid with (more) numOfGridCoordsToTest',isMainPeakDigital)

    # another fit
    out <-  fitGridPeriod(gridHeights,numOfGridCoordsToTest, minPeriod)
    maxShift <- out$maxShift; maxGridCoord <- out$maxGridCoord; maxPeriod <- out$maxPeriod; maxSum <- out$maxSum

    outPeaks  <-  translateBinsToPeaks()
    dPeaks <- outPeaks$dPeaks;   nCopyPeaks_dig <- outPeaks$nCopyPeaks_dig;  nCopyMaxDigPeakOnGrid <- outPeaks$nCopyMaxDigPeakOnGrid

    iterationStats <- plyr::rbind.fill(iterationStats,
                                       data.frame(iteration=gridIteration, numOfGridCoordsToTest=numOfGridCoordsToTest,
                                                  minPeriod=minPeriod,period=maxPeriod, numPeaks=numPeaks, numDigPeaks=length(nCopyPeaks_dig),
                                                  gridScore=round(maxSum/length(nCopyPeaks_dig),3)))
    plotDigitalGrid(omitAnnotations = omitAnnotations)
    gridHeights[gridModes]
    gridHeights[dPeaks]
  }

  return(list(dPeaks=dPeaks,
              nCopyPeaks_dig=nCopyPeaks_dig,
              period=maxPeriod,
              maxSum=maxSum,
              gridIteration=gridIteration,
              numOfGridCoordsToTest=numOfGridCoordsToTest,
              iterationStats=iterationStats))

}

#' Determines ploidy of a given sample
#'
#' Fit the peaks from the read depth distributions to a digital grid.
#' Evaluate the heterozygosity score to determine if first digital peak is 1N or 2N.
#' Then find the expected number of reads in the 2N peak and normalize that value to one bp. Tumor percent is calcuated from the two biggest digital peaks.
#'
#' @param centroArray array with the positions of the centromeres for each chromosome
#' @param hetScoreData heterozygosity scores determined per 30 kb bin over a 1 Mb region
#' @param numChroms number of chromosomes in the reference genome to consider
#' @param minGridHeight minimum value that can be assigned to the gridHeights
#' @param grabDataPercentManual portion of main peak data to grab, other peaks will be scaled based on read depth (x location), set to -1 to base off of mainPeak width
#' @param origMaxPercentCutoffManual peaks smaller than this portion of the max peak are not considered; set to -1 to use default value
#' @param pause pause execution until user prompts to continue, available interactively only, useful during testing
#' @param skipExtras logical to turn on/off plots used for testing and debugging
#' @param minPeriodManual manually set \code{minPeriod} within \code{calculatePloidy}
#' @param maxPeriodManual manually set \code{maxPeriod} within \code{calculatePloidy}
#' @param forceFirstDigPeakCopyNum value to force copy number of first digital peak, use only when ploidy calculation is wrong
#' @param minReasonableSegmentSize initial smallest segment size to include in ploidy test segments; want to keep as large as possible to avoid 0N segments, but will decrease size if not enough segments are found
#'
#' @inheritParams commonParameters
#'
#' @example inst/examples/calculatePloidyExample.R
#'
#' @return expReadsIn2NPeak_1bp, percentTumor, peakInfo, hetScoreQuantiles
#'
calculatePloidy <- function(sampleId, outputDir,alternateId=NULL,
                            readDepthPer30kbBin=NULL, readDepthPer100kbBin=NULL,
                            segmentation, segmentationBinSize=30000,
                            hetScoreData,
                            numChroms=24,centroArray,
                            dPeaksCutoff=0.01,    penaltyCoefForAddingGrids=0.49, minGridHeight=0.2, minPeriodManual=-1, maxPeriodManual=-1,    # digital peaks
                            grabDataPercentManual= -1,  origMaxPercentCutoffManual=-1,  #  peaksByDensity
                            pause=FALSE, noPdf=FALSE,skipExtras=FALSE,forceFirstDigPeakCopyNum=-1,
                            minReasonableSegmentSize=5.5e6,
                            omitAnnotations = FALSE,
                            heterozygosityScoreThreshold=0.98,  # If segment hetScore is more than this, the segment is heterozygous
                            allowedTumorPercent = 106,
                            hsNormMat=NULL
){
  ### defaults
  # dPeaksCutoff=0.01; penaltyCoefForAddingGrids=0.49; minGridHeight=0.2; minPeriodManual=-1;maxPeriodManual=-1;grabDataPercentManual= -1; origMaxPercentCutoffManual=-1; pause=FALSE; skipExtras=FALSE; heterozygosityScoreThreshold=0.98

  # TODO: consider moving this out of function, as part of set up etc.
  if (!noPdf) {
    ploidyPdfFile <- getTypedFile('ploidyReportPdf', dir=outputDir, values=list(sampleId=sampleId))
    # TODO: check that this directory is created by the pipeline
    if(!dir.exists(dirname(ploidyPdfFile@path))) {
      loginfo("will create output directory: %s", dirname(ploidyPdfFile@path))
      dir.create(path=file.path(dirname(ploidyPdfFile@path)),mode = "0775")
    }
    pdf(file = ploidyPdfFile@path, paper="a4r", width=8, height=10, title=paste0('Ploidy_',sampleId))
    on.exit(dev.off(),add = TRUE)
  }


  xind <-23 # index of chrX
  coords <- getLinearCoordinates(chromosomes = 1:numChroms)
  maxcn <- numChroms


  ################################'
  ### peaks By density------
  ################################'
  #   normalized so x values represent read depth at 1 bp, peak magnitudes are normalized so the max peak = 1

  ### Use density on read depth to determine peaks, peaksByDensity() ----------------------------
  # peaks are returned in order of read depth with most common (highest frequency) read depth peak listed first
  qualityPostNorm <- NULL

  resultPBD <- peaksByDensity(sampleId, readDepthPer100kbBin=readDepthPer100kbBin,
                              segmentation=segmentation, segmentationBinSize=segmentationBinSize,
                              grabDataPercentManual= grabDataPercentManual, origMaxPercentCutoffManual=origMaxPercentCutoffManual,
                              addAreaLinesToPlot=!omitAnnotations,qualityPostNorm=qualityPostNorm,pause=FALSE, omitAnnotations=omitAnnotations,
                              alternateId=alternateId)

  #
  ####################'
  # create peakInfo
  if(TRUE){
    numPeaks  <-  length(resultPBD[[1]])
    # ranked By Freq/Magnitude/height/Size of each peak - bin height
    peakReadDepthList_per1bp_rankByHeight <- resultPBD$peakReadDepthList_per1bp   # x  very small changes are not significant and can mess with the ordering
    peakHeightList_rankByHeight <- resultPBD$peakHeightList_maxIsOne              # y
    reorderIndex <- order(peakReadDepthList_per1bp_rankByHeight)

    cnColors <- getCNcolors()

    peakInfoTemp <- data.frame(peakCol       = cnColors[1:numPeaks],
                               scaledGrabDataPercent = round(resultPBD$scaledGrabDataPercentPerPeak[, 'scaledGrabDataPercent'],4),
                               peakReadDepth_1bp = peakReadDepthList_per1bp_rankByHeight,         # normalized so x values represent readDepth at 1 bp
                               peakHeight    = peakHeightList_rankByHeight,     # normalized so biggest peak = 1
                               rankByHeight  = 1:length(peakReadDepthList_per1bp_rankByHeight)
    )

    ## peak positions and heights are ranked by x position (read depth) - for digital grid purposes
    peakInfo <- peakInfoTemp[reorderIndex,]
    peakInfo[,'peakReadDepth_normX'] <- (peakInfo$peakReadDepth_1bp/min(peakInfo$peakReadDepth_1bp)) # normalize so peak with left-most x peak position sits at position x=1
    peakInfo[,'rankByX'] <- 1:nrow(peakInfo)

    peakInfo
  }

  ## assign bonuses to the peaks -------------------
  if(TRUE){
    ## tallest peak gets the biggest bonus, smallest peak gets the smallest bonus
    ## these values must be in the same order as gridModes,   'rankByX'

    ## multiple options for generating the bonus values:
    ## 1) bonusLinear..... range: max peak = 2, min peak = 1, equal intervals inbetween
    if(numPeaks==1){
      bonusLinear <- 1
    }else{
      bonusLinear <- rev(1 + seqFwd(0, numPeaks - 1) / (numPeaks - 1)) # bonusLinear go from 2 to 1, most prominent peak should have the biggest bonus/ highest dot
    }
    bonusLinear <- bonusLinear[reorderIndex]

    ## 2)  normalized to the max peak which is 1
    bonusFactor <- 8  # heuristic, if it changes how should penalty change?
    bonusRelative <- peakInfo$peakHeight*bonusFactor

    ## 3)  both options
    bonusOther <- 1 # bonusRelative* bonusLinear

    # now add this column with the value that was choosen
    bonus <- bonusRelative
    peakInfo$bonus <- round(bonus, 4)
    print(peakInfo)
  }


  ### load data ----------
  if(TRUE){
    wsz <- readDepthPer30kbBin$windowSize

    ### get frequency array ---------------
    frq00 <- readDepthPer30kbBin$readDepthArray       # y axis
    wdnsMSK00 <- readDepthPer30kbBin$goodWindowArray      # x axis

    getSegmentsWeWant <- function(segmentation,minReasonableSegmentSize, xind){
      segmentsWeWant <- which( segmentation[,'chr'] < xind &
                                 (segmentation[,'end'] - segmentation[,'start']) > minReasonableSegmentSize
                               # & segmentation[,'expectedCount'] > 0  # expectedCount should not be in the file, and this should already be taken into account when loading the segmentation file
      )
      return(segmentsWeWant)
    }


    ###  grab segment data from segmentation
    getAndMakeSegmentData <- function(segmentation, minReasonableSegmentSize, xind=23){
      # minReasonableSegmentSize=5.5e6
      # what is the correct minReasonableSegmentSize? choosing 5.5e6  helps avoid 0N segments, and thus avoids a 0N peak.

      # keys=which(segmentation$chr==10); cbind(segmentation[keys,],size=segmentation[keys,'end']-segmentation[keys,'start'])
      # segmentsWeWant <- NULL; minReasonableSegmentSize <- 5.5e6
      minReasonableSegmentSizeStart <- minReasonableSegmentSize



      ### Get only segments that are longer than minReasonableSegmentSize
      segmentsWeWant <- getSegmentsWeWant(segmentation,minReasonableSegmentSize,xind)


      # TODO V1: remove/skip the centromeres and genomeGaps (heterochromatin etc.), they have noisy cnv and loh...might be skipped in wdnsMSK00, but not in segmentation
      # gapArray    <- getGenomeGapPositions(ideogram=exampleIdeogram() )
      # centroArray <- getCentromerePositions(ideogram = ideogram )

      # TODO: what is an appropriate min number of segments?  just enough so it doesn't crash?
      #       ie. PT626 has only 48 segments (none of which are in the first digital peak) for minReasonableSegmentSize <- 5.5e6, but 65 at minReasonableSegmentSize <- 5.0e6

      minNumberOfSegments <- 109  # 5-24-22 reduced from 120 to 109
      while(length(segmentsWeWant) < minNumberOfSegments ){
        minReasonableSegmentSize <- minReasonableSegmentSize-0.5e6
        if(minReasonableSegmentSize < 2e6){
          logerror('not enough segments (%i) and minReasonableSegmentSize is too small: %s',length(segmentsWeWant), minReasonableSegmentSize)
          break
        }
        loginfo('only %i fragments, reducing minReasonableSegmentSize to %s',length(segmentsWeWant), minReasonableSegmentSize)
        segmentsWeWant <- getSegmentsWeWant(segmentation,minReasonableSegmentSize,xind=xind)
      }
      minReasonableSegmentSizeFinal <- minReasonableSegmentSize

      # segment further into smaller chucks? ----------------
      # this step is important for providing enough data to create a density score
      # also hetScores may vary within a segment, and don't want to average various levels together, want distinct values.
      segmentFurther <- TRUE
      if(!segmentFurther){
        segmentData <- data.frame(
          segmentation[segmentsWeWant, c('chr', 'start', 'end')],
          lohScoreMean=NA,
          lohScoreMedian=NA,
          pkmod=NA)
        segmentData[,'size'] <- segmentData$end - segmentData$start
        segmentData[,'rd'] <- segmentation[segmentsWeWant,'rd']

        loginfo(' %i segments of size %i Mb or bigger.', nrow(segmentData),minReasonableSegmentSize/1000000)

      }else{
        # split up the really big segments into smaller sections
        chunkSize <- 3e6  # TODO: the chunk size could be a lot smaller than minReasonableSegmentSizeFinal, like 1MB?
        #       this is for getting a mean/median hetScore across the CNV segment so why not, it would provide more data for the hetScore density which would be great!

        segmentTemp <- data.frame( segmentation[segmentsWeWant,  c('chr', 'start', 'end')])
        segmentTemp[,'size'] <- segmentTemp$end - segmentTemp$start
        segmentTemp[,'rd'] <- segmentation[segmentsWeWant,'rd']

        segmentData <- data.frame()
        for(jRow in 1:nrow(segmentTemp)){
          # jRow=2
          if(segmentTemp[jRow, 'size'] / chunkSize > 2){
            newRows <- segmentTemp[FALSE,] # initialize
            numOfNewSegs <- floor(segmentTemp[jRow, 'size'] / chunkSize)

            startStops <- seq(from=segmentTemp[jRow,'start'], to= segmentTemp[jRow,'end'], len=numOfNewSegs+1)

            for(k in 1:numOfNewSegs){
              newRows[k,'start'] <- startStops[k]
              newRows[k,'end'] <- startStops[k+1]
            }
            newRows[,'chr'] <- segmentTemp[jRow,'chr']
            newRows[,'size'] <- newRows$end - newRows$start
            newRows[,'rd'] <- segmentTemp[jRow,'rd'] # TODO: make sure this is ok to just copy to the new segments

          }else{
            newRows <- segmentTemp[jRow,]
          }
          segmentData <- rbind(segmentData,newRows)
        }

        segmentData[,'size'] <- segmentData$end - segmentData$start
        loginfo('segmentData was broken into smaller chunks: from %i to %i segments', nrow(segmentTemp), nrow(segmentData) )
        if(minReasonableSegmentSizeStart!=minReasonableSegmentSizeFinal){
          loginfo('start and final segment size min: %i to  %i ', minReasonableSegmentSizeStart, minReasonableSegmentSizeFinal )
        }

      }
      return(list(segmentData=segmentData, minReasonableSegmentSizeFinal=minReasonableSegmentSizeFinal))
    }

    segmentDataReturned <- getAndMakeSegmentData(segmentation,minReasonableSegmentSize=minReasonableSegmentSize, xind=xind)
    segmentData <- segmentDataReturned$segmentData
    minReasonableSegmentSizeFinal <- segmentDataReturned$minReasonableSegmentSizeFinal

    minSegment <- round(min(segmentData$size)/1000000,1)
    maxSegment <- round(max(segmentData$size)/1000000,1)
    meanSegment <- round(mean(segmentData$size)/1000000,1)
    logdebug('segmentData from segmentation: min = %s Mb  max = %s Mb  mean = %s Mb', minSegment, maxSegment, meanSegment)

    segmentDataSizes=list(minReasonableSegmentSizeFinal=minReasonableSegmentSizeFinal,
                          minSegment=minSegment,
                          maxSegment=maxSegment,
                          meanSegment=meanSegment)
  }


  #' wrap steps to repeat into a function:
  #' add 'nCopy' to peaksInfo, construct hetScoreQuantiles, calculate Tumor
  repeatSteps <- function(peakInfo,keyHetScoresPerPeak,nCopyPeaks_step, allowed=allowedTumorPercent){

    ## 1 ## which peakRD_1bp (normPeak) is which of the dPeaks?
    peakInfo[,'nCopy'] <- NA_integer_
    for(i in 1:max(c(nCopyPeaks_step, length(nCopyPeaks_step)))){  # need length() in the event of a 0N level, ie 0,1,2
      key <- which(peakInfo[,'dPeaks']==dPeaks[i])
      peakInfo[key,'nCopy'] <- nCopyPeaks_step[i]
    }

    ## 2 ## HetScore quantile scores for each digital peak
    hetScoreQuantiles <-  data.frame(sampleId,
                                     alternateId=ifelse(is.null(alternateId),NA,alternateId),
                                     'nCopy'=peakInfo$nCopy,
                                     keyHetScoresPerPeak )

    ## 3 ## Calculate tumor percentage
    # can only calculate Tumor percentage if there is more than one nCopyPeaks_step
    numDigitalPeaks <- length(which(!is.na(peakInfo$dPeaks)))
    calcTumor <- ifelse( numDigitalPeaks > 1, TRUE, FALSE)
    if(calcTumor){
      percentTumor <- calcTumorFromPloidyPeaks(peakCopyNum = peakInfo$nCopy,
                                               peakHeight  = peakInfo$peakHeight,
                                               peakReadDepth_1bp=peakInfo$peakReadDepth_1bp,
                                               dPeaks = peakInfo$dPeaks)
      if(percentTumor > allowed){
        logwarn('percentTumor: %s more than %i algorithm error?', round(percentTumor,2),allowed)
      }else if(percentTumor <40){
        logwarn('low percentTumor: %s ', round(percentTumor,2))
      }
    }else{
      percentTumor <- NA
      loginfo('cannot calculate Tumor Percentage, only one nCopyPeak found')
    }

    return(list(peakInfo=peakInfo, percentTumor=percentTumor,hetScoreQuantiles=hetScoreQuantiles))
  }

  ## layout Grid coordinates (red dots) ------
  ## heights (bonus) and positions (at grid modes)
  n00 <- 100 # Precision with which the "digital grid" is created
  gridHeights <- assignGridHeights(peakInfo, n00, minGridHeight )
  nonZeroGridCoords <- which(gridHeights > 0)


  ################################'
  ## While loop-Digital grid ----------------------------------
  ################################'
  ### TODO:  needs to be more robust yet flexible, most problems occur when it can't detect subClones from actual CN levels, or the first peak is a subclone
  # penaltyCoefForAddingGrids TODO: what is the correct value?

  iterationStatsAll <- list(iterationStats=data.frame(),
                            digitalPeakZone=NULL)
  gridIteration <- 1; period <- NULL; numOfGridCoordsToTest <- NULL
  # pause <- FALSE
  while(gridIteration < numPeaks | gridIteration==1){  # if there is only one peak, need to run the while loop once

    loginfo('gridIteration: %.02f',gridIteration)

    digGridResult <- digitalGrid(peakInfo, gridHeights,
                                 dPeaksCutoff,penaltyCoefForAddingGrids,
                                 n00=n00, bonusFactor = bonusFactor, sampleId=sampleId, alternateId=alternateId,
                                 gridIteration=gridIteration, minPeriodManual=minPeriodManual, maxPeriodManual=maxPeriodManual,
                                 numOfGridCoordsToTest=numOfGridCoordsToTest,
                                 previousPeriod = period,minGridHeight=minGridHeight,
                                 omitAnnotations=omitAnnotations,
                                 iterationStats=iterationStatsAll$interationStats)
    if(is.na(digGridResult$nCopyPeaks_dig)[1]){
      gridIteration <- gridIteration+10000 # high value signal

      # use previous digitalGrid() result
      break
    }

    dPeaks <- digGridResult$dPeaks
    nCopyPeaks_dig <- digGridResult$nCopyPeaks_dig
    gridIteration <- digGridResult$gridIteration
    numOfGridCoordsToTest <-digGridResult$numOfGridCoordsToTest
    period <- digGridResult$period
    numDigPeaks     <- length(dPeaks)
    iterationStats <- digGridResult$iterationStats

    ### How far can a data point be (from first digital peak) to still count as part of that peak
    # TODO: what is the correct value?  hard-coded constant for all peaks? based on the spread of data? grabDataPercent (see below)?
    if(period <= 12){
      digitalPeakZone <- 0.02
    }else if(period <= 15){
      digitalPeakZone <- 0.035
    }else{
      digitalPeakZone <- 0.05
    }

    iterationStatsAll <- list(iterationStats=iterationStats,
                              digitalPeakZone=digitalPeakZone)

    ### Add digital peak results to peakInfo ----------------------------------
    # Now that we know how the peaks fall on a grid
    if(TRUE){

      ### match up the peaks in peakRD_1bp, peakInfo$peakReadDepth_normX (peakRD_normX) to the peaks in dPeaks
      # remaining NAs are subclonal or noise?
      peakInfo[,'dPeaks'] <- NA_integer_
      for(i in 1:length(dPeaks)){
        key <- which.min(abs(dPeaks[i] - (peakInfo$peakReadDepth_normX * n00)))
        peakInfo[key,'dPeaks'] <- dPeaks[i]
      }

      dPeak_rankByHeight <- peakInfo[!is.na(peakInfo$dPeaks), 'rankByHeight']

      mainPeakIndex <- which(peakInfo$rankByHeight==1) # out of all the peaks
      if(is.na(peakInfo[mainPeakIndex, 'dPeaks'])){
        logwarn('mainPeak is not a digital peak, this is a problem, look closer at the digital peak algorithm and cutoffs/coefficients')
      }

      firstDigPeakIndex <- which.min(peakInfo$dPeaks) # may not be 1 if the grid is not forced to start at the first peak
      expReadsInFirstDigitalPeak_wsz <- peakInfo[firstDigPeakIndex, 'peakReadDepth_1bp']* wsz   # peakRD_1bp is value for 1bp window

      ### number the digital peaks within the total peaks
      # NA=not a digital peak
      digitalPeakIndex <- rep(NA, numPeaks)
      increment <- 1
      for(i in 1:numPeaks){
        if(!is.na(peakInfo[i,'dPeaks'])){
          digitalPeakIndex[i] <- increment
          increment <- increment+1
        }
      }

      nCopyPeaks_final <- NULL # reinitialize - to be sure to avoid stupid mistakes during debugging
    } # CLOSE run this chunk



    ############################'
    ## by Segment ------
    ############################'

    bySegment <- TRUE
    if(bySegment){

      ### get segment het scores, pkmod, and plot
      # this cannot be moved out/above while loop, pkmod will depend on expReadsInFirstDigitalPeak_wsz

      # TODO only have to do some of this once, can skip parts for multiple iterations of while loop
      if(TRUE){
        ###############'
        ### Plot frequency (and segments and het score by segment) -------'
        plotFreqVsRd <- FALSE
        if(plotFreqVsRd){
          yMax <- max(frq00  / expReadsInFirstDigitalPeak_wsz)
          yMaxPlot <- min(6, yMax) # do not extend y axis too far
          # plot on normalized y axis...first digital peak will be y=1
          plot(wdnsMSK00,frq00/expReadsInFirstDigitalPeak_wsz,"l",col='gray', main="",
               ylim = c(0,  yMaxPlot),
               xlim = c(0,  temp00chrEnd[maxcn]), # temp00chrEnd[maxcn] will extend the axis to the end of chrY rather than the end of the Y data
               xaxt="n", xaxs='i')
          if(!is.null(alternateId)){
            mtext(3, text=c(sampleId, alternateId),adj=c(0,1))
          }else{
            mtext(3, text=c(sampleId),adj=c(0))
          }
          markChromEdges(chromStarts = temp00chrStart,maxcn = maxcn, vCol='gray30')
        }

        # mean and median Het. score for each segment
        # pkmod,  read depth normalized to the read depth of the first digital peak, if the first digital peak changes, pkmod will too
        # valid (1 means both avg and median loh != 0)

        ### initialize "valid" as 0
        # 0 = don't use this segment, hetScore data is not valid avg and/or median = 0
        # 1 = valid, will be filled in later after validity is confirmed
        segmentData[,'lohScoreMean'] <- NA
        segmentData[,'lohScoreMedian'] <- NA
        segmentData[,'pkmod'] <- NA
        segmentData[,'valid'] <- 0

        for(segmentId in seq_len(nrow(segmentData))) {
          # keys=which(segmentData$chr==1); segmentData[keys,];  keys=which(segmentData$lohScoreMean>1.01); segmentData[keys,]
          # segmentId=1+segmentId
          segment <- segmentData[segmentId,]
          cn1 <- segment[['chr']] # Chromosome

          chrStart <- binnedPosStart(coords@chromStart[cn1], binSize = wsz)
          chrEnd   <- binnedPosEnd(coords@chromEnd[cn1], binSize = wsz)

          ### heterozygosity scores
          lohSeqname <- convertChromToCharacter(cn1, withChrPrefix=TRUE)
          intersectingLoh <- which(hetScoreData[,'seqnames']==lohSeqname &
                                     hetScoreData[,'start'] >= segment[['start']] &
                                     hetScoreData[,'start'] < segment[['end']]) # TODO: should this be hetScoreData[,'end']?????

          # 'NoMask' not masked or checked with the hsNormMat
          meanLohNoMask    <- mean(hetScoreData[intersectingLoh, 'score'])
          medianLohNoMask  <- median(hetScoreData[intersectingLoh, 'score'])

          # check the hetScores in the normal samples
          # for each row of data, how many samples have a value outside of the normal range?
          if(!is.null(hsNormMat)){
            lohChrMed <- apply(hsNormMat[intersectingLoh,,drop = FALSE],1,function(x) (sum(x<0.975 | x> 2)))
            numLohRefSamples<-ncol(hsNormMat) #  23, as in the 23 TCGA normals
            hsNormMatCutoff <-  numLohRefSamples/2
            whichToKeep <- which(lohChrMed < hsNormMatCutoff)
            passedMaskedIntersectingLoh <- (intersectingLoh)[whichToKeep]

            # plot hsNormMat
            if(FALSE){
              plot(hsNormMat[,1],type = 'p', col=1)
              for(i in 2:23){
                points(hsNormMat[,i],type = 'p', col=i)
              }
              temp00chrStart <- binnedPosStart(coords@chromStart, binSize = wsz) # start positions must be same bin size as the plotted data
              markChromEdges(chromStarts = temp00chrStart,maxcn = maxcn, vCol='gray30')
            }

            if(length(passedMaskedIntersectingLoh) > 0){
              meanLoh   <- mean(hetScoreData[passedMaskedIntersectingLoh, 'score'])
              medianLoh <- median(hetScoreData[passedMaskedIntersectingLoh, 'score'])
            }else{
              meanLoh   <- 0
              medianLoh <- 0
            }

            # plot hetScores masked and unmasked for segmentId
            if(FALSE){
              op <- par(mfrow=c(2,1),mar=c(2.5, 3.5, 1.5, 1),mgp=c(1.5, 0.5,0))
              hist(lohChrMed, breaks = 23, col='gray', main=paste('chrom: ', cn1, 'segmentId:', segmentId))
              abline(v=hsNormMatCutoff)
              if(!is.null(alternateId)){
                mtext(3, text=c(sampleId, alternateId),adj=c(0,1))
              }else{
                mtext(3, text=c(sampleId),adj=c(0))
              }
              plot(intersectingLoh,hetScoreData[intersectingLoh, 'score'])
              abline(h=meanLoh, col='red')
              abline(h=meanLohNoMask, col='black')
              abline(h=0.98, lty='dashed')
              points(intersectingLoh,hetScoreData[intersectingLoh, 'score'],)
              points(passedMaskedIntersectingLoh,hetScoreData[passedMaskedIntersectingLoh, 'score'],pch=19, col='red',cex=.3)
              par(op)
            }

          }else{
            # the mask is not provided so use the unmasked data
            meanLoh   <- meanLohNoMask
            medianLoh <- medianLohNoMask
          }

          # reads per segment
          chromStart <- coords@chromStart[cn1]
          segmentStartWindowIndex <- which.min(abs(binnedPosStart(chromStart + segment[['start']], binSize=wsz) - wdnsMSK00))
          segmentEndWindowIndex   <- which.min(abs(binnedPosEnd(chromStart + segment[['end']], binSize=wsz) - wdnsMSK00))
          meanReadsPerSeg <- mean(frq00[segmentStartWindowIndex:segmentEndWindowIndex],na.rm = TRUE)   # ME26210 has an outlier on chr 4:33116154 36132923, maybe median is better?
          medianReadsPerSeg <- median(frq00[segmentStartWindowIndex:segmentEndWindowIndex],na.rm = TRUE)
          pkmod <- medianReadsPerSeg / expReadsInFirstDigitalPeak_wsz        # normalize to the read depth of the first digital peak
          distanceFromFirstDigitalPeak    <- abs(1 - pkmod)                # subtract from 1, the y axis is frq/expReadsInFirstDigitalPeak_wsz

          # update array if heterozygosity scores are valid, and color code plot accordingly
          if(meanLoh>0 & medianLoh>0){
            segmentData[segmentId, 'valid'] <- 1
            meanCol <- 'cyan'
            medianCol <- 'magenta'
            segmentCol <- 'black'
          }else{
            meanCol <- 'orange'
            medianCol <- 'orange'
            segmentCol <- 'green3'
          }

          # add segments to data.frame
          segmentData[segmentId, 'pkmod']       <- pkmod
          segmentData[segmentId, 'lohScoreMean'] <- meanLoh
          segmentData[segmentId, 'lohScoreMedian'] <- medianLoh

          if(plotFreqVsRd){
            # add segments to plot and data.frame
            lines(x=wdnsMSK00[c(segmentStartWindowIndex, segmentEndWindowIndex)], y=rep(pkmod, 2),
                  col=ifelse(distanceFromFirstDigitalPeak < digitalPeakZone, 'red', segmentCol), lwd=3)

            # add LOH score to plot
            points(x=mean(wdnsMSK00[c(segmentStartWindowIndex, segmentEndWindowIndex)]), y=meanLoh, col=meanCol, bg='black', pch=21, cex=.4)
            points(x=mean(wdnsMSK00[c(segmentStartWindowIndex, segmentEndWindowIndex)]), y=medianLoh, col=medianCol, pch=1)
          }
        } # CLOSE: for(segmentId in seq_len(nrow(segmentData)))

        numSegmentsValid <- length(which(segmentData$valid==1))
        percentSegmentsValid <- numSegmentsValid/nrow(segmentData)*100
        if(numSegmentsValid<150){
          logwarn('numSegmentsValid is low, may want to consider smaller segment chunk size? %i < %i',numSegmentsValid,150)
        }
        segmentCountText <- paste0( numSegmentsValid, '/', nrow(segmentData),' (',round(percentSegmentsValid,1),'%) ', 'segments with valid het. score')
        if(plotFreqVsRd){
          mtext(1, text= segmentCountText,adj=0)

          abline(h=heterozygosityScoreThreshold, col='gray', lty='dotted', lwd=0.65)
          legend("topright", legend=c('1st Digital', 'Others', 'mean score', 'median score ','invalid score'), col=c('red', 'black', 'cyan', 'magenta','orange'),
                 lwd=c(2,2,NA,NA,NA), pch=c(NA,NA,21,1,1), pt.cex = .6, pt.bg = 'black',
                 inset = c(0, .05),cex=.8)
        }
      } # CLOSE: TRUE chunk


      # These are indices of segments close to First Digital Peak
      firstDigPeakIndex <- which.min(peakInfo$dPeaks) # may not be 1 if the grid is not forced to start at the first peak
      normXofFirstDigPeak <- peakInfo$peakReadDepth_normX[firstDigPeakIndex] # 1 if first digital peak is first peak, but get this value in case it isn't
      # translate in case the first digital peak is not the first peak
      peakModes <- peakInfo[, 'peakReadDepth_normX']/normXofFirstDigPeak  # normXofFirstDigPeak will be 1 if first peak is first digital peak (dPeak)


      ### key het score quantiles for each Peak (digital and non-digital) -----------------------

      # initialize
      segmentsCloseToPeak <- matrix(nrow=nrow(segmentData),ncol=numPeaks)
      probs <- c(0, 0.05, .5, 0.75, 0.85, 0.95, 0.99,1)
      keyHetScoresPerPeakAccum <- matrix(nrow=numPeaks, ncol = length(probs))
      colnames(keyHetScoresPerPeakAccum) <- paste0('q',probs)

      ### fill in
      # qProbs will be NA if there are no valid segments in that peak
      for(iPeak in seqFwd(1,numPeaks)){
        segmentsCloseToPeak[,iPeak] <- abs(peakModes[iPeak]-segmentData$pkmod) < digitalPeakZone # TRUE or FALSE
        validPeakIndexes <- which(segmentsCloseToPeak[,iPeak] &
                                    segmentData$valid==1)

        lohScores <- segmentData[validPeakIndexes,'lohScoreMedian']
        qProbs   <- quantile(lohScores,na.rm =TRUE,probs=probs)
        keyHetScoresPerPeakAccum[iPeak,] <- qProbs

      }


      # add some peak info
      keyHetScoresPerPeak <- cbind('rankByHeight' = peakInfo$rankByHeight,
                                   'dPeak'       = digitalPeakIndex,
                                   'dPeakPos'    = peakInfo$dPeaks,
                                   keyHetScoresPerPeakAccum)

      ### make sure first digital peak is labeled as 1
      if(nCopyPeaks_dig[1]!=1){
        stop(sprintf('nCopyPeaks_dig[1]= %i and not 1 like expected, what happened in digital peaks?',nCopyPeaks_dig[1]))
      }



      ### assign nCopyPeaks_seg by checking het score of the first digital peak ------------------
      if(TRUE){
        ## multiple methods:

        hetTestScoreFor1stDigPeak <- NA

        ## instead of doing mean over all, find a counter-example - a long enough region that has het. score > .98
        ##    That alone proves this peak is 2N and not 1N ?
        ## should only have to check the first digital peak, but there are cases where the first digital peak does not have valid LOH scores,
        ##    then go to the next digital peak i.e. 58024,82012 ? or just assign this peak as 1N?


        ## method 1) for looking at multiple digital peaks to find a valid digital peak
        if(FALSE){
          dPeak_forMaxLOHscore <- 0
          while(is.na(hetTestScoreFor1stDigPeak)){
            dPeak_forMaxLOHscore <- dPeak_forMaxLOHscore+1  # the digital peak used to get the maxLOHscore, look at the first and then move up as needed
            # check for runaway while loop
            if(dPeak_forMaxLOHscore>numDigPeaks){
              stop("something went screwy with getting the max LOH score from a non-digital peak")
            }
            hetTestScoreFor1stDigPeak <- keyHetScoresPerPeak[which(digitalPeakIndex==dPeak_forMaxLOHscore) ,'q1'] # use which() to eliminate NAs
          }

          if(dPeak_forMaxLOHscore==1){
            ## STANDARD, most of the time the first digital peak will be the evaluated peak
            if(sum(segmentsCloseToPeak[,dPeak_forMaxLOHscore], na.rm=TRUE) > 0 && hetTestScoreFor1stDigPeak >= heterozygosityScoreThreshold) {
              # if TRUE, then the digital peak does not have LOH and cannot correspond to deletion, therefore it must be a 2N peak (or 4N or 6N...)
              nCopyPeaks_seg <- nCopyPeaks_dig + 1
            }else{
              nCopyPeaks_seg <- nCopyPeaks_dig
            }
          }else{
            ## could not use first digital peak because there were no valid het. scores.
            logwarn('using het. scores of digital peak number %i to determine copy number', dPeak_forMaxLOHscore)
            if(sum(segmentsCloseToPeak[,dPeak_forMaxLOHscore], na.rm=TRUE) > 0 && hetTestScoreFor1stDigPeak >= heterozygosityScoreThreshold) {
              # if TRUE, then the digital peak does not have LOH and cannot correspond to deletion, therefore it must be a 2N peak (or 4N or 6N...)

              if((nCopyPeaks_dig[dPeak_forMaxLOHscore]%%2)==0){ # if it was previously labeled 2N, 4N, 6N etc.
                nCopyPeaks_seg <- nCopyPeaks_dig
              }else {
                nCopyPeaks_seg <- nCopyPeaks_dig-(dPeak_forMaxLOHscore-1)
              }
            }else{
              # if FALSE, then the digital peak does have LOH or is a deletion, therefore it must be a 1N peak
              if((nCopyPeaks_dig[dPeak_forMaxLOHscore]%%2)==0){ # if it was previously labeled 2N, 4N, 6N etc.
                nCopyPeaks_seg <- nCopyPeaks_dig-(dPeak_forMaxLOHscore-1)
              }else {
                nCopyPeaks_seg <- nCopyPeaks_dig
              }
            }
          }
        }else if(FALSE){
          ## method 2: for looking at only the first digital peak, if it doesn't have valid scores, assume it is LOH or AOH.
          ## using max scoredPeak_maxLOHscore <- 1  # the digital peak used to get the maxLOHscore, look at the first only
          hetTestScoreFor1stDigPeak <- keyHetScoresPerPeak[which(digitalPeakIndex==dPeak_forMaxLOHscore) ,'q1'] # use which() to eliminate NAs

          if(is.na(hetTestScoreFor1stDigPeak)){
            nCopyPeaks_seg <- nCopyPeaks_dig
            ## STANDARD, most of the time the first digital peak will be the evaluated peak
          }else if(sum(segmentsCloseToPeak[,dPeak_forMaxLOHscore], na.rm=TRUE) > 0 && hetTestScoreFor1stDigPeak >= heterozygosityScoreThreshold) {
            # if TRUE, then the digital peak does not have LOH and cannot correspond to deletion, therefore it must be a 2N peak (or 4N or 6N...)
            nCopyPeaks_seg <- nCopyPeaks_dig + 1
          }else{
            nCopyPeaks_seg <- nCopyPeaks_dig
          }
        }else if(TRUE){
          ### method 3: check the max MODE from het score density of the first digital peak ------------------
          op <- par(mfrow=c(2,1),mar=c(2.75, 3.5, 2, 1.5),mgp=c(1.5, 0.5,0))

          densityFirstDigPeak <- hetScoreDensity(segmentsCloseToPeak,segmentData, index=firstDigPeakIndex,sampleId,alternateId, skipPlot = skipExtras, plotTextPrefix='1st digital peak:',
                                                 heterozygosityScoreThreshold=heterozygosityScoreThreshold)
          hetTestScoreFor1stDigPeak <- densityFirstDigPeak$testScore
          loginfo('het test score for first digital peak: %.3f',densityFirstDigPeak$testScore)

          # not used right away but is used later
          if(all(peakInfo$dPeaks != min(peakInfo$dPeaks, na.rm = TRUE))){
            stop('all NAs need to inspect and debug this step %s %s', sampleId, alternateId)
          }

          if(numDigPeaks>1){
            secondDigPeak  <-  min(peakInfo$dPeaks[peakInfo$dPeaks != min(peakInfo$dPeaks, na.rm = TRUE)], na.rm = TRUE)
            secondDigPeakIndex <-     which(peakInfo$dPeaks== secondDigPeak)
            densitySecondDigPeak <- hetScoreDensity(segmentsCloseToPeak,segmentData, index=secondDigPeakIndex,sampleId,alternateId, skipPlot = skipExtras, plotTextPrefix='2nd digital peak:',
                                                    heterozygosityScoreThreshold=heterozygosityScoreThreshold)
            loginfo('het test score for second digital peak: %.3f',densitySecondDigPeak$testScore)
          }else{
            secondDigPeakIndex <- NA
            densitySecondDigPeak <- list(testScore=0,numHetClusters=NA,observ=0)
          }



          ### first digital peak logic -------------
          ## determine nCopyPeaks_seg based on density of het scores from the segments of the first digital peak
          if(is.na(densityFirstDigPeak$testScore)){
            stop(sprintf('densityFirstDigPeak$testScore is NA, check what happened? why?'))
          }else if(sum(segmentsCloseToPeak[,firstDigPeakIndex], na.rm=TRUE) > 0 && # TODO: require more than 1 segment for this test?
                   densityFirstDigPeak$testScore >= heterozygosityScoreThreshold) {
            # if TRUE, then the digital peak has heterozygosity and can NOT be 1N (deletion), therefore it must be at least a 2N peak
            # TODO: is this assuming nCopyPeaks_dig[1] ==1
            nCopyPeaks_seg <- nCopyPeaks_dig + 1
            nCopyAltered <- TRUE
          }else{
            # otherwise there is LOH, this could be 1N, or there are no valid segments and can't judge so assume 1N,
            #   or a higher CN level where the first peak is entirely loh, or even 3N with a very low tumor percent or low coverage, but that scenario is addressed later
            nCopyPeaks_seg <- nCopyPeaks_dig
            nCopyAltered <- FALSE
          }
        }

      } # CLOSE: if(TRUE) run this chunk -- check het score of the first digital peak

    } # CLOSE bySegment


    ### define nCopyPeaks_while
    if(forceFirstDigPeakCopyNum >=0){
      nCopyIncrement <- forceFirstDigPeakCopyNum - nCopyPeaks_seg[1]
      nCopyPeaks_while <- nCopyPeaks_seg + nCopyIncrement
      loginfo('ploidy results determined by segment, first peak is %iN', nCopyPeaks_seg[1])
      loginfo(' but forceFirstDigPeakCopyNum is enabled so now first peak is %iN', nCopyPeaks_while[1])
    }else{
      nCopyPeaks_while <- nCopyPeaks_seg
      loginfo('ploidy results determined by segment, first peak is %iN', nCopyPeaks_seg[1])
    }





    whileLoopStep <- repeatSteps(peakInfo,keyHetScoresPerPeak,nCopyPeaks_step=nCopyPeaks_while)
    peakInfo <- whileLoopStep$peakInfo
    percentTumor <- whileLoopStep$percentTumor
    hetScoreQuantiles <- whileLoopStep$hetScoreQuantiles


    ###############################################'
    # redo digital peaks? ----------------------
    ###############################################'

    if(TRUE){
      # TODO: what is the correct value for 'allowedTumorPercent'
      #     PDXs will be high ie 78014=108.4 78024=103.4, perhaps they need a different limit than tumor?
      #     need a bit of a tolerance; 49052 is 106%, 26295 was finding 106% for 2N, 3N-main, which is wrong but digital peaks isn't having this problem anymore

      ####check tumor percent -----'
      checkTumorPercent <- function(percentTumor, allowed=allowedTumorPercent){
        if(is.na(percentTumor)){
          pTumorPasses <- NA
        }else if(percentTumor > allowed){
          pTumorPasses <- FALSE
          logerror('percentTumor: %s NOT VALID more than allowed limit of %i', round(percentTumor,2),allowed)
        }else{
          pTumorPasses <- TRUE
        }
        return(pTumorPasses)
      }

      pTumorPasses <- checkTumorPercent(percentTumor)

      #### main peak is 4N or less -----'
      if( !is.na(peakInfo[mainPeakIndex,'nCopy']<=4) && peakInfo[mainPeakIndex,'nCopy']<=4 ){
        mainMaxCopyPasses <-TRUE
      }else{
        mainMaxCopyPasses <-FALSE
        logerror('mainPeak is %s; -too- high (or NA); redo digital peaks', peakInfo[mainPeakIndex,'nCopy'] )
      }


      #### 2N level present -----'
      # twoNPresentPasses <- FALSE 58081
      if(length(which(peakInfo$nCopy==2)) > 0){
        twoNPresentPasses <- TRUE
      }else{
        twoNPresentPasses <- FALSE
        logerror('no diploid level, redo digital peaks')
      }



      ### test each metric, set gridIteration
      if(forceFirstDigPeakCopyNum >= 0){
        # even if stuff fails, we are forcing this situation for a reason so ...
        # Done with while loop
        gridIteration <- gridIteration+100 # high value will exit out of while loop
      }else if(all(c(mainMaxCopyPasses,pTumorPasses,twoNPresentPasses), na.rm=TRUE)){
        # Done with while loop
        gridIteration <- gridIteration+100 # high value will exit out of while loop
      }else{
        if(gridIteration+1 < numPeaks){
          gridIteration <- gridIteration+1
          logwarn("rerun digital grid: %s", gridIteration)
        }else{
          gridIteration <- gridIteration+1000 # high value will exit out of while loop
          logwarn('digital grid not quite right,(not TRUE: %s), but too many gridIterations; exiting while loop to try something different',
                  paste(c('mainMaxCopyPasses','pTumorPasses','twoNPresentPasses')[-which(c(mainMaxCopyPasses,pTumorPasses,twoNPresentPasses))]))
        }

      }
    } ### CLOSE: TRUE chunk

  } ### end while loop





  #####'
  ## check hetScores are valid -----------
  #####'
  if(densityFirstDigPeak$testScore >= heterozygosityScoreThreshold &
     densitySecondDigPeak$testScore >= heterozygosityScoreThreshold){
    logwarn('both of the first two digital peaks have a het score above the threshold: %.3f & %.3f >= %.2f',
            densityFirstDigPeak$testScore, densitySecondDigPeak$testScore, heterozygosityScoreThreshold )
  }

  ############'
  ### two final adjustments to nCopyPeaks_while before we move on to testing
  ############'

  firstDigPeakIndex <- which.min(peakInfo$dPeaks)                              # will not be 1 if first peak is not the first digital peak

  ## invalid tumor percent check ---------
  # lets say even though invalid, the best possible digital peaks was been found, so lets trying adjusting the cn calls
  if(!is.na(pTumorPasses)){
    if(!pTumorPasses){

      # assume the first peak is correct and all the other peak CNs have to be bumped up
      # ie 26306,58130, 43010
      peakInfo['nCopyOrg'] <- peakInfo['nCopy']
      peakInfo['nCopy'] <- peakInfo$nCopyOrg+1
      peakInfo[firstDigPeakIndex, 'nCopy']  <-  peakInfo[firstDigPeakIndex, 'nCopyOrg']

      newTuPercent <- calcTumorFromPloidyPeaks(peakCopyNum = peakInfo$nCopy, peakHeight  = peakInfo$peakHeight, peakReadDepth_1bp=peakInfo$peakReadDepth_1bp,dPeaks = peakInfo$dPeaks)
      newTuPercentValid <-  checkTumorPercent(newTuPercent)
      if(!newTuPercentValid){
        stop("new tumor percent is STILL not valid")
      }else{
        loginfo('lets roll with this new Tumer Percent: %.1f and see how it goes', newTuPercent)
        percentTumor <- newTuPercent
        nCopyPeaks_while <- peakInfo[!is.na(peakInfo$nCopy),'nCopy']
      }

    }
  }


  minObservations  <-  20 # what is the correct min value? at least 15 see 80235
  ## TODO: consider making this a standalone test since this does not require a 3N peak?
  firstPeakMustBeTwoNorMore <-  hasNorMoreClusters(hetScoreDensityResult = densityFirstDigPeak,N=2, heterozygosityScoreThreshold, minObservations )  # true if it has multiple het score clusters
  if(nCopyPeaks_while[1]==1 & firstPeakMustBeTwoNorMore){
    nCopyPeaks_test <- nCopyPeaks_while+1
    testFirstPeakStep <- repeatSteps(peakInfo,keyHetScoresPerPeak,nCopyPeaks_step=nCopyPeaks_test)

    if(checkTumorPercent(testFirstPeakStep$percentTumor)){
      logwarn('must increase copy number by 1 because first peak has two hetScore clusters but was only 1N')
      nCopyPeaks_while <- nCopyPeaks_while+1
      peakInfo <- testFirstPeakStep$peakInfo
      percentTumor <- testFirstPeakStep$percentTumor
    }else{
      logwarn('will not increase copy number by 1 even though first peak has two hetScore clusters and is 1N because the tumor percent will not be valid')
      firstPeakMustBeTwoNorMore <- FALSE
    }

  }



  # add hetScore to peak info
  peakInfo[,'hetScore'] <- NA_real_
  peakInfo[,'hetScoreObserv'] <- NA_integer_

  peakInfo[firstDigPeakIndex,'hetScore'] <- densityFirstDigPeak$testScore
  peakInfo[firstDigPeakIndex,'hetScoreObserv'] <- densityFirstDigPeak$observ
  if(!is.na(secondDigPeakIndex)){
    peakInfo[secondDigPeakIndex,'hetScore'] <- densitySecondDigPeak$testScore
    peakInfo[secondDigPeakIndex,'hetScoreObserv'] <- densitySecondDigPeak$observ
  }



  ################################'
  ## testing and finalizing results ----------------------------------
  ################################'

  # forceFirstDigPeakCopyNum  - use this if ploidy calc is wrong and you need to specify it manually to get it correct
  # otherwise see if both of the first two peaks have high het Scores (happens in low tumor), if so, find which one is bigger and more reliable and assign it to the 2N peak
  # otherwise run through options for checking the copy number of first digital peak by testing the 3N peak

  test3Npeaks <- TRUE # but might change to FALSE below...

  # set 'test3Npeaks' to FALSE if necessary
  if(forceFirstDigPeakCopyNum >= 0){
    nCopyIncrement <- forceFirstDigPeakCopyNum- nCopyPeaks_while[1]
    nCopyPeaks_final <- nCopyPeaks_while+nCopyIncrement
    forceStep <- repeatSteps(peakInfo,keyHetScoresPerPeak,nCopyPeaks_step=nCopyPeaks_final)
    peakInfo <- forceStep$peakInfo
    percentTumor <- forceStep$percentTumor
    hetScoreQuantiles <- forceStep$hetScoreQuantiles
    hetTestScoreFor3NPeak <- NA
    test3Npeaks <- FALSE
    loginfo("first digital peak forced to %iN with %.1f tumor percent", nCopyPeaks_final[1], percentTumor)
  }else if(length(nCopyPeaks_while)==1 &
           # numPeaks==1 &  # test length of nCopyPeaks_while rather than numPeaks
           nCopyPeaks_while[1]==1){
    # there is only one peak and it is below threshold (hence it is called 1N) suspect depressed hetScores so bump up the copy number
    nCopyPeaks_final <- nCopyPeaks_while+1
    onePeakStep <- repeatSteps(peakInfo,keyHetScoresPerPeak,nCopyPeaks_step=nCopyPeaks_final)
    peakInfo <- onePeakStep$peakInfo
    percentTumor <- onePeakStep$percentTumor
    hetScoreQuantiles <- onePeakStep$hetScoreQuantiles
    hetTestScoreFor3NPeak <- NA
    test3Npeaks <- FALSE
    loginfo("first digital peak changed from 1N to 2N even though hetscore was below the threshold (%.3f) because it is the only peak", densityFirstDigPeak$testScore)

  }else if(densityFirstDigPeak$testScore >= heterozygosityScoreThreshold &
           densitySecondDigPeak$testScore >=heterozygosityScoreThreshold ){
    loginfo('the first two peaks both have het scores above the threshold')

    # see if both of the first two peaks have high het Scores (happens in low tumor), if so, find which is bigger and more reliable and assign it to the 2N peak.
    # ie. 58047, 80235, 58091
    if( which.max(c(densityFirstDigPeak$testScore, densitySecondDigPeak$testScore))==2 &
        which.max(c(densityFirstDigPeak$observ,    densitySecondDigPeak$observ   ))==2){
      #  strong support that of the two peaks, the 2nd dig peak is more likely the 2N peak

      if(nCopyPeaks_while[2]==2){
        test3Npeaks <- FALSE

        loginfo('the first two peaks both have high het scores but no change to copy number needed')
        nCopyPeaks_final <- nCopyPeaks_while
        bothHighHetStep <- repeatSteps(peakInfo,keyHetScoresPerPeak,nCopyPeaks_step=nCopyPeaks_final)
        peakInfo <- bothHighHetStep$peakInfo
        percentTumor <- bothHighHetStep$percentTumor
        hetScoreQuantiles <- bothHighHetStep$hetScoreQuantiles
        hetTestScoreFor3NPeak <- NA

      }else if(nCopyPeaks_while[2]==3 & !firstPeakMustBeTwoNorMore){
        # if there is a third peak, check it to see if it could be 4N(2:2), if not then ok to reduce copy number.
        # if so, do the other checks
        if(length(nCopyPeaks_while)>=3){
          thirdDigPeak  <-  min(peakInfo$dPeaks[!peakInfo$dPeaks %in% c(peakInfo$dPeaks[firstDigPeakIndex], secondDigPeak) ], na.rm = TRUE)
          thirdDigPeakIndex <-     which(peakInfo$dPeaks== thirdDigPeak)
          densityThirdDigPeak <- hetScoreDensity(segmentsCloseToPeak,segmentData, index=thirdDigPeakIndex,sampleId,alternateId, skipPlot = skipExtras, plotTextPrefix='1st digital peak:',
                                                 heterozygosityScoreThreshold=heterozygosityScoreThreshold)
          if(densityThirdDigPeak$testScore < heterozygosityScoreThreshold &
             densityThirdDigPeak$observ > 15){
            logwarn('reducing copy number by 1 (because it is safe to do so: 1st peak will not be 0N and does not have multiple het score clusters ) and
            3rd peak is not 4N (2;2) and
            because it is more likely that the first two peaks are 1N 2N, rather than 2N 3N')

            nCopyPeaks_final <- nCopyPeaks_while-1
            test3Npeaks <- FALSE
            bothHighHetStep <- repeatSteps(peakInfo,keyHetScoresPerPeak,nCopyPeaks_step=nCopyPeaks_final)
            peakInfo <- bothHighHetStep$peakInfo
            percentTumor <- bothHighHetStep$percentTumor
            hetScoreQuantiles <- bothHighHetStep$hetScoreQuantiles
            hetTestScoreFor3NPeak <- NA
          }

        }else{
          # no third peak so go ahead
          logwarn('reducing copy number by 1 (because it is safe to do so: 1st peak will not be 0N and does not have multiple het score clusters ) and
            because it is more likely that the first two peaks are 1N 2N, rather than 2N 3N')

          nCopyPeaks_final <- nCopyPeaks_while-1
          test3Npeaks <- FALSE
          bothHighHetStep <- repeatSteps(peakInfo,keyHetScoresPerPeak,nCopyPeaks_step=nCopyPeaks_final)
          peakInfo <- bothHighHetStep$peakInfo
          percentTumor <- bothHighHetStep$percentTumor
          hetScoreQuantiles <- bothHighHetStep$hetScoreQuantiles
          hetTestScoreFor3NPeak <- NA
        }
      }


    }else if( which.max(c(densityFirstDigPeak$testScore, densitySecondDigPeak$testScore))==1 &
              which.max(c(densityFirstDigPeak$observ,    densitySecondDigPeak$observ   ))==1){
      #  strong support that of the two peaks, the 1st dig peak is more likely the 2N peak
      if(nCopyPeaks_while[1]==2){
        loginfo('the first two peaks both have high het scores but no change to copy number needed')
        nCopyPeaks_final <- nCopyPeaks_while
      }else if(nCopyPeaks_while[1]==1){
        logwarn('increasing copy number by 1 because the first two peaks both have high het scores and
                  it is more likely that the first two peaks are 2N 3N, rather than 1N 2N,
                  because the first peak hetScore is bigger AND has more observations')
        nCopyPeaks_final <- nCopyPeaks_while+1
      }else if(nCopyPeaks_while[1]>=3){
        stop(sprintf('nCopyPeaks_while[1] = %i, this is not possible, is it? inspect further!',nCopyPeaks_while[1]))
      }
      test3Npeaks <- FALSE
      bothHighHetStep <- repeatSteps(peakInfo,keyHetScoresPerPeak,nCopyPeaks_step=nCopyPeaks_final)
      peakInfo <- bothHighHetStep$peakInfo
      percentTumor <- bothHighHetStep$percentTumor
      hetScoreQuantiles <- bothHighHetStep$hetScoreQuantiles
      hetTestScoreFor3NPeak <- NA

    }else{
      logwarn('the first two peaks both have high het scores but it is not clear which one should be the 2N peak')
    }

  }






  ################################'
  ## test 3N peak ----------------------------------
  ################################'
  if(test3Npeaks){

    # TODO: add exit status values so I know which test the sample exited on. Track which tests never get used or seem to fail
    ### density of mean loh scores in the 3N peak: get mode, see if mode fits the expected side of cutoff for 3N, adjust accordingly

    testsForIncreasingCNof4Npeak <- function(){

      if(!is.na(densityResult4$testScore)){
        # testing if it is ok to make the tested 3N peak the 4N peak, this test could fail if the testScore of the 3N peak (proposed 4N) has LOH, uncommon but maybe common enough to cause problems?
        if(densityResult4$testScore < densityResult3$testScore){   # proposed 5N peak < proposed 4N peak ?
          logwarn('increase copy number by 1 because it is safe to do so: het score of 4N peak (proposed 5N) is less than 3N peak (proposed 4N): %.3f < %.3f',  densityResult4$testScore, densityResult3$testScore)
          nCopyPeaks_final <- nCopyPeaks_while+1
        }else{
          loginfo('NOT increasing copy number by 1 because it is NOT safe to do so: het score of 4N peak (proposed 5N) would be more than 3N peak (proposed 4N): %.3f !< %.3f',  densityResult4$testScore, densityResult3$testScore)
          nCopyPeaks_final <- nCopyPeaks_while
        }
      }else{
        loginfo('NOT increasing copy number by 1 because it is NOT safe to do so: no het score from 4N peak (proposed 5N) to test')
        nCopyPeaks_final <- nCopyPeaks_while
      }
      return(nCopyPeaks_final)
    }



    ### does 2N, 3N, 4N exist?
    cn2indexTemp <- which(peakInfo$nCopy==2); cn2existsTemp <- ifelse(length(cn2indexTemp)!=0, TRUE, FALSE)
    cn3indexTemp <- which(peakInfo$nCopy==3); cn3existsTemp <- ifelse(length(cn3indexTemp)!=0, TRUE, FALSE)
    cn4indexTemp <- which(peakInfo$nCopy==4); cn4existsTemp <- ifelse(length(cn4indexTemp)!=0, TRUE, FALSE)

    if(cn2existsTemp){
      if(cn2indexTemp==firstDigPeakIndex){
        # already have this, don't need to redo
        densityResult2 <- densityFirstDigPeak
      }else if(cn2indexTemp==secondDigPeakIndex){
        # already have this, don't need to redo
        densityResult2 <- densitySecondDigPeak
      }else{
        densityResult2 <- hetScoreDensity(segmentsCloseToPeak,segmentData, index=cn2indexTemp,sampleId,alternateId, skipPlot = skipExtras,plotTextPrefix='initial 2N peak:',
                                          heterozygosityScoreThreshold=heterozygosityScoreThreshold)
      }
    }else{
      densityResult2 <- NA
    }


    if(!cn3existsTemp){
      loginfo('no 3N peak to test')
      nCopyPeaks_final <- nCopyPeaks_while
      hetTestScoreFor3NPeak <- NA
      densityResult3 <- NA
    }else{
      if(cn3existsTemp==secondDigPeakIndex){
        # already have this, don't need to redo
        densityResult3 <- densitySecondDigPeak
      }else {
        densityResult3 <- hetScoreDensity(segmentsCloseToPeak,segmentData, index=cn3indexTemp,sampleId,alternateId, skipPlot = skipExtras,plotTextPrefix='initial 3N peak:',
                                          heterozygosityScoreThreshold=heterozygosityScoreThreshold)
      }

      if(densityResult3$observ < 15){
        logwarn('only %i segments in 3N peak, not enough to test reliably, not changing copy number of 3N peak', densityResult3$observ)
        nCopyPeaks_final <- nCopyPeaks_while
        hetTestScoreFor3NPeak <- NA
      }else{
        hetTestScoreFor3NPeak <- densityResult3$testScore # hetTestScoreFor3NPeak does not get updated if copy number is increased or reduced during test3Npeaks


        if(cn4existsTemp){
          densityResult4 <- hetScoreDensity(segmentsCloseToPeak,segmentData, index=cn4indexTemp,sampleId,alternateId, skipPlot = skipExtras,plotTextPrefix='initial 4N peak:',
                                            heterozygosityScoreThreshold=heterozygosityScoreThreshold)
        }else{
          densityResult4 <- NA
        }

        # can only calculate Tumor percentage if there is more than one nCopyPeaks_step
        calcTumorStatus <- ifelse( numDigPeaks > 1, TRUE, FALSE) #   numDigPeaks = length(which(!is.na(peakInfo$dPeaks)))
        if(calcTumorStatus){
          loginfo('before testing 3N peak...')
          percentTumorTemp <- calcTumorFromPloidyPeaks(peakCopyNum = peakInfo$nCopy,
                                                       peakHeight  = peakInfo$peakHeight,
                                                       peakReadDepth_1bp=peakInfo$peakReadDepth_1bp,
                                                       dPeaks = peakInfo$dPeaks)
        }else{
          percentTumorTemp <- NA_real_
        }

        # TODO: 26278 do not want to decrease copy number because then the new 3N peak will have the highest het score, even bigger than the new 4N peak
        #       test for this occurrence before deciding to increase/decrease copy number.
        # TODO: 58061 do not want to increase copy number because then the new 5N peak will have the highest het score, even bigger than the new 4N peak
        #       test for this occurrence before deciding to increase/decrease copy number.

        # to compare to the 3N peak het scores

        ### test 3N peak for three het score clusters
        # TODO: test lowering the heterozygosityScoreThreshold for this test? ie. LU19102 & LU19103 right cluster=0.978
        peakHas3HetScoreClusters <- hasNorMoreClusters(densityResult3,N=3,heterozygosityScoreThreshold = heterozygosityScoreThreshold-.005 )

        if(peakHas3HetScoreClusters){
          # incrementBy=
          logwarn('temp3N peak has 3 or more het score clusters, must be 4N or more, increasing copy number by 1')
          nCopyPeaks_final <- nCopyPeaks_while+1

        }else if(hetTestScoreFor3NPeak >= heterozygosityScoreThreshold){
          logwarn('temp 3N Peak het score is above threshold: %.3f >= %.2f',  hetTestScoreFor3NPeak,heterozygosityScoreThreshold )

          # make sure percentTumor isn't too low, otherwise checking the 3N peak isn't a good test.
          if(!is.na(percentTumorTemp) & percentTumorTemp >= 15){ # TODO: what is the correct value? This percent tumor isn't the final value though, will change if copy number is changed
            # alter copy number in order of confidence (confidence is not the right word here, safety?)
            # reduce copy number except for a couple of situations when it is not safe to reduce copy number, in which case increase copy number

            ## can't decrease copy number if....
            # if 1st peak has more than 1 hetScore cluster, and the copy number is only 2N, otherwise increase
            # that is if NumClusters =2 and nCopy=2, can not reduce nCopy to 1 because 1N can't have two clusters

            # PT58126 - 1N peak 1 cluster, 2N peak 2 clusters
            # PT58147 -  fails because it is increased when it shouldn't, the first peak has 2 distinct peaks and this should be incorporated.

            ##  First see if we --MUST-- increase copy number?
            if(firstPeakMustBeTwoNorMore &&
               nCopyPeaks_while[1]==1){  # if numHetClusters=2, copyNum must be at least 2, so increasing is ok, in fact it is necessary
              logwarn('increasing copy number by 1 because first peak needs it too')
              nCopyPeaks_final <- nCopyPeaks_while+1

              # do some checks before decreasing copy number such that   3-->2 and 2-->1
            }else if(nCopyPeaks_while[1] > 2 |                                # confirm first digital peak will not become 0N, we are assuming that 1st digital peak will never be 0N  # TODO is ">2" correct or should it be ">=2"?
                     (!firstPeakMustBeTwoNorMore & nCopyPeaks_while[1]==2) ){  # confirm first digital peak will not become 1N if we know it must be 2N

              if(cn4existsTemp){
                if(!is.na(densityResult4$testScore)){
                  # if(densityResult4$testScore < heterozygosityScoreThreshold){ # check proposed 3N peak -- this is not as reliable in low tumor samples
                  if(densityResult4$testScore < densityResult3$testScore ){       # compare to proposed 3N peak -- only works if there is a higher peak to test
                    logwarn('reducing copy number by 1 because it is safe to do so: no 0N peak and het score of 4N peak (proposed 3N) is less than 3N peak (proposed 2N): %.3f < %.3f',  densityResult4$testScore, densityResult3$testScore)
                    nCopyPeaks_final <- nCopyPeaks_while-1
                  }else{
                    # TODO: what about when all the het Scores are above threshold (due to low tumor % ie 80180)
                    #       for example could reduce copy number when densityResult1$testScore & densityResult2$testScore & densityResult3$testScore are all above .98
                    loginfo('NOT reducing copy number by 1 because it is NOT safe to do so: het score of 4N peak (proposed 3N) would be more than 3N peak (proposed 2N): %.3f !< %.3f',  densityResult4$testScore, densityResult3$testScore)
                    nCopyPeaks_final <- nCopyPeaks_while
                  }
                }else{
                  loginfo('NOT reducing copy number by 1 because it is NOT safe to do so: no het score of 4N peak (proposed 3N) to test')
                  nCopyPeaks_final <- nCopyPeaks_while
                }
              }else{
                # check het score of the temp2N peak to make sure it is ok to be 1N, but only do this if it is a reliable test
                if(densityResult2$observ >= minObservations &  densityResult2$testScore >= heterozygosityScoreThreshold){
                  # this peak needs to stay 2N do not change
                  logwarn('not reducing copy number by 1 because it is not safe to do so: would result in a 1N peak with high het score. No new 3N peak to test')
                  nCopyPeaks_final <- nCopyPeaks_while
                }else{
                  logwarn('reducing copy number by 1 because it is safe to do so: no 0N peak will result and no new 3N peak to test')
                  nCopyPeaks_final <- nCopyPeaks_while-1
                }
              }

              # do some checks before increasing copy number....
            }else{
              if(nCopyPeaks_dig[1] < nCopyPeaks_while[1]){
                logwarn('NOT changing copy number, it was already increased in while loop') # increasing again would make first peak 3N, which would more likely be an error than doing nothing
                nCopyPeaks_final <- nCopyPeaks_while
              }else if(cn2existsTemp){
                if( densityResult2$observ >= 15 &
                    !is.na(densityResult2$testScore)){
                  # testing if it is ok to make the tested 3N peak the 4N peak,
                  # this test could fail if the testScore of the tested 3N peak (proposed 4N) has LOH, uncommon but maybe common enough to cause problems?
                  if((densityResult2$testScore+0.004) <= densityResult3$testScore){   # proposed 3N peak <= proposed 4N peak ?
                    # TODO: changed from 0.005 to 0.004 because 0.005 is too much for 26258, but I could fix it by forcingFirst digital peak if 0.005 is required for some other sample
                    logwarn('increase copy number by 1 because it is safe to do so: het score of 2N peak (proposed 3N) is less than (or equal to) tested 3N peak (proposed 4N): %.3f+0.004 <= %.3f',  densityResult2$testScore, densityResult3$testScore)
                    nCopyPeaks_final <- nCopyPeaks_while+1
                  }else{
                    loginfo('NOT increasing copy number by 1 because it is NOT safe to do so: het score of 2N peak (proposed 3N) is not less than (or equal to) tested 3N peak (proposed 4N): %.3f+0.004 !<= %.3f',  densityResult2$testScore+0.005, densityResult3$testScore)
                    nCopyPeaks_final <- nCopyPeaks_while
                  }
                }else{
                  if(cn4existsTemp){
                    if(densityResult4$observ >= 15){  # TODO what is a good number of observ for the temp cn4 peak?
                      loginfo('not enough data in 2N peak (proposed 3N) so check 4N peak (proposed 5N)')
                      nCopyPeaks_final <- testsForIncreasingCNof4Npeak()
                    }else{
                      loginfo('not enough data in 2N peak (proposed 3N) and not enough data (%i observ) in 4N peak (proposed 5N) to test',densityResult4$observ)
                      loginfo('thus NOT increasing copy number by 1 because it is NOT safe to do so: no het scores to test')
                      nCopyPeaks_final <- nCopyPeaks_while
                    }
                  }else if(!cn4existsTemp){
                    loginfo('not enough data in 2N peak (proposed 3N) and no 4N peak (proposed 5N) to test')
                    loginfo('thus NOT increasing copy number by 1 because it is NOT safe to do so: no het scores to test')
                    nCopyPeaks_final <- nCopyPeaks_while
                  }
                }
              }else if(cn4existsTemp){
                loginfo('no 2N peak (proposed 3N) so check 4N peak (proposed 5N)')
                # TODO: do i need to check length of valid 4N peak indexes?
                nCopyPeaks_final <- testsForIncreasingCNof4Npeak()
              }else{
                loginfo('NOT increasing copy number by 1, because it is NOT safe to do so: no proposed 3N or proposed 5N peak to test')
                nCopyPeaks_final <- nCopyPeaks_while
              }
            }
          }else{
            # when tumor percent is low, or coverage is low, het scores for odd peaks will creep past the cut off of 0.98,
            # in which case we will allow the 3N peak to go over the threshold and will not change copy number
            logwarn('tumor percent with 3N as main peak is less than 15%, too low to test so will not change copy number')
            nCopyPeaks_final <- nCopyPeaks_while
          }
        }else{
          #hetTestScoreFor3NPeak passes
          nCopyPeaks_final <- nCopyPeaks_while
        }
      }


      test3NStep <- repeatSteps(peakInfo,keyHetScoresPerPeak,nCopyPeaks_step=nCopyPeaks_final)
      peakInfo <- test3NStep$peakInfo
      percentTumor <- test3NStep$percentTumor
      hetScoreQuantiles <- test3NStep$hetScoreQuantiles
      if(percentTumor > allowedTumorPercent){
        if(!all(nCopyPeaks_final == nCopyPeaks_while)){
          logerror('undo change to copy number, percentTumor cannot be above 100')
          nCopyPeaks_final <- nCopyPeaks_while
          peakInfo <- whileLoopStep$peakInfo
          percentTumor <- whileLoopStep$percentTumor
          hetScoreQuantiles <- whileLoopStep$hetScoreQuantiles
        }
      }else{
        logdebug('percent tumor after testing 3N peak passes')
      }

    }


  }

  ## calc read depth values for 2N peak -------
  #   expReadsIn2NPeak_1bp
  #   rdNormX_2Npeak
  if(TRUE){
    diploidIndex <- which(nCopyPeaks_final==2) # find this index in dPeaks, digitalPeaks, which could skip subClonal peaks that are included in peakRD_1bp and peakRD_normX

    if(length(diploidIndex)==0){
      # when we specified/forced the first digital peak, or bumped up the peaks, we ended up with no 2N peak, so we need to extrapolate it
      loginfo('extrapolating 2N peak')

      bestPeaksResult <- getTwoBestPeakIndexes(peakCopyNum = peakInfo$nCopy, peakHeight = peakInfo$peakHeight)

      expReadsIn2NPeak_1bp <- extrapRD(cn = 2,
                                       cn1 = peakInfo[bestPeaksResult$bestPeakIndex, 'nCopy'],             cn2 = peakInfo[bestPeaksResult$secondBestPeakIndex, 'nCopy'],
                                       rd1 = peakInfo[bestPeaksResult$bestPeakIndex, 'peakReadDepth_1bp'], rd2 =peakInfo[bestPeaksResult$secondBestPeakIndex, 'peakReadDepth_1bp'])
      expReadsIn2NPeak_wsz <- expReadsIn2NPeak_1bp* wsz  # peakRD_1bp is value for 1bp window, need to convert for given wsz

      # extrapolate read depth for 2N level
      rdNormX_2Npeak <- extrapRD(cn = 2,
                                 cn1 = peakInfo[bestPeaksResult$bestPeakIndex, 'nCopy'],             cn2 = peakInfo[bestPeaksResult$secondBestPeakIndex, 'nCopy'],
                                 rd1 = peakInfo[bestPeaksResult$bestPeakIndex,'peakReadDepth_normX'],rd2 = peakInfo[bestPeaksResult$secondBestPeakIndex,'peakReadDepth_normX']
      )

    }else{
      # which peakRD_1bp, peakRD_normX is diploid peak
      normPeaks2Nindex <- which.min(abs(dPeaks[diploidIndex] - (peakInfo$peakReadDepth_normX * n00)))
      if(which(peakInfo$nCopy==2) != normPeaks2Nindex){
        logerror('peakInfo$peakReadDepth_normX and peakInfo are out of sync, normPeaks2Nindex = %i',normPeaks2Nindex)
      }
      expReadsIn2NPeak_1bp <- peakInfo[normPeaks2Nindex, 'peakReadDepth_1bp']

      if(diploidIndex==1){
        #we already have this value, good to go
        expReadsIn2NPeak_wsz <- expReadsInFirstDigitalPeak_wsz
      }else{
        # compute
        expReadsIn2NPeak_wsz <- peakInfo[normPeaks2Nindex, 'peakReadDepth_1bp']* wsz # peakRD_1bp is value for 1bp window, need to convert for given wsz
      }

      ### for adding copy number level, NRD
      rdNormX_2Npeak      <- peakInfo[ which(peakInfo$nCopy==2), 'peakReadDepth_normX']
    }
  }


  ### add NRD and copy number level ----------
  nrd <- pkmodToNRD(segmentData, peakInfo, rdNormX_2Npeak)
  segmentData[,'nrd'] <- nrd
  segmentData[,'cnLevel'] <- calcCopyNumber(NRD=nrd, # read depth normalized to 2N level
                                            tau= percentTumor/100)



  # these are the final indexes
  cn1index <- which(peakInfo$nCopy==1)
  cn2index <- which(peakInfo$nCopy==2)
  cn3index <- which(peakInfo$nCopy==3)
  cn4index <- which(peakInfo$nCopy==4)
  cn5index <- which(peakInfo$nCopy==5)



  ##########################'
  ## confirmation tests ----------------
  ##########################'
  # not using for tests just for printing warning messages
  # and adding additional hetScore results to peakInfo
  if(TRUE){

    ### het score clusters in first digital peak --------
    # test 1) if there are 2 or more loh clusters(modes) in first digital peak and max het score is above threshold, we expect this peak to be 2N or greater
    # test 2) if test score (max of the two clusters) is above threshold, we expect this peak to also be an even numbered copy level
    # NOTE: this was originally part of a test in the while loop but was redundant (more specific) than the current test

    ### test 1)
    firstPeakMustBeTwoNorMore <- hasNorMoreClusters(hetScoreDensityResult=densityFirstDigPeak,N=2,heterozygosityScoreThreshold=heterozygosityScoreThreshold)
    # first peak is allowed to be 1 no matter what 'firstPeakMustBeTwoNorMore' says, if forceFirstDigPeakCopyNum was set to "1"
    if(firstPeakMustBeTwoNorMore & forceFirstDigPeakCopyNum!=1){
      testFirstDigPeakPasses_1 <- ifelse(nCopyPeaks_final[1]>=2, TRUE, FALSE)
    }else{
      testFirstDigPeakPasses_1 <- TRUE
    }
    if(!testFirstDigPeakPasses_1){
      logwarn(sprintf('Confirmation failed: %i clusters in first digital peak, but it is not 2N or higher', densityFirstDigPeak$numHetClusters))
    }

    ### test 2)
    if(densityFirstDigPeak$testScore >= heterozygosityScoreThreshold){ # if the test score (the max of the two scores) is above cutoff
      testFirstDigPeakPasses_2 <- ifelse(nCopyPeaks_final[1]%%2==0, TRUE, FALSE)
    }else{
      testFirstDigPeakPasses_2 <- TRUE
    }
    if(!testFirstDigPeakPasses_2){
      logwarn('test score of first digital peak is above threshold, %.3f > %.2f,  but copy number is %i',densityFirstDigPeak$testScore, heterozygosityScoreThreshold, nCopyPeaks_final[1])
    }



    ### threeNpasses:  --------
    #  1) should have 2 or fewer het score clusters and
    #  2) het. score should be lower than 0.98, but...
    #     low tumor samples will have elevated het scores due to the presence of normal DNA,
    #    therefore test the right-mode mode from the density() rather than the max het score
    threeNpasses <- NA
    if(length(cn3index) > 0){
      # add hetScore to peak info
      testScoreForPeak <- hetScoreDensity(segmentsCloseToPeak,segmentData, index=cn3index,sampleId,alternateId, skipPlot = skipExtras, plotTextPrefix= paste('double check 3N peak:'),
                                          heterozygosityScoreThreshold=heterozygosityScoreThreshold)
      peakInfo[cn3index,'hetScore'] <- testScoreForPeak$testScore
      peakInfo[cn3index,'hetScoreObserv'] <- testScoreForPeak$observ

      # threeNpasses <- ifelse(hetScoreQuantiles[cn3index,'q0.85'] < 1.0, TRUE, FALSE)
      validPeakIndexes <- which(segmentsCloseToPeak[,cn3index] &
                                  segmentData$valid==1)
      if(length(validPeakIndexes) >= 15){
        threeNpasses <- TRUE
        peakHas3HetScoreClustersDC <- hasNorMoreClusters(testScoreForPeak,N=3,heterozygosityScoreThreshold )

        if(peakHas3HetScoreClustersDC){  # 3N peaks can not have more than 2 het score clusters
          threeNpasses <- FALSE
          logwarn("3N het score has too many clusters: %.3f > 2", testScoreForPeak$numHetClusters)
        }
        # check testScore
        if( testScoreForPeak$testScore >= heterozygosityScoreThreshold){
          threeNpasses <- FALSE
          logwarn("3N het score is high: %.3f >= %.2f (percentTumor: %.1f)",testScoreForPeak$testScore, heterozygosityScoreThreshold,percentTumor)
        }
      }
    }


    ### twoNpasses: --------
    # 2N peak should be one of the top three biggest peaks  # FAILS: 58124, 74001 but ok
    twoNpasses <- NA
    if(length(cn2index) > 0){
      twoNpasses <- ifelse(peakInfo[cn2index,'rankByHeight'] <=4, TRUE, FALSE)
      if(!twoNpasses){
        logwarn("twoNpasses=FALSE: 2N peak is not one of the top four prominent peaks: %i of %i",peakInfo[cn2index,'rankByHeight'],max(peakInfo[,'rankByHeight']))
      }
    }else{
      twoNpasses <- FALSE
      logwarn("twoNpasses=FALSE: 2N peak does not exist")
    }

    ### fourNpasses:  --------
    ### test 4N het scores
    ###     4N 2:2 should have higher het scores than 1N, 3N, 5N, but might be 3:1 segments rather than 2:2 i.e. 58109
    ###     therefore this is just for informational purposes.  FYI: 2N could be 2NcnLOH so don't check it,
    # TODO: what is the correct quantile to check?  58124 requires 85th, 75th will fail
    # TODO: instead of quantiles, check het score densities, assuming there are enough observations
    fourNpasses <- NA
    if(length(cn4index) > 0){

      # add hetScore to peak info
      testScoreForPeak <- hetScoreDensity(segmentsCloseToPeak,segmentData, index=cn4index,sampleId,alternateId, skipPlot = TRUE,
                                          heterozygosityScoreThreshold=heterozygosityScoreThreshold)
      peakInfo[cn4index,'hetScore'] <- testScoreForPeak$testScore
      peakInfo[cn4index,'hetScoreObserv'] <- testScoreForPeak$observ



      if( length(cn1index) > 0){
        fourNgt1N <- ifelse(hetScoreQuantiles[cn4index,'q0.85'] > hetScoreQuantiles[cn1index,'q0.85'], TRUE, FALSE)
      }else{fourNgt1N <- NA}

      if( length(cn3index) > 0){
        fourNgt3N <- ifelse(hetScoreQuantiles[cn4index,'q0.85'] > hetScoreQuantiles[cn3index,'q0.85'], TRUE, FALSE)
      }else{fourNgt3N <- NA}

      if( length(cn5index) > 0){
        # add hetScore to peak info
        testScoreForPeak <- hetScoreDensity(segmentsCloseToPeak,segmentData, index=cn5index,sampleId,alternateId, skipPlot = TRUE,
                                            heterozygosityScoreThreshold=heterozygosityScoreThreshold)
        peakInfo[cn5index,'hetScore'] <- testScoreForPeak$testScore
        peakInfo[cn5index,'hetScoreObserv'] <- testScoreForPeak$observ

        fourNgt5N <- ifelse(hetScoreQuantiles[cn4index,'q0.85'] > hetScoreQuantiles[cn5index,'q0.85'], TRUE, FALSE)
      }else{fourNgt5N <- NA}

      ### if any are false then fourNpasses <- FALSE
      if(all(c(fourNgt1N, fourNgt3N,fourNgt5N), na.rm=TRUE)){
        fourNpasses <- TRUE
      }else{
        logwarn('4N q85th het. score: %s is lower than 1 or more -odd- CN levels (1N, 3N and/or 5N)', round(hetScoreQuantiles[cn4index,'q0.85'],4)  )
        fourNpasses <- FALSE
      }
    }

  }

  # if(!skipExtras){
  par(op) # done with hetScore density plots
  # }

  if(FALSE){
    # add NRD value- this will replace the value that was there if it was loaded from a pre-run pipeline version), the previous value will now be labeled 'nrd_org'
    segmentDataAccum <- addNewNRD(wszNormalPeak = wsz, ploidyBasedNormalBin=expReadsIn2NPeak_1bp ,segmentation=segmentData)

    for(segmentId in seq_len(nrow(segmentDataAccum))) {  # keys=which(segmentDataAccum$chr==10); segmentDataAccum[keys,];  keys=which(segmentDataAccum$lohScoreMedian>1.05); segmentDataAccum[keys,]
      # segmentId =455
      segment <- segmentDataAccum[segmentId,]
      # TODO: get the cutoffs from pipeline
      if(segment[['nrd']] >= 1.85 & segment[['nrd']] <= 2.15){
        cnvLevel <- 2
        cnvState <- 2
      }else if(segment[['nrd']] < 1.85){
        cnvLevel  <-  max(0,round((segment[['nrd']]-2)/(percentTumor/100) + 2)) # can not be smaller than 0 # TODO:floor?
        cnvState <- 1
      }else{
        cnvLevel  <-  round((segment[['nrd']]-2)/(percentTumor/100) + 2) # TODO: ceiling?
        cnvState <- 3
      }
      # add segments to data.frame
      segmentDataAccum[segmentId, 'cnvState'] <- cnvState
      segmentDataAccum[segmentId, 'cnvLevel'] <- cnvLevel
    }

    table(segmentDataAccum$cnvLevel)
    table(segmentDataAccum$cnvState)

    considerPeakCutoff <- 0.02
    maxLevel <- max(segmentDataAccum$cnvLevel)

    maxHetDensityPeak <- data.frame()
    for(iLevel in 1:maxLevel){
      # iLevel <- 2
      iKeys <- which(segmentDataAccum$cnvLevel==iLevel &
                       segmentDataAccum$valid==1)
      if(length(iKeys) >= 20){ # TODO: what is an appropriate min number of observations
        hetToUse <- segmentDataAccum[iKeys, 'lohScoreMedian']

        bandwidth <- bw.nrd(hetToUse)
        hetDen <- density(hetToUse,bw <- bandwidth)
        plot(hetDen, main=paste('density of mean het scores for segments in cnvLevel:',iLevel))
        polygon(hetDen$x, hetDen$y, col='gray92',border='gray92')
        dx <- hetDen$x
        dy <- hetDen$y
        maxPeakY <- dy[which.max(dy)]
        maxPeakX <- dx[which.max(dy)]

        #Find the peaks
        tp <- pastecs::turnpoints(dy)

        #Get all of them that seem interesting (i.e. big enough)
        topConsidX <- dx[tp$pos[which(tp$peaks & tp$points > considerPeakCutoff*maxPeakY)]]
        topConsidY <- dy[tp$pos[which(tp$peaks & tp$points > considerPeakCutoff*maxPeakY)]]
        abline(v=topConsidX, lty='dotted', col='red')

        # second max
        max2  <-  function(x) max( x[-which.max(x)]) #// return a second max even if it happens to equal the first max, while the following code will not:  max( x[x!=max(x)])
        max2PeakY <-  max2(topConsidY)
        max2PeakX <- dx[which(dy==max2PeakY)]

        maxHetDensityPeak <- rbind(maxHetDensityPeak,
                                   data.frame(cnLevel=iLevel, n=length(iKeys), maxPeakX,maxPeakY, max2PeakX,max2PeakY))
      }
    }
  }


  ### 2D plot -------
  # do last, after all the density plots
  plotHetScoreVsPeakMode(peakInfo,segmentData,sampleId, alternateId, percentTumor,
                         segmentCountText,digitalPeakZone,keyHetScoresPerPeak,segmentDataSizes,
                         heterozygosityScoreThreshold,
                         origMaxPercentCutoffManual,grabDataPercentManual ,minPeriodManual,maxPeriodManual)




  if(interactive() & pause){BBmisc::pause()}

  if (!noPdf) {
    # dev.off() # now executed by on.exit()
    loginfo("wrote file: %s", ploidyPdfFile@path)
  }

  #############################################'
  ## check hetScores are valid -----------
  #############################################'
  testScores <- NULL
  dPeakIndexes <- c(1:nrow(peakInfo))[!is.na(peakInfo$dPeaks)]
  for( i in dPeakIndexes){
    # i=1
    hsdResult <- hetScoreDensity(segmentsCloseToPeak,segmentData, index=i,sampleId,alternateId, skipPlot = T,
                                 heterozygosityScoreThreshold=heterozygosityScoreThreshold)
    testScores <- c(testScores,hsdResult$testScore)

    # TODO: replace this code below which is scattered about, using the stuff above
    # # add hetScore to peak info
    # testScoreForPeak <- hetScoreDensity(segmentsCloseToPeak,segmentData, index=cn3index,sampleId,alternateId, skipPlot = skipExtras, plotTextPrefix= paste('double check? 3N peak'),
    #                                     heterozygosityScoreThreshold=heterozygosityScoreThreshold)
    # peakInfo[cn3index,'hetScore'] <- testScoreForPeak$testScore
    # peakInfo[cn3index,'hetScoreObserv'] <- testScoreForPeak$observ


  }
  if(!any(testScores>=heterozygosityScoreThreshold)){
    logerror('none of the digital peaks have a test score above the threshold (%s): %s;
             suspect something wrong with the raw data or a critial peak was skipped',
             heterozygosityScoreThreshold, paste( testScores, collapse = ", "))
  }



  mainPeakIndex <- which(peakInfo$rankByHeight==1)
  loginfo('Main peak is %sN',peakInfo[mainPeakIndex,'nCopy'])

  ## return output -------------
  ploidyOutput <- list(expReadsIn2NPeak_1bp=expReadsIn2NPeak_1bp,
                       percentTumor=round(percentTumor,1),
                       peakInfo=peakInfo,
                       segmentedData=segmentData,
                       # hetScoreQuantiles=hetScoreQuantiles,
                       iterationStatsAll=iterationStatsAll)

  loginfo('END OF FUNCTION:calculatePloidy')
  return(ploidyOutput)
}

## Plot segments-Het. score vs peak mode --------------------------
plotHetScoreVsPeakMode <- function(peakInfo,segmentData,sampleId, alternateId, percentTumor,
                                   segmentCountText,digitalPeakZone,keyHetScoresPerPeak,segmentDataSizes,
                                   heterozygosityScoreThreshold,
                                   origMaxPercentCutoffManual,grabDataPercentManual ,minPeriodManual,maxPeriodManual){

  # segmentData=result$segmentData; percentTumor=result$percentTumor;  peakInfo= result$peakInfo

  mainPeakIndex <- which(peakInfo$rankByHeight==1)
  rdNormX_2Npeak      <- peakInfo[ which(peakInfo$nCopy==2), 'peakReadDepth_normX']
  firstDigPeakIndex <- which.min(peakInfo$dPeaks) # may not be 1 if the grid is not forced to start at the first peak
  normXofFirstDigPeak <- peakInfo$peakReadDepth_normX[firstDigPeakIndex] # 1 if first digital peak is first peak, but get this value in case it isn't
  # translate in case the first digital peak is not the first peak
  peakModes <- peakInfo[, 'peakReadDepth_normX']/normXofFirstDigPeak  # normXofFirstDigPeak will be 1 if first peak is first digital peak (dPeak)





  firstDigPeakKey <- which.min(peakInfo$dPeaks) # may not be 1 if the grid is not forced to start at the first peak
  logdebug('firstDigPeakKey %i ',firstDigPeakKey)

  rdNormX_firstDigPeak <- peakInfo[firstDigPeakKey, 'peakReadDepth_normX']  # index = firstDigPeakIndex
  numDigPeaks = sum(!is.na(peakInfo$dPeaks))

  ### one color for each chromosome
  segColors = getSegmentColors()

  rowsToView <- 1:nrow(segmentData) # used for debugging, to look at a few rows at a time, ie 1:17 for chr1 in sample...
  maxXplot <- max(segmentData$lohScoreMedian[rowsToView],1)

  #sometime the coverage is so low, there is little to no hetScore data
  nValidHetScores <- sum(segmentData$lohScoreMedian>0)
  # TODO: 100 is a heuristic, maybe should be based on a percentage of total segments?
  # need to avoid failing when there are very few data points.
  if(nValidHetScores<100){
    minXplot <- 0
    legendPlacement <- "topright"
    minYplot <- quantile(segmentData$pkmod,probs = .005) # range(segmentData$pkmod)[1]
    maxYplot <- quantile(segmentData$pkmod,probs = .995) # range(segmentData$pkmod)[2]
  }else{
    minXplot <- 0.08
    legendPlacement <- "topleft"

    # use some data to figure out how big to make the y axis, but be careful not to make it too big
    if(numDigPeaks == 1){
      # do not know tumor percentage if there is only one peak
      atMaxCN <- NULL
    }else{
      maxCN <- max(peakInfo$nCopy,na.rm = TRUE)+2 # expand for plotting
      nrdMaxCN <- 2+((maxCN-2)* percentTumor/100)  # these values are in NRD coordinates, need to convert 1st dig Peak read depth
      atMaxCN <- nrdMaxCN/2*(rdNormX_2Npeak/rdNormX_firstDigPeak) # rdNormX_2Npeak and rdNormX_firstDigPeak are calculated in calculatePloidy, pkmodToNRD
    }

    yQuantile98 <- quantile(segmentData$pkmod[segmentData$lohScoreMedian>0],probs = .98)
    maxYplot <- min(c(5,                         # y axis max will be between 1.4 - 5 pkmod
                      max( c(1.4, yQuantile98*1.2, atMaxCN), na.rm = TRUE)) )
    minYplot <- min(c(0.7, min(segmentData$pkmod[segmentData$lohScoreMedian>0])))  # y axis min will be at least 0.7
  }

  segmentPloidy <- sum(segmentData$size * segmentData$cnLevel)/sum(segmentData$size) # weighted by length of segment

  graphics::layout(1)
  op <- par(mar=c(5.5, 4, 3.5, 3.5),mgp=c(1.5, 0.5,0))
  plot(segmentData$lohScoreMedian[rowsToView], segmentData$pkmod[rowsToView], type='n',
       xlim=c(minXplot,maxXplot),
       ylim= c(minYplot,maxYplot),
       xlab="Heterozygosity Score",
       ylab="read depth per segment normalized by 1st digital peak",
       main=segmentCountText)
  if(!is.null(alternateId)){
    mtext(3, text=c(sampleId, alternateId),adj=c(0,1))
  }else{
    mtext(3, text=c(sampleId),adj=c(0))
  }

  # annotate normalized peak locations
  # add boxes to show the zones around the peaks
  rect(xleft  = 0, xright = par('usr')[2],
       ybottom = peakInfo[!is.na(peakInfo$dPeaks), 'peakReadDepth_normX']/normXofFirstDigPeak - digitalPeakZone,
       ytop    = peakInfo[!is.na(peakInfo$dPeaks), 'peakReadDepth_normX']/normXofFirstDigPeak + digitalPeakZone,
       border = NA,  col="#f9edfa")
  box()
  mtext(1, text=paste('segment ploidy:', round(segmentPloidy,1)), adj = 0, line=-2, cex=.9)


  # add points overtop the zone boxes
  points(segmentData$lohScoreMedian[rowsToView], segmentData$pkmod[rowsToView],
         pch=segmentData$chr[rowsToView],
         col=scales::alpha(segColors[segmentData$chr[rowsToView]],.5))


  # identify scores that are not valid
  notValidRows <- segmentData$valid==0
  points(segmentData$lohScoreMedian[notValidRows], segmentData$pkmod[notValidRows],
         pch=1,col='orange', cex=1.75)

  ### add lines -----------------'
  ## vertical lines to mark cutoff
  abline(v=heterozygosityScoreThreshold, col='red', lty='dotted')
  mtext(heterozygosityScoreThreshold, side=1, at = heterozygosityScoreThreshold, cex=.65, line = -1, col='red')

  ## horizontal lines to mark peaks
  # abline(h=peakInfo[!is.na(peakInfo$dPeaks), 'peakReadDepth_normX']/normXofFirstDigPeak, col='gray', lty='dotted')     #  on the digital grid
  abline(h=peakInfo[ is.na(peakInfo$dPeaks), 'peakReadDepth_normX']/normXofFirstDigPeak, col='gray', lty='dotted')    #  NOT on the digital grid

  ## horizontal lines to mark boundary of first digital peak
  abline(h= peakInfo[firstDigPeakKey, 'peakReadDepth_normX']/normXofFirstDigPeak - digitalPeakZone ,col="red", lwd=1.3)
  # abline(h= peakInfo[which(digitalPeakIndex==firstDigPeakKey) , 'peakReadDepth_normX']/normXofFirstDigPeak - digitalPeakZone ,col="red", lwd=1.5)
  abline(h= peakInfo[firstDigPeakKey, 'peakReadDepth_normX']/normXofFirstDigPeak + digitalPeakZone ,col="red", lwd=1.3)

  ## vertical line to mark 3N peak hetScore mode
  # NOTE: do Not use 'hetTestScoreFor3NPeak' from calculatePloidy, does not get updated if there is a decrease or increase to copy number during test3Npeaks
  threeNPeakKey <- which(peakInfo$nCopy==3)
  if(length(threeNPeakKey)==1){
    hetScore3NPeak <- peakInfo[threeNPeakKey, 'hetScore']
    abline(v=hetScore3NPeak, col='gray40', lty='dashed')
    mtext(paste('3N:',round(hetScore3NPeak,3)), side=1, at = hetScore3NPeak, cex=.65, line = -.25, col='gray40')
  }



  ### add annotations -----------------'
  # peak rank labels
  ymin <- par('usr')[3]
  ymax <- par('usr')[4]
  colorCodeRank <- ifelse(is.na(peakInfo$dPeaks),'gray40', 'black')
  mtext(text=peakInfo$rankByHeight, side=4, at=peakInfo$peakReadDepth_normX/normXofFirstDigPeak, las=2, line=-.25, col=colorCodeRank, adj=1, cex=.7)
  # mtext(text='Rank', side=4, at=ymin*1.1, col='black', las=1, line=-1.5, cex=.8, col='gray') # bottom
  mtext(text='Peak', side=4, at=ymax*.98, col='gray40', las=1, line=-1.7, cex=.75)
  mtext(text='Rank', side=4, at=ymax*.96, col='gray40', las=1, line=-1.7, cex=.75) # top

  legend(legendPlacement, legend=c(1:22), col=segColors, pch=1:22, cex=1.1) # do legend after lines so the legend is on top
  mtext(paste('pink zone: +/-',digitalPeakZone), side=1, adj=0, line=.75, cex=.8, col='hotpink')
  # mtext(paste('cyan: 75th'), side=1, adj=0, line=0, cex=.8, col='cyan')
  # mtext(paste('magenta: 85th'), side=1, adj=0, line=.75, cex=.8, col='magenta')

  ## segment Data stats
  mtext(text=paste('Initial segment min:', round(segmentDataSizes$minReasonableSegmentSizeFinal/1e6, 1),'Mb'), side=1, adj=0, line=1.5, cex=.8) # cnv segment size before breaking up into smaller chunks
  mtext(text=paste('min:', segmentDataSizes$minSegment,'Mb'), side=1, adj=0, line=2.25, cex=.8)
  mtext(text=paste('max:', segmentDataSizes$maxSegment,'Mb'), side=1, adj=0, line=3,   cex=.8)
  mtext(text=paste('mean:',segmentDataSizes$meanSegment,'Mb'),side=1, adj=0, line=3.75, cex=.8)

  ## label peak CN
  if(numDigPeaks == 1){
    mtext(text='CN', side=4, at=ymax*.99, col='black', las=1, line=1)
    colorCodeCN <- ifelse(is.na(peakInfo$dPeaks),'gray', 'black')
    axis(side=4, at=peakInfo[!is.na(peakInfo$dPeaks), 'peakReadDepth_normX']/normXofFirstDigPeak, labels = FALSE, tick=TRUE) # tick marks only
    mtext(text=paste0(peakInfo$nCopy,"N"), side=4, at=peakInfo$peakReadDepth_normX/normXofFirstDigPeak, las=1, line=1, col=colorCodeCN) #, adj=1
  }else{
    # CN level annotations based on Tumor percent
    maxCN <- max(20,max(peakInfo$nCopy,na.rm = TRUE)) # expand for plotting
    tau <- percentTumor/100
    iCN <-0:maxCN
    nrdCN <- 2+((iCN-2)*tau)  # these values are in NRD coordinates, need to convert 1st dig Peak read depth
    cnAt <- nrdCN/2*(rdNormX_2Npeak/rdNormX_firstDigPeak) # rdNormX_2Npeak and rdNormX_firstDigPeak are calculated in calculatePloidy, pkmodToNRD

    axis(side=4, at=cnAt, labels=paste0(iCN,'N'), tick=TRUE, las=1)
  }

  mtext(paste0('tumor: ', round(percentTumor,0), '%'),side = 1, adj=1, line=1.9)
  if(forceFirstDigPeakCopyNum >=0){
    mtext(1, text=paste0('manual 1st digital Peak: ', forceFirstDigPeakCopyNum, 'N'),adj=1, line=2.5, cex=.7, col=2)
  }
  if(origMaxPercentCutoffManual > 0){
    mtext(1, text=paste0('manual peak height % cutoff: ', origMaxPercentCutoffManual),adj=1, line=3.0, cex=.7, col=2)
  }
  if(grabDataPercentManual > 0){
    mtext(1, text=paste0('manual main peak width factor: ', grabDataPercentManual),        adj=1, line=3.5, cex=.7, col=2)
  }
  if(minPeriodManual > 0){
    mtext(1, text=paste0('manual minPeriod: ', minPeriodManual),               adj=1, line=4.0, cex=.7, col=2)
  }
  if(maxPeriodManual > 0){
    mtext(1, text=paste0('manual maxPeriod: ', maxPeriodManual),               adj=1, line=4.0, cex=.7, col=2)
  }

  par(op)
} # CLOSE: plotHetScoreVsPeakMode()
