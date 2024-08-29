#' Find the NRD of the main Peak in the ---cnvBinned data---
#' Will equal 2 if main peak is the normal peak
#' @param calcPloidyResult object returned from \code{calculatePloidy}
#' @keywords internal
getMainPeakNRD=function(calcPloidyResult){
  expReadsIn2NPeak_1bp= calcPloidyResult$expReadsIn2NPeak_1bp
  mainPeakKey=which(calcPloidyResult$peakInfo$rankByHeight==1)
  mainPeakReadDepth_1bp = calcPloidyResult$peakInfo[mainPeakKey,'peakReadDepth_1bp']
  mainPeakNRD = 2*(mainPeakReadDepth_1bp / expReadsIn2NPeak_1bp)
  # loginfo('mainPeakNRD = %.3f', mainPeakNRD)
  return(mainPeakNRD)
}

#' Find the NRD of the diploid peak from the binned read depth data
#' @param calcPloidyResult object returned from \code{calculatePloidy}
#' @keywords internal
getDiploidPeakNRD=function(calcPloidyResult){
  mainPeakKey=which(calcPloidyResult$peakInfo$rankByHeight==1)
  rdNormX_Mainpeak = calcPloidyResult$peakInfo[mainPeakKey,'peakReadDepth_normX']

  diploidPeakKey=which(calcPloidyResult$peakInfo$nCopy==2)
  rdNormX_2Npeak = calcPloidyResult$peakInfo[diploidPeakKey,'peakReadDepth_normX']

  mainPeakNRD=getMainPeakNRD(calcPloidyResult)
  diploidPeakNRD <-round( mainPeakNRD*rdNormX_2Npeak/rdNormX_Mainpeak, 3)
  # loginfo('diploidPeakNRD: %.3f',diploidPeakNRD)  #dipVal => diploidPeakNRD
  return(diploidPeakNRD)
}

#' determine if the peak has a minimum number of clusters
#'
#' @param hetScoreDensityResult hetScore values returned from density function
#' @param N minimum number of clusters
#' @param minObservations required number of input values to determine number of clusters present
#'
#' @inheritParams commonParameters
#' @keywords internal
hasNorMoreClusters <- function(hetScoreDensityResult, N, heterozygosityScoreThreshold, minObservations=20){
  # hetScoreDensityResult=densityFirstDigPeak
  if(hetScoreDensityResult$observ > minObservations){
    if(!is.na(hetScoreDensityResult$numHetClusters) && hetScoreDensityResult$numHetClusters >= N){
      if(N==2){
        if(hetScoreDensityResult$testScore < .45){  #19063 =.508 because there is a clone at .508  TODO: what is the 2:0 value when tu% is 99% - it is close to .45 but is .51 too much?
          test <- FALSE # the max testScore is too low, if there are multiple clusters, they are subclones
        }else{
          test <- TRUE
        }
      }else if(N==3){
        if(hetScoreDensityResult$testScore >= heterozygosityScoreThreshold){  # not a requirement for N=2 because 3N could have two peaks but won't be above the threshold
          test <- TRUE
        }else{
          test <- FALSE # testScore is too low
        }
      }else{
        test <- TRUE
      }
    }else{
      test <- FALSE # incorrect # clusters
    }
  }else{
    test  <-  FALSE # too few observations to know for sure
  }
  return(test)
}



hetScoreDensity <- function(segmentsCloseToPeak,segmentData, index,sampleId,folderId, plotTextPrefix=NULL,skipPlot=FALSE, minObserv=15, heterozygosityScoreThreshold){
  # plotTextPrefix=NULL;skipPlot=FALSE; minObserv=15
  validPeakIndexes <- which(segmentsCloseToPeak[,index] &
                              segmentData$valid==1)
  lohScores <- segmentData[validPeakIndexes,'lohScoreMedian']
  if(is.null(plotTextPrefix)){
    plotTextPrefix  <-  paste('peak',index)
  }else{
    plotTextPrefix <- paste(plotTextPrefix, ' peak',index)
  }
  # use max score if not enough data points, otherwise use mode # TODO: what is 'enough'
  densityResult <- getModeOrMaxScore(dataIn = lohScores,plotTextPrefix = plotTextPrefix, sampleId=sampleId,folderId=folderId,skipPlot=skipPlot,
                                     minObserv=minObserv, heterozygosityScoreThreshold=heterozygosityScoreThreshold  )
  return(densityResult)
}


#' Calculate density of loh scores
#'
#' If enough data ( >minObserv ) return testScore=right-most mode that is big enough (> consider PeakCutoff * max peak) and numHetClusters > 20% max peak
#' If not enough data return max loh score and numHetClusters =NA
#' @keywords internal
#' @param dataIn hetScores mean hetScores for selected copy number segments
#' @param considerPeakCutoff peaks must be bigger than this percentage of maxPeak to be considered # TODO: should this be adjusted higher when there are more data points, ie >100 or 125? 26297 requires >=0.06 not 0.05 which is what I've done all my testing on
#' @param countPeakCutoff peaks must be bigger than this percentage of maxPeak to be counted for numClusters
#' @param plotTextPrefix text to add to front of plot title
#' @param default default return value
#' @param skipPlot logical to control making density plot
#' @param minObserv minimum number of hetScores to compute the mode of a peak, otherwise the max value is used
#' @param folderId alternate id
#' @inheritParams commonParameters
#'
getModeOrMaxScore <- function(dataIn, considerPeakCutoff=0.06, countPeakCutoff=0.25, plotTextPrefix=NULL,
                              folderId=NULL, sampleId=NULL, skipPlot=FALSE, minObserv=15,default=0, heterozygosityScoreThreshold) {
  # dataIn=lohScores; considerPeakCutoff=0.06; countPeakCutoff=0.25; minObserv=15;skipPlot=FALSE;default=0

  # TODO: what is the minimum for minObserv 10? 15? 20?
  data <- dataIn[dataIn <= 1] # do not take density with data greater than 1 aka noisy data
  observ <- length(data)

  if (observ >= 2) {
    xlimMin  <-  min( min(dataIn)*0.98, 0.6) # make sure the scale isn't too zoomed in which can be distracting and misleading
    xlimits  <-  c(xlimMin, max(1,dataIn)) # specify max to make sure the plot includes the data that was removed

    #  bw.SJ(data) : need at least 2 data points

    # which bandwidth to use?
    # bandwidth <- bw.ucv(data); bwMethod <- 'bw.ucv'                  # this option seems to do the least smoothing
    # bandwidth <- bw.nrd(data); bwMethod <- 'bw.nrd' # PT58126 with nrd first peak only has 1 density peak but with nrd0 first peak has 2 which is not good
    # bandwidth <- bw.nrd0(data); bwMethod <- 'bw.nrd0'
    # bandwidth <- bw.bcv(data); bwMethod <- 'bw.bcv'                    # this option seems to do the MOST smoothing
    # bandwidth <- bw.SJ(data); bwMethod <- 'bw.SJ' # defaults to 'bw.SJ-ste'
    # bandwidth <- bw.SJ(data,method='ste'); bwMethod <- 'bw.SJ-ste' # this option seems to do the 2nd least smoothing
    bandwidth <- bw.SJ(data,method='dpi'); bwMethod <- 'bw.SJ-dpi'

    hetDen <- density(data,bw=bandwidth)
    maxPeakX <- hetDen$x[which.max(hetDen$y)]
    maxPeakY <- hetDen$y[which.max(hetDen$y)]

    #Find the peaks
    tp <- pastecs::turnpoints(hetDen$y)

    if(FALSE){
      summary(tp)
      tp$points # =hetDen$y[tp$pos]; data points minus the exaequos positions
      tp$pos    # =index into hetDen$y
      tp$tppos  # index of turning points in original data (hetDen$y)
      # plot using two different x axis scales
      plot(hetDen$y, type='l');
      points(tp$tppos,hetDen$y[tp$tppos],pch='*')

      plot(hetDen,xlim=xlimits, main=paste(plotTextPrefix,'het score density via', bwMethod),type='l');rug(dataIn)
      points(hetDen$x[tp$tppos],hetDen$y[tp$tppos],pch='*' )
    }


    # some stuff for plotting the densities with different bandwidth methods to visualize its affect on number of peaks
    # see 58089 lohScores of peak 4 which has 3 clusters, "/research/labs/experpath/vasm/shared/NextGen/johnsonsh/Routput/58089_peak4HetScores.Rdata"
    # nrd: only finds 2 peaks, the max being between the right 2 clusters,
    # SJ-dpi: only finds two but the max is the far right of the right 2 clusters
    # SJ-ste: finds three clusters
    if(FALSE){
      op <- par(mfrow=c(3,1),mar=c(3, 3.5, 1.5, 1),mgp=c(1.5, 0.5,0))
      bwMethodsAll <-  c( "ucv", "nrd","nrd0","bcv", "SJ-dpi",  "SJ-ste")
      for(bwMethod in bwMethodsAll ){
        # bwMethod <- bwMethodsAll[2]
        hetDen <- density(data,bw=bwMethod)
        maxPeakY <- hetDen$y[which.max(hetDen$y)]
        # Calculate turning points for this series
        hetDenY.tp <- turnpoints(hetDen$y)
        summary(hetDenY.tp)
        # plot(hetDenY.tp)
        #
        ## Add envelope and median line to original data
        # plot(hetDen$y, type = "l", main=bwMethod)
        # lines(hetDenY.tp)
        ## Note that lines() applies to the graph of original dataset

        ## plot turning point object
        # plot(hetDenY.tp)

        plot(hetDen,main=bwMethod ); rug(dataIn)
        # Get scores for peaks only
        #peakScores <- hetDen$x[pastecs::extract(hetDenY.tp, no.tp = FALSE, peak = TRUE, pit = FALSE)]

        ## Get all the scores that seem interesting (i.e. big enough)
        peakConsidX <- hetDen$x[hetDenY.tp$pos[which(hetDenY.tp$peaks & hetDenY.tp$points > considerPeakCutoff*maxPeakY)]]
        pitConsidX   <- hetDen$x[hetDenY.tp$pos[which(hetDenY.tp$pit  & hetDenY.tp$points > considerPeakCutoff*maxPeakY)]]

        abline(v=peakConsidX, lty='dotted', col='darkgray');
        abline(v=pitConsidX, lty='dotted', col='darkgray');

        abline(h=considerPeakCutoff*maxPeakY,col='darkgrey')
        abline(h=countPeakCutoff*maxPeakY,col='darkgrey')
        mtext(round(peakConsidX,3), side=3, at=peakConsidX, cex=0.7, line=-1)
      }
      par(op)
    }

    if(FALSE){
      data <- firstPeakLohScores
      # data <- lohScores # third peak

      plot(density(data),xlim=xlimits);       rug(dataIn)
      lines(density(data, bw = "nrd"), col = 2)
      lines(density(data, bw = "ucv"), col = 3)
      lines(density(data, bw = "bcv"), col = 4)
      lines(density(data, bw = "SJ-ste"), col = 5)
      lines(density(data, bw = "SJ-dpi"), col = 6)
      legend('topleft',legend = c("nrd0", "nrd", "ucv", "bcv", "SJ-ste", "SJ-dpi"),
             col = 1:6, lty = 1)
    }

    if(FALSE){
      bwMethodsAll <-  c("nrd0", "nrd", "ucv", "bcv", "SJ-ste", "SJ-dpi")
      plot(density(data,bw=bwMethodsAll[1]))
      rug(dataIn)
      for(i in 1:length(bwMethodsAll)){
        bwMethod <- bwMethodsAll[i]
        print(bwMethod)
        hetDen <- density(data,bw=bwMethod)
        hetDenY.tp <- pastecs::turnpoints(hetDen$y)
        lines(hetDen, col = i)
        peakScores <- hetDen$x[pastecs::extract(hetDenY.tp, no.tp = FALSE, peak = TRUE, pit = FALSE)]
        print(peakScores)
        abline(v=peakScores, col=i);
        # mtext(round(peakScores,3), side=3, at=peakScores, cex=0.5)
      }

      legend('topleft',legend = bwMethodsAll,col = 1:6, lty = 1, cex=0.9)
    }


    #Get all peaks that seem interesting (i.e. big enough)
    topConsidX <- hetDen$x[tp$pos[which(tp$peaks & tp$points > considerPeakCutoff*maxPeakY)]]
    topConsidY <- hetDen$y[tp$pos[which(tp$peaks & tp$points > considerPeakCutoff*maxPeakY)]]

    (peakScoresX  <- hetDen$x[pastecs::extract(tp, no.tp = FALSE, peak = TRUE, pit = FALSE)])
    (peakScoresY  <- hetDen$y[pastecs::extract(tp, no.tp = FALSE, peak = TRUE, pit = FALSE)])

    # density height cut off for counting peaks
    # countPeakCutoff=0.25
    topCountX <- hetDen$x[tp$pos[which(tp$peaks & tp$points > countPeakCutoff*maxPeakY)]] # 20 was too low for 58047 while using bw.SJ
    topCountY <- hetDen$y[tp$pos[which(tp$peaks & tp$points > countPeakCutoff*maxPeakY)]]

    (pitScoresX  <- hetDen$x[pastecs::extract(tp, no.tp = FALSE, peak = FALSE, pit = TRUE)])
    (pitScoresY  <- hetDen$y[pastecs::extract(tp, no.tp = FALSE, peak = FALSE, pit = TRUE)])

    # instead of counting the number of peaks to get the number of clusters, check that there is a deep enough pit(saddle) between the peaks
    # valid pits are those that 1) are between two peaks that are far enough apart to be legit
    #                           2) are high enough (countPeakCutoff*maxPeakY) i.e. 25% of the maxPeak,
    #                           3) but have a deep enough saddle, aka dip low enough between its adjacent peaks
    #
    validPits <- NULL
    if(length(topCountX)>1){
      for(i in 2:length(topCountX)){
        peakToPeakDistance <- topCountX[i]-topCountX[i-1]
        if(peakToPeakDistance>0.028){  # 26243 needs >0.0177; 19019>=0.023; 19033>0.027 # TODO: what is the correct value? don't go smaller than 3:1 and 2:2 clusters for 4N at low tumor which is probably around .04
          # find all the pits between two top peaks
          pitsBetween <- pitScoresX[data.table::between(pitScoresX,topCountX[i-1],topCountX[i])]
          if(length(pitsBetween)>0){
            # get the min (y) pit
            testPitY <- min(pitScoresY[(pitScoresX %in% pitsBetween)])
            # test the y value
            if(abs(topCountY[i-1] - testPitY) > topCountY[i-1]*0.27 &  # TODO: what is the correct value, 26% too small: 66291,58149;58172
               abs(topCountY[i]   - testPitY) > topCountY[i]*0.27){
              iValidPits <- pitScoresX[which(pitScoresY== testPitY)][1] # take just the first one if there are multiple # ie LU19091 peak 3
              validPits <- c(validPits,iValidPits)
            }
          }else{
            next
          }
        }else{
          next
        }

      }
    }


    ## second max
    # max2 = function(x) max( x[-which.max(x)]) #// return a second max even if it happens to equal the first max, while the following code will not:  max( x[x!=max(x)])
    # max2PeakY = max2(topConsidY)
    # max2PeakX <- hetDen$x[which(hetDen$y==max2PeakY)]

    # define test score
    if (observ >= minObserv) {
      # Enough for density
      # the right most mode of the top peaks
      testScore <- max(topConsidX);
      numHetClusters <- length(validPits) + 1 ; # is length(validPits) + 1    better than    numHetClusters = length(topCountX)
      plotCol <- 'red'
    }else{
      # first Max less than 1 must suffice
      testScore <- max( data[which(data < 1)] );
      # loginfo('%i observations for %s, too few for density, returning the first max<1: %.3f',observ,plotTextPrefix,testScore)
      numHetClusters <- NA_integer_; plotCol <- 'orange'
    }



    #plot density, modes
    if(!skipPlot){
      plot(hetDen,xlim=xlimits, main=paste(plotTextPrefix,'het score density via', bwMethod),type='l')
      # points(hetDen$x[tp$tppos],hetDen$y[tp$tppos],pch='*', col=ifelse(hetDen$x[tp$tppos] %in% pitScoresX, 'red','blue') )
      # polygon(hetDen$x, hetDen$y, col='gray92',border='gray92')
      rug(dataIn)
      abline(v=heterozygosityScoreThreshold, col='gray', lwd=0.5)

      if(!is.null(sampleId)){mtext(text=sampleId,3, adj=0)};
      if(!is.null(folderId)){mtext(text=folderId,3, adj=1)}

      abline(h=considerPeakCutoff*maxPeakY, col='gray')
      mtext(text = paste0(considerPeakCutoff*100,'%'), side=2, line=-2, at = considerPeakCutoff*maxPeakY, cex=.75, col='gray', las=2,padj = -1 )
      abline(h=countPeakCutoff*maxPeakY, col='gray')
      mtext(text = paste0(countPeakCutoff*100,'%'), side=2, line=-2, at = countPeakCutoff*maxPeakY, cex=.75, col='gray', las=2,padj = -1 )

      abline(v=topConsidX, lty='dotted')                # top peaks (modes)
      abline(v=testScore, lty='dotted',  col=plotCol)   # testScore to use/output
      if(!is.null(validPits))mtext('*', side=3, at=validPits, line=-1, col='green3')         # mark valid pits
      mtext(round(topConsidX,3), side=3, at=topConsidX, line=-.75, cex=.7)
      mtext(text=round(testScore,3), side=1, at=testScore, cex=.75, col=plotCol, line=-0.2)
    }
  }else{
    testScore <- default;
    numHetClusters <- NA_integer_;
  }

  # limit test score to 3 sig figs more is not relevant. See ME26301: peak2 vs peak3... 0.9907994 vs 0.9909797 these should be considered the same value
  return(list(testScore=round(testScore, 3), numHetClusters=numHetClusters, observ=observ))

}





#' find tumor fraction from the ratio between two copy number peaks and their corresponding read depth
#'
#' @param rd1,rd2  read depth of two peaks
#' @param cn1,cn2  copy number  of two peaks
#' @examples calcTumorRatio(rd1=1.2,rd2=2,cn1=2,cn2=4)
#' @examples calcTumorRatio(rd1=2, rd2=3.5, cn1=1,cn2=2) # 58163 = .155
#' @family copyNumberCalcs
#' @export
#' @keywords internal
calcTumorRatio=function(rd1,rd2,cn1,cn2){
  D <- abs(rd1-rd2)/abs(cn1-cn2) # expected difference between two peaks  copy number, i.e. 1N and 2N
  T <- 2*D/(rd1+D*(2-cn1))

  # loginfo('tumorRatio: %.2f', T)
  return(T)
}


#' find tumor fraction from copy number and normalized read depth
#'
#' @param cn copy number
#' @param NRD normalized read depth
#' @examples calcTau(NRD=2.9,cn=3)
#' @family copyNumberCalcs
#'
#' @export
#' @keywords internal
calcTau <- function(cn, NRD){
  tau  <-  (NRD-2)/ (cn-2)
  return(tau)
}


#' find copy number from read depth and tumor fraction
#'
#' @param NRD normalized read depth
#' @param tau tumor fraction
#' @examples calcCopyNumber(NRD=23,tau=0.85)
#' @examples calcCopyNumber(NRD=1.4,tau=0.33)
#' @examples calcCopyNumber(NRD=4,tau=0.75)
#' @family copyNumberCalcs
#'
#' @export
#' @keywords internal
calcCopyNumber <- function(NRD, tau){
  #nrd=2+((iCN-2)*tau)
  cn <- 2+ ((NRD-2)/tau)
  return(round(cn,1))
}


#' find normalized read depth from copy number and tumor ratio
#'
#' @param cn copy number
#' @param tau tumor fraction
#' @examples calcNrd(cn=3,tau=.7)
#' @family copyNumberCalcs
#' @export
#' @keywords internal
calcNrd <- function(cn, tau){
  nrd <- 2+((cn-2)*tau)
  return(round(nrd,3))
}

#' find read depth per window size from NRD
#' @param nrd normalized read depth
#' @param wsz window size
#' @param expReadsIn2NPeak_1bp expected number of reads in a 1 bp bin for the diploid peak
calcRD=function(nrd, wsz, expReadsIn2NPeak_1bp){
  rd=(nrd/2)*wsz*expReadsIn2NPeak_1bp
  return(round(rd,3))
}

#' converts pkmod rather than rd to nrd, normalized read depth
#' @param rdNormX_2Npeak read depth of the diploid peak normalized by the x axis
#' @inheritParams commonParameters
#'
pkmodToNRD <- function(segmentData, peakInfo, rdNormX_2Npeak){
  firstDigPeakIndex <- which.min(peakInfo$dPeaks) # may not be 1 if the grid is not forced to start at the first peak
  rdNormX_firstDigPeak <- peakInfo[firstDigPeakIndex, 'peakReadDepth_normX']   #
  nrdNormCoef <- 2/rdNormX_2Npeak*rdNormX_firstDigPeak   # will be 2 if the first digital peak is the 2N peak
  nrd <-  segmentData[,'pkmod'] * nrdNormCoef
  return(nrd)

}


#' extrapolate or interpolate read depth for a give copy number
#'
#' given read depth and copy number for two other positions
#' @examples extrapRD(cn=2,cn1=4,cn2=5,rd1=2,rd2=2.376)
#'
#' @export
#' @keywords internal
extrapRD <- function(cn,cn1,cn2,rd1,rd2){
  rd  <-  rd1 + ((cn - cn1) / (cn2 - cn1)) * (rd2 - rd1)
  return(rd)
}

#' find the two best digital peaks for calculating tumor percentage or extrapolating other peaks
#'
#' Best peak will be the 2N peak if it exists and is suitable, otherwise use main peak
#' Second best peak will be the taller of the two digital peaks adjacent to the best peak
#' Don't rely on 2N and then the 1N/3N because they may be small and/or missing
#' 2N peak is considered suitable if it is one of the top 3 biggest peaks
#' Do not default to two biggest peaks, we want to use lower copy number peaks which will typically align to the digital grid better due to less (sub) clonal mixing
#' @param peakCopyNum from peakInfo, copy number for the peak
#' @param peakHeight from peakInfo, height of each peak normalized so max peak = 1
#' @keywords internal
getTwoBestPeakIndexes <- function(peakCopyNum, peakHeight){
  # peakCopyNum = peakInfo$nCopy; peakHeight  = peakInfo$peakHeight;
  # peakHeight=c(0.568, NA, 1.0, 0.501, 0.584, 0.086, 0.039); peakCopyNum=c(3, NA,  4,  5,  7,  8, 12)

  twoNexists <- ifelse(any(peakCopyNum==2, na.rm = TRUE), TRUE, FALSE)
  peakMagnitudes <- peakHeight[!is.na(peakCopyNum)]
  peakMagRankings  <-  order(peakMagnitudes, decreasing = TRUE)

  # is 2N suitable?
  if(twoNexists){
    twoNPeakCNIndex <- which(peakCopyNum[!is.na(peakCopyNum)]==2)
    # use 2N peak if it is one of the top 3 biggest peaks
    twoNSuitable <- ifelse(which(peakMagRankings==twoNPeakCNIndex) <= 5, TRUE, FALSE) # use 2N peak as much as possible! i.e. 74002

  }else{
    twoNSuitable <- FALSE
  }

  # get best peak index for peakMagnitudes and its two neighbor peak indexes
  # use 2N peak (if suitable) otherwise the tallest peak
  # BEWARE peakMagnitudes and peakMagRankings do not include the NAs while peakHeight does! so indexes into the two arrays may differ
  if(twoNSuitable){
    twoNPeakIndex  <- which(peakCopyNum[!is.na(peakCopyNum)]==2)
    twoNpeak       <-  peakMagnitudes[twoNPeakIndex]
    bestPeakIndex  <-  which(peakHeight==twoNpeak)

    nextPeakIndex  <-  twoNPeakIndex+1 # these are digital peaks
    prevPeakIndex  <-  twoNPeakIndex-1
  }else{
    # NO- use tallest peak and second tallest peak - one of these peaks could be a really high copy number and the bigger the copy number the more likely it is a mixture and not reliable
    # YES-use tallest peak and its tallest of the two neighbor peak as these peaks will probably give you the most reliable value

    # BEWARE peakMagnitudes and peakMagRankings do not include the NAs while peakHeight does! so indexes into the two arrays may differ
    tallestPeak       <- peakMagnitudes[peakMagRankings[1]]  # max(peakMagnitudes)
    #secondTallestPeak <- peakMagnitudes[peakMagRankings[2]]  # max(peakMagnitudes[peakMagnitudes!=max(peakMagnitudes)])

    tallestPeakIndex       <- which(peakMagnitudes==tallestPeak)
    bestPeakIndex        <-   which(peakHeight==tallestPeak)

    nextPeakIndex  <-  tallestPeakIndex+1 # these are digital peaks
    prevPeakIndex  <-  tallestPeakIndex-1
  }

  if(FALSE){
    # use the left peak unless it is too small ( less than 10% of the bestPeak)
    # this mimics the choice of peaks pre-ploidy in tumorEstFunc() but has not been tested for use in ploidy
    if(!is.na(prevPeakIndex)){
      if(  peakMagnitudes[prevPeakIndex]/ peakMagnitudes[bestPeakIndex]>.1){
        secondBestPeakIndex <- prevPeakIndex
      }else{
        secondBestPeakIndex <- nextPeakIndex
      }
    }else{
      secondBestPeakIndex <- nextPeakIndex
    }

  }else{
    # a slight preference is given to the lower CN of the two peaks
    # for example, 1N vs 3N, if the 1N peak is less than 30% smaller, use it to find the tumor % as higher peaks are not as accurate (mixing, subclones etc.)
    tallestNeighborPeakIndex  <- which.max(c(peakMagnitudes[prevPeakIndex]*1.3, peakMagnitudes[nextPeakIndex ]))
    tallestNeighborPeak <-                 c(peakMagnitudes[prevPeakIndex],     peakMagnitudes[nextPeakIndex ])[tallestNeighborPeakIndex]
    secondBestPeakIndex <- which(peakHeight==tallestNeighborPeak)
  }

  # loginfo(' bestPeakIndex = %i; secondBestPeakIndex = %i',bestPeakIndex,secondBestPeakIndex)

  if(bestPeakIndex==secondBestPeakIndex){
    logerror('something is wrong, bestPeakIndex = secondBestPeakIndex %i = %i, this should not happen', bestPeakIndex,secondBestPeakIndex)
    return(NA)
  }
  if(is.na( peakCopyNum[bestPeakIndex]) | is.na( peakCopyNum[secondBestPeakIndex])){
    logerror('something is wrong, nCopy for one of the peak indexes is NA, this should not happen')
    return(NA)
  }
  return(list(bestPeakIndex=bestPeakIndex, secondBestPeakIndex=secondBestPeakIndex))

}



#' calculate tumor from ploidy read depth and copy number calls (peak positions)
#'
#' use \code{getTwoBestPeakIndexes} to find best peak and the biggest of its neighbor peaks.
#'
#' @param peakCopyNum from peakInfo, copy number for the peak
#' @param peakHeight from peakInfo, height of each peak normalized so max peak = 1, to determine two best peaks
#' @param peakReadDepth_1bp from peakInfo, read depth per peak normalized to 1bp window
#'
#' @keywords internal
calcTumorFromPloidyPeaks <- function(peakCopyNum, peakHeight,peakReadDepth_1bp,dPeaks){
  # peakInfo=calcPloidyResult$peakInfo
  # peakCopyNum = peakInfo$nCopy; peakHeight  = peakInfo$peakHeight; peakReadDepth_1bp=peakInfo$peakReadDepth_1bp; dPeaks=peakInfo$dPeaks

  bestPeaksResult <- getTwoBestPeakIndexes(peakCopyNum, peakHeight)
  rd1a <- peakReadDepth_1bp[bestPeaksResult$bestPeakIndex]
  rd1b <- dPeaks[bestPeaksResult$bestPeakIndex]
  rd2a <- peakReadDepth_1bp[bestPeaksResult$secondBestPeakIndex]
  rd2b <- dPeaks[bestPeaksResult$secondBestPeakIndex]
  cn1 <- peakCopyNum[bestPeaksResult$bestPeakIndex]
  cn2 <- peakCopyNum[bestPeaksResult$secondBestPeakIndex]


  percentTumor_a <- calcTumorRatio(rd1a,rd2a,cn1,cn2)*100 # using the readDepth at the max y position of each peak - will depend on which peaks you use, not spaced equal distances apart
  percentTumor_b <- calcTumorRatio(rd1b,rd2b,cn1,cn2)*100 # using the readDepth at the digital peak position each peak - does not depend on which peaks you use, all spaced equal distance apart

  percentTumorMean <- mean(c(percentTumor_a, percentTumor_b))
  # loginfo('mean percent tumor calc: %.1f',percentTumorMean)
  return(percentTumor=percentTumor_a)
}



#' add vertical lines to linear genome plot
#'
#' @param chromStarts must be same window size scale as the plotted frequency data
#' @param maxcn max chromosome number
#' @keywords internal
markChromEdges <- function(chromStarts,maxcn,vCol='gray90'){
  chroms <- convertChromToCharacter(1:maxcn)
  abline(v=chromStarts[1:maxcn], col=vCol)
  for(ic in 1:maxcn) {
    mtext(side = 1, text = chroms[ic], line=-1, at = (chromStarts[ic]+chromStarts[ic+1])/2, cex=0.7)
    mtext(side = 3, text = chroms[ic], line=-1, at = (chromStarts[ic]+chromStarts[ic+1])/2, cex=0.7)
  }
}


#' list of colors to use for the peaks in \code{peaksByDensity}
#' @keywords internal
getCNcolors <- function(){
  # cnColors <- c(palette()[-1], "orange", 'white', palette()[-1]) # remove black, add orange, white, and repeat to make sure we have enough colors

  # R 4.0 changed palette() to different shades of the original colors, use palette_R3.0 to get back to those original shades, which had standard names, back
  palette_R3.0 <- c("black","red","green3","blue","cyan","magenta","yellow","gray")
  cnColors <- c(palette_R3.0[-1], "#E69F00", palette_R3.0[1], palette_R3.0[-1]) # add orange, put black last, and repeat to make sure we have enough colors
  return(cnColors)
}



#' chromosome segment colors
#'
#' one color for each chromosome 1-
#' @keywords internal
getSegmentColors=function(){
  ### one color for each chromosome 1-22
  segColors <- c(
    #E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#A65628","#F781BF",
    "red","blue","cyan","gray45","magenta",
    "#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#A6D854","#FFD92F","#E5C494","gray75",
    "purple","#00FF92FF")


  # if(FALSE){
  #   # how the colors above were created
  #   segColors <- c(RColorBrewer::brewer.pal(9, 'Set1')[-c(6,9)] ,                     # remove the yellow and gray
  #                  'red', 'blue', 'cyan', 'gray45','magenta',                         # to get to 22
  #                  RColorBrewer::brewer.pal(8, 'Set2')[-c(8)], 'gray75',              # replace the gray with a slightly lighter color
  #                  'purple', "#00FF92FF"
  #   )
  #
  #   # past color scheme that was not ideal
  #   segColorsOld <- c(
  #     'black',
  #     RColorBrewer::brewer.pal(9, 'Set1')[-c(6)] ,                     # remove the yellow and gray
  #     palette()[-7],
  #     '#FF0000FF', '#00FF92FF', '#FFDB00FF',  '#0092FFFF',  '#4900FFFF',  '#FF00DBFF'
  #   )
  #   plot(1:22, type='n');    #abline(v=7.5);    abline(v=13.5)
  #   points(x=1:22, y=rep(5,22), pch=1:22, cex=2,col=segColors)
  #   points(x=1:22, y=rep(3,22), pch=1:22, cex=2,col=segColorsOld)
  # }

  return(segColors)
}



#' parameters available for non-default config inputs to calculatePloidy() as determined through testing
#' @export
#' @keywords internal
getAdjustablePloidyConfigParams <- function(){
  columns <- c( "origMaxPercentCutoffManual","grabDataPercentManual","forceFirstDigPeakCopyNum","minPeriodManual","allowedTumorPercent","heterozygosityScoreThreshold",
                "continueOnPloidyFail", "maxPeriodManual","minReasonableSegmentSize")
  return(columns)
}



