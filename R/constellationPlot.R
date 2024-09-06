#' load all the data necessary to make the constellation plot
#'
#' object created is a required input to \link{plotStarsInTheClouds} and be can be save to file
#' for future reference
#'
#' @param mainPeakNRD the NRD of the main Peak in the cnvBinned data.  Will equal 2 if main peak is the normal peak
#' @param expReadsIn2NPeak_1bp expected number of reads in a 1 bp bin for the diploid peak
#' @param hetScorePerBinFile full path to object output from \link{calculateHetScore}
#' @inheritParams commonParameters
#' @example inst/examples/constellationPlotExample.R
#' @export
loadStarsInTheClouds <- function(sampleId, inputDir, readDepthPer30kbBin,hetScorePerBinFile, hsNormMat, testVals, readDepthBinSize=30000,
                                 mainPeakNRD, expReadsIn2NPeak_1bp){
  # testVals - not used, just passed out with the other data needed for the plot

  yind=23
  mainChromsNoY <- 1:23
  mainChroms <- 1:24
  autosomes= 1:22
  coords <- getLinearCoordinates(chromosomes = mainChroms)

  ## load hetScore wig file  ---------------
  lohIn <- loadHetScoreFromWig(hetScorePerBinFile)
  loginfo('loading ref and alt counts from dir: %s', inputDir)

  lohTot <- list()
  covTot <- list()
  snpTot <- list()
  ## load LOH Count and LOH SNP data  ---------------
  for (chrNum in mainChromsNoY) {
    ichrChar=convertChromToCharacter(chrNum)
    loginfo('chrom %s',ichrChar)

    # snpFull: SNP information for the matching lohCountBpFull for sva
    # countBPFull: Counts of ref/alt occurrences for SNP position for each sva
    # loading BACDAC .Rds inputs
    iFile=file.path(inputDir, paste0(sampleId,'_','refAltCount_', ichrChar,'.Rds'))
    # logdebug('loading %s',iFile)

    iRefAltCount = readRDS(file=iFile )
    countBPFull=iRefAltCount[,c('ref', 'alt')]
    snpFull=iRefAltCount[,'pos']

    lohTot[[chrNum]] <- countBPFull #lohData
    covTot[[chrNum]] <- apply(countBPFull,1,sum) # coverage of all the SNPS, in each chrom
    snpTot[[chrNum]] <- snpFull #snpData
  }


  ### get frequency array ---------------
  # need two arrays - count(freq count) and array indexes (after masking)

  ## make sure copy number chroms are the same as the hetScore chroms... remove Y from copy number data
  # otherwise you get errors when it gets to a Y chrom position;
  chrYStartWsz <- binnedPosStart(coords@chromStart[yind], binSize = readDepthBinSize)
  w23   <- which(readDepthPer30kbBin$goodWindowArray < chrYStartWsz)
  frq23 <- readDepthPer30kbBin$readDepthArray[w23]  # number of fragments, in a 1kb window converted to readDepthBinSize (wsz), for chrs 1-23
  wdnsMSK23 <- readDepthPer30kbBin$goodWindowArray[w23]

  chrStartKey <- array(0, coords@maxcn)
  for(ichr in coords@chroms){
    # Round the start coordinate to the new bin size
    chrStartKey[ichr] <- 1+as.integer((coords@chromStart[ichr]-1)/readDepthBinSize)
  }



  ## masking and processing of some sort --------------
  # normalize read depth (rd) to NRD with normal peak = 2
  # hetScore also on a 30K bin, masked.
  cnvListFull <- 2*frq23/(readDepthBinSize*expReadsIn2NPeak_1bp) # all the rd data (except for chrY), normalized to normal rd=2, NRD=2
  posListFull <- wdnsMSK23
  chrStart <- c(chrStartKey,max(wdnsMSK23))
  chrList <- as.character(lohIn$seqnames)
  lohChrOutFull <- NULL
  cnvListChrFull <- NULL
  chrEnd <- NULL

  for(chrNum in autosomes) {
    # logdebug('process chrom %i',chrNum)
    lohSeqname <- convertChromToCharacter(chrNum, withChrPrefix=TRUE)
    lohChr <- lohIn$score[which(chrList==lohSeqname)]
    posListChr <- posListFull[which(posListFull%in%c(chrStart[chrNum]:(chrStart[chrNum+1]-1)))]
    cnvListChr <- cnvListFull[which(posListFull%in%c(chrStart[chrNum]:(chrStart[chrNum+1]-1)))]
    lohChrOut <- lohChr[posListChr-chrStart[chrNum]+1]
    #Use the mask to find places where the hetScore, in more than half of 23 TCGA normals, dropped below the 0.975 cutoff.
    lohMean <-0.9875 #
    oneSD <-0.0125
    cutoff <-lohMean-oneSD

    lohChrMed <- apply(hsNormMat[which(chrList==lohSeqname),],1,
                       function(x) length(which(x< cutoff))) # a value between 0-23 (0=none below, 23=all below)
    numLohRefSamples<-ncol(hsNormMat) # currently 23, as in the 23 TCGA normals
    whichPlot <- which(lohChrMed[posListChr-chrStart[chrNum]+1] < numLohRefSamples/2) # keep the data where less than half the samples were below..
    lohChrOutFull <- c(lohChrOutFull,lohChrOut[whichPlot])
    cnvListChrFull <- c(cnvListChrFull,cnvListChr[whichPlot])
    chrEnd <- c(chrEnd,length(lohChrOutFull))
  }
  chrEnd <- c(0,chrEnd)

  # logdebug('calculating theoretical lambda, the poisson that best fits the coverage') crux of the
  # algorithm, where it determines the theoretical lambda based on the assumption of a poisson for
  # extrapolation to higher or lower copy number
  #
  # Sample data from the mainPeak (not the 2N peak) to be sure to get a good sampling, in case the
  # 2N peak is small. Then find the poisson that best matches the sample data and save the
  # lambda that corresponds to that poisson.

  ### sample data in main peak
  ### find lambdaMain, the poisson that best fits the coverage
  ### mainPeakNRD is NRD of main peak in cnvBinnedData

  sampLoc <- sample(x =  posListFull[(which(cnvListFull>= (mainPeakNRD-0.025) & cnvListFull<(mainPeakNRD+0.025)))],
                    size = min(length(which(cnvListFull>= (mainPeakNRD-0.025) & cnvListFull<(mainPeakNRD+0.025))),1000)) # take 1000 unless there aren't 1000 data points available
  ## get the chrom number of the sampled positions
  chrNumList <-linearToBima(globalPos = sampLoc, binSize =readDepthBinSize )$svaNumber
  posValTemp <- (sampLoc-chrStart[chrNumList]+1)*30000-15000

  # hist of the coverage you've sampled
  tabTemp <- hist(unlist(sapply(1:length(chrNumList),
                                function(x) ((covTot[[chrNumList[x]]])[which(snpTot[[chrNumList[x]]]>posValTemp[x]-500000&snpTot[[chrNumList[x]]]<posValTemp[x]+500000)]))),
                  breaks=c(seq(-0.5,1000.5,1),1000000),plot=FALSE)$counts[1:1001]

  midVal <- which.max(tabTemp)-1
  pSim <- sapply(seq(max(0.5,floor(midVal/2)),min(1000,2*midVal),0.01),
                 function(y) sapply(0:1000,
                                    function(x) dpois(x,y)))

  seqTo <-max(0.5,floor(midVal/2))
  seqFrom <-min(1000,2*midVal)
  lambdaMainOrig <- seq(seqTo,seqFrom,0.01)[which.min(apply(pSim,2,
                                                            function(x) sum(abs(tabTemp/sum(tabTemp)-x))))]

  cnvListChrFullOrig <- cnvListChrFull

  starCloudPlotInputs <-(list(testVals=testVals,
                              lambdaMainOrig=lambdaMainOrig,
                              cnvListChrFullOrig=cnvListChrFullOrig,
                              lohChrOutFull=lohChrOutFull,
                              chrEnd=chrEnd))
  loginfo('END OF FUNCTION: loadStarsInTheClouds')
  return(starCloudPlotInputs)
}


#' plot heterozygosity (actual and theoretical) vs NRD for a given tumor purity
#'
#' annotations in top outer margin include values for ploidy, tumor purity, 2N+LOH and ploidy state
#'
#' @param diploidPeakNRD the NRD of the diploid peak, don't assume it is 2, may be
#'    choosing a different peak than from a previous calculation
#' @param mainPeakNRD the NRD of the main Peak
#' @param tau tumor ratio
#' @param digitalPeakZone size of the shaded pink boxes, hetScores in this areas were used to
#'  compute the hetScore mode of the peak
#' @param cleanMode logical if TRUE do not include extra annotations
#' @param starCloudPlotInputs object returned from \link{loadStarsInTheClouds}
#' @param plotEachChrom logical if TRUE show each chromosome in a separate plot
#' @param plotCex factor to alter the cex for plots
#' @param addSegmentLegend logical if TRUE color code segments by chromosome and add corresponding legend to the plot
#' @inheritParams commonParameters
#' @example inst/examples/constellationPlotExample.R
#' @export
plotStarsInTheClouds <- function(sampleId, alternateId=NULL, starCloudPlotInputs, plotEachChrom=FALSE,
                                 mainPeakNRD,  diploidPeakNRD, tau, segmentData=NULL, peakInfo=NULL,
                                 forceFirstDigPeakCopyNum=-1,grabDataPercentManual=-1, origMaxPercentCutoffManual=-1,minPeriodManual=-1,maxPeriodManual=-1,
                                 minReasonableSegmentSize=5.5e6,
                                 digitalPeakZone = 0.05, heterozygosityScoreThreshold = 0.98,
                                 cleanMode=FALSE,
                                 plotCex=1,addSegmentLegend=FALSE
){

  autosomes=1:22
  testValsSt <- starCloudPlotInputs$testVals;
  lambdaMainOrig <- starCloudPlotInputs$lambdaMainOrig;
  cnvListChrFullOrig <- starCloudPlotInputs$cnvListChrFullOrig;
  lohChrOutFull <- starCloudPlotInputs$lohChrOutFull
  chrEnd <- starCloudPlotInputs$chrEnd

  # find the normalized read depth for the diploid peak (do not assume it is 2)
  if(!is.null(peakInfo)){

    # get mainPeak info,  mainPeak may not be the normal peak
    rdNormX_Mainpeak    <- peakInfo[which(peakInfo$rankByHeight==1),        'peakReadDepth_normX']

    twoNexists <-any(peakInfo$nCopy==2,na.rm = TRUE) # in case it doesn't exist
    if(twoNexists){
      rdNormX_2Npeak      <- peakInfo[which(peakInfo$nCopy==2),               'peakReadDepth_normX']
    }else{
      # find the best two indexes for doing the extrapolation (or tumor percent)
      bestPeaksResult <- getTwoBestPeakIndexes(peakCopyNum=peakInfo$nCopy, peakHeight=peakInfo$peakHeight)
      # extrapolate read depth for 2N level
      rdNormX_2Npeak <- extrapRD(cn = 2,
                                 cn1 = peakInfo[bestPeaksResult$bestPeakIndex, 'nCopy'],              cn2 = peakInfo[bestPeaksResult$secondBestPeakIndex, 'nCopy'],
                                 rd1 = peakInfo[bestPeaksResult$bestPeakIndex,'peakReadDepth_normX'], rd2 = peakInfo[bestPeaksResult$secondBestPeakIndex,'peakReadDepth_normX']
      )
    }

    if(is.null(diploidPeakNRD)){
      diploidPeakNRD <-round( mainPeakNRD*rdNormX_2Npeak/rdNormX_Mainpeak, 3)
      # logdebug('derived diploidPeakNRD: %.3f',diploidPeakNRD)
    }

  }else{
    if(is.null(diploidPeakNRD)){
      logerror('must provide either peakInfo or diploidPeakNRD')
    }
  }


  cnvListChrFull <- cnvListChrFullOrig*(2/diploidPeakNRD) # scale cnvBinned NRD to be around the new 2N peak.
  lambdaMain <- lambdaMainOrig*(diploidPeakNRD/mainPeakNRD)

  if(is.null(segmentData) || is.null(peakInfo)){
    plotMaxCN  <- round(min(12, quantile(cnvListChrFull,probs = .99)+2)) # add 2 to expand the plot just a bit but not too much
    ylimMinTemp <- calcNrd(cn=0, tau) # min(cnvListChrFull)
    loginfo('max(cnvListChrFull):%.2f',max(cnvListChrFull))
  }else{
    if(!'nrd' %in% names(segmentData)){
      stop('nrd is missing from segmentData')
    }

    nearMaxNRD <- quantile(segmentData$nrd,probs = .99)
    nearMaxCN <- 2+(nearMaxNRD-2)/tau
    plotMaxCN  <- round(min(12, nearMaxCN+1.5)) # add 1.5 to expand the plot just a bit
    ylimMinTemp <- quantile(segmentData$nrd,probs = .01)
  }


  seqRange <- seq(0,plotMaxCN,0.01)
  pPlot <- unlist(sapply(seqRange,function(x) ((1-tau))/(x*tau+2*(1-tau)))) # prob for the binomial
  plotRange <- 2+tau*(seqRange-2);

  # renormalize all of the probabilities by the sum of the probabilities
  shiftVals <- sapply(1:length(plotRange),
                      function(x) sum(((dpois(4:1000,lambdaMain*(plotRange[x]/2)))/sum(dpois(4:1000,lambdaMain*(plotRange[x]/2))))*testValsSt[(4:1000)-3,which.min(abs(seq(0.001,0.5,0.001)-pPlot[x]))]))
  medVals   <- sapply(1:length(plotRange),
                      function(x) sum(((dpois(4:1000,lambdaMain*(plotRange[x]/2)))/sum(dpois(4:1000,lambdaMain*(plotRange[x]/2))))*testValsSt[(4:1000)-3,which.min(abs(seq(0.001,0.5,0.001)-0.5))]))

  seqRange <- seq(1,plotMaxCN,1) # draw stars just for the CN levels
  plotStarRange <- unlist(sapply(seqRange,function(x) sapply(1:(1+floor(x/2)),function(y) x*tau+2*(1-tau))))
  plotStarRange <- c(2*(1-tau),plotStarRange)  # NRD for each star, y values
  pPlot <- unlist(sapply(seqRange,
                         function(x) sapply(1:(1+floor(x/2)),
                                            function(y) ((y-1)*tau+(1-tau))/(x*tau+2*(1-tau)))))
  pPlot <- c(0.5,pPlot)

  # 3 or less reads is not informative (possibility of mis-sequencing, we don't know if there is an actual snp, etc), so start with 4
  # het score value for each allele fraction, for each coverage level
  starVals <- sapply(1:length(plotStarRange),
                     function(x) sum(((dpois(4:1000,lambdaMain*(plotStarRange[x]/2)) )/sum(dpois(4:1000,lambdaMain*(plotStarRange[x]/2))))*testValsSt[(4:1000)-3,which.min(abs(seq(0.001,0.5,0.001)-pPlot[x]))]))
  # perfect het score value, for each coverage level
  medStarVals <- sapply(1:length(plotStarRange),
                        function(x) sum(((dpois(4:1000,lambdaMain*(plotStarRange[x]/2)))/sum(dpois(4:1000,lambdaMain*(plotStarRange[x]/2))))*testValsSt[(4:1000)-3,which.min(abs(seq(0.001,0.5,0.001)-0.5))]))

  # define "ploidy" for the sample
  ploidyNRD <- mean(cnvListChrFull)
  ploidyCN <- calcCopyNumber(NRD=ploidyNRD, tau)


  xlimTemp <- min(0.6,0.9*min(starVals/medStarVals,na.rm = TRUE))
  ylimMinTemp <- min(ylimMinTemp, min(plotStarRange))

  ### one color for each chromosome 1-22
  if(addSegmentLegend){
    segColors = getSegmentColors()
  }else{
    segColors = rep('gray30', length(autosomes) )
  }


  if(!is.null(segmentData)){
    starInfoTemp <- list(starVals=starVals,medStarVals=medStarVals,plotStarRange=plotStarRange)
    starLookUp <- makeStarLookUpTable(starCloudResult=starInfoTemp,percentTumor=tau*100)
    allelicSegments <- allelicCNV(starLookUp, segmentData=segmentData)
    lohContent <- getLohContent(allelicSegments)  ## getLohContent() is in ploidy.R
    lohContentA_maj2_min0 <- lohContent$lohContentA_maj2_min0
    loginfo('2N+LOH content: %s',round(lohContentA_maj2_min0,3))
  }else{
    lohContentA_maj2_min0 <- NA
  }

  ploidyState = getPloidyState(twoNpLOH=lohContentA_maj2_min0, ploidy=ploidyCN)

  # logdebug('poisson Cov Fit= %s\t diploidPeakNRD= %s\t mainPeakNRD=%s',lambdaMainOrig,diploidPeakNRD, round(mainPeakNRD,3))

  # plot once with all chromosome data -------
  if(TRUE){

    # set up plot
    xlimits <- c(xlimTemp, 1.1)
    ylimits <- c(ylimMinTemp, max(plotRange))
    plot(lohChrOutFull,cnvListChrFull,type='n',
         xlim=xlimits,
         ylim=ylimits,
         xlab="Heterozygosity Score",
         ylab="NRD",
         cex.lab=plotCex,cex.axis=plotCex
    )

    if(!cleanMode){
      # annotate discretized peak locations with boxes to show the zones around the peaks
      # draw first so everything else can be layered on top
      if(!is.null(segmentData) && !is.null(peakInfo)){
        xMax <-par('usr')[2]
        rect(xleft  = 0, xright = xMax,
             ybottom = peakInfo[!is.na(peakInfo$dPeaks), 'peakReadDepth_normX']*2/rdNormX_2Npeak - digitalPeakZone,
             ytop    = peakInfo[!is.na(peakInfo$dPeaks), 'peakReadDepth_normX']*2/rdNormX_2Npeak + digitalPeakZone,

             border = NA,  col="#f9edfa")
      }
      box() # redraw border

      ### vertical line to mark cutoff
      abline(v=heterozygosityScoreThreshold, col='red',lty='dashed')
      mtext(text=heterozygosityScoreThreshold, side=1, at = heterozygosityScoreThreshold, cex=.8, col='red', line= -1)
    }

    # CN level annotations
    iCN <-0:plotMaxCN
    abline(h=2+((iCN-2)*tau), lty=2, col='gray85')
    # mtext(text='CN', side=4, col='black', las=1, line = 1.5)
    axis(side=4, at=2+((iCN-2)*tau), labels=paste0(iCN,'N'), tick=TRUE, las=1, cex.axis=plotCex)


    # add -actual- hetScore binned data
    points(lohChrOutFull, cnvListChrFull,
           pch='.',
           col=scales::alpha('black', 0.25))


    # add -median- hetScore PER segment
    if(!is.null(segmentData) && !is.null(peakInfo)){
      points(segmentData$lohScoreMedian, segmentData$nrd,
             pch=segmentData$chr,cex=0.9,
             col=scales::alpha(segColors[segmentData$chr],0.75))


      # peak rank labels
      if(!cleanMode){
        ymax <- par('usr')[4]
        colorCodeRank <- ifelse(is.na(peakInfo$dPeaks),'gray80', 'gray50')
        mtext(text=peakInfo$rankByHeight, side=4, at=peakInfo$peakReadDepth_normX*2/rdNormX_2Npeak, las=2, line=-.5, col=colorCodeRank, adj=1, cex=.7)
        mtext(text='Peak', side=4, at=ymax*.98, col='gray50', las=1, line=-2, cex=.7)
        mtext(text='Rank', side=4, at=ymax*.96, col='gray50', las=1, line=-2, cex=.7)
      }
    }

    # LOH line
    lines(shiftVals/medVals,plotRange,type='l',xlim=c(0,1),col='purple')
    # stars for each heterozygosity ratio for each cn level - make sure this layers on top the segment data
    points(starVals/medStarVals,plotStarRange,pch="*",cex=2,xlim=c(0,1),col='green2')

    legend('bottomleft',legend = c('LOH', 'theoretical', 'actual'),col = c('purple',  'green2', 'black'),
           lty=c(1,NA,NA), lwd=c(1.5,NA,NA), pch=c(NA, "*", "."), cex=.9*plotCex, pt.cex=2,
           bty='n')


    # outer margin plot annotations
    mtext(sampleId,    side=3, adj=0.05, line=1, outer = TRUE,cex = plotCex)
    mtext(alternateId, side=3, adj=0.05, line=0, outer = TRUE,cex = plotCex)
    mtext(side=3, text=paste('ploidy: ',round(ploidyCN,1)), adj=0.5, line=1, outer = TRUE,cex = plotCex)
    mtext(side=3, paste0('tumor: ',round(tau*100), '%'),    adj=0.5, line=0, outer = TRUE,cex = plotCex)
    mtext(side=3, paste0('ploidy state: ',ploidyState),    adj=0.95, line=1, outer = TRUE,cex = plotCex)
    if(!is.na(lohContentA_maj2_min0)){
      mtext(3, text=paste("2N+LOH:", round(lohContentA_maj2_min0,3)), adj = 0.95, line=0, outer = TRUE,cex = plotCex)
    }


    if(addSegmentLegend){
      # chrom symbol legends
      legend("topleft", legend=c(1:22), col=segColors, pch=1:22, cex=.95*plotCex)
    }


    if(!cleanMode){
      # manual input annotations
      if(forceFirstDigPeakCopyNum > 0){
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
    }

  }

  # plot each chromosome individually -----
  if(plotEachChrom){
    op=par(mfrow=c(3,3), mgp=c(2, 0.5, 0),mar=c(3.5, 3.5,2, 2.5) ,oma=c(0,0,2,0))
    for(chrNum in autosomes) {
      plot(lohChrOutFull[(chrEnd[chrNum]+1):chrEnd[chrNum+1]],cnvListChrFull[(chrEnd[chrNum]+1):chrEnd[chrNum+1]],
           pch='.',
           ylim=c(0,max(plotRange)),
           xlim=c(xlimTemp,1.1),
           xlab="Heterozygosity Score",
           ylab="NRD",
           col=scales::alpha('black',0.3),
           cex.lab=1.2)
      lines(shiftVals/medVals,plotRange,type='l',xlim=c(0,1),col='purple')
      points(starVals/medStarVals,plotStarRange,pch="*",cex=2,xlim=c(0,1),col='green2')
      title(paste0("chr=",chrNum))


      # CN level annotations
      iCN <-0:plotMaxCN
      abline(h=2+((iCN-2)*tau), lty=2, col='gray85')
      # mtext(text='CN', side=4, col='black', las=1, line = 1.5)
      axis(side=4, at=2+((iCN-2)*tau), labels=paste0(iCN,'N'), tick=TRUE, las=1)

      # print once for each page
      if(chrNum %in% c(1, 10, 19)){
        # mtext(side=3, paste0('tumor: ',round(tau*100), '%'),    adj=0.05, outer = TRUE, line=-2)
        # mtext(side=3, text=paste('ploidy: ',round(ploidyCN,1)), adj=0.05,  outer = TRUE,line=-1)
        # mtext(c(sampleId, alternateId), side=3, line=0, adj=c(.05,.95), outer = TRUE)
        mtext(paste(sampleId, alternateId), side=3, line=0, adj=c(.05), outer = TRUE)
        mtext(paste(' tumor:',round(tau*100), '% ',' ploidy:',round(ploidyCN,1) ), side=3, line=0, adj=c(.95), outer = TRUE)
      }
    }
    par(op)
  }

  return(list(starVals=starVals,
              medStarVals=medStarVals,
              plotStarRange=plotStarRange,
              shiftVals=shiftVals,
              medVals=medVals,
              diploidPeakNRD=diploidPeakNRD,
              ploidyCN=ploidyCN,
              lohContent=lohContent,
              plotAxisLimits = list(hetScoreAxisLims=xlimits,
                                    nrdAxisLims=ylimits),
              allelicSegments=allelicSegments
  ))
}



#' draw constellation plot left of the linear genome plot
#'
#' @param starCloudPlotInputs object output from \link{loadStarsInTheClouds}
#' @param calcPloidyResult object output from \link{calculatePloidy}
#' @inheritParams commonParameters
#' @export
twoPanelReport=function(starCloudPlotInputs, calcPloidyResult, readDepthPer30kbBin, segmentation,
                        sampleId=NULL, alternateId=NULL, gainColor='blue', lossColor= 'red',
                        noPdf=FALSE, outputDir=NULL){

  if (!noPdf) {
    if(!is.null(outputDir)){
      ploidyPdfFile <-file.path(dir=outputDir, paste0(sampleId, 'BACDAC_ploidyReport.pdf'))
      loginfo('will plot to pdf: %s', ploidyPdfFile)
      if(!dir.exists(dirname(ploidyPdfFile))) {
        loginfo("creating output directory: %s", dirname(ploidyPdfFile))
        dir.create(path=file.path(dirname(ploidyPdfFile)),mode = "0775")
      }
      pdf(file = ploidyPdfFile,width=12, height=4, title=paste0('BACDAC_ploidy_',sampleId)) # paper="a4r",
      on.exit(dev.off(),add = TRUE)
    }else{
      logwarn('can not output pdf, outputDir=NULL, will print to plot window')
    }
  }

  # 1234 x 445
  graphics::layout( matrix(c(1,2),nrow=1),
                    heights= c(1),
                    widths = c(.5,1))   # Widths of the two columns
  labelCex=1.5
  plotCex= ifelse(noPdf,1.2,0.9)
  leftFigLabel=NULL;rightFigLabel=NULL
  diploidPeakNRD=getDiploidPeakNRD(calcPloidyResult)

  # left panel
  op <- par(mar=c(3.5,3.5,1,2),mgp=c(1.5,0.5,0),oma=c(0,0,2,0)+.2)
  starCloudResult=plotStarsInTheClouds(sampleId=sampleId, starCloudPlotInputs=starCloudPlotInputs, diploidPeakNRD=diploidPeakNRD,
                                       tau=min(1,calcPloidyResult$percentTumor/100),
                                       plotEachChrom=FALSE, mainPeakNRD=getMainPeakNRD(calcPloidyResult),
                                       segmentData=calcPloidyResult$segmentData, peakInfo=calcPloidyResult$peakInfo,
                                       digitalPeakZone =calcPloidyResult[['iterationStatsAll']][['digitalPeakZone']],
                                       plotCex=plotCex)
  myAt=starCloudResult$plotAxisLimits$nrdAxisLims[2]
  mtext(leftFigLabel, side=2, at=myAt,cex=labelCex,las=1,line=1.5)

  # right panel
  # convert the nrd axis limits in the constellation plot to rd, so the linear genome plot can be on the same scale
  rdAxisLimits= calcRD(nrd=starCloudResult$plotAxisLimits$nrdAxisLims, wsz=readDepthPer30kbBin$windowSize, calcPloidyResult$expReadsIn2NPeak_1bp)
  linearGenomePlot( readDepthPer30kbBin=readDepthPer30kbBin, readDepthBinSize=readDepthPer30kbBin$windowSize, segmentation=segmentation,
                    allelicSegments=starCloudResult$allelicSegments,
                    gainColor = gainColor, lossColor= lossColor, yAxisLimits = rdAxisLimits,plotCex=plotCex)
  myAt=rdAxisLimits[2]
  mtext(rightFigLabel, side=2, cex=labelCex,at =myAt,las=1,line=1.5)
  par(op)

  return(starCloudResult)
}



# make a blank plot, no data, used as a stand in when combining lots of samples into one pdf
blankPlotStarsInTheClouds=function(
    heterozygosityScoreThreshold=0.98,
    tau=0.8
    ){


  # set up plot
  plot(x=1:10,y=1:10,type='n',
       ylim=c(0, 5), yaxt='n',
       xlim=c(.4, 1.1),
       xlab="Heterozygosity Score",
       ylab="")

  ### vertical line to mark cutoff
  abline(v=heterozygosityScoreThreshold, col='red',lty='dashed')
  mtext(text=heterozygosityScoreThreshold, side=1, at = heterozygosityScoreThreshold, cex=.75, col='red', line=-0.2)

  ### one color for each chromosome 1-22
  segColors = getSegmentColors()

  # legend
  legend("topleft", legend=c(1:22), col=segColors, pch=1:22, cex=.95)


  # CN level annotations
  iCN <-0:8
  abline(h=2+((iCN-2)*tau), lty=2, col='gray85')
  # mtext(text='CN', side=4, col='black', las=1, line = 1.5)
  axis(side=4, at=2+((iCN-2)*tau), labels=paste0(iCN,'N'), tick=TRUE, las=1)

}
