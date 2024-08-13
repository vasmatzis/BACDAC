# code and example to run Jamie's stars in the clouds plot - the actual and theoretical heterozygosity scores vs NRD for a given tumor ratio
# source this file first, to make the functions below available


#' load all the data necessary to make the plot
#' @param mainPeakcnvBinnedNRD the NRD of the main Peak in the ---cnvBinned data---, if mainPeak based, will be 2, but if ploidybased, may be different
loadStarsInTheClouds <- function(sampleId, postProcessingDir, rgdObject, cnvBinnedData,  wsz, mainPeakcnvBinnedNRD){

  # make sure we have a couple of necessary values first
  if(is.na(cnvBinnedData$expectedNormalBin)){
    # note this is a function, but this will probably bypass some of the purpose of testing for this value:     expectedNormalBin <- getExpectedNormalBin(cnvBinnedData)
    # cnvBinnedData$expectedNormalBin=ploidyBasedNormalBin
    stop('cnvBinnedData$expectedNormalBin is NA')
  }
  if(is.na(mainPeakcnvBinnedNRD)){
    stop('mainPeakcnvBinnedNRD is NA')
  }


  logdebug('mainPeakcnvBinnedNRD=%.3f',mainPeakcnvBinnedNRD)
  bimaVersion <-rgdObject$processInformation$bimaVersion

  lohChroms <- c(svaAutosomes(rgdObject), svaAllosomes(rgdObject, chromosome = "X"))
  standardChroms <- unique(sort(c(svaAutosomes(rgdObject), svaAllosomes(rgdObject))))
  coords <- cnvBinnedData[['coords']] # another way: getLinearCoordinates(rgdObject, standardChroms)

  lohTot <- list()
  covTot <- list()
  snpTot <- list()
  ## load LOH Count and LOH SNP data  ---------------
  for(chrNum in lohChroms) {
    loginfo('load loh data for chr %i',chrNum)
    # Counts of ref/alt occurrences for SNP position for each sva
    lohCountBpFile <- getTypedFile("lohCountBpFull",postProcessingDir,values=list(sampleId=sampleId, svaNumber=chrNum), legacy = TRUE)@path
    lohData <- loadRdata(lohCountBpFile)
    # SNP information for the matching lohCountBpFull for sva
    lohSnpFile     <- getTypedFile("lohSnpFull",    postProcessingDir,values=list(sampleId=sampleId, svaNumber=chrNum), legacy = TRUE)@path
    snpData <- loadRdata(lohSnpFile)
    lohTot[[chrNum]] <- lohData
    covTot[[chrNum]] <- apply(lohData,1,sum) # coverage of all the SNPS, in each chrom
    snpTot[[chrNum]] <- snpData
  }



  ## load two reference files  ---------------
  # lohMat: LOH analysis mask, used to look for places in the 23 TCGA normals where more than half dropped below the a (i.e. 0.975) cutoff.
  # testVals: used to find each possible heterozygosity value for each copy number level (find the right spots for the stars)
  lohMat   <- loadRdata(file.path(mainDir, 'NextGen/Misc/pipelineInputs/hetScoreAnalysis/lohMat.Rdata'))
  testVals <- loadRdata(file.path(mainDir, 'NextGen/Misc/pipelineInputs/hetScoreAnalysis/testVals.Rdata'))


  ## load LOH wig file  ---------------
  lohAnalysisWigFile <- getTypedFile("lohAnalysisWig",postProcessingDir,values=list(sampleId=sampleId), legacy = TRUE)@path
  lohIn <- rtracklayer::import.wig(lohAnalysisWigFile)


  # load copy number data:   ---------------
  # need two arrays - count(freq count) and array indexes (after masking)
  shrunkData  <- bmdSvPipeline:::shrinkCnvBinned(data=cnvBinnedData, wsz = wsz)
  temp1 <- bmdSvPipeline:::cnvBinnedToLegacy(shrunkData)

  ## make sure copy number chroms are the same as the loh chroms... remove Y from copy number data
  # otherwise you get errors in tabTemp when it gets to a Y chrom position; not likely to happen unless Y is gained as it is in SA43002
  yind  <- svaAllosomes(rgdObject, chromosome = "Y")
  chrYStartWsz <- binnedPosStart(cnvBinnedData[['coords']]@chromStart[yind], binSize = wsz)
  w23   <- which(temp1$wdnsMSK < chrYStartWsz)
  frq23 <- temp1$frq[w23]  # number of fragments, in a 1kb window converted to wsz, for chrs 1-23
  wdnsMSK23 <- temp1$wdnsMSK[w23]

  chrStartKey <- array(0, shrunkData$coords@maxcn)
  for(ichr in shrunkData$coords@chroms){
    # Round the start coordinate to the new bin size
    chrStartKey[ichr] <- 1+as.integer((shrunkData$coords@chromStart[ichr]-1)/wsz)
  }



  ## processing of some sort??? --------------
  # normalize cnv to NRD with normal peak = 2
  # LOH (hetScore) also on a 30K bin, masked.

  cnvListFull <- 2*frq23/(wsz*cnvBinnedData$expectedNormalBin) # all the CNV data (except for chrY), normalized to normal CNV=2, NRD=2, normal might be ploidy based or main peak based
  posListFull <- wdnsMSK23
  chrStart <- c(chrStartKey,max(wdnsMSK23))
  chrList <- as.character(lohIn@seqnames)
  lohChrOutFull <- NULL
  cnvListChrFull <- NULL
  chrEnd <- NULL
  # colList <- colors()[seq(1,length(colors()),length.out=24)[2:23]]
  for(chrNum in svaAutosomes(rgdObject)) {
    logdebug('process chr %i',chrNum)
    lohSeqname <- convertChromToCharacter(chrNum, rgdObject = rgdObject, withChrPrefix=TRUE)
    lohChr <- lohIn$score[which(chrList==lohSeqname)]
    posListChr <- posListFull[which(posListFull%in%c(chrStart[chrNum]:(chrStart[chrNum+1]-1)))]
    cnvListChr <- cnvListFull[which(posListFull%in%c(chrStart[chrNum]:(chrStart[chrNum+1]-1)))]
    lohChrOut <- lohChr[posListChr-chrStart[chrNum]+1]
    #This is where the mask is used. Looks for places where in the 23 TCGA normals more than half dropped below the 0.975 cutoff.
    #Not optimized or anything.
    lohMean <-0.9875 # as used in the pipeline
    oneSD <-0.0125
    cutoff <-lohMean-oneSD

    lohChrMed <- apply(lohMat[which(chrList==lohSeqname),],1,
                       function(x) length(which(x< cutoff))) # a value between 0-23 (0=none below, 23=all below)
    numLohRefSamples<-ncol(lohMat) # currently 23, as in the 23 TCGA normals
    whichPlot <- which(lohChrMed[posListChr-chrStart[chrNum]+1] < numLohRefSamples/2) # keep the data where less than half the samples were below..
    lohChrOutFull <- c(lohChrOutFull,lohChrOut[whichPlot])
    cnvListChrFull <- c(cnvListChrFull,cnvListChr[whichPlot])
    chrEnd <- c(chrEnd,length(lohChrOutFull))
  }
  chrEnd <- c(0,chrEnd)

  loginfo('calc theoretical lambda')
  # crux of the algorithm, where it determines the theoretical lambda based on the assumption of a poisson for extrapolation to higher or lower copy number
  #
  # Sample data from the mainPeak (not the 2N peak) to be sure to get a good sampling, in case the 2N peak is small ie. 26273
  # then find the poisson that best matches the sample data and save the lambda that corresponds to that poisson.
  # TODO: could this ever break for really low coverage samples? possibly so there might be modifications needed

  ### sample data in main peak
  ### find lambdaMain, the poisson that best fits the coverage
  ### mainPeakcnvBinnedNRD is NRD of main peak in cnvBinnedData

  sampLoc <- sample(x =  posListFull[(which(cnvListFull>= (mainPeakcnvBinnedNRD-0.025) & cnvListFull<(mainPeakcnvBinnedNRD+0.025)))],
                    size = min(length(which(cnvListFull>= (mainPeakcnvBinnedNRD-0.025) & cnvListFull<(mainPeakcnvBinnedNRD+0.025))),1000)) # take 1000 unless there aren't 1000 data points available
  ## get the chrom number of the sampled positions
  # replaces Jamie's... sapply(sampLoc,function(x) which(order(c(x,chrStart))==1)-1)
  chrNumList <-linearToBima(rgdObject,globalPos = sampLoc, binSize =wsz )$svaNumber
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
                              bimaVersion=bimaVersion))
  loginfo('finished getting inputs')
  return(starCloudPlotInputs)
}


#' plot heterozygosity (actual and theoretical) vs NRD for a given tumor ratio
#'
#' @param dipVal               the NRD of the diploid peak in  ---cnvBinned data---, may not be 2 if you choose a different peak from the pipeline to normalize by
#' @param mainPeakcnvBinnedNRD the NRD of the main Peak in the ---cnvBinned data---, if mainPeak based, will be 2, but if ploidybased, may be different
#' @param tau tumor ratio
plotStarsInTheClouds <- function(sampleId, folderId, starCloudPlotInputs, dipVal, tau, plotEachChrom=FALSE, mainPeakcnvBinnedNRD,
                                 segmentData=NULL, peakInfo=NULL, # add the ploidy inputs
                                 bimaVersion=NULL, forceFirstDigPeakCopyNum=-1,grabDataPercentManual=-1, origMaxPercentCutoffManual=-1,minPeriodManual=-1,maxPeriodManual=-1,minReasonableSegmentSize=5.5e6, # plot annotations
                                 digitalPeakZone = 0.05,heterozygosityScoreThreshold = 0.98,
                                 paperMode=FALSE  # without all the extra fluff, crisp and clean for papers and presentations
){
  #  paperMode=TRUE; digitalPeakZone = 0.05; heterozygosityScoreThreshold = 0.98; dipVal=NULL; tau=min(1,result$percentTumor/100); plotEachChrom=FALSE; segmentData=result$segmentData; peakInfo=result$peakInfo


  testVals <- starCloudPlotInputs$testVals;
  lambdaMainOrig <- starCloudPlotInputs$lambdaMainOrig;
  cnvListChrFullOrig <- starCloudPlotInputs$cnvListChrFullOrig;
  lohChrOutFull <- starCloudPlotInputs$lohChrOutFull

  # do not assume the diploid peak NRD for cnvBinned = 2, so get it

  if(!is.null(peakInfo)){

    # use mainPeakcnvBinnedNRD then you don't have to assume  mainPeak was normalized to a read depth of 2
    # cnvBinned, if run in pipeline with other ploidy output may not be the same as the ploidy output here.
    rdNormX_Mainpeak    <- peakInfo[which(peakInfo$rankByHeight==1),        'peakReadDepth_normX']

    twoNexists <-any(peakInfo$nCopy==2,na.rm = TRUE) # just in case it doesn't exist, i.e. ME26265
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

    if(is.null(dipVal)){
      dipVal <-round( mainPeakcnvBinnedNRD*rdNormX_2Npeak/rdNormX_Mainpeak, 3)
      loginfo('dipVal: %.3f',dipVal)
    }
  }
  if(is.null(dipVal) | is.null(peakInfo)){
    logerror('must provide peakInfo or dipVal')
  }



  cnvListChrFull <- cnvListChrFullOrig*(2/dipVal) # scale cnvBinned NRD to be around the new 2N peak.
  lambdaMain <- lambdaMainOrig*(dipVal/mainPeakcnvBinnedNRD)

  if(is.null(segmentData) || is.null(peakInfo)){
    plotMaxCN  <- round(min(12, quantile(cnvListChrFull,probs = .99)+2)) # add 2 to expand the plot just a bit but not too much
    ylimMinTemp <- calcNrd(cn=0, tau) # min(cnvListChrFull)
    loginfo('max(cnvListChrFull):%.2f',max(cnvListChrFull))
  }else{
    if(!'nrd' %in% segmentData){
      NRD <- bmdSvPipeline:::pkmodToNRD(segmentData, peakInfo, rdNormX_2Npeak)
      segmentData[,'nrd'] <-  NRD
    }




    nearMaxNRD <- quantile(segmentData$nrd,probs = .99)
    nearMaxCN <- 2+(nearMaxNRD-2)/tau
    plotMaxCN  <- round(min(12, nearMaxCN+2)) # add 2 to expand the plot just a bit but not too much
    ylimMinTemp <- quantile(segmentData$nrd,probs = .01)
  }


  seqRange <- seq(0,plotMaxCN,0.01)
  pPlot <- unlist(sapply(seqRange,function(x) ((1-tau))/(x*tau+2*(1-tau)))) # prob for the binomial
  plotRange <- 2+tau*(seqRange-2);

  # renormalize all of the probabilities by the sum of the probabilities
  shiftVals <- sapply(1:length(plotRange),
                      function(x) sum(((dpois(4:1000,lambdaMain*(plotRange[x]/2)))/sum(dpois(4:1000,lambdaMain*(plotRange[x]/2))))*testVals[(4:1000)-3,which.min(abs(seq(0.001,0.5,0.001)-pPlot[x]))]))
  medVals   <- sapply(1:length(plotRange),
                      function(x) sum(((dpois(4:1000,lambdaMain*(plotRange[x]/2)))/sum(dpois(4:1000,lambdaMain*(plotRange[x]/2))))*testVals[(4:1000)-3,which.min(abs(seq(0.001,0.5,0.001)-0.5))]))

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
                     function(x) sum(((dpois(4:1000,lambdaMain*(plotStarRange[x]/2)) )/sum(dpois(4:1000,lambdaMain*(plotStarRange[x]/2))))*testVals[(4:1000)-3,which.min(abs(seq(0.001,0.5,0.001)-pPlot[x]))]))
  # perfect het score value, for each coverage level
  medStarVals <- sapply(1:length(plotStarRange),
                        function(x) sum(((dpois(4:1000,lambdaMain*(plotStarRange[x]/2)))/sum(dpois(4:1000,lambdaMain*(plotStarRange[x]/2))))*testVals[(4:1000)-3,which.min(abs(seq(0.001,0.5,0.001)-0.5))]))

  # define "ploidy" for the sample
  ploidyNRD <- mean(cnvListChrFull)
  ploidyCN <- calcCopyNumber(NRD=ploidyNRD, tau)

  # this doesn't seem to be used for anything?
  if(FALSE){
    whichScore <- rep(0,length(cnvListChrFull))
    for(x in 1:length(starVals)) {
      whichScore[which(cnvListChrFull >  plotStarRange[x]*0.9 &
                         cnvListChrFull <  plotStarRange[x]*1.1 &
                         lohChrOutFull  > (starVals[x]/medStarVals[x])*0.9 &
                         lohChrOutFull  < (starVals[x]/medStarVals[x])*1.1)] <- 1
    }
    optScore <- length(which(whichScore==1))
  }

  xlimTemp <- min(0.6,0.9*min(starVals/medStarVals,na.rm = TRUE))
  ylimMinTemp <- min(ylimMinTemp, min(plotStarRange))

  ### one color for each chromosome 1-22
  segColors <- c(RColorBrewer::brewer.pal(9, 'Set1')[-c(6,9)] ,                     # remove the yellow and gray
                 'red', 'blue', 'cyan', 'gray45','magenta',                         # to get to 22
                 RColorBrewer::brewer.pal(8, 'Set2')[-c(8)], 'gray75',              # replace the gray with a slightly lighter color
                 'purple', "#00FF92FF"
  )

  # cnLohRatioA = lohContentA_maj2_min0
  # cnLohRatioB = lohContentB_maj1_min0
  # mcnGtEq2Ratio = lohContentC_maj2



  if(!is.null(segmentData)){
    starInfoTemp <- list(starVals=starVals,medStarVals=medStarVals,plotStarRange=plotStarRange)
    starLookUp <- makeStarLookUpTable(starInfo=starInfoTemp,percentTumor=tau*100)
    allelicSegments <- allelicCNV(starLookUp, segmentDataIn=segmentData)
    lohContent <- bmdSvPipeline:::getLohContent(allelicSegments)  ## getLohContent() is in ploidy.R
    lohContentA_maj2_min0 <- lohContent$lohContentA_maj2_min0
    loginfo('2N+loh content (lohContentA_maj2_min0): %s',round(lohContentA_maj2_min0,3))
  }else{
    lohContentA_maj2_min0 <- NA
  }

  # plot once with all chromosome data -------
  if(TRUE){
    # set up plot
    xlimits <- c(xlimTemp, 1.1)
    ylimits <- c(ylimMinTemp, max(plotRange))
    plot(lohChrOutFull,cnvListChrFull,type='n',
         xlim=xlimits,
         ylim=ylimits,
         xlab="Heterozygosity Score",
         ylab="NRD")

    if(!paperMode){
      # annotate peak locations with boxes to show the zones around the peaks
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
      abline(v=heterozygosityScoreThreshold, col='red')
      mtext(text=heterozygosityScoreThreshold, side=1, at = heterozygosityScoreThreshold, cex=.75, col='red', line=-0.2)
    }

    # add -actual- hetScore data
    points(lohChrOutFull, cnvListChrFull,
           pch='.',
           col=scales::alpha('black', 0.25))


    # add -median- hetScore PER test segment
    if(!is.null(segmentData) && !is.null(peakInfo)){
      # add median loh per test segment
      points(segmentData$lohScoreMedian, segmentData$nrd,
             pch=segmentData$chr,cex=0.9,
             col=scales::alpha(segColors[segmentData$chr],0.75))

      # legend
      legend("topleft", legend=c(1:22), col=segColors, pch=1:22, cex=.95)

      # peak rank labels
      if(!paperMode){
        ymax <- par('usr')[4]
        colorCodeRank <- ifelse(is.na(peakInfo$dPeaks),'gray80', 'gray50')
        mtext(text=peakInfo$rankByHeight, side=4, at=peakInfo$peakReadDepth_normX*2/rdNormX_2Npeak, las=2, line=-.5, col=colorCodeRank, adj=1, cex=.8)
        mtext(text='Peak', side=4, at=ymax*.98, col='gray50', las=1, line=-1.7, cex=.8)
        mtext(text='Rank', side=4, at=ymax*.96, col='gray50', las=1, line=-1.7, cex=.8)
      }
    }

    # LOH line
    lines(shiftVals/medVals,plotRange,type='l',xlim=c(0,1),col='purple')
    # stars for each heterozygosity ratio for each cn level - make sure this layers on top the segment data
    points(starVals/medStarVals,plotStarRange,pch="*",cex=2,xlim=c(0,1),col='green')

    # plot annotation
    mtext(paste0('lambdaMainOrig=',lambdaMainOrig,'   dipVal=',dipVal,'   mainPeakcnvBinnedNRD=',round(mainPeakcnvBinnedNRD,3)), side=3, line=2)

    mtext(c(sampleId, folderId), side=3, adj=c(0,1))
    mtext(side=1, text=paste('ploidy: ',round(ploidyCN,1)), adj=0, line=1.7)
    if(!is.null(bimaVersion)){
      mtext(1, text=bimaVersion, adj = 0, line=3.0, cex=.9, col='gray')
    }
    if(!is.na(lohContentA_maj2_min0)){
      mtext(1, text=round(lohContentA_maj2_min0,3), adj = 0, line=3.9, cex=.9, col='gray')
    }

    mtext(side=1, paste0('tumor: ',round(tau*100), '%'),    adj=1, line=1.7)

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

    # CN level annotations
    iCN <-0:plotMaxCN
    abline(h=2+((iCN-2)*tau), lty=2, col='gray85')
    # mtext(text='CN', side=4, col='black', las=1, line = 1.5)
    axis(side=4, at=2+((iCN-2)*tau), labels=paste0(iCN,'N'), tick=TRUE, las=1)

    # legend
    legend('bottomleft',legend = c('LOH', 'theoretical', 'actual'),col = c('purple',  'green', 'black'),lty=c(1,NA,NA), pch=c(NA, "*", "."),cex=.85,pt.cex=2)
  }

  # plot each chromosome individually -----
  if(plotEachChrom){
    for(chrNum in 1:22) {
      plot(lohChrOutFull[(chrEnd[chrNum]+1):chrEnd[chrNum+1]],cnvListChrFull[(chrEnd[chrNum]+1):chrEnd[chrNum+1]],
           pch='.',
           ylim=c(0,max(plotRange)),
           xlim=c(xlimTemp,1.1),
           xlab="Heterozygosity Score",
           ylab="NRD",
           col=scales::alpha('black',0.3),
           cex=2)
      lines(shiftVals/medVals,plotRange,type='l',xlim=c(0,1),col='purple')
      points(starVals/medStarVals,plotStarRange,pch="*",cex=2,xlim=c(0,1),col='green')
      title(paste0("tau=",tau,", dipVal=",dipVal,", chr=",chrNum))
      mtext(c(sampleId, folderId), side=3, adj=c(0,1))
      # CN level annotations
      for(i in 1:plotMaxCN) {
        # lines(c(0,1.1),c(2+((i-2)*tau),2+((i-2)*tau)),lty=2,col='green')
        abline(h= 2+((i-2)*tau), lty=2, col='green')
        text(x=1.1,y=2+((i-2)*tau),labels=as.character(i))
      }
    }
  }

  return(list(starVals=starVals,
              medStarVals=medStarVals,
              plotStarRange=plotStarRange,
              shiftVals=shiftVals,
              medVals=medVals,
              dipVal=dipVal,
              ploidyCN=ploidyCN,
              lohContent=lohContent,
              plotAxisLimits = list(hetScoreAxisLims=xlimits,
                                    nrdAxisLims=ylimits)
  ))
}


#' make a blank plot, no data, used as a stand in when combining lots of samples into one pdf
blankPlotStarsInTheClouds=function(
    noPdf=FALSE,
    heterozygosityScoreThreshold=0.98,
    tau=0.8,
    sampleId = 'BLANK',
    outputDir='/research/labs/experpath/vasm/shared/NextGen/johnsonsh/Routput/ploidy'){


  ### make the pdf
  if (!noPdf) {
    constellationPdfFile <- getTypedFile('constellationPlot',dir = outputDir, list(sampleId = sampleId))
    pdf(file = constellationPdfFile@path, paper="a4r", width=8, height=10, title=paste0('hetScVsCN_BLANK'))
    loginfo('writing pdf: %s', constellationPdfFile@path)
  }

  # set up plot
  plot(x=1:10,y=1:10,type='n',
       ylim=c(0, 5),
       xlim=c(.4, 1.1),
       xlab="Heterozygosity Score",
       ylab="NRD")

  ### vertical line to mark cutoff
  abline(v=heterozygosityScoreThreshold, col='red')
  mtext(text=heterozygosityScoreThreshold, side=1, at = heterozygosityScoreThreshold, cex=.75, col='red', line=-0.2)

  ### one color for each chromosome 1-22
  segColors <- c(RColorBrewer::brewer.pal(9, 'Set1')[-c(6,9)] ,                     # remove the yellow and gray
                 'red', 'blue', 'cyan', 'gray45','magenta',                         # to get to 22
                 RColorBrewer::brewer.pal(8, 'Set2')[-c(8)], 'gray75',              # replace the gray with a slightly lighter color
                 'purple', "#00FF92FF"
  )

  # legend
  legend("topleft", legend=c(1:22), col=segColors, pch=1:22, cex=.95)


  # CN level annotations
  iCN <-0:8
  abline(h=2+((iCN-2)*tau), lty=2, col='gray85')
  # mtext(text='CN', side=4, col='black', las=1, line = 1.5)
  axis(side=4, at=2+((iCN-2)*tau), labels=paste0(iCN,'N'), tick=TRUE, las=1)

  if (!noPdf) {
    dev.off()
  }

}
