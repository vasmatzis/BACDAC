#' Combines data from lohDetection
#'
#' Outputs heterozygosity score by bin (.wig) and by arm (.csv) for each chromosome
#' To plot see \link{plotLohAnalysis} and the three graph plot function \link{makeLohFullReportPdf}
#'
#' @param minSnpsToCalculateStatistic  Minimum SNPs required in the window (1Mb) to calculate the statistic
#' @param samplingStep How frequently to try to summarize the data (bp), to produce overlapping windows
#' @param extraWindow Size of the window, how many bps to look at during each sampling step
#' @param maximumCoverage Do not process SNPs covered more than this
#' @param trimFromAlt Bins to discard from the 'alt' side of the distribution
#' @param trimFromRef Bins to discard from the 'ref' side of the distribution
#' @param trimExtraPerCoverage Fraction of bins to trim per each extra coverage
#' @inheritParams commonParameters
#'
#' @example inst/examples/lohAnalysisExample.R
#'
#' @export
lohAnalysis <- function(
    sampleId,
    rgd,
    cytoBandFile,
    outputDir,
    postProcessingDir = outputDir,
    noPdf = FALSE,
    maximumCoverage = 1000,
    trimFromAlt = 2,
    trimFromRef = 1,
    trimExtraPerCoverage = 0.1,
    minSnpsToCalculateStatistic = 20,
    samplingStep = 30000,
    extraWindow = 1000000
) {
  logCall()

  # maximumCoverage = 1000;  trimFromAlt = 2;  trimFromRef = 1;  trimExtraPerCoverage = 0.1;  minSnpsToCalculateStatistic = 20;  samplingStep = 30000;  extraWindow = 1000000

  cytoBands  <- loadIdeogram(path = rcfPrefix(cytoBandFile))
  rgdObject <- loadRgd(rgd)
  mainChroms <- unique(sort(c(svaAutosomes(rgdObject), svaAllosomes(rgdObject)))) # Allosomes

  chrYRefNumber <- refNameToNumber(rgdObject, "chrY")
  chrYSvaNumber <- genomeToBima(rgdObject, chrYRefNumber, 1)$svaNumber
  # We skip Y chromosome because LOH does not make much sense there
  mainChromsNoY <- mainChroms[mainChroms != chrYSvaNumber]

  coords <- getLinearCoordinates(rgdObject, mainChroms)

  # We will be writing these files out, let's check that we can
  lohAnalysisWigFileInfo <- getTypedFile("lohAnalysisWig", dir = outputDir, values=list(sampleId=sampleId))
  validateFileBeforeWrite(lohAnalysisWigFileInfo)

  lohPerArmFileInfo <- getTypedFile('lohPerArm', dir=outputDir, values = list(sampleId=sampleId))
  validateFileBeforeWrite(lohPerArmFileInfo)

  het_score_wigFileInfo <- getTypedFile("heterozygosity_score",dir = outputDir, values=list(sampleId=sampleId))  # genViz
  validateFileBeforeWrite(het_score_wigFileInfo)

  # We want to calculate our test statistic for a wide range of coverages
  # we can encounter in practice.
  valSave <- calculateLohTestStatisticPerCoverage(maximumCoverage, trimFromAlt, trimFromRef, trimExtraPerCoverage)

  seqListTotal <- list()
  seqValsTotal <- list()



  for (i in mainChromsNoY) {
    snpFull <- loadRdata(getTypedFile("lohSnpFull", dir = postProcessingDir, values=list(sampleId=sampleId, svaNumber=i), legacy = TRUE))
    countBPFull <- loadRdata(getTypedFile("lohCountBpFull", dir = postProcessingDir, values=list(sampleId=sampleId, svaNumber=i), legacy = TRUE))

    # Determine total coverage (how many times we see the ref/alt alleles)
    covVals <- countBPFull[['ref']] + countBPFull[['alt']]

    # this filtering does match the test statistic
    # We will gradually remove bins as coverage grows because the spilling
    # effect seems to touch extra bin each extra 10 coverages
    snpsToUse <- which(covVals >= 4 &
                         countBPFull[,'alt']>=trimFromAlt+floor(trimExtraPerCoverage*covVals) &
                         countBPFull[,'ref']>=trimFromRef+floor(trimExtraPerCoverage*covVals)
    )

    countBPFour <- countBPFull[snpsToUse,]
    snpFour <- snpFull[snpsToUse]
    covVals2 <- countBPFour[['ref']] + countBPFour[['alt']]

    seqList <- seq(1, coords@chromEnd[i]-coords@chromStart[i]+1, samplingStep)
    seqListTotal[[i]] <- seqList
    seqVals <- rep_len(0, length(seqList))
    for (seqI in seq_len(length(seqList))) {
      # Only retain SNPs of interest
      extraWindowShift <- -floor(extraWindow / 2) # Center the extra window around position of interest
      whichSeg <- which(
        snpFour>=(seqI-1)*samplingStep+1+extraWindowShift &
          snpFour<(seqI-1)*samplingStep+extraWindowShift+extraWindow
      )

      # Drop those SNPs that have coverage out of bounds
      covgOutOfBounds <- covVals2[whichSeg] > maximumCoverage
      if (sum(covgOutOfBounds)>0) {
        whichSeg <- whichSeg[-covgOutOfBounds]
      }

      # If we have at least this many data points, calculate the statistic
      if (length(whichSeg)>=minSnpsToCalculateStatistic) {
        minList <- pmin(countBPFour[whichSeg,'ref'], countBPFour[whichSeg,'alt'])

        # Sum up what we observed, divide by expected value based on coverage
        # Note that na.rm - the covVals2 can be > maximum supported coverage, which would
        # produce an NA when we try to get the valSave for those values
        seqVals[seqI] <- sum(minList) / sum(valSave[covVals2[whichSeg]], na.rm=TRUE)
      }
    }
    seqValsTotal[[i]] <- seqVals
  }

  # Save our data in a .wig format

  saveLohAnalysisToWig(
    wigFile = lohAnalysisWigFileInfo,
    seqListTotal = seqListTotal,
    seqValsTotal = seqValsTotal,
    rgdObject = rgdObject,
    chromsToSave = mainChromsNoY,
    samplingStep = samplingStep)

  # Summarize the data and write out
  centroArray <- getCentromerePositions(ideogram = cytoBands, rgd = rgdObject)
  lohPerArm <- lohSummary(seqValsTotal = seqValsTotal, centroArray = centroArray, coords = coords, chromosomes=mainChromsNoY)
  bmd.write.csv(lohPerArm, file=lohPerArmFileInfo)

  # plotLohAnalysis plotting code was copied to a function makeLohFullReportPdf(), and is now called in genomePlot.R

  # write GenViz output file
  hetScoreToWig(postProcessingDir_SV=postProcessingDir,sampleId_SV=sampleId,outputDir=outputDir,wsz = samplingStep)
  # loh_analysis_wigFileInfo=lohToGenViz(postProcessingDir,sampleId,outputDir=outputDir) # not needed

  loginfo("END OF SCRIPT")
}

#' Save the results of LOH analysis as a wig file
#'
#' @param wigFile \link{typedFile-class} for the wig to save
#' @param samplingStep Only uniformly sampled data can be used, use this sampling step
#' @inheritParams plotLohAnalysis
#'
#' @export
#' @family loh
saveLohAnalysisToWig <- function(wigFile, seqListTotal, seqValsTotal, rgdObject, chromsToSave,
                                 samplingStep) {
  # We go through GRanges object which is a bit of an overkill
  # but it allows us to use different formats than just wig if we wanted to
  refNumbers <- bimaToGenome(rgdObject, chromsToSave, rep(1, length(chromsToSave)))[['refNumber']]
  seqNames <- refNumberToName(rgdObject, refNumbers)
  data <- NULL
  for (i in chromsToSave) {
    part <- data.frame(seqname=seqNames[i], start=seqListTotal[[i]], value=seqValsTotal[[i]])
    data <- rbind(data, part)
  }

  grange <- GenomicRanges::GRanges(
    seqnames=data[['seqname']],
    ranges=IRanges::IRanges(start=data[['start']], width=samplingStep))
  BiocGenerics::score(grange) <- data[['value']]

  ensureDirExists(dirname(wigFile@path))
  rtracklayer::export.wig(object = grange, con = wigFile@path)
}

#' Load data from the wig file
#'
#' @return data.frame with seqnames, start, end, strand, score columns. Note that the strand column is not set (all values will be '*')
#'
#' @export
loadLohAnalysisFromWig <- function(postProcessingDir, sampleId, typeId = 'lohAnalysisWig') {
  data <- as.data.frame(
    rtracklayer::import.wig(
      getTypedFile(typeId = typeId, dir=postProcessingDir, values=list(sampleId=sampleId), legacy=TRUE)@path))
  return(data)
}

#' Make 3-subplot LOH analysis PDF
#'
#' Will show 3 separate plots contrasting CNV, heterozygosity score per 30K and heterozygosity score per arm.
#'
#' @inheritParams commonParameters
#' @export
makeLohFullReportPdf <- function(postProcessingDir, sampleId, outputDir,
                                 rgdObject,
                                 noPdf) {

  # Skip the Y chromosome
  mainChroms <- unique(sort(c(svaAutosomes(rgdObject), svaAllosomes(rgdObject)))) # Allosomes
  coords <- getLinearCoordinates(rgd = rgdObject, chromosomes = mainChroms)
  chrYRefNumber <- refNameToNumber(rgdObject, "chrY")
  chrYSvaNumber <- genomeToBima(rgdObject, chrYRefNumber, 1)$svaNumber
  # We skip Y chromosome because LOH does not make much sense there
  mainChromsNoY <- mainChroms[mainChroms != chrYSvaNumber]

  # load cnvIntervals
  cnvIntervalsFile <-getTypedFile("cnvIntervals",dir = postProcessingDir,values = list(sampleId = sampleId),legacy = TRUE )
  if(file.exists(cnvIntervalsFile@path)) {
    cnvIntervals <- bmd.read.csv(cnvIntervalsFile)
    cnvMetadata <- readMetadata(cnvIntervalsFile)
    if(!is.null(cnvMetadata[['normalPostProcessingDir']])) {
      normalPeakMethod <- cnvMetadata[['normalPeakMethod']]
    }else{
      logwarn('missing normalPeakMethod from cnvIntervals metadata')

    }
  }else{
    logerror('missing cnvIntervals file: %s',cnvIntervalsFile@path)
  }

  # load allelicSegData if possible
  allelicSegData <- NULL
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

  # TODO: does the graphicsDeviceOpen stuff have to be done here, like it is done for other plots? not sure what this feature does for us.
  if(!noPdf) {
    LOHPlotFile <- getTypedFile("lohFullReportPdf", dir=outputDir, values=list(sampleId=sampleId))
    ensureDirExists(dirname(LOHPlotFile@path))
    pdf(file=LOHPlotFile@path, width=11, height=8,  paper="a4r", title=paste0('heteroScore_',sampleId))
    loginfo('loh Full Report Pdf  file: %s',LOHPlotFile@path)
  }

  # georgeTempFile
  biallelicAnalysisFile <- getTypedFile(typeId = 'biallelicAnalysisWig', dir=postProcessingDir, values=list(sampleId=sampleId), legacy=TRUE)@path
  if(file.exists(biallelicAnalysisFile)){
    op <- par(mfrow=c(4,1),oma=c(0, 1, 3, 1), mar=c(2, 4, 0.5, 0))  # define an outer margin for placing a title using mtext
  }else{
    op <- par(mfrow=c(3,1),oma=c(0, 1, 3, 1), mar=c(2, 4, 0.5, 0))  # define an outer margin for placing a title using mtext
  }


  # Row 1: linear genome plot
  linearGenomePlot(
    postProcessingDir = postProcessingDir,
    rgd = rgdObject,
    sampleId=sampleId, # must provide in order to load other files (allelic, hetScore stuff)
    cnvIntervals=cnvIntervals,
    allelicSegData=allelicSegData)

  # annotate with sampleId and folderId in the upper left and upper right respectively
  mtext(sampleId,           side = 3, line= 1, outer=TRUE, cex= 1, adj=0)
  mtext(rgdObject$folderId, side = 3, line= 1, outer=TRUE, cex= 1, adj=1)

  # Row 2a: heterozygosity scores -per bin-
  lohdata <- loadLohAnalysisFromWig(postProcessingDir=postProcessingDir, sampleId=sampleId)
  plotLohAnalysis(lohdata, coords=coords, chromsToPlot = mainChromsNoY, rgdObject=rgdObject,allelicSegData=allelicSegData )

  if(file.exists(biallelicAnalysisFile)){
    # Row 2b: biallelic scores -per bin-
    biallelicData <- loadLohAnalysisFromWig(postProcessingDir=postProcessingDir, sampleId=sampleId, typeId='biallelicAnalysisWig')
    plotLohAnalysis(biallelicData, coords=coords, chromsToPlot = mainChromsNoY, rgdObject=rgdObject,allelicSegData=allelicSegData, ylab='biallelic score')

  }



  # Row 3: heterozygosity scores -per arm-
  lohPerArmFile <- getTypedFile('lohPerArm', dir=postProcessingDir, values = list(sampleId=sampleId), legacy=TRUE)
  lohPerArm <- bmd.read.csv(file=lohPerArmFile)
  lohSummaryPlot(lohPerArm=lohPerArm, rgdObject = rgdObject, coords=coords)

  # add indicator for BIMA version and use of BIMA indels
  bimaIndels <- extractSvaIndelColumnEnabled(commandLine=rgdObject$commandLine, file=rgdObject$file)
  bimaVersion <- rgdObject$processInformation$bimaVersion
  bimaInfo <- paste0(bimaVersion, ' indels:', bimaIndels)
  mtext(bimaInfo, side = 1, line= 1, outer=F, cex= .7, adj = 0)

  par(op)
  if (!noPdf) {
    dev.off()
  }
}

#' Create multiple plots for the LOH analysis using heterozygosity scores.
#'
#' These plots show the heterozygosity scores per bin, not per arm (summaries).
#' For per arm values, see \link{lohSummaryPlot}. \code{sampleId} and \code{sampleAlias} are optional,
#' for plot annotation only.
#'
#' @param lohdata heterozygosity scores as loaded from the loh wig file
#' @param coords The coordinate system description from RGD
#' @param chromsToPlot Vector of chromosome numbers to plot
#' @param addIndividualChrPlots option to plot each chromosome individually
#' @param rgdObject To convert chromosome number to text
#' @param yMap How to transform y coordinates for drawing purposes
#' @param allelicSegData ploidy segments with mean and median het scores and allelic ratios
#' @param ylab label for y axis default "Heterozygosity Score" refers to the low coverage method,
#' 'biallelic Score' refers to high coverage method
#'
#' @inheritParams commonParameters
#'
#' @export
#' @family loh
#' @example inst/examples/plotLohAnalysisExample.R
plotLohAnalysis <- function(lohdata, coords, chromsToPlot, rgdObject, sampleId=NULL, sampleAlias=NULL, addIndividualChrPlots=FALSE,
                            yMap=function(y) { y },ylab="Heterozygosity Score",allelicSegData=NULL) {
  # Make an overview plot
  # Get a plot started
  maxX <- max(coords@chromEnd[coords@maxcn])
  dataYRange <- c(0, 1)
  yRange <- yMap(dataYRange)
  yRange[2] <- yRange[1] + (yRange[2]-yRange[1]) * 1.1 # Add 10 % on the top
  plot(x=0, y=0, type="n",  xaxs="i", xlim=c(1,maxX), ylim=yRange,
       xaxt="n", yaxt="n", ylab=ylab, xlab='')

  # Custom Y axis
  ticksPosition <- seq(dataYRange[1], dataYRange[2], 0.2)
  axis(side = 2,
       at=yMap(ticksPosition),
       labels= ticksPosition)

  title(xlab='chromosome', line=0)

  ## option for Zebra bars for the chromosomes. Draw first so the dots can go over
  for(i in seq_len(length(chromsToPlot))) {
    if(i%%2==1){
      rect(xleft = coords@chromStart[i],
           ybottom= yRange[1],
           xright = coords@chromEnd[i],
           ytop   = yRange[2]*1.04,
           col=rgb(0.5,0.5,0.5,alpha=0.15), border=NA)
    }
  }

  for (i in chromsToPlot) {
    chromName <- convertChromToCharacter(i, rgdObject, withChrPrefix=TRUE)

    lohDataForChrom <- lohdata[
      lohdata[['seqnames']]==chromName,
      c('start', 'score')
    ]

    points(lohDataForChrom[['start']]+coords@chromStart[i]-1,
           yMap(lohDataForChrom[['score']]),
           col=ifelse(i%%2==0, 'black', 'black'),pch=".") # was alternating black, orange, but got rid of the orange and instead alternate shading the background
  }

  # TODO: add allelicSegData somehow
  # add purple lines for LOH segments
  if(!is.null(allelicSegData)){
    # we don't want to include the 1N segments
    minorZeroSegments <- allelicSegData[which(allelicSegData$minor==0 & allelicSegData$major>=2),]
    if(nrow(minorZeroSegments)>0){
      # convert to linear coordinates
      linPosStart <- abs(bimaToLinear(rgd=rgdObject,  svaNumber=minorZeroSegments$chr, svaPos=minorZeroSegments$start) )
      linPosEnd   <- abs(bimaToLinear(rgd=rgdObject,  svaNumber=minorZeroSegments$chr, svaPos=minorZeroSegments$end) )

      yloc <- par('usr')[3]/2 # to put it just below 0 but still on the plot
      # add allele=0 segments to plot and data.frame
      segments(x0 = linPosStart, y0 = yloc, x1 = linPosEnd, y1 = yloc, col = 'purple', lwd=3)
    }

  }


  ## annotations
  chrCharacters <- convertChromToCharacter(chromsToPlot,rgdObject = rgdObject) # required for X aka 23
  abline(h=yMap(1.0), col='green')
  axis(side=3, at=coords@chromEnd[chromsToPlot], labels=NA, lwd=0, lwd.ticks = 1, tck = 1, col='lightgray') # Draw ticks
  axis(side=3, at=(coords@chromEnd[chromsToPlot]+coords@chromStart[chromsToPlot])/2, line = -2, labels = chrCharacters, cex.axis=0.85, lwd=0, padj=0) # Draw labels
  # add if sampleId is provided
  if(!is.null(sampleId)) {
    title(main= sampleId,    adj= 0)
    title(main= sampleAlias, adj= 1)
  }

  ## include a plot per each chromosome
  if(addIndividualChrPlots){
    for(i in chromsToPlot) {
      chromName <- convertChromToCharacter(i, rgdObject, withChrPrefix=TRUE)
      lohDataForChrom <- lohdata[lohdata[['seqnames']]==chromName, c('start', 'score')]

      plot(lohDataForChrom[['start']]+coords@chromStart[i]-1, yMap(lohDataForChrom[['score']]),
           ylim=yRange, col='black',type='l', ylab=ylab,xlab='linear genome (top) and chromosome (bot) position (bp)', yaxt="n")
      abline(h=yMap(1.0), col='green3')
      mtext(text=paste('chr',chrCharacters[i]), side = 3)
      if(!is.null(sampleId)) {
        title(main= sampleId,    adj= 0)
        title(main= sampleAlias, adj= 1)
      }
      # add chromosome positions to x axis
      xlimMax <- bmdSvPipeline:::convertpos(coords@chromEnd[i], coords)[2]
      chrPosLabels <- round(seq(1,xlimMax, length.out = length(axTicks(3))),-6)
      axis(side = 1, at = axTicks(3), labels = chrPosLabels,line = 1,tick = FALSE)
    }
  }
}

#' Calculate the value of test statistic for LOH analysis
#'
#' The statistic is expected number of times we observe reference/alternate
#' allele (whichever count is smaller) for given coverage, assuming there is no
#' loss of heterozygosity, so the ref/alt observations are of equal probability.
#'
#' Let's say coverage is 4. We can observe 16 different combinations of
#' ref/alt on these 4 reads. These will contain:
#' 1 times - 0 ref alleles << problem
#' 4 times - 1 ref allele  << problem
#' 6 times - 2 ref alleles = 2 alt alleles
#' 4 times - 3 ref alleles
#' 1 times - 4 ref alleles << problem
#'
#' Now, there is a problem with 0/1 ref allele bins. For each SNP location we look at,
#' we have no idea if the person in question is actually heterozygous. So many times
#' we would encounter 0 alt / 0 ref situations. Since we do not know how many
#' SNPs the sample is heterozygous four, we cannot use these numbers.
#' Moreover, it appears that the bin #1 (1 ref allele) tends to get 'tainted'
#' because just by sheer luck sometimes you measure a wrong base pair even though
#' the sample is homozygous. So the values "spill" from bin #0 into bin #1.
#' Similar situation is on the other end (not symmetrical- due to mismappings on the left side)
#' Thus we define the \code{trimFromAlt} and \code{trimFromRef} values that will
#' let us skip calculating the test statistic if the bin is a problematic one.
#'
#' @param maximumCoverage We return an array with test statistics for coverage going
#' from 1 to maximumCoverage
#' @param trimFromAlt How many bins to ignore at the "alt" side
#' @param trimFromRef How many bins to ignore at the "ref" side
#' @param trimExtraPerCoverage Fraction of bins to trim for coverage.
#' If coverage == 10 and \code{trimExtraPerCoverage} is 0.1, we trim one extra bin
#'
#' @return Array with test statistic for each coverage from 1 to maximumCoverage
#'
#' @export
#' @family loh
calculateLohTestStatisticPerCoverage <- function(maximumCoverage, trimFromAlt=2, trimFromRef=1, trimExtraPerCoverage=0.1) {
  # This is binomial distribution (initially for coverage of 1)
  # We will use "Pascal's triangle" method of computing this distribution
  # for each coverage.
  # The binomial distribution tells us how many different ways can we observe
  # x amount of alt alleles. So for instance if coverage is 4, the array
  # should be 1, 4, 6, 4, 1 (the first value corresponds to 0 alt alleles)
  binomDist <- c(1, 1) # 1 times 0 alt, 1 times 1 alt

  # This is our resulting array
  valSave <- rep_len(0, maximumCoverage)

  for (i in seqFwd(1, maximumCoverage)) {
    # "Blank" the bins that are not to be considered
    actualDist <- binomDist
    nBins <- length(binomDist)
    actualDist[seqFwd(1, min(trimFromAlt + floor(trimExtraPerCoverage*i), nBins))] <- 0
    actualDist[seqFwd(max(1, nBins-trimFromRef + 1 + floor(trimExtraPerCoverage*i)), nBins)] <- 0

    # Calculate probability that we observe a given bin value
    testSum <- sum(actualDist)
    if (testSum == 0) {
      probOfBin <- actualDist # It is all zeroes
    } else {
      probOfBin <- actualDist/testSum
    }

    # testMin would be the reported value - minimum of #ref and #alt, for each bin
    startMin1 <- 0:i
    startMin2 <- i:0
    binValue <- pmin(startMin1, startMin2)

    # This is the expected value to be reported
    valSave[i] <- sum(binValue*probOfBin)

    # Determine the next level of the binomial distribution, Pascal's triangle way
    binomDist <- c(0, binomDist) + c(binomDist, 0)
  }

  return(valSave)
}

#' Summarize LOH levels per chromosome arm
#'
#' @param seqValsTotal Calculated values, 30K binned, list of one array per chromosome
#' @param centroArray 2D array, first dimension is chromosome number, second is 1=start, 2=end of centromere
#' @param coords Coordinate system from rgd
#' @param chromosomes List of chromosomes to calculate the value for (by default all but Y)
#' @param noPArm Do not return calculation for these p arms
#'
#' @return data.frame with chr, arm, val columns
#'
#' @export
#'
#' @family loh
lohSummary <- function(seqValsTotal, centroArray, coords, chromosomes=1:23, noPArm = c(13, 14, 15, 21, 22)) {
  numChromosomes <- length(chromosomes)

  pVals <- chromosomes * 0
  qVals <- chromosomes * 0

  binSize <- 30000

  densityRange <- c(0.01, 1.2) # This is where we smooth our values to look for peak
  densityN <- 5000 # How smoothly to estimate density

  for (i in chromosomes) {
    # Separate values for bins for p and q arms
    centromereStartBin <- floor(centroArray[i,1]/binSize)
    centromereEndBin <- floor(centroArray[i,2]/binSize)
    chromosomeEndBin <- floor((coords@chromEnd[i]-coords@chromStart[i]+1)/binSize)

    lohTempP <- (seqValsTotal[[i]])[seqFwd(1, centromereStartBin)]
    lohTempQ <- (seqValsTotal[[i]])[seqFwd(centromereEndBin, chromosomeEndBin)]

    denseP <- density(lohTempP,from=densityRange[1],to=densityRange[2],n=densityN)
    denseQ <- density(lohTempQ,from=densityRange[1],to=densityRange[2],n=densityN)

    pVals[i] <- denseP$x[which.max(denseP$y)]
    qVals[i] <- denseQ$x[which.max(denseQ$y)]

    # print(paste(i,denseP$x[which.max(denseP$y)],denseQ$x[which.max(denseQ$y)]))
  }

  # Serialize the p and q arm values into one long vector (1p, 1q, 2p, 2q, ... 23q)
  # The as.vector trick turns a data.frame into a long list by reading out values by column
  chrOut <- as.vector(rbind(chromosomes, chromosomes))
  armOut <- as.vector(rbind(rep("p", numChromosomes), rep("q", numChromosomes)))
  valOut <- as.vector(rbind(pVals,qVals))
  lohPerArm <- data.frame(chr=chrOut, arm=armOut, val=valOut)
  # Drop the missing p arms
  lohPerArm <- lohPerArm[lohPerArm[,"arm"]!="p" | (!(lohPerArm[,"chr"] %in% noPArm)),,drop=FALSE]
  rownames(lohPerArm) <- NULL # Otherwise indexing fails

  return(lohPerArm)
}

#' Plot the heterozygosity scores for evaluating LOH, summary values
#'
#' @param lohPerArm A matrix with chr, arm and val columns
#' @param rgdObject To convert chromosome number to text
#' @param lohMean Mean expected value
#' @param lohStDev standard deviation of the loh values
#' @param yMap A function that turns the actual y value into a position on screen
#'
#' @export
#'
#' @family loh
#'
#' @example inst/examples/lohSummaryPlotExample.R
lohSummaryPlot <- function(lohPerArm, rgdObject, coords, lohMean=0.9875, lohStDev=0.0125, yMap=function(y) { 2 ^ (10*y) }) {
  dataYRange <- c(0, 1) # Where the data falls

  extraPercentOnTop <- 0.16 # Take the range, add this much on top for labels

  yRange <- yMap(dataYRange) # Range of our plot
  yRangeWithLabel <- c(yRange[1], yRange[1] + (yRange[2]-yRange[1]) * (1+extraPercentOnTop))

  chromsToPlot <- unique(sort(as.integer(lohPerArm[,'chr'])))
  chrCharacters <- convertChromToCharacter(chromsToPlot,rgdObject = rgdObject) # required for X aka 23
  maxX <- max(coords@chromEnd[coords@maxcn])
  plot(c(0,0),c(0,0), type="n",  xaxs="i",
       xlim=c(1,maxX),
       ylim=yRangeWithLabel,
       xaxt="n",
       yaxt="n",
       ylab="Heterozygosity Score", xlab='')
  title(xlab='chromosome', line=0)
  ticksPosition <- seq(dataYRange[1], dataYRange[2], 0.05)
  axis(side = 2,
       at=yMap(ticksPosition),
       labels= ticksPosition)

  ## reference lines for LOH analysis
  abline(h=yMap(dataYRange[2]),      col='green',lwd=1) # This is where we should be in theory
  abline(h=yMap(lohMean - lohStDev),            col='blue', lwd=1) # 1 stdev from mean
  abline(h=yMap(lohMean - 2*lohStDev), col='red',  lwd=1) # 2 stdevs from mean


  ## option for vertical lines/ticks for the chromosomes
  # axis(side=3, at= coords@chromEnd[chromsToPlot], labels=NA, lwd=0, lwd.ticks = 1, tck = 1, col='darkgray') # Draw ticks

  ## option for Zebra bars for the chromosomes. Draw first so the dots can go over
  for(i in seq_len(length(chromsToPlot))) {
    if(i%%2==1){
      rect(xleft = coords@chromStart[i],
           ybottom= yRangeWithLabel[1],
           xright = coords@chromEnd[i],
           ytop   = yRangeWithLabel[2]*1.04,
           col=rgb(0.5,0.5,0.5,alpha=0.15), border=NA)
    }
  }

  # chr labels
  chromLabelPos <- (coords@chromEnd[chromsToPlot]+coords@chromStart[chromsToPlot])/2
  names(chromLabelPos) <- chromsToPlot

  pqLabelPos <- coords@chromStart[lohPerArm[,'chr']] +
    (coords@chromEnd[lohPerArm[,'chr']] - coords@chromStart[lohPerArm[,'chr']]) *
    ifelse(lohPerArm[,'arm']=='p', 0.25, 0.75)

  axis(side=3, at=(coords@chromEnd[chromsToPlot]+coords@chromStart[chromsToPlot])/2, line = -2, labels = chrCharacters, cex.axis=0.85, lwd=0, padj=0) # Draw labels
  axis(side=3, at=pqLabelPos, line=-2.6, labels = lohPerArm[,'arm'], cex.axis=0.66, tick = FALSE) # Draw labels
  # pq labels and points
  for(i in seq_len(nrow(lohPerArm))){
    xArm <- pqLabelPos[i]
    # Plot the values that are out of range as arrow pointing up/down
    points(x=xArm,
           y=yMap(pmin(dataYRange[2], pmax(dataYRange[1], lohPerArm[i,'val']))),
           ylim=yRange,
           pch=ifelse(lohPerArm[i,'val'] < dataYRange[1], 25,
                      ifelse(lohPerArm[i,'val'] > dataYRange[2], 24, 16)), col="darkmagenta", bg="magenta")
  }

}


#' Run LOH analysis
#'
#' @family command
#'
#' @export
lohAnalysisCmd <- function(args = commandArgs(TRUE)) {
  bmdCommand("lohAnalysis", args)
}
