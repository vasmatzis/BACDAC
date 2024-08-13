#' Combines ref/alt SNP data from each chromosome 1-22,X
#'
#' Outputs heterozygosity score by bin (.wig) and by arm (.csv) for each chromosome
#' To plot see \link{plotHetScorePerBin} and the three graph plot function \link{makeHetScoreReportPdf}
#'
#' @param sampleId sample identifier, used as prefix for all input and output file names \code{<sampleId>_snpVals_<chromosome>.Rdata} and \code{<sampleId>_countBP_<chromosome>.Rdata}
#' @param inputDir full path to snpVal and countBP Rdata
#' @param outputDir full path for output files
#' @param segmentation read depth data.frame with required columns: chr, start, end, rd; optional: cnvState for color coded linear genome plot
#' @param noPdf When TRUE, pdf files will not be generated, instead plots are drawn on default device
#' @param maximumCoverage Do not process SNPs covered more than this
#' @param trimFromAlt Bins to discard from the 'alt' side of the distribution
#' @param trimFromRef Bins to discard from the 'ref' side of the distribution
#' @param trimExtraPerCoverage Fraction of bins to trim per each extra coverage
#' @param minSnpsToCalculateStatistic  Minimum SNPs required in the window (1Mb) to calculate the statistic
#' @param samplingStep How frequently to try to summarize the data (bp), to produce overlapping windows
#' @param extraWindow Size of the window, how many bps to look at during each sampling step
#' @param readDepthBinnedData list of two arrays: one with (normalized) read count per each window, second with linear genome position of each window (masked windows have been removed)
#' @param readDepthBinSize size of each window
#' @param gainColor color to use for gains in linear genome plot, default is blue
#' @param lossColor color to use for losses in linear geneome plot, default is red
#' @param noDelAmpDetection do not color code deletions and gains in genome plot
#'
#' @examples
#' sampleId='TCGA-14-1402-02A_ds'
#' inputDir='/research/labs/experpath/vasm/shared/NextGen/Projects/MethodDev/MD66301/GRCh38/svar-1/loh'
#' outputDir='/research/labs/experpath/vasm/shared/NextGen/johnsonsh/Routput/BACDAC'
#'
#'
#'  sampleId='TCGA-14-1402-02A_ds'
#'  inputDir <- system.file('extdata', package = "BACDAC")
#'  outputDir <- tempdir()
#'  segmentationFile <- file.path(inputDir, paste0(sampleId, '_segmentation.csv'))
#'  segmentation <- read.csv(segmentationFile,comment.char = '#', header = TRUE)
#'  # check for required columns
#'  requiredColumns=c('chr', 'start', 'end','rd')
#'  missingColumnKey=which(!requiredColumns %in% names(segmentation))
#'  if(length(missingColumnKey)>0){
#'    logerror('missing required column: %s',requiredColumns[missingColumnKey])
#'  }
#'
#' @example inst/examples/calculateHetScoreExample.R
#'
#' @export
calculateHetScore <- function(
    sampleId,
    inputDir,
    outputDir=inputDir,
    segmentation=NULL,
    noPdf = FALSE,
    maximumCoverage = 1000,
    trimFromAlt = 2,
    trimFromRef = 1,
    trimExtraPerCoverage = 0.1,
    minSnpsToCalculateStatistic = 20,
    samplingStep = 30000,
    extraWindow = 1000000,

    # inputs just for makeHetScoreReportPdf
    readDepthBinnedData=NULL,  # if null will make an empty plot
    readDepthBinSize=30000,
    gainColor='blue',
    lossColor='red',
    noDelAmpDetection=FALSE  # FALSE = gray plot, TRUE= use gain and loss info to color plot
) {
  # inst/examples/lohAnalysisExample.R

  # sysdata loaded automatically: rgdObject, ideogram

  # maximumCoverage = 1000;  trimFromAlt = 2;  trimFromRef = 1;  trimExtraPerCoverage = 0.1;  minSnpsToCalculateStatistic = 20;  samplingStep = 30000;  extraWindow = 1000000

  # sampleId='TCGA-14-1402-02A_ds'; inputDir='/research/labs/experpath/vasm/shared/NextGen/Projects/MethodDev/MD66301/GRCh38/svar-1/loh'; outputDir='/research/labs/experpath/vasm/shared/NextGen/johnsonsh/Routput/BACDAC'
  # sampleId='TCGA-14-1402-02A_ds'; inputDir='/research/labs/experpath/vasm/shared/NextGen/johnsonsh/Routput/BACDAC/data'; outputDir='/research/labs/experpath/vasm/shared/NextGen/johnsonsh/Routput/BACDAC'

  mainChroms <- 1:24
  # We skip Y chromosome because hetScore does not make much sense there
  mainChromsNoY <- 1:23
  coords <- getLinearCoordinates(rgdObject, mainChroms)

  # we will be writing to this path, make sure it exists # TODO: do we need to check that the path is writable?
  if(!dir.exists(file.path(outputDir, 'reports'))){
    dir.create(path = file.path(outputDir, 'reports'))
    logdebug('creating output directory: \n\t%s:', file.path(outputDir, 'reports'))
  }
  # specify output file names
  hetScorePerArmFile <- file.path(outputDir, 'reports', paste0(sampleId, '_hetScorePerArm.csv'))
  hetScorePerBinWigFile <- file.path(outputDir, 'reports', paste0(sampleId, '_hetScorePerBin.wig.gz'))

  # We want to calculate our test statistic for a wide range of coverages
  # we can encounter in practice.
  valSave <- calculateHetScoreTestStatisticPerCoverage(maximumCoverage, trimFromAlt, trimFromRef, trimExtraPerCoverage)

  seqListTotal <- list()
  seqValsTotal <- list()


  # check to see if we are loading inputs from internal bmdSvPipeline users or from external BACDAC users
  inputDirIsNextGenProjects=ifelse(grepl('shared/NextGen/Projects', x=inputDir), TRUE, FALSE)
  loginfo('loading ref and alt counts from dir: %s', inputDir)

  for (i in mainChromsNoY) {
    loginfo('loading chrom %i',i)
    ichrChar=convertChromToCharacter(i)


    if(inputDirIsNextGenProjects){
      # loading .Rdata output from bmdSvPipeline using bmdSvPipeline functions
      snpFull= bmdTools::loadRdata(file.path(inputDir,    paste0(sampleId, '_snpVals_',i,'.Rdata')), verbose = TRUE) # snpFull
      countBPFull = bmdTools::loadRdata(file.path(inputDir,paste0(sampleId, '_countBP_',i,'.Rdata')), verbose = TRUE) # countBPFull

      iRefAltCount=data.frame('chr'=ichrChar, 'pos'=snpFull, 'ref'=countBPFull$ref, 'alt'=countBPFull$alt)
    }else{
      # loading BACDAC .Rds inputs
      iFile=file.path(inputDir, paste0(sampleId,'_','refAltCount_', ichrChar,'.Rds'))
      loginfo('%i loading %s',i, iFile)

      iRefAltCount = readRDS(file=iFile )
      countBPFull=iRefAltCount[,c('ref', 'alt')]
      snpFull=iRefAltCount[,'pos']
    }

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

  # Save hetScore binned data in a .wig format
  saveHetScoreToWig(
    wigFile = hetScorePerBinWigFile,
    seqListTotal = seqListTotal,
    seqValsTotal = seqValsTotal,
    chromsToSave = mainChromsNoY,
    samplingStep = samplingStep
  )

  # make and save hetScore arm data in a .csv format
  hetScorePerArm=makeAndSaveHetScorePerArm(
    hetScorePerArmFile,
    seqValsTotal,
    chromsToSave=mainChromsNoY,
    noPArm = c(13, 14, 15, 21, 22)
  )

  # make plot with segmentation, hetScore by bin, hetScore by arm
  makeHetScoreReportPdf(
    hetScorePerBinWigFile=hetScorePerBinWigFile,
    hetScorePerArmFile=hetScorePerArmFile,
    segmentation=segmentation,
    readDepthBinnedData=readDepthBinnedData,
    readDepthBinSize=readDepthBinSize,
    sampleId=sampleId,
    outputDir=outputDir,
    gainColor = gainColor, lossColor= lossColor,noDelAmpDetection=noDelAmpDetection,
    noPdf=noPdf)


  loginfo("END OF SCRIPT")
}

#' Save the binned results of \code{calculateHetScore} as a wig file
#'
#' @param wigFile full Path to wig output file containing the Heterozygosity Score per bin as calculated by \code{calculateHetScore}
#' @param seqListTotal Start position for each 30K bin, 1 based, list of one array per chromosome
#' @param seqValsTotal Calculated values, 30K binned, list of one array per chromosome
#' @param chromsToSave List of chromosomes to calculate the value for (by default all but Y)
#' @param samplingStep Only uniformly sampled data can be used, use this sampling step
#'
#' @family hetScore
saveHetScoreToWig <- function(wigFile, seqListTotal, seqValsTotal, chromsToSave,
                                 samplingStep) {
  # We go through GRanges object which is a bit of an overkill
  # but it allows us to use different formats than just wig if we wanted to
  seqNames = convertChromToCharacter(chromsToSave, rgdObject, withChrPrefix = TRUE)
  data <- NULL
  for (i in chromsToSave) {
    part <- data.frame(seqname=seqNames[i], start=seqListTotal[[i]], value=seqValsTotal[[i]])
    data <- rbind(data, part)
  }

  grange <- GenomicRanges::GRanges(
    seqnames=data[['seqname']],
    ranges=IRanges::IRanges(start=data[['start']], width=samplingStep))
  BiocGenerics::score(grange) <- data[['value']]

  if(!dir.exists(dirname(wigFile))){
    dir.create(path = dirname(wigFile))
    logdebug('creating output directory: \n\t%s:', dirname(wigFile))
  }
  rtracklayer::export.wig(object = grange, con = wigFile)
  logdebug('wrote hetScore per 30kb bin to wig file: \n\t%s', wigFile)
}


#' summarize the hetScore from bins to chromosome arm
#'
#' @param hetScorePerArmFile full path to Heterozygosity Score per arm csv file as created by \code{calculateHetScore}
#' @param seqValsTotal Calculated values, 30K binned, list of one array per chromosome
#' @param chromsToSave List of chromosomes to calculate the value for (by default all but Y)
#' @param noPArm Do not return calculation for these p arms
#'
#' @return data.frame with chr, arm, hetScore columns
#'
#' @family hetScore
makeAndSaveHetScorePerArm <- function(hetScorePerArmFile, seqValsTotal, chromsToSave=1:23, noPArm = c(13, 14, 15, 21, 22)) {
  # previously lohSummary and without the write to file

  #centroArray 2D array, first dimension is chromosome number, second is 1=start, 2=end of centromere
  centroArray <- getCentromerePositions(ideogram = ideogram)

  coords <- getLinearCoordinates(rgdObject, chromosomes = 1:24)
  numChromosomes <- length(chromsToSave)

  pVals <- chromsToSave * 0
  qVals <- chromsToSave * 0

  binSize <- 30000

  densityRange <- c(0.01, 1.2) # This is where we smooth our values to look for peak
  densityN <- 5000 # How smoothly to estimate density

  for (i in chromsToSave) {
    # Separate values for bins for p and q arms
    centromereStartBin <- floor(centroArray[i,1]/binSize)
    centromereEndBin <- floor(centroArray[i,2]/binSize)
    chromosomeEndBin <- floor((coords@chromEnd[i]-coords@chromStart[i]+1)/binSize)

    lohTempP <- (seqValsTotal[[i]])[seqFwd(1, centromereStartBin)]
    lohTempQ <- (seqValsTotal[[i]])[seqFwd(centromereEndBin, chromosomeEndBin)]

    denseP <- BiocGenerics::density(lohTempP,from=densityRange[1],to=densityRange[2],n=densityN)
    denseQ <- BiocGenerics::density(lohTempQ,from=densityRange[1],to=densityRange[2],n=densityN)

    pVals[i] <- denseP$x[which.max(denseP$y)]
    qVals[i] <- denseQ$x[which.max(denseQ$y)]

    # print(paste(i,denseP$x[which.max(denseP$y)],denseQ$x[which.max(denseQ$y)]))
  }

  # Serialize the p and q arm values into one long vector (1p, 1q, 2p, 2q, ... 23q)
  # The as.vector trick turns a data.frame into a long list by reading out values by column
  chrOut <- as.vector(rbind(chromsToSave, chromsToSave))
  armOut <- as.vector(rbind(rep("p", numChromosomes), rep("q", numChromosomes)))
  valOut <- as.vector(rbind(pVals,qVals))
  hetScorePerArm <- data.frame('chr'=chrOut, 'arm'=armOut, 'hetScore'=round(valOut,3))
  # Drop the missing p arms
  hetScorePerArm <- hetScorePerArm[hetScorePerArm[,"arm"]!="p" | (!(hetScorePerArm[,"chr"] %in% noPArm)),,drop=FALSE]
  rownames(hetScorePerArm) <- NULL # Otherwise indexing fails
  hetScorePerArm$chr=convertChromToCharacter(hetScorePerArm$chr)

  if(!dir.exists(dirname(hetScorePerArmFile))){
    dir.create(path = dirname(hetScorePerArmFile))
    loginfo('creating output directory: \n\t%s:', dirname(hetScorePerArmFile))
  }

  write.csv(hetScorePerArm, file=hetScorePerArmFile)
  loginfo('wrote hetScore per arm to csv file: \n\t%s', hetScorePerArmFile)

  return(hetScorePerArm)
}



#' Load hetScore as data.frame from the binned wig file
#' @param wigFile full path to Heterozygosity Score per bin wig file as created by \code{calculateHetScore}
#' @return data.frame with seqnames, start, end, strand, score columns. Note that the strand column is not set (all values will be '*')
#'
#' @export
loadHetScoreFromWig <- function(wigFile) {
  data <- as.data.frame(rtracklayer::import.wig(wigFile))
  return(data)
}

#' Make 3-subplot hetScore analysis PDF
#'
#' Will show 3 separate plots contrasting CNV, heterozygosity score per 30K and heterozygosity score per arm.
#'
#' @param segmentation read depth data.frame with required columns: chr, start, end, rd
#' @param hetScorePerBinWigFile full path to Heterozygosity Score per bin wig file as created by \code{calculateHetScore}
#' @param hetScorePerArmFile full path to Heterozygosity Score per arm csv file as created by \code{calculateHetScore}
#' @param segmentation read depth data.frame with required columns: chr, start, end, rd; optional: cnvState for color coded linear genome plot
#' @param sampleId sample id, will be used as the prefix for all input and output file names
#' @param outputDir full path for output files
#' @param readDepthBinnedData list of two arrays: one with (normalized) read count per each window, second with linear genome position of each window (masked windows have been removed)
#' @param readDepthBinSize size of each window
#' @param gainColor color to use for gains in linear genome plot, default is blue
#' @param lossColor color to use for losses in linear geneome plot, default is red
#' @param noDelAmpDetection do not color code deletions and gains in genome plot
#' @param noPdf When TRUE, pdf files will not be generated, instead plots are drawn on default device
#'
#' @examples
#' sampleId='TCGA-14-1402-02A_ds'
#' inputDir='/research/labs/experpath/vasm/shared/NextGen/johnsonsh/Rprojects/BACDAC/inst/extdata'
#' hetScoreDir='/research/labs/experpath/vasm/shared/NextGen/johnsonsh/Routput/BACDAC/reports'
#' hetScorePerArmFile <- file.path(hetScoreDir, paste0(sampleId, '_hetScorePerArm.csv'))
#' hetScorePerBinWigFile <- file.path(hetScoreDir, paste0(sampleId, '_hetScorePerBin.wig.gz'))
#' segmentationFile = file.path(inputDir, paste0(sampleId, '_segmentation.csv'))
#' segmentation <- read.csv(segmentationFile,comment.char = '#', header = TRUE); dim(segmentation)
#' thirtyKbFile=file.path(inputDir, paste0(sampleId,'_','readDepthPer30kbBin.Rds'))
#' readDepthBinnedData = readRDS(file=thirtyKbFile )

#'
#'  op <- par(mfrow=c(3,1),mai=c(.25,0.5, 0.3,0.25), mgp=c(2, .5, 0))
#' # default cnv coloring (by default) and annotations
#' makeHetScoreReportPdf(hetScorePerBinWigFile=hetScorePerBinWigFile,hetScorePerArmFile=hetScorePerArmFile,
#'                       segmentation=segmentation,readDepthBinnedData=readDepthBinnedData,
#'                       readDepthBinSize=readDepthBinnedData$windowSize,sampleId=sampleId,noPdf=TRUE)
#' @export
makeHetScoreReportPdf <- function(hetScorePerBinWigFile,
                                  hetScorePerArmFile,
                                  segmentation=NULL,
                                  readDepthBinnedData=NULL,
                                  readDepthBinSize=30000,
                                  sampleId,
                                  outputDir,
                                  gainColor = 'blue', lossColor= 'red', noDelAmpDetection=FALSE,
                                  noPdf=FALSE) {
  # previously named makeLohFullReportPdf and called from genomePlot in svaTools pipeline
  # gainColor = 'blue'; lossColor= 'red'; noDelAmpDetection=FALSE;  noPdf=TRUE

  mainChroms <- 1:24
  # We skip Y chromosome because hetScore does not make much sense there
  mainChromsNoY <- 1:23
  coords <- getLinearCoordinates(rgdObject, mainChroms)

  # check segmentation data (cnvIntervals)
  if(is.null(segmentation)) {
    logwarn('missing required input - segmentation - data frame, will not plot')
  }else{
    # check for required columns
    requiredColumns=c('chr', 'start', 'end','rd','cnvState')
    missingColumnKey=which(!requiredColumns %in% names(segmentation))
    if(length(missingColumnKey)>0){
      logerror('segmentation data is missing required column: %s',requiredColumns[missingColumnKey])
    }
  }

  # set up/open pdf
  if(!noPdf) {
    # we will be writing to this path, make sure it exists # TODO: do we need to check that the path is writable?
    if(!dir.exists(file.path(outputDir, 'reports'))){
      dir.create(path = file.path(outputDir, 'reports'))
      loginfo('creating output directory for hetScoreReport PDF: \n\t%s:', file.path(outputDir, 'reports'))
    }
    hetScoreReportPdf <- file.path(outputDir, 'reports', paste0(sampleId, '_hetScoreReport.pdf'))
    pdf(file=hetScoreReportPdf, width=11, height=8,  paper="a4r", title=paste0('hetScoreReport_',sampleId))
    logdebug('writing hetScore report to PDF: \n\t%s:', hetScoreReportPdf)
  }

  op <- par(mfrow=c(3,1),oma=c(0, 1, 3, 1), mar=c(2, 4, 0.5, 0))  # define an outer margin for placing a title using mtext

  # Row 1: linear genome plot ----
  if(!is.null(readDepthBinnedData)){
    linearGenomePlot( readDepthBinnedData, wsz=readDepthBinSize, segmentation=segmentation,
                      sampleId=NULL, noDelAmpDetection = noDelAmpDetection,
                      gainColor = gainColor, lossColor= lossColor)
  }else{
    plotEmptyLinearGenomePlot(chromsToPlot=mainChromsNoY,coords=coords)
  }

  # annotate with title and sampleId in the upper left
  title(main= 'Heterozygosity Score Report', outer = TRUE)
  mtext(sampleId,side=3, adj=0)
  # mtext(sampleId, side = 3, line= 1, outer=TRUE, cex= 1, adj=0)


  # Row 2: heterozygosity scores -per bin- ----
  hetScore <- loadHetScoreFromWig(hetScorePerBinWigFile)
  plotHetScorePerBin(hetScore,
                     chromsToPlot = mainChromsNoY,
                     ylab="Heterozygosity Score by bin" ) # aka plotLohAnalysis()


  # Row 3: heterozygosity scores -per arm- ----
  hetScorePerArm <- read.csv(file=hetScorePerArmFile, header = TRUE, comment.char = '#')
  plotHetScorePerArm(hetScorePerArm=hetScorePerArm,
                     chromsToPlot=mainChromsNoY,
                     ylab="Heterozygosity Score by arm") # aka lohSummaryPlot()

  par(op)
  if (!noPdf) {
    dev.off()
  }
}




#' Plot the binned heterozygosity score
#'
#' These plots show the heterozygosity scores per bin.
#' For per arm values, see \link{plotHetScorePerArm}. \code{sampleId} is optional,
#' for plot annotation only.
#'
#' @param hetScore binned heterozygosity scores as loaded from the hetScore wig file
#' @param chromsToPlot Vector of chromosome numbers to plot
#' @param sampleId sample id, will be used as the prefix for all input and output file names
#' @param yMap A function that turns the actual y value into a position on screen, transform y coordinates for drawing purposes
#' @param segmentation segments with major and minor allelic specific copy number
#' @param ylab label for y axis default "Heterozygosity Score"
#'
#' @export
#' @family hetScore

plotHetScorePerBin <- function(hetScore, chromsToPlot, sampleId=NULL,
                            yMap=function(y) { y },
                            ylab="Heterozygosity Score",
                            segmentation=NULL) {
  # @example inst/examples/plotLohAnalysisExample.R

  mainChroms <- 1:24
  coords <- getLinearCoordinates(rgdObject, mainChroms)

  # Make an overview plot
  # Get a plot started
  maxX <- max(coords@chromEnd[coords@maxcn])
  dataYRange <- c(0, 1)
  yRange <- yMap(dataYRange)
  yRange[2] <- yRange[1] + (yRange[2]-yRange[1]) * 1.1 # Add 10 % on the top
  plot(x=0, y=0, type="n",
       xaxs="i",
       xlim=c(1,maxX), ylim=yRange,
       xaxt="n", yaxt="n",
       xlab='', ylab=ylab)

  # Custom Y axis
  ticksPosition <- seq(dataYRange[1], dataYRange[2], 0.2)
  axis(side = 2,
       at=yMap(ticksPosition),
       labels= ticksPosition)
  title(xlab='chromosome', line=0, cex.lab=1.5)


  ## option for Zebra bars for the chromosomes. Draw first so the dots can go over top
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

    hetScoreForChrom <- hetScore[
      hetScore[['seqnames']]==chromName,
      c('start', 'score')
    ]

    points(hetScoreForChrom[['start']]+coords@chromStart[i]-1,
           yMap(hetScoreForChrom[['score']]),
           col=ifelse(i%%2==0, 'black', 'black'),pch=".") # was alternating black, orange, but got rid of the orange and instead alternate shading the background
  }

  # TODO: test for minor, major columns segmentation
  # add purple lines for LOH segments
  if(!is.null(segmentation)){
    # we don't want to include the 1N segments
    minorZeroSegments <- segmentation[which(segmentation$minor==0 & segmentation$major>=2),]
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
  axis(side=3, at=coords@chromEnd[chromsToPlot], labels=NA, lwd=0, lwd.ticks = 1, tck = 1, col='lightgray',cex.axis=1.2) # Draw ticks/lines
  axis(side=3, at=(coords@chromEnd[chromsToPlot]+coords@chromStart[chromsToPlot])/2, line = -2, labels = chrCharacters, cex.axis=1.2, lwd=0, padj=0) # Draw labels
  # add if sampleId is provided
  if(!is.null(sampleId)) {
    title(main= sampleId,    adj= 0)
  }

}


#' Plot the heterozygosity scores for evaluating LOH, summary values
#'
#' @param hetScorePerArm A matrix with chr, arm and hetScore columns
#' @param chromsToPlot Vector of chromosome numbers to plot
#' @param sampleId sample id, will be used as the prefix for all input and output file names
#' @param hetScoreMean Mean expected hetScore
#' @param hetScoreStDev standard deviation of the hetScore values
#' @param yMap A function that turns the actual y value into a position on screen, transform y coordinates for drawing purposes
#' @param ylab label for y axis
#'
#' @export
#' @family hetScore
plotHetScorePerArm <- function(hetScorePerArm, chromsToPlot,sampleId=NULL,
                               hetScoreMean=0.9875, hetScoreStDev=0.0125,
                               yMap=function(y) { 2 ^ (10*y) },
                               ylab="Heterozygosity Score"
) {
  # @example inst/examples/hetScoreSummaryPlotExample.R
  mainChroms <- 1:24
  coords <- getLinearCoordinates(rgdObject, mainChroms)


  # Make an overview plot
  # Get a plot started
  maxX <- max(coords@chromEnd[coords@maxcn])
  dataYRange <- c(0, 1) # Where the data falls
  yRange <- yMap(dataYRange) # Range of our plot
  extraPercentOnTop <- 0.16 # Take the range, add this much on top for labels
  yRangeWithLabel <- c(yRange[1], yRange[1] + (yRange[2]-yRange[1]) * (1+extraPercentOnTop))
  plot(x=0, y=0, type="n",
       xaxs="i",
       xlim=c(1,maxX), ylim=yRangeWithLabel,
       xaxt="n", yaxt="n",
       xlab='', ylab=ylab)

  # Custom Y axis
  ticksPosition <- seq(dataYRange[1], dataYRange[2], 0.05)
  axis(side = 2,
       at=yMap(ticksPosition),
       labels= ticksPosition)
  title(xlab='chromosome', line=0, cex.lab=1.5)


  ## reference lines for hetScore analysis
  abline(h=yMap(dataYRange[2]),      col='green',lwd=1)                      # This is where we should be in theory
  abline(h=yMap(hetScoreMean - hetScoreStDev),            col='blue', lwd=1) # 1 stdev from mean
  abline(h=yMap(hetScoreMean - 2*hetScoreStDev), col='red',  lwd=1)          # 2 stdevs from mean


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
  ## option for vertical lines/ticks for the chromosomes
 axis(side=3, at= coords@chromEnd[chromsToPlot], labels=NA, lwd=0, lwd.ticks = 1, tck = 1, col='darkgray') # Draw ticks/lines

  pqLabelPos <- coords@chromStart[hetScorePerArm[,'chr']] +
    (coords@chromEnd[hetScorePerArm[,'chr']] - coords@chromStart[hetScorePerArm[,'chr']]) *
    ifelse(hetScorePerArm[,'arm']=='p', 0.25, 0.75)

  ## annotations
  chrCharacters <- convertChromToCharacter(chromsToPlot) # required for X aka 23
  axis(side=3, at=(coords@chromEnd[chromsToPlot]+coords@chromStart[chromsToPlot])/2, line = -2, labels = chrCharacters, cex.axis=1.2, lwd=0, padj=0) # Draw labels
  axis(side=3, at=pqLabelPos, line=-2.6, labels = hetScorePerArm[,'arm'], cex.axis=0.85, tick = FALSE) # Draw labels
  # add if sampleId is provided
  if(!is.null(sampleId)) {
    title(main= sampleId,    adj= 0)
  }

  # pq labels and points
  for(i in seq_len(nrow(hetScorePerArm))){
    xArm <- pqLabelPos[i]
    # Plot the values that are out of range as arrow pointing up/down
    points(x=xArm,
           y=yMap(pmin(dataYRange[2], pmax(dataYRange[1], hetScorePerArm[i,'hetScore']))),
           ylim=yRange,
           pch=ifelse(hetScorePerArm[i,'hetScore'] < dataYRange[1], 25,
                      ifelse(hetScorePerArm[i,'hetScore'] > dataYRange[2], 24, 16)), col="darkmagenta", bg="magenta")
  }

}



#' Calculate the value of test statistic for HetScore analysis
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
#' @family hetScore
calculateHetScoreTestStatisticPerCoverage <- function(maximumCoverage, trimFromAlt=2, trimFromRef=1, trimExtraPerCoverage=0.1) {
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

