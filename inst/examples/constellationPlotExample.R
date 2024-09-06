\dontrun{
  library(BACDAC)
  library(logging)

  # This example first loads the inputs for plotting, then shows three options for
  # creating the constellation plot. Lastly, output allelic segments to file.

  noPdf <- TRUE      # TRUE= print to screen, FALSE=print to pdf in outputDir

  # directory for output files:
  outputDir <- tempdir()

  sampleId <- 'TCGA-14-1402-02A_ds';
  alternateId <- 66301

  ## LOAD DATA ------------------
  ### 1 - load two reference files
  # NOTE/WARNING: if these files do not exist, an attempt will be made to download from Zenodo.
  hsNormMatFile <- "./referenceFiles/hetScoreNormMat.Rds"
  hsNormMat <- loadHsNormMat(hsNormMatFile)

  testValsFile <-  "./referenceFiles/testVals.Rds"
  testVals <- loadTestVals(testValsFile)


  # directory with input files:
  inputDir <- system.file('extdata', package = "BACDAC")
  ### 2 - load read depth data
  readDepthPer30kbBin <-readRDS(file.path(inputDir, paste0(sampleId,'_','readDepthPer30kbBin.Rds')))
  readDepthBinSize  <- 30000


  # initialize
  starCloudPlotInputs <- NULL

  ### load hetScore data
  # hetScorePerBinWigFile is the result from calculateHetScore, but will load from example data here
  hetScorePerBinWigFile <- file.path(inputDir, paste0(sampleId, '_hetScorePerBin.wig.gz'))

  # result from calculatePloidy, but will use example data here
  calcPloidyResultOutputFile <- file.path(inputDir, paste0(sampleId, '_calculatePloidyResult.Rds'))
  calcPloidyResult <- readRDS(calcPloidyResultOutputFile)
  mainPeakNRD <- getMainPeakNRD(calcPloidyResult)
  expReadsIn2NPeak_1bp <- calcPloidyResult$expReadsIn2NPeak_1bp

  ### get the needed input values for the plot
  if(is.null(starCloudPlotInputs)){     # run once only... it takes about 3-5 minutes
    starCloudPlotInputs <- loadStarsInTheClouds(
      sampleId=sampleId, inputDir=inputDir, readDepthPer30kbBin=readDepthPer30kbBin,
      hetScorePerBinFile=hetScorePerBinWigFile, readDepthBinSize=readDepthBinSize,
      hsNormMat, testVals=testVals, mainPeakNRD=mainPeakNRD,
      expReadsIn2NPeak_1bp=expReadsIn2NPeak_1bp)
  }

  # MULTIPLE PLOTTING OPTIONS -----------

  if (!noPdf) {
    # we will be writing to this path, make sure it exists
    if(!dir.exists(file.path(outputDir))){
      dir.create(path = file.path(outputDir))
      loginfo('creating output directory for pdf: \n\t%s:', file.path(outputDir))
    }
    constellationPdfFile <- file.path(outputDir, paste0(sampleId, '_BACDAC_ConstellationPlot.pdf'))
    pdf(file = constellationPdfFile, paper="a4r", width=8, height=10,
        title=paste0('constellationPlot_',sampleId))
    loginfo('opening constellation plot pdf: %s', constellationPdfFile)
  }

  ### OPTION 1
  ### One Panel --------------------
  ### constellation plot only

  op <- par(mfrow=c(1,1),mar=c(3.5,3.5,1,2),mgp=c(1.5,0.5,0),oma=c(0,0,2,0)+.2)
  starCloudResult <- plotStarsInTheClouds(
    sampleId, alternateId,starCloudPlotInputs, diploidPeakNRD=NULL,
    tau=min(1,calcPloidyResult$percentTumor/100),
    plotEachChrom=FALSE, mainPeakNRD=mainPeakNRD,
    segmentData=calcPloidyResult$segmentData, peakInfo=calcPloidyResult$peakInfo,
    digitalPeakZone =calcPloidyResult[['iterationStatsAll']][['digitalPeakZone']],
    addSegmentLegend = TRUE)
  par(op)


  ### OPTION 2
  ### One Panel plus individual chromosome plot    --------------------
  ### constellation plot then chromosomes 1-22 plotted in a 4x4 matrix of individual plots
  op <- par(mfrow=c(1,1),mar=c(3.5,3.5,1,2),mgp=c(1.5,0.5,0),oma=c(0,0,2,0)+.2)

  starCloudResult <- plotStarsInTheClouds(
    sampleId, alternateId,starCloudPlotInputs, diploidPeakNRD=NULL,
    tau=min(1,calcPloidyResult$percentTumor/100),
    plotEachChrom=TRUE, mainPeakNRD=mainPeakNRD,
    segmentData=calcPloidyResult$segmentData, peakInfo=calcPloidyResult$peakInfo,
    digitalPeakZone =calcPloidyResult[['iterationStatsAll']][['digitalPeakZone']],
    addSegmentLegend = TRUE)
  par(op)

  if (!noPdf) {
    loginfo('closing pdf: %s', constellationPdfFile)
    dev.off()  # CLOSE the pdf file for option 1 and 2, option 3 will have a separate file
  }

  ### OPTION 3
  ### Two Panel --------------------
  ### constellation plot in left panel, linear genome plot in right panel

  # load segmentation data
  segmentationFile <- file.path(inputDir, paste0(sampleId, '_segmentation.csv'))
  segmentation <- read.csv(segmentationFile, comment.char = '#')
  # check for required columns: # chr, start, end, rd and optionally cnvState
  segmentation <- checkSegmentation(segmentation)

  starCloudResult <- twoPanelReport(
    starCloudPlotInputs=starCloudPlotInputs, calcPloidyResult=calcPloidyResult,
    readDepthPer30kbBin=readDepthPer30kbBin,segmentation=segmentation,
    sampleId=sampleId, alternateId=alternateId, gainColor='blue', lossColor= 'red',
    noPdf=noPdf,outputDir=outputDir)




  # WRITE TO FILE --------
  # write segments in constellation plot to file
  fileToWrite=file.path(outputDir, paste0(sampleId, '_BACDAC_allelic_segments.csv'))
  loginfo('writing allelic segments to file: \n\t%s', fileToWrite)
  write.csv(starCloudResult$allelicSegments,file =fileToWrite,row.names = FALSE )

}

