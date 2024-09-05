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


  ## load two reference files  ---------------
  # NOTE/WARNING: if these files do not exist, an attempt will be made to download from Zenodo.
  hsNormMatFile <- "../reference/hetScoreNormMat.Rds"
  hsNormMat <- loadHsNormMat(hsNormMatFile)

  testValsFile <-  "../reference/testVals.Rds"
  testVals <- loadTestVals(testValsFile)


  # directory with input files:
  inputDir <- system.file('extdata', package = "BACDAC")
  readDepthPer30kbBin <-readRDS(file.path(inputDir, paste0(sampleId,'_','readDepthPer30kbBin.Rds')))
  readDepthBinSize  <- 30000


  # initialize
  starCloudPlotInputs <- NULL

  # result from calculateHetScore, but will use example data here
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

  ### OPTION 1
  ##### One Panel OPTION    --------------------
  ### constellation plot only
  if (!noPdf) {
    # we will be writing to this path, make sure it exists
    if(!dir.exists(file.path(outputDir))){
      dir.create(path = file.path(outputDir))
      loginfo('creating output directory for pdf: \n\t%s:', file.path(outputDir))
    }
    constellationPdfFile <- file.path(outputDir, paste0(sampleId, '_constellationPlot.pdf'))
    pdf(file = constellationPdfFile, paper="a4r", width=8, height=10,
        title=paste0('constellationPlot_',sampleId))
    loginfo('writing pdf: %s', constellationPdfFile)
  }

  op <- par(mfrow=c(1,1),mar=c(5,4,3.5,3.5),mgp=c(1.5, 0.5,0))
  starCloudResult <- plotStarsInTheClouds(
    sampleId, alternateId,starCloudPlotInputs, diploidPeakNRD=NULL,
    tau=min(1,calcPloidyResult$percentTumor/100),
    plotEachChrom=FALSE, mainPeakNRD=mainPeakNRD,
    segmentData=calcPloidyResult$segmentData, peakInfo=calcPloidyResult$peakInfo,
    digitalPeakZone =calcPloidyResult[['iterationStatsAll']][['digitalPeakZone']],
    addAnnotations = TRUE)
  par(op)


  ### OPTION 2
  ##### One Panel OPTION plus individual chromosome plot    --------------------
  ### constellation plot then chromosomes 1-22 plotted individually
  op <- par(mfrow=c(1,1),mar=c(5,4,3.5,3.5),mgp=c(1.5, 0.5,0))
  starCloudResult <- plotStarsInTheClouds(
    sampleId, alternateId,starCloudPlotInputs, diploidPeakNRD=NULL,
    tau=min(1,calcPloidyResult$percentTumor/100),
    plotEachChrom=TRUE, mainPeakNRD=mainPeakNRD,
    segmentData=calcPloidyResult$segmentData, peakInfo=calcPloidyResult$peakInfo,
    digitalPeakZone =calcPloidyResult[['iterationStatsAll']][['digitalPeakZone']],
    addAnnotations = TRUE)
  par(op)


  ### OPTION 3
  ##### Two Panel OPTION   --------------------
  ### constellation plot with linear genome plot to the right

  # load segmentation data
  segmentationFile <- file.path(inputDir, paste0(sampleId, '_segmentation.csv'))
  segmentation <- read.csv(segmentationFile, comment.char = '#')
  # check for required columns: # chr, start, end, rd and optionally cnvState
  segmentation <- checkSegmentation(segmentation)

# added the below because rgdObject couldn't be found
#  rgdObject = BACDAC:::rgdObject

  ### draw constellation plot left of the linear genome plot ----
  starCloudResult <- twoPanelReport(
    starCloudPlotInputs=starCloudPlotInputs, calcPloidyResult=calcPloidyResult,
    readDepthPer30kbBin=readDepthPer30kbBin,segmentation=segmentation,
    sampleId=sampleId, gainColor='blue', lossColor= 'red')


  # write segments to file
  fileToWrite=file.path(outputDir, paste0(sampleId, '_BACDAC_allelic_segments.csv'))
  loginfo('writing allelic segments to file: \n\t%s', fileToWrite)
  write.csv(starCloudResult$allelicSegments,file =fileToWrite,row.names = FALSE )

  if (!noPdf) {
    dev.off()  # CLOSE the pdf file
  }
}

