\dontrun{
  library(BACDAC)
  library(logging)

  # initialize
  starCloudPlotInputs=NULL

  noPdf=TRUE

  ## load two reference files  ---------------
  # hsNormMat/lohMat: LOH analysis mask, used to look for places in 23 TCGA normals where more than half dropped below the a (i.e. 0.975) cutoff.
  # testVals: used to find each possible heterozygosity value for each copy number level (find the right spots for the stars)

  ## load two reference files  ---------------
  # if these files do not exist, an attempt will be made to download them from Zenodo.
  hsNormMatFile <- "../reference/hetScoreNormMat.Rds"
  hsNormMat = loadHsNormMat(hsNormMatFile, destfile=hsNormMatFile)

  testValsFile <-  "../reference/testVals.Rds"
  myTestVals = loadTestVals(testValsFile,destfile=testValsFile)

  sampleId='TCGA-14-1402-02A_ds'; alternateId=66301

  # directory with input files:
  inputDir <- system.file('extdata', package = "BACDAC")
  readDepthPer30kbBin=  readRDS(file.path(inputDir, paste0(sampleId,'_','readDepthPer30kbBin.Rds')))
  readDepthBinSize=30000

  # directory for output files:
  outputDir=tempdir()

  # result from calculateHetScore, but will use example data here
  hetScorePerBinWigFile <- file.path(inputDir, paste0(sampleId, '_hetScorePerBin.wig.gz'))

  # result from calculatePloidy, but will use example data here
  calcPloidyResultOutputFile=file.path(inputDir, paste0(sampleId, '_calculatePloidyResult.Rds'))
  calcPloidyResult = readRDS(calcPloidyResultOutputFile)
  mainPeakNRD=getMainPeakNRD(calcPloidyResult)
  expReadsIn2NPeak_1bp=calcPloidyResult$expReadsIn2NPeak_1bp

  ### get the needed input values for the plot
  if(is.null(starCloudPlotInputs)){     # run once only... it takes about 3-5 minutes
    starCloudPlotInputs=loadStarsInTheClouds(
      sampleId=sampleId, inputDir=inputDir, readDepthPer30kbBin=readDepthPer30kbBin,
      hetScorePerBinWigFile=hetScorePerBinWigFile, readDepthBinSize=readDepthBinSize,
      hsNormMat, testVals, mainPeakNRD=mainPeakNRD, expReadsIn2NPeak_1bp=expReadsIn2NPeak_1bp)
  }

  ##### One Panel OPTION --------------------
  if (!noPdf) {
    # we will be writing to this path, make sure it exists
    if(!dir.exists(file.path(outputDir))){
      dir.create(path = file.path(outputDir))
      loginfo('creating output directory for pdf: \n\t%s:', file.path(outputDir))
    }
    constellationPdfFile <- file.path(outputDir, paste0(sampleId, '_constellationPlot.pdf'))
    pdf(file = constellationPdfFile, paper="a4r", width=8, height=10, title=paste0('constellationPlot_',sampleId))
    loginfo('writing pdf: %s', constellationPdfFile)
  }

  op <- par(mfrow=c(1,1),mar=c(5,4,3.5,3.5),mgp=c(1.5, 0.5,0))
  starCloudResult=plotStarsInTheClouds(
    sampleId, alternateId,starCloudPlotInputs, diploidPeakNRD=NULL, tau=min(1,calcPloidyResult$percentTumor/100),
    plotEachChrom=FALSE, mainPeakNRD=mainPeakNRD,
    segmentData=calcPloidyResult$segmentData, peakInfo=calcPloidyResult$peakInfo,
    digitalPeakZone =calcPloidyResult[['iterationStatsAll']][['digitalPeakZone']],
    addAnnotations = TRUE)
  par(op)

  if (!noPdf) {
    dev.off()
  }

  ##### One Panel OPTION with individual chromosome plots   --------------------
  op <- par(mfrow=c(1,1),mar=c(5,4,3.5,3.5),mgp=c(1.5, 0.5,0))
  starCloudResult=plotStarsInTheClouds(
    sampleId, alternateId,starCloudPlotInputs, diploidPeakNRD=NULL, tau=min(1,calcPloidyResult$percentTumor/100),
    plotEachChrom=TRUE, mainPeakNRD=mainPeakNRD,
    segmentData=calcPloidyResult$segmentData, peakInfo=calcPloidyResult$peakInfo,
    digitalPeakZone =calcPloidyResult[['iterationStatsAll']][['digitalPeakZone']],
    addAnnotations = TRUE)
  par(op)


  ##### Two Panel OPTION   --------------------
  ### plot constellation plot with linear plot

  # load segmentation data
  segmentationFile <- file.path(inputDir, paste0(sampleId, '_segmentation.csv'))
  segmentation <- read.csv(segmentationFile, comment.char = '#') # chr, start, end, rd
  segmentation <- checkSegmentation(segmentation)

  ### draw constellation plot left of the linear genome plot ----
  starCloudResult = twoPanelReport(
    starCloudPlotInputs=starCloudPlotInputs, calcPloidyResult=calcPloidyResult,
    readDepthPer30kbBin=readDepthPer30kbBin,segmentation=segmentation,
    sampleId=sampleId, gainColor='blue', lossColor= 'red')

  # write segments to file
  fileToWrite=file.path(outputDir, paste0(sampleId, '_bacdacAllelicSegments.csv'))
  write.csv(starCloudResult$allelicSegments,file =fileToWrite,row.names = FALSE )

}

