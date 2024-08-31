
  library(BACDAC)
  library(logging)
  # calculatePloidyExample.R
  basicConfig("DEBUG")
  noPdf=TRUE          # TRUE= print to screen, FALSE=print to pdf (i e. outputDir/dev/ploidy)
  outputDir = NULL    # output folder for pdfs, only needed if noPdf=FALSE

  sampleId='TCGA-14-1402-02A_ds'; alternateId=66301

  ### load data ###

  # if this file does not exist, an attempt will be made to download from Zenodo.
  hsNormMatFile <- "../reference/hetScoreNormMat.Rds"
  hsNormMat = loadHsNormMat(hsNormMatFile)

  exampleDataDir <- system.file('extdata', package = "BACDAC")
  inputDir <- exampleDataDir


  # segmentation data
  segmentationFile <- file.path(inputDir, paste0(sampleId, '_segmentation.csv'))
  segmentation <- read.csv(segmentationFile, comment.char = '#') # chr, start, end, rd per
  segmentation <- checkSegmentation(segmentation)

  # read depth data
  thirtyKbFile=file.path(inputDir, paste0(sampleId,'_','readDepthPer30kbBin.Rds'))
  readDepthPer30kbBin = readRDS(file=thirtyKbFile )
  hundredKbFile=file.path(inputDir, paste0(sampleId,'_','readDepthPer100kbBin.Rds'))
  readDepthPer100kbBin = readRDS(file=hundredKbFile )

  # hetScore data - the output from calculateHetScore()
  # hetScoreDir is typically the outputDir specified in `calculateHetScore` but here we will load
  # from the package example data.
  hetScoreDir <- exampleDataDir
  hetScorePerBinWigFile <- file.path(hetScoreDir, paste0(sampleId, '_hetScorePerBin.wig.gz'))
  hetScoreData <- loadHetScoreFromWig(hetScorePerBinWigFile)

  ### defaults ###
  segmentationBinSize=30000; numChroms=24;
  omitAnnotations = FALSE;
  dPeaksCutoff=0.01;    penaltyCoefForAddingGrids=0.49; minGridHeight=0.2; minPeriodManual=-1;
  maxPeriodManual=-1; forceFirstDigPeakCopyNum=-1;
  grabDataPercentManual= -1; origMaxPercentCutoffManual=-1;
  minReasonableSegmentSize=5.5e6;
  heterozygosityScoreThreshold=0.98;
  allowedTumorPercent = 106


  ### call calculatePloidy, the function to do all the ploidy work ---------
  loginfo('calculate ploidy for %s ', sampleId)
  calcPloidyResult=calculatePloidy(sampleId=sampleId, outputDir = outputDir, noPdf=noPdf, alternateId=alternateId,
                         readDepthPer30kbBin = readDepthPer30kbBin, readDepthPer100kbBin= readDepthPer100kbBin,
                         segmentation=segmentation, hetScoreData = hetScoreData,
                         segmentationBinSize=segmentationBinSize,

                         dPeaksCutoff=dPeaksCutoff,    penaltyCoefForAddingGrids=penaltyCoefForAddingGrids, minGridHeight=minGridHeight, minPeriodManual=minPeriodManual, maxPeriodManual=maxPeriodManual, forceFirstDigPeakCopyNum=forceFirstDigPeakCopyNum,   # digital peaks
                         grabDataPercentManual= grabDataPercentManual,  origMaxPercentCutoffManual=origMaxPercentCutoffManual,  #  peaksByDensity

                         minReasonableSegmentSize=minReasonableSegmentSize,
                         heterozygosityScoreThreshold=heterozygosityScoreThreshold,  # If segment hetScore is more than this, the segment is heterozygous
                         allowedTumorPercent = allowedTumorPercent,
                         hsNormMat=hsNormMat
  )

  if(!is.null(outputDir)){
    calcPloidyResultOutputFile=file.path(outputDir, paste0(sampleId, '_calculatePloidyResult.Rds'))
    loginfo('saving result to: %s',calcPloidyResultOutputFile)
    saveRDS(calcPloidyResult, file=calcPloidyResultOutputFile)
  }

  loginfo('tumor percentage: %s ',round(calcPloidyResult$percentTumor) )
  segmentPloidy <- sum(calcPloidyResult$segmentData$size * calcPloidyResult$segmentData$cnLevel)/sum(calcPloidyResult$segmentData$size) # weighted by length of segment
  loginfo('approximate ploidy: %s ',round(segmentPloidy,1) )
  logwarn('now go create the contellation plot to confirm this result')

