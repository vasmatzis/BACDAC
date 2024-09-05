
  library(BACDAC)
  library(logging)
  # calculatePloidyExample.R

  noPdf <- TRUE          # TRUE= print to screen, FALSE=print to pdf in outputDir
  outputDir  <-  NULL    # output folder for pdfs, only needed if noPdf=FALSE

  sampleId <- 'TCGA-14-1402-02A_ds';
  alternateId <- 66301

  ### load data  ---------------
  ### 1 - load two reference files
  # NOTE/WARNING: if these files do not exist, an attempt will be made to download from Zenodo.
  hsNormMatFile <- "./referenceFiles/hetScoreNormMat.Rds"
  hsNormMat  <-  loadHsNormMat(hsNormMatFile)


  exampleDataDir <- system.file('extdata', package = "BACDAC")
  inputDir <- exampleDataDir

  ### 2 - load read depth data
  thirtyKbFile <- file.path(inputDir, paste0(sampleId,'_','readDepthPer30kbBin.Rds'))
  readDepthPer30kbBin <- readRDS(file=thirtyKbFile )
  hundredKbFile <- file.path(inputDir, paste0(sampleId,'_','readDepthPer100kbBin.Rds'))
  readDepthPer100kbBin <- readRDS(file=hundredKbFile )


  ### 3 - load segmentation data
  segmentationFile <- file.path(inputDir, paste0(sampleId, '_segmentation.csv'))
  segmentation <- read.csv(segmentationFile, comment.char = '#')
  # check for required columns: # chr, start, end, rd and optionally cnvState
  segmentation <- checkSegmentation(segmentation)
  segmentationBinSize <- 30000;


  ### 4 - load heterozygosity score data
  # hetScore data - the output from calculateHetScore()
  # hetScoreDir will typically be the outputDir specified in `calculateHetScore` but here we will
  # load from the BACDAC package example data.
  hetScoreDir <- exampleDataDir
  hetScorePerBinWigFile <- file.path(hetScoreDir, paste0(sampleId, '_hetScorePerBin.wig.gz'))
  hetScoreData <- loadHetScoreFromWig(hetScorePerBinWigFile)

  ### defaults --------
  omitAnnotations <- FALSE;
  dPeaksCutoff <- 0.01; penaltyCoefForAddingGrids <- 0.49; minGridHeight <- 0.2;
  minPeriodManual <- -1;
  maxPeriodManual <- -1; forceFirstDigPeakCopyNum <- -1;
  grabDataPercentManual <- -1; origMaxPercentCutoffManual <- -1;
  minReasonableSegmentSize <- 5.5e6;
  heterozygosityScoreThreshold <- 0.98;
  allowedTumorPercent <- 106


  ### call calculatePloidy, the function to do the initial ploidy configurations ---------
  loginfo('calculate ploidy for %s ', sampleId)
  calcPloidyResult <- calculatePloidy(
    sampleId=sampleId, outputDir=outputDir, noPdf=noPdf, alternateId=alternateId,
    readDepthPer30kbBin=readDepthPer30kbBin, readDepthPer100kbBin=readDepthPer100kbBin,
    segmentation=segmentation, hetScoreData=hetScoreData,
    segmentationBinSize=segmentationBinSize,
    # digital peaks
    dPeaksCutoff=dPeaksCutoff, penaltyCoefForAddingGrids=penaltyCoefForAddingGrids,
    minGridHeight=minGridHeight, minPeriodManual=minPeriodManual, maxPeriodManual=maxPeriodManual,
    forceFirstDigPeakCopyNum=forceFirstDigPeakCopyNum,
    #  peaksByDensity
    grabDataPercentManual=grabDataPercentManual,
    origMaxPercentCutoffManual=origMaxPercentCutoffManual,

    minReasonableSegmentSize=minReasonableSegmentSize,
    heterozygosityScoreThreshold=heterozygosityScoreThreshold,
    allowedTumorPercent=allowedTumorPercent,
    hsNormMat=hsNormMat
  )

  if(!is.null(outputDir)){
    calcPloidyResultOutputFile <- file.path(
      outputDir, paste0(sampleId, '_calculatePloidyResult.Rds'))
    loginfo('saving result to: %s',calcPloidyResultOutputFile)
    saveRDS(calcPloidyResult, file=calcPloidyResultOutputFile)
  }

  loginfo('tumor percentage: %s ',round(calcPloidyResult$percentTumor) )
  segmentPloidy <-  # weighted by length of segment
    sum(calcPloidyResult$segmentData$size *
          calcPloidyResult$segmentData$cnLevel)/sum(calcPloidyResult$segmentData$size)
  loginfo('approximate ploidy: %s ',round(segmentPloidy,1) )
  logwarn('now go create the contellation plot to confirm this result')

