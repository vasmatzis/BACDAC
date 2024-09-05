\dontrun{
  library(BACDAC)
  library(logging)

  noPdf <- TRUE          # TRUE= print to screen, FALSE=print to pdf in outputDir
  myOutputDir  <-  tempdir()

  mySampleId <- 'TCGA-14-1402-02A_ds';
  myAlternateId <- 66301

  ### load data  ---------------

  ### 1 - load two reference files
  # NOTE/WARNING: if these files do not exist, an attempt will be made to download from Zenodo.
  hsNormMatFile <- "./referenceFiles/hetScoreNormMat.Rds"
  myHsNormMat <- loadHsNormMat(hsNormMatFile)

  testValsFile <-  "./referenceFiles/testVals.Rds"
  myTestVals <- loadTestVals(testValsFile)


  # inputDir is the path to the loaded package data
  myInputDir <- system.file('extdata', package = "BACDAC")

  ### 2 - read depth data
  thirtyKbFile <- file.path(myInputDir, paste0(mySampleId,'_','readDepthPer30kbBin.Rds'))
  readDepthPer30kbBin <- readRDS(file=thirtyKbFile )
  hundredKbFile <- file.path(myInputDir, paste0(mySampleId,'_','readDepthPer100kbBin.Rds'))
  readDepthPer100kbBin <- readRDS(file=hundredKbFile )

  ### 3 - segmentation data
  segmentationFile <- file.path(myInputDir, paste0(mySampleId, '_segmentation.csv'))
  mySegmentation <- read.csv(segmentationFile, comment.char = '#')
  # check for required columns: # chr, start, end, rd and optionally cnvState
  mySegmentation <- checkSegmentation(mySegmentation)

  mySegmentationBinSize <- 30000


  result <- runBACDAC(sampleId=mySampleId,
                      alternateId=myAlternateId,
                      outputDir=myOutputDir,
                      inputDir=myInputDir,
                      noPdf,
                      readDepthPer30kbBin, readDepthPer100kbBin,
                      segmentation=mySegmentation,
                      segmentationBinSize=mySegmentationBinSize,
                      hsNormMat=myHsNormMat,
                      testVals=myTestVals)

}
