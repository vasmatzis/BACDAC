\dontrun{

  library(BACDAC)
  library(logging)
  # calculateHetScoreExample.R
  # generate heterozygosity score by bin and by chromosome arm and save to the output directory

  # make a 3-panel plot, top panel is the segmentation data and the next two panels are the
  # resulting hetscore by bin and by arm

  basicConfig("DEBUG")
  sampleId='TCGA-14-1402-02A_ds'; alternateId=66301
  outputDir <- tempdir()

  # inputDir is the path to the load package data
  inputDir <- system.file('extdata', package = "BACDAC")
  segmentationFile <- file.path(inputDir, paste0(sampleId, '_segmentation.csv'))
  segmentation <- read.csv(segmentationFile, comment.char = '#') # chr, start, end, rd per
  segmentation <- checkSegmentation(segmentation)

  thirtyKbFile=file.path(inputDir, paste0(sampleId,'_','readDepthPer30kbBin.Rds'))
  readDepthPer30kbBin = readRDS(file=thirtyKbFile )

  calculateHetScore(
    sampleId=sampleId,
    inputDir=inputDir,
    outputDir=outputDir,
    segmentation=segmentation,
    noPdf = TRUE,
    #optional - used for the top panel of the three panel plot
    readDepthPer30kbBin=readDepthPer30kbBin,
    readDepthBinSize=readDepthPer30kbBin$windowSize
   )
}
