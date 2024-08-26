\dontrun{

  library(BACDAC)
  # calculateHetScoreExample.R
  # generate heterozygosity score by bin and by arm
  # make a 3-panel plot, top panel the segmentation data and the next two panels are the resulting hetscore by bin and by arm

  basicConfig("DEBUG")
  sampleId='TCGA-14-1402-02A_ds'
  outputDir <- tempdir()

  # inputDir is the path to the load package data
  inputDir <- system.file('extdata', package = "BACDAC")
  segmentationFile <- file.path(inputDir, paste0(sampleId, '_segmentation.csv'))
  segmentation= loadSegmentationFile(segmentationFile) # chr, start, end, rd per segment

  calculateHetScore(
    sampleId,
    inputDir,
    outputDir=outputDir,
    segmentation=segmentation,
    noPdf = TRUE)

}
