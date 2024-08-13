\dontrun{

  library(BACDAC)
  # calculateHetScoreExample.R
  basicConfig("DEBUG")
  sampleId='TCGA-14-1402-02A_ds'
  outputDir <- tempdir()

  inputDir <- system.file('extdata', package = "BACDAC")
  segmentationFile <- file.path(inputDir, paste0(sampleId, '_segmentation.csv'))
  segmentation= loadSegmentationFile(segmentationFile)

  calculateHetScore(
    sampleId,
    inputDir,
    outputDir=outputDir,
    segmentation=segmentation,
    noPdf = TRUE)
}
