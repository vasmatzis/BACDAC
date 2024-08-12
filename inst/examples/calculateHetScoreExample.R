\dontrun{

  library(BACDAC)
  # calculateHetScoreExample.R
  basicConfig("DEBUG")
  sampleId='TCGA-14-1402-02A_ds'
  inputDir <- system.file('extdata', package = "BACDAC")
  outputDir <- tempdir()

  segmentationFile <- file.path(inputDir, paste0(sampleId, '_segmentation.csv'))
  if(file.exists(segmentationFile)){
    segmentation <- read.csv(segmentationFile,comment.char = '#', header = TRUE)

    # check for required columns
    requiredColumns=c('chr', 'start', 'end','rd')
    missingColumnKey=which(!requiredColumns %in% names(segmentation))
    if(length(missingColumnKey)>0){
      logerror('missing required column: %s',requiredColumns[missingColumnKey])
    }
  }

  calculateHetScore(
    sampleId,
    inputDir,
    outputDir=outputDir,
    segmentation=segmentation,
    noPdf = TRUE)
}
