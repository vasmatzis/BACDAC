\dontrun{
  library(BACDAC)
  library(logging)
  # 1. load needed data
  # 2. call peaksByDensity
  # 3. save the returned object to outputDir

  noPdf <- TRUE          # TRUE= print to screen, FALSE=print to pdf in outputDir
  outputDir  <-  NULL    # output folder for pdfs, only needed if noPdf=FALSE

  sampleId <- 'TCGA-14-1402-02A_ds';
  alternateId <- 66301

  ### load data  ---------------
  exampleDataDir <- system.file('extdata', package = "BACDAC")
  inputDir <- exampleDataDir

  ### 1 - load read depth data
  hundredKbFile <- file.path(inputDir, paste0(sampleId,'_','readDepthPer100kbBin.Rds'))
  readDepthPer100kbBin <- readRDS(file=hundredKbFile )

  ### 2 - load segmentation data
  segmentationFile <- file.path(inputDir, paste0(sampleId, '_segmentation.csv'))
  segmentation <- read.csv(segmentationFile, comment.char = '#')
  # check for required columns: # chr, start, end, rd and optionally cnvState
  segmentation <- checkSegmentation(segmentation)
  segmentationBinSize <- 30000;


  ### call peaksByDensity  ---------
  ### the function to
  loginfo('peaksByDensity %s ', sampleId)

  resultPBD <- peaksByDensity(sampleId,readDepthPer100kbBin, segmentation, segmentationBinSize=30000, wszPeaks = 100000, grabDataPercentManual= -1, origMaxPercentCutoffManual=-1,
                           addAreaLinesToPlot=FALSE, omitAnnotations=FALSE,alternateId=NULL)

  ### resultPBD
  # $peakReadDepthList_per1bp
  # [1] 0.02725890 0.01569796 0.02342289 0.03178275 0.03866112
  #
  # $peakHeightList_maxIsOne
  # [1] 1.00000000 0.12749274 0.04343474 0.03383797 0.03103924
  #
  # $scaledGrabDataPercentPerPeak
  # peakRank scaledGrabDataPercent
  # 1        1            0.08300000
  # 2        2            0.10937314
  # 3        3            0.08953895
  # 4        4            0.07686640
  # 5        5            0.06969395
}
