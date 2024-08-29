

#' run all the steps in BACDAC, will output pdfs and intermediate files if outputDir is provided
#' @inheritParams commonParameters
#' @example inst/examples/runBACDAC_Example.R
#' @export
runBACDAC=function(sampleId, alternateId,
                   inputDir,
                   outputDir,
                   noPdf=FALSE,
                   readDepthPer30kbBin, readDepthPer100kbBin,segmentation,segmentationBinSize,
                   hsNormMat,testVals,
                   dPeaksCutoff=0.01,penaltyCoefForAddingGrids=0.49, minGridHeight=0.2, minPeriodManual=-1, maxPeriodManual=-1, forceFirstDigPeakCopyNum=-1,   # digital peaks
                   grabDataPercentManual= -1, origMaxPercentCutoffManual=-1,  #  peaksByDensity
                   minReasonableSegmentSize=5.5e6,
                   heterozygosityScoreThreshold=0.98,
                   allowedTumorPercent = 106){

  ### call calculateHetScore ---------
  # generate heterozygosity score by bin and by arm
  # make a 3-panel plot, top panel the segmentation data and the next two panels are the resulting hetscore by bin and by arm
  loginfo('calculate heterozygosity score for %s ', sampleId)

  listOfFiles=calculateHetScore(
    sampleId=sampleId,
    inputDir=inputDir,
    outputDir=outputDir,
    segmentation=segmentation,
    noPdf = noPdf,
    #optional
    readDepthBinnedData=readDepthPer30kbBin,
    readDepthBinSize=readDepthPer30kbBin$windowSize
  )

  # hetScore data - the output from calculateHetScore()
  hetScoreData <- as.data.frame(rtracklayer::import.wig(listOfFiles$hetScorePerBinFile))


  ### call calculatePloidy ---------
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

  calcPloidyResultOutputFile=file.path(outputDir, paste0(sampleId, '_calculatePloidyResult.Rds'))


  loginfo('tumor Percent: %s',round(calcPloidyResult$percentTumor))
  segmentPloidy <- sum(calcPloidyResult$segmentData$size * calcPloidyResult$segmentData$cnLevel)/sum(calcPloidyResult$segmentData$size) # weighted by length of segment
  loginfo('approximate ploidy: %s ',round(segmentPloidy,1) )


  mainPeakNRD=getMainPeakNRD(calcPloidyResult)
  expReadsIn2NPeak_1bp=calcPloidyResult$expReadsIn2NPeak_1bp

  # initialize
  starCloudResult=NULL
  starCloudPlotInputs=NULL


  ### load and make input values for the constellation plot ----
  if(is.null(starCloudPlotInputs)){     # run once only... it takes about 3-5 minutes
    starCloudPlotInputs=loadStarsInTheClouds(sampleId=sampleId, inputDir=inputDir, readDepthPer30kbBin=readDepthPer30kbBin,hetScorePerBinFile=listOfFiles$hetScorePerBinFile,
                                             hsNormMat=hsNormMat, testVals=testVals, wsz=30000, mainPeakNRD=mainPeakNRD, expReadsIn2NPeak_1bp=expReadsIn2NPeak_1bp)
  }

  ### draw constellation plot with linear plot ----
  twoPanelReport(starCloudPlotInputs, calcPloidyResult,readDepthPer30kbBin,segmentation)


}
