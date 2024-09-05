#' run all the steps in BACDAC, will output pdfs, intermediate and allele specific copy number data if outputDir is provided
#' @inheritParams commonParameters
#' @example inst/examples/runBACDAC_Example.R
#' @export
runBACDAC=function(sampleId, alternateId,
                   inputDir,
                   outputDir,
                   noPdf=FALSE,
                   readDepthPer30kbBin, readDepthPer100kbBin,segmentation,segmentationBinSize,
                   hsNormMat,testVals,
                   gainColor='blue', lossColor= 'red',
                   dPeaksCutoff=0.01,penaltyCoefForAddingGrids=0.49, minGridHeight=0.2, minPeriodManual=-1,
                   maxPeriodManual=-1, forceFirstDigPeakCopyNum=-1,   # digital peaks
                   grabDataPercentManual= -1, origMaxPercentCutoffManual=-1,  #  peaksByDensity
                   minReasonableSegmentSize=5.5e6,
                   heterozygosityScoreThreshold=0.98,
                   allowedTumorPercent = 106){
  # initialize
  starCloudPlotInputs=NULL

  ### call calculateHetScore ---------
  # generate heterozygosity score by bin and by arm
  # make a 3-panel plot, top panel the segmentation data and the next two panels are the resulting hetscore by bin and by arm
  loginfo('calculate heterozygosity score for %s ', sampleId)

  listOfHetScoreFiles=calculateHetScore(
    sampleId=sampleId,
    inputDir=inputDir,
    outputDir=outputDir,
    segmentation=segmentation,
    noPdf = noPdf,
    #optional
    readDepthPer30kbBin=readDepthPer30kbBin,
    readDepthBinSize=readDepthPer30kbBin$windowSize
  )

  # hetScore data - the output from calculateHetScore()
  hetScoreData <- as.data.frame(rtracklayer::import.wig(listOfHetScoreFiles$hetScorePerBinFile))


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

  loginfo('tumor Percent: %s',round(calcPloidyResult$percentTumor))
  segmentPloidy <- sum(calcPloidyResult$segmentData$size * calcPloidyResult$segmentData$cnLevel)/sum(calcPloidyResult$segmentData$size) # weighted by length of segment
  loginfo('approximate ploidy: %s ',round(segmentPloidy,1) )
  loginfo('will now create the contellation plot to confirm this result')

  mainPeakNRD=getMainPeakNRD(calcPloidyResult)
  expReadsIn2NPeak_1bp=calcPloidyResult$expReadsIn2NPeak_1bp


  loginfo("loadStarsInTheClouds")

  ### load and make input values for the constellation plot ----
  if(is.null(starCloudPlotInputs)){     #   takes about 3-5 minutes
    starCloudPlotInputs=loadStarsInTheClouds(sampleId=sampleId, inputDir=inputDir, readDepthPer30kbBin=readDepthPer30kbBin,
                                             hetScorePerBinFile=listOfHetScoreFiles$hetScorePerBinFile,
                                             hsNormMat=hsNormMat, testVals=testVals, readDepthBinSize=30000,
                                             mainPeakNRD=mainPeakNRD, expReadsIn2NPeak_1bp=expReadsIn2NPeak_1bp)
  }

  ### draw constellation plot left of the linear genome plot ----
  loginfo('twoPanelReport')
  starCloudResult= twoPanelReport(starCloudPlotInputs=starCloudPlotInputs, calcPloidyResult=calcPloidyResult,
                                  readDepthPer30kbBin=readDepthPer30kbBin,segmentation=segmentation,
                                  sampleId=sampleId,gainColor=gainColor, lossColor= lossColor)

  loginfo('%s ploidy: %s ',sampleId, round(starCloudResult$ploidyCN,1) )
  loginfo('%s tumor Percent: %s',sampleId, round(calcPloidyResult$percentTumor))

  # write segments to file
  fileToWrite=file.path(outputDir, paste0(sampleId, '_BACDAC_allelic_segments.csv'))
  loginfo('writing allelic segments to file: \n\t%s', fileToWrite)
  write.csv(starCloudResult$allelicSegments,file =fileToWrite,row.names = FALSE )

  return(starCloudResult)

}
