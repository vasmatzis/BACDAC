\dontrun{
  # make figures


  library(BACDAC)
  library(logging)
  writeFileStatus=FALSE
  ### frequency array 1Kb, 30Kb, 100Kb ---------------
  #'TCGA-14-1402-02A_ds'; 19070
  folder=19070
  sampleId=bmdSvPipeline::getSampleId(folder); alternateId=folder
  postProcessingDir=bmdSvPipeline::getPostProcessingDir(folder)
  # directory for output files:
  outputDir='/research/labs/experpath/vasm/shared/NextGen/johnsonsh/Routput/BACDAC'


  cnvBinnedFile <- file.path(postProcessingDir, 'cnv/cnvBinned.Rdata')
  if(file.exists(cnvBinnedFile)){
    cnvBinnedData <- bmdTools::loadRdata(cnvBinnedFile)
  }else{
    logerror('cant find cnvBinned file: %s', cnvBinnedFile)
  }
  outputWsz <- 1000
  readDepthPer1kbBin <- bmdSvPipeline:::getFreqArrayFromCnvBinned(cnvBinnedData, newWindowSize=outputWsz, maxChrom=24)
  wszLinear <- 30000
  readDepthPer30kbBin <- bmdSvPipeline:::getFreqArrayFromCnvBinned(cnvBinnedData, newWindowSize=wszLinear, maxChrom=24)
  wszPeaks <- 100000
  readDepthPer100kbBin <- bmdSvPipeline:::getFreqArrayFromCnvBinned(cnvBinnedData, newWindowSize=wszPeaks, maxChrom=22)

  if(writeFileStatus){
    # oneKbFile=file.path(outputDir, paste0(sampleId,'_','readDepthPer1kbBin.Rds'))
    # loginfo('writing %s',oneKbFile)
    # saveRDS(readDepthPer1kbBin, file=oneKbFile )

    thirtyKbFile=file.path(outputDir, paste0(sampleId,'_','readDepthPer30kbBin.Rds'))
    loginfo('writing %s',thirtyKbFile)
    saveRDS(readDepthPer30kbBin, file=thirtyKbFile )

    hundredKbFile=file.path(outputDir, paste0(sampleId,'_','readDepthPer100kbBin.Rds'))
    loginfo('writing %s',hundredKbFile)
    saveRDS(readDepthPer100kbBin, file=hundredKbFile )
  }

  segmentationFile <- file.path(postProcessingDir,'cnv', paste0(sampleId, '_cnvIntervals.csv'))
  segmentation= loadSegmentationFile(segmentationFile)

  # result from heterozgygosityScore
  hetScorePerBinWigFile <- file.path(postProcessingDir, 'reports', paste0(sampleId, '_loh.wig.gz'))
  hetScoreData <- loadHetScoreFromWig(hetScorePerBinWigFile)

  # intialize
  calcPloidyResult=NULL
  starCloudPlotInputs=NULL
  starCloudResult=NULL;
  saveHetScoreRdata=FALSE
  noPdf=TRUE

  ## load two reference files  ---------------
  # hsNormMat/lohMat: LOH analysis mask, used to look for places in 23 TCGA normals where more than half dropped below the a (i.e. 0.975) cutoff.
  # testVals: used to find each possible heterozygosity value for each copy number level (find the right spots for the stars)
  hsNormMat <- bmdTools::loadRdata('/research/labs/experpath/vasm/shared/NextGen/Misc/pipelineInputs/hetScoreAnalysis/lohMat.Rdata') # aka lohMat
  testVals <-  bmdTools::loadRdata(file.path('/research/labs/experpath/vasm/shared/NextGen/Misc/pipelineInputs/hetScoreAnalysis/testVals.Rdata'))


  # defaults
  segmentationBinSize=30000; numChroms=24;
  pause=FALSE; skipExtras=FALSE; omitAnnotations = FALSE;
  dPeaksCutoff=0.01;    penaltyCoefForAddingGrids=0.49; minGridHeight=0.2; minPeriodManual=-1; maxPeriodManual=-1; forceFirstDigPeakCopyNum=-1;   # digital peaks
  grabDataPercentManual= -1; origMaxPercentCutoffManual=-1;  #  peaksByDensity
  minReasonableSegmentSize=5.5e6;
  heterozygosityScoreThreshold=0.98;  # If segment hetScore is more than this, the segment is heterozygous
  allowedTumorPercent = 106


  # result from calculatePloidy
  loginfo('calculate ploidy for %s ', sampleId)
  calcPloidyResult=calculatePloidy(sampleId=sampleId, outputDir = outputDir, noPdf=noPdf, alternateId=alternateId,
                                   readDepthPer30kbBin = readDepthPer30kbBin, readDepthPer100kbBin= readDepthPer100kbBin,
                                   segmentation=segmentation, hetScoreData = hetScoreData,

                                   segmentationBinSize=segmentationBinSize, numChroms=numChroms,
                                   pause=pause, skipExtras=skipExtras, omitAnnotations = omitAnnotations,

                                   dPeaksCutoff=dPeaksCutoff,    penaltyCoefForAddingGrids=penaltyCoefForAddingGrids, minGridHeight=minGridHeight, minPeriodManual=minPeriodManual, maxPeriodManual=maxPeriodManual, forceFirstDigPeakCopyNum=forceFirstDigPeakCopyNum,   # digital peaks
                                   grabDataPercentManual= grabDataPercentManual,  origMaxPercentCutoffManual=origMaxPercentCutoffManual,  #  peaksByDensity

                                   minReasonableSegmentSize=minReasonableSegmentSize,
                                   heterozygosityScoreThreshold=heterozygosityScoreThreshold,  # If segment hetScore is more than this, the segment is heterozygous
                                   allowedTumorPercent = allowedTumorPercent,
                                   hsNormMat=hsNormMat
  )
  mainPeakNRD=getMainPeakNRD(calcPloidyResult)
  expReadsIn2NPeak_1bp=calcPloidyResult$expReadsIn2NPeak_1bp

  ### get the needed input values for the plot
  if(is.null(starCloudPlotInputs)){     # run once only... it takes about 3-5 minutes
    starCloudPlotInputs=loadStarsInTheClouds(sampleId, inputDir=postProcessingDir, readDepthPer30kbBin,hetScorePerBinWigFile, hsNormMat, testVals, mainPeakNRD=mainPeakNRD, expReadsIn2NPeak_1bp=expReadsIn2NPeak_1bp)
  }

  ### make the plot
  if (!noPdf) {
    constellationPdfFile <- getTypedFile('constellationPlot',dir = outputDir, list(sampleId = sampleId))
    pdf(file = constellationPdfFile@path, paper="a4r", width=8, height=10, title=paste0('hetScVsCN_',sampleId))
    loginfo('writing pdf: %s', constellationPdfFile@path)
    on.exit(dev.off(), add=TRUE)
  }
  op <- par(mfrow=c(1,1),mar=c(5,4,3.5,3.5),mgp=c(1.5, 0.5,0))
  starCloudResult=plotStarsInTheClouds(sampleId, alternateId,starCloudPlotInputs, diploidPeakNRD=NULL, tau=min(1,calcPloidyResult$percentTumor/100),
                                       plotEachChrom=FALSE, mainPeakNRD=mainPeakNRD,
                                       segmentData=calcPloidyResult$segmentData, peakInfo=calcPloidyResult$peakInfo,
                                       digitalPeakZone =calcPloidyResult[['iterationStatsAll']][['digitalPeakZone']],
                                       paperMode=TRUE)
  par(op)



  ##### TWO PLOT OPTION ######
  ### plot constellation plot with linear plot ----





  layout( matrix(c(1,2,2),
                 nrow=1),
          heights= c(1),
          widths = c(1.5,2))   # Widths of the two columns
  layout.show(2)

  # figure 1
  op <- par(mar=c(5,3, 2,2.5),mgp=c(1.5, 0.5,0))
  starCloudResult=plotStarsInTheClouds(sampleId=NULL, alternateId=NULL,starCloudPlotInputs, diploidPeakNRD=NULL, tau=min(1,calcPloidyResult$percentTumor/100),
                                       plotEachChrom=FALSE, mainPeakNRD=mainPeakNRD,
                                       segmentData=calcPloidyResult$segmentData, peakInfo=calcPloidyResult$peakInfo,
                                       digitalPeakZone =calcPloidyResult[['iterationStatsAll']][['digitalPeakZone']],
                                       paperMode=FALSE,plotCex=1.3)
  par(op)

  # figure 2
  op <- par(mar=c(5, 3, 2,1),mgp=c(1.5, 0.5,0))
  yAxisLimits=convertYlimitsToNrd(starCloudResult, wsz=readDepthPer30kbBin$windowSize, expReadsIn2NPeak_1bp)
    linearGenomePlot( readDepthBinnedData=readDepthPer30kbBin, wsz=readDepthPer30kbBin$windowSize, segmentation=segmentation,
                    allelicSegments=starCloudResult$allelicSegments,
                    gainColor = 'blue', lossColor= 'red', yAxisLimits = yAxisLimits)
  par(op)


  ### save the results for future inspection ----
  if(saveHetScoreRdata){
    starLookUp=makeStarLookUpTable(starCloudResult,percentTumor=calcPloidyResult$percentTumor)

    # TODO: maybe this should be moved inside plotStarsInTheClouds?
    starCloudData <- list(
      #outputs:
      starLookUp=starLookUp,
      expReadsIn2NPeak_1bp=expReadsIn2NPeak_1bp,
      percentTumor=calcPloidyResult$percentTumor,
      ploidyCN=starCloudResult$ploidyCN,
      diploidPeakNRD=starCloudResult$diploidPeakNRD,
      mainPeakNRD=mainPeakNRD,
      ploidySegments=calcPloidyResult$segmentData,
      lohContent=starCloudResult$lohContent,
      plotAxisLimits=starCloudResult$plotAxisLimits,
      # inputs:
      forceFirstDigPeakCopyNum=forceFirstDigPeakCopyNum,
      grabDataPercentManual=grabDataPercentManual,
      minPeriodManual=minPeriodManual,
      origMaxPercentCutoffManual=origMaxPercentCutoffManual,
      allowedTumorPercent=allowedTumorPercent
    )

    starCloudResultFile=file.path(outputDir, paste0(sampleId, '_StarCloudResult.Rds'))
    loginfo('saving result to: %s',starCloudResultFile)
    saveRDS(starCloudData, file=starCloudResultFile)
  }


}
