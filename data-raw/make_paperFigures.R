\dontrun{
  # make figures
  # Figure 1 26214  LU-1416305913
  # Figure 2 58184  PT644
  #          58197  PT651
  #          74002  S1-OM-S
  # Figure 3 66301  TCGA-14-1402-02A_ds
  #          19070  TCGA-14-1402-02A

  library(BACDAC)
  library(logging)

  ## load two reference files  ---------------
  # hsNormMat/lohMat: LOH analysis mask, used to look for places in 23 TCGA normals where more than half dropped below the a (i.e. 0.975) cutoff.
  # testVals: used to find each possible heterozygosity value for each copy number level (find the right spots for the stars)
  hsNormMat <- bmdTools::loadRdata('/research/labs/experpath/vasm/shared/NextGen/Misc/pipelineInputs/hetScoreAnalysis/lohMat.Rdata') # aka lohMat
  testVals <-  bmdTools::loadRdata(file.path('/research/labs/experpath/vasm/shared/NextGen/Misc/pipelineInputs/hetScoreAnalysis/testVals.Rdata'))


  ### frequency array 1Kb, 30Kb, 100Kb ---------------
 if(TRUE){
  folder=74002
  sampleId=bmdSvPipeline::getSampleId(folder); alternateId=folder
  postProcessingDir=bmdSvPipeline::getPostProcessingDir(folder)

  # directory for output files:
  outputDir=file.path('/research/labs/experpath/vasm/shared/NextGen/johnsonsh/Routput/BACDAC',sampleId)
  # we will be writing to this path, make sure it exists
  if(!dir.exists(outputDir)){
    dir.create(path = outputDir)
    loginfo('creating output directory: \n\t%s', outputDir)
  }

  cnvBinnedFile <- file.path(postProcessingDir, 'cnv/cnvBinned.Rdata')
  if(file.exists(cnvBinnedFile)){
    cnvBinnedData <- bmdTools::loadRdata(cnvBinnedFile)
  }else{
    logerror('cant find cnvBinned file: %s', cnvBinnedFile)
  }

  # oneKbFile=file.path(outputDir, paste0(sampleId,'_','readDepthPer1kbBin.Rds'))
  # if(!file.exists(oneKbFile)){
  #   outputWsz <- 1000
  #   readDepthPer1kbBin <- bmdSvPipeline:::getFreqArrayFromCnvBinned(cnvBinnedData, newWindowSize=outputWsz, maxChrom=24)
  #   loginfo('writing %s',oneKbFile)
  #   saveRDS(readDepthPer1kbBin, file=oneKbFile )
  # }

  thirtyKbFile=file.path(outputDir, paste0(sampleId,'_','readDepthPer30kbBin.Rds'))
  if(!file.exists(thirtyKbFile)){
    loginfo('writing %s',thirtyKbFile)
    wszLinear <- 30000
    readDepthPer30kbBin <- bmdSvPipeline:::getFreqArrayFromCnvBinned(cnvBinnedData, newWindowSize=wszLinear, maxChrom=24)
    saveRDS(readDepthPer30kbBin, file=thirtyKbFile )
  }else{
    readDepthPer30kbBin=readRDS(readDepthPer30kbBin, file=thirtyKbFile )
  }

  hundredKbFile=file.path(outputDir, paste0(sampleId,'_','readDepthPer100kbBin.Rds'))
  if(!file.exists(hundredKbFile)){
    loginfo('writing %s',hundredKbFile)
    wszPeaks <- 100000
    readDepthPer100kbBin <- bmdSvPipeline:::getFreqArrayFromCnvBinned(cnvBinnedData, newWindowSize=wszPeaks, maxChrom=22)
    saveRDS(readDepthPer100kbBin, file=hundredKbFile )
  }else{
    readDepthPer100kbBin=readRDS(readDepthPer100kbBin, file=hundredKbFile )
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
  saveStarCloudResult=FALSE
  noPdf=TRUE

  # defaults
  segmentationBinSize=30000; omitAnnotations = FALSE;
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

                                   segmentationBinSize=segmentationBinSize,

                                   dPeaksCutoff=dPeaksCutoff,    penaltyCoefForAddingGrids=penaltyCoefForAddingGrids, minGridHeight=minGridHeight, minPeriodManual=minPeriodManual, maxPeriodManual=maxPeriodManual, forceFirstDigPeakCopyNum=forceFirstDigPeakCopyNum,   # digital peaks
                                   grabDataPercentManual= grabDataPercentManual,  origMaxPercentCutoffManual=origMaxPercentCutoffManual,  #  peaksByDensity

                                   minReasonableSegmentSize=minReasonableSegmentSize,
                                   heterozygosityScoreThreshold=heterozygosityScoreThreshold,  # If segment hetScore is more than this, the segment is heterozygous
                                   allowedTumorPercent = allowedTumorPercent,
                                   hsNormMat=hsNormMat
  )

  calcPloidyResultOutputFile=file.path(outputDir, paste0(sampleId, '_calculatePloidyResult.Rds'))
  if(!file.exists(calcPloidyResultOutputFile)){
    loginfo('saving result to: %s',calcPloidyResultOutputFile)
    saveRDS(calcPloidyResult, file=calcPloidyResultOutputFile)
  }


  ### get the needed input values for the plot
  starCloudPlotInputsFile=file.path(outputDir, paste0(sampleId, '_starCloudPlotInputs.Rds'))
  if(is.null(starCloudPlotInputs) & !file.exists(starCloudPlotInputsFile)){     # run once only... it takes about 3-5 minutes
      starCloudPlotInputs=loadStarsInTheClouds(sampleId, inputDir=postProcessingDir, readDepthPer30kbBin,hetScorePerBinWigFile, hsNormMat,
                                               testVals, mainPeakNRD=getMainPeakNRD(calcPloidyResult),
                                               expReadsIn2NPeak_1bp=calcPloidyResult$expReadsIn2NPeak_1bp)

      loginfo('saving result to: %s',starCloudPlotInputsFile)
      saveRDS(starCloudPlotInputs, file=starCloudPlotInputsFile)
  }

  ### make the plot
  if (!noPdf) {
    constellationPdfFile <- file.path(outputDir, paste0(sampleId, 'constellationPlot'))
    pdf(file = constellationPdfFile, paper="a4r", width=8, height=10, title=paste0('constellationPlot',sampleId))
    loginfo('writing pdf: %s', constellationPdfFile)
    on.exit(dev.off(), add=TRUE)
  }

  ##### ONE PLOT OPTION ######
  op <- par(mfrow=c(1,1),mar=c(5,4,3.5,3.5),mgp=c(1.5, 0.5,0))
  starCloudResult=plotStarsInTheClouds(sampleId, alternateId, starCloudPlotInputs, diploidPeakNRD=NULL, tau=min(1,calcPloidyResult$percentTumor/100),
                                       plotEachChrom=FALSE, mainPeakNRD=getMainPeakNRD(calcPloidyResult),
                                       segmentData=calcPloidyResult$segmentData, peakInfo=calcPloidyResult$peakInfo,
                                       digitalPeakZone =calcPloidyResult[['iterationStatsAll']][['digitalPeakZone']],
                                       addAnnotations=TRUE)
  par(op)



  ##### TWO PLOT OPTION ######
  ### plot constellation plot with linear plot ----

  layout( matrix(c(1,2,2),nrow=1),
          heights= c(1),
          widths = c(1.5,2))   # Widths of the two columns
  layout.show(2)
  labelCex=1.5
  # left figure
  op <- par(mar=c(5,3,2,3),mgp=c(1.5, 0.5,0))
  figLabel='E'
  starCloudResult=plotStarsInTheClouds(sampleId=NULL, alternateId=NULL,starCloudPlotInputs, diploidPeakNRD=NULL, tau=min(1,calcPloidyResult$percentTumor/100),
                                       plotEachChrom=FALSE, mainPeakNRD=getMainPeakNRD(calcPloidyResult),
                                       segmentData=calcPloidyResult$segmentData, peakInfo=calcPloidyResult$peakInfo,
                                       digitalPeakZone =calcPloidyResult[['iterationStatsAll']][['digitalPeakZone']],
                                       paperMode=FALSE,plotCex=1.3)
  myAt=starCloudResult$plotAxisLimits$nrdAxisLims[2]
  mtext(figLabel, side=2, at=myAt,cex=labelCex,las=1,line=1.5)
  par(op)

  # right figure
  op <- par(mar=c(5,3,2,1),mgp=c(1.5, 0.5,0))
  figLabel='F'
  yAxisLimits=convertYlimitsToRD(starCloudResult, wsz=readDepthPer30kbBin$windowSize, calcPloidyResult$expReadsIn2NPeak_1bp)
  linearGenomePlot( readDepthBinnedData=readDepthPer30kbBin, wsz=readDepthPer30kbBin$windowSize, segmentation=segmentation,
                    allelicSegments=starCloudResult$allelicSegments,
                    gainColor = 'blue', lossColor= 'red', yAxisLimits = yAxisLimits)
  myAt=yAxisLimits[2]
  mtext(figLabel, side=2, cex=labelCex,at =myAt,las=1,line=1.5)
  par(op)

}

  # left figure
  op <- par(mar=c(5,3,2,3),mgp=c(1.5, 0.5,0))
  blankPlotStarsInTheClouds()
  par(op)

  # right figure
  op <- par(mar=c(5,3,2,1),mgp=c(1.5, 0.5,0))
  figLabel='C'
  plotHetScorePerBin(hetScoreData,
                     ylab="Heterozygosity Score by bin",
                     allelicSegments=starCloudResult$allelicSegments,)
  myAt=1*1.1 # par('usr')[4]
  mtext(figLabel, side=2, cex=labelCex,at =myAt,las=1,line=1.5)
  par(op)


  ### save the results for future inspection ----
  if(saveStarCloudResult){
    starLookUp=makeStarLookUpTable(starCloudResult,percentTumor=calcPloidyResult$percentTumor)

    starCloudData <- list(
      #outputs:
      starLookUp=starLookUp,
      expReadsIn2NPeak_1bp=calcPloidyResult$expReadsIn2NPeak_1bp,
      percentTumor=calcPloidyResult$percentTumor,
      ploidyCN=starCloudResult$ploidyCN,
      diploidPeakNRD=starCloudResult$diploidPeakNRD,
      mainPeakNRD=getMainPeakNRD(calcPloidyResult),
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
