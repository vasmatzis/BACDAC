\dontrun{

  library(BACDAC)

  # intialize
  starCloudResult=NULL;
  saveHetScoreRdata=FALSE
  starCloudPlotInputs=NULL
  noPdf=TRUE

  ## load two reference files  ---------------
  # hsMat/lohMat: LOH analysis mask, used to look for places in the 23 TCGA normals where more than half dropped below the a (i.e. 0.975) cutoff.
  # testVals: used to find each possible heterozygosity value for each copy number level (find the right spots for the stars)
  hsMat   <- bmdTools::loadRdata(file.path('/research/labs/experpath/vasm/shared/NextGen/Misc/pipelineInputs/hetScoreAnalysis/lohMat.Rdata'))
  testVals <-  bmdTools::loadRdata(file.path('/research/labs/experpath/vasm/shared/NextGen/Misc/pipelineInputs/hetScoreAnalysis/testVals.Rdata'))

  sampleId='TCGA-14-1402-02A_ds'; alternateId=66301

  # directory with input files:
  inputDir <- system.file('extdata', package = "BACDAC") # or '/research/labs/experpath/vasm/shared/NextGen/johnsonsh/Rprojects/BACDAC/inst/extdata'
  readDepthPer30kbBin=  readRDS(file.path(inputDir, paste0(sampleId,'_','readDepthPer30kbBin.Rds')))

  # directory for output files:
  outputDir='/research/labs/experpath/vasm/shared/NextGen/johnsonsh/Routput/BACDAC'

  # result from heterozgygosityScore
  hetScorePerBinWigFile <- file.path(outputDir, 'reports', paste0(sampleId, '_hetScorePerBin.wig.gz'))

  # result from calculatePloidy
  calcPloidyResultOutputFile=file.path(outputDir, paste0(sampleId, '_calculatePloidyResult.Rds'))
  calcPloidyResult = readRDS(calcPloidyResultOutputFile)
  mainPeakNRD=getMainPeakNRD(calcPloidyResult)
  expReadsIn2NPeak_1bp=calcPloidyResult$expReadsIn2NPeak_1bp

  ### get the needed input values for the plot
  if(is.null(starCloudPlotInputs)){ # to make sure I don't accidentally rerun this... it takes about 3-5 minutes
    starCloudPlotInputs=loadStarsInTheClouds(sampleId, inputDir, readDepthPer30kbBin,hetScorePerBinWigFile, hsMat, testVals, mainPeakNRD=mainPeakNRD, expReadsIn2NPeak_1bp=expReadsIn2NPeak_1bp)
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


    ### plot constellation plot with linear plot ----
    # load read depth data
    thirtyKbFile=file.path(inputDir, paste0(sampleId,'_','readDepthPer30kbBin.Rds'))
    readDepthPer30kbBin = readRDS(file=thirtyKbFile )
    # load segmentation data
    segmentationFile <- file.path(inputDir, paste0(sampleId, '_segmentation.csv'))
    segmentation= loadSegmentationFile(segmentationFile)
    hetScorePerBinWigFile <- file.path(outputDir, 'reports', paste0(sampleId, '_hetScorePerBin.wig.gz'))
    hetScore <- loadHetScoreFromWig(hetScorePerBinWigFile)


    ##### THREE PLOT OPTION ######
    layout( matrix(c(1,1,2,3),
                   byrow=FALSE,
                   nrow=2),
                   heights= c(1,1),
                   widths = c(1.5,2))   # Widths of the two columns
    layout.show(3)

    # figure 1
    op <- par(mar=c(5,3, 2,2.5),mgp=c(1.5, 0.5,0))
    starCloudResult=plotStarsInTheClouds(sampleId=NULL, alternateId=NULL,starCloudPlotInputs, diploidPeakNRD=NULL, tau=min(1,calcPloidyResult$percentTumor/100),
                                  plotEachChrom=FALSE, mainPeakNRD=mainPeakNRD,
                                  segmentData=calcPloidyResult$segmentData, peakInfo=calcPloidyResult$peakInfo,
                                  digitalPeakZone =calcPloidyResult[['iterationStatsAll']][['digitalPeakZone']],
                                  paperMode=FALSE)
    par(op)

    # figure 2
    op <- par(mar=c(1,3,2,1),mgp=c(1.5, 0.5,0))
    plotHetScorePerBin(hetScore,  allelicSegments=starCloudResult$allelicSegments, sampleId=sampleId)

    # figure 3
    op <- par(mar=c(2, 3, 1,1),mgp=c(1.5, 0.5,0))
    linearGenomePlot( readDepthBinnedData=readDepthPer30kbBin, wsz=readDepthPer30kbBin$windowSize, segmentation=segmentation,
                      allelicSegments=starCloudResult$allelicSegments,
                      gainColor = 'blue', lossColor= 'red')
    par(op)


    ##### TWO PLOT OPTION ######
    layout( matrix(c(1,1,1,0,2,0),
                   byrow=FALSE,
                   nrow=3),
            heights= c(.5,1,.5),
            widths = c(1.5,2))   # Widths of the two columns
    layout.show(2)

    # figure 1
    op <- par(mar=c(5,3, 2,2.5),mgp=c(1.5, 0.5,0))
    starCloudResult=plotStarsInTheClouds(sampleId=NULL, alternateId=NULL,starCloudPlotInputs, diploidPeakNRD=NULL, tau=min(1,calcPloidyResult$percentTumor/100),
                                  plotEachChrom=FALSE, mainPeakNRD=mainPeakNRD,
                                  segmentData=calcPloidyResult$segmentData, peakInfo=calcPloidyResult$peakInfo,
                                  digitalPeakZone =calcPloidyResult[['iterationStatsAll']][['digitalPeakZone']],
                                  paperMode=FALSE)
    par(op)

    # figure 2
    op <- par(mar=c(2, 3, 1,1),mgp=c(1.5, 0.5,0))
    linearGenomePlot( readDepthBinnedData=readDepthPer30kbBin, wsz=readDepthPer30kbBin$windowSize, segmentation=segmentation,
                      allelicSegments=starCloudResult$allelicSegments,
                      gainColor = 'blue', lossColor= 'red')
    par(op)


    ### save some data ----

    starLookUp=makeStarLookUpTable(starCloudResult,percentTumor=calcPloidyResult$percentTumor)

    # save the data for future use
    # WARNING: this must match exactly what is in cnvDetect2 where it is typically saved
    # TODO: maybe this (and the code in cnvDetct should be moved inside plotStarsInTheClouds?)
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
      allowedTumorPercent=allowedTumorPercent,
      continueOnPloidyFail=continueOnPloidyFail
    )

    if(saveHetScoreRdata){
      hetScores_StarCloudDataInfo <- getTypedFile("hetScores_StarCloudData", dir = outputDir, values = list(sampleId = sampleId))
      saveRdata(file=hetScores_StarCloudDataInfo,
                data=starCloudData,
                fileLabel = "theoretical het scores per major and minor alleles", deleteIfExists = TRUE)
    }
  }

}
