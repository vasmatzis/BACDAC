starInfo=NA;
saveHetScoreRdata=FALSE

  if(is.na(result$percentTumor)){
    logerror('percent tumor is undefined, can not make star-cloud plot')
  }else{
    expReadsIn2NPeak_1bp= result$expReadsIn2NPeak_1bp
    # need mainPeakNRD for the cnvBinned Data, can't assume it is 2; pipeline ploidy output may not be the same as the ploidy output here.
    # cnvBinned, if run in pipeline with other ploidy output may not be the same as the ploidy output here.

    # source('/research/labs/experpath/vasm/shared/NextGen/johnsonsh/Rprojects/bmdSvPipeline_merge/R/hetScores_starsInTheClouds.R')

    cnvBinnedNormalBin=cnvBinnedData$expectedNormalBin  # cnvBinnedNormalBin is the same as expReadsIn2NPeak_1bp, but is the value for the data on the cluster
    mainPeakcnvBinnedNRD=2*(result$peakInfo[which(result$peakInfo$rankByHeight==1),'peakReadDepth_1bp'] / cnvBinnedNormalBin)
    loginfo('mainPeakcnvBinnedNRD = %.3f', mainPeakcnvBinnedNRD)

    ### get the needed input values for the plot
    if(is.null(starCloudPlotInputs)){ # to make sure I don't accidentally rerun this... it takes about 3-5 minutes
      starCloudPlotInputs <- loadStarsInTheClouds(sampleId, postProcessingDir, rgdObject, cnvBinnedData,  wsz=30000, mainPeakcnvBinnedNRD=mainPeakcnvBinnedNRD)
    }

    ### make the plot
    if (!noPdf) {
      constellationPdfFile <- getTypedFile('constellationPlot',dir = outputDir, list(sampleId = sampleId))
      pdf(file = constellationPdfFile@path, paper="a4r", width=8, height=10, title=paste0('hetScVsCN_',sampleId))
      loginfo('writing pdf: %s', constellationPdfFile@path)
      on.exit(dev.off(), add=TRUE)
    }
    op <- par(mfrow=c(1,1),mar=c(5,4,3.5,3.5),mgp=c(1.5, 0.5,0))
    starInfo=plotStarsInTheClouds(sampleId, folderId,starCloudPlotInputs, dipVal=NULL, tau=min(1,result$percentTumor/100),
                                  plotEachChrom=FALSE, mainPeakcnvBinnedNRD=mainPeakcnvBinnedNRD,
                                  segmentData=result$segmentData, peakInfo=result$peakInfo,
                                  bimaVersion=bimaVersion,
                                  forceFirstDigPeakCopyNum=forceFirstDigPeakCopyNum,grabDataPercentManual=grabDataPercentManual,
                                  origMaxPercentCutoffManual=origMaxPercentCutoffManual,minPeriodManual = minPeriodManual,
                                  digitalPeakZone =result[['iterationStatsAll']][['digitalPeakZone']],heterozygosityScoreThreshold=heterozygosityScoreThreshold,
                                  paperMode=FALSE)


    par(op)

    starLookUp=makeStarLookUpTable(starInfo,percentTumor=result$percentTumor)

    # save the data for future use
    # WARNING: this must match exactly what is in cnvDetect2 where it is typically saved
    # TODO: maybe this (and the code in cnvDetct should be moved inside plotStarsInTheClouds?)
    starCloudData <- list(
      #outputs:
      starLookUp=starLookUp,
      expReadsIn2NPeak_1bp=expReadsIn2NPeak_1bp,
      percentTumor=result$percentTumor,
      ploidyCN=starInfo$ploidyCN,
      dipVal=starInfo$dipVal,
      mainPeakcnvBinnedNRD=mainPeakcnvBinnedNRD,
      ploidySegments=result$segmentData,
      lohContent=starInfo$lohContent,
      plotAxisLimits=starInfo$plotAxisLimits,
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
