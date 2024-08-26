\dontrun{
  library(BACDAC)
  # calculateHetScoreExample.R
  basicConfig("DEBUG")
  noPdf=TRUE                          # TRUE= print to screen, FALSE=print to pdf (i e. outputDir/dev/ploidy)
  # outputDir = tempdir();              # output folder for pdfs etc.
  outputDir='/research/labs/experpath/vasm/shared/NextGen/johnsonsh/Routput/BACDAC'
  sampleId='TCGA-14-1402-02A_ds'; alternateId=66301

  ### load data ###
  inputDir <- system.file('extdata', package = "BACDAC") # or '/research/labs/experpath/vasm/shared/NextGen/johnsonsh/Rprojects/BACDAC/inst/extdata'

  # hetScores from 23 Normals,  101046 x 23 rows, one row for each 30kb segment of the genome, 1-22, X and a part of Y. Columns are values for each of the 23 Normals for each segment
  # TODO: can I also make this example data? or should it be downloaded separately? file size= 18.6M
  hsNormMat <- bmdTools::loadRdata('/research/labs/experpath/vasm/shared/NextGen/Misc/pipelineInputs/hetScoreAnalysis/lohMat.Rdata') # aka lohMat

  # segmentation data
  segmentationFile <- file.path(inputDir, paste0(sampleId, '_segmentation.csv'))
  segmentation= loadSegmentationFile(segmentationFile)

  # hetScore data - the output from calculateHetScore()
  hetScoreDir='/research/labs/experpath/vasm/shared/NextGen/johnsonsh/Routput/BACDAC/reports'
  hetScorePerBinWigFile <- file.path(hetScoreDir, paste0(sampleId, '_hetScorePerBin.wig.gz'))
  hetScoreData <- as.data.frame(rtracklayer::import.wig(hetScorePerBinWigFile))

  # read depth data
  thirtyKbFile=file.path(inputDir, paste0(sampleId,'_','readDepthPer30kbBin.Rds'))
  readDepthPer30kbBin = readRDS(file=thirtyKbFile )
  hundredKbFile=file.path(inputDir, paste0(sampleId,'_','readDepthPer100kbBin.Rds'))
  readDepthPer100kbBin = readRDS(file=hundredKbFile )


  # defaults
  segmentationBinSize=30000; numChroms=24;
  pause=FALSE; skipExtras=FALSE; omitAnnotations = FALSE;
  dPeaksCutoff=0.01;    penaltyCoefForAddingGrids=0.49; minGridHeight=0.2; minPeriodManual=-1; maxPeriodManual=-1; forceFirstDigPeakCopyNum=-1;   # digital peaks
  grabDataPercentManual= -1; origMaxPercentCutoffManual=-1;  #  peaksByDensity
  minReasonableSegmentSize=5.5e6;
  heterozygosityScoreThreshold=0.98;  # If segment hetScore is more than this, the segment is heterozygous
  allowedTumorPercent = 106


  ### call calculatePloidy, the function to do all the ploidy work ---------
  loginfo('calculate ploidy for %s ', sampleId)
  result=calculatePloidy(sampleId=sampleId, outputDir = outputDir, noPdf=noPdf, alternateId=alternateId,
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

  calcPloidyResultOutputFile=file.path(outputDir, paste0(sampleId, '_calculatePloidyResult.Rds'))
  loginfo('saving result to: %s',calcPloidyResultOutputFile)
  saveRDS(result, file=calcPloidyResultOutputFile)

  print(result)

  mainPeakIndex = which(result$peakInfo$rankByHeight==1)
  loginfo('Main peak is %sN',result$peakInfo[mainPeakIndex,'nCopy'])
  loginfo('tumor percentage: %s ',round(result$percentTumor) )
  loginfo('approximate ploidy: %s ',round( mean(result$segmentData$cnLevel),1) ) # based on segment copy number and not adjusted for size


}
