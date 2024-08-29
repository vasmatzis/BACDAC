if(FALSE){

  # convert from .Rdata to .Rds for better file management. To be loaded in to Zenodo
  myHsNormMat <- bmdTools::loadRdata(file.path(bmdTools::mainDir, 'NextGen/Misc/pipelineInputs/hetScoreAnalysis/lohMat.Rdata')) # aka lohMat
  myTestVals <-  bmdTools::loadRdata(file.path(bmdTools::mainDir, 'NextGen/Misc/pipelineInputs/hetScoreAnalysis/testVals.Rdata'))
  saveRDS(myHsNormMat, file=file.path(bmdTools::mainDir, 'NextGen/Misc/pipelineInputs/hetScoreAnalysis/hetScoreNormMat.Rds'))
  saveRDS(myTestVals,  file=file.path(bmdTools::mainDir, 'NextGen/Misc/pipelineInputs/hetScoreAnalysis/testVals.Rds') )



  folder=66301
  sampleId=bmdSvPipeline::getSampleId(folder)
  postProcessingDir=bmdSvPipeline::getPostProcessingDir(folder)
  outputDir=file.path(bmdTools::mainDir, 'NextGen/johnsonsh/Rprojects/BACDAC/inst/extdata')

  # we will be writing to this path, make sure it exists # TODO: do we need to check that the path is writable?
  if(!dir.exists(file.path(outputDir))){
    dir.create(file.path(outputDir))
  }
  if(!dir.exists(file.path(outputDir))){
    dir.create(file.path(outputDir))
  }

  copyFileToNewLocation=function(oldFileNamePath, newFileNamePath){
    # cpCmd = paste('cp', oldFileNamePath, newFileNamePath)
    if(file.exists(oldFileNamePath)){
      params <- c(oldFileNamePath, newFileNamePath)
      params <- shQuote(params)
      loginfo("Executing '%s %s'", 'cp', paste(params, collapse = " "))
      ret <- system2('cp', params)
      if(ret!=0) {
        stop(sprintf("cp '%s' returned nonzero value %d\nPlease review ", paste(params, collapse = " "), ret))
      }
    }
  }


  ### segmentation data ----
  # convert <sampleId>_cnvIntervals.csv to <sampleId>_segmentation.csv
  #      file.path(bmdTools::mainDir, 'NextGen/Projects/MethodDev/MD66301/GRCh38/svar-1/cnv/TCGA-14-1402-02A_ds_cnvIntervals.csv')
  outputDir = system.file('extdata', package = "BACDAC")
  outputDir=file.path(bmdTools::mainDir, 'NextGen/johnsonsh/Rprojects/BACDAC/inst/extdata')
  newFileNamePath=file.path(outputDir,  paste0(sampleId,'_segmentation.csv'))
  oldFileNamePath=file.path(postProcessingDir, 'cnv', paste0(sampleId,'_cnvIntervals.csv'))

  copyFileToNewLocation(oldFileNamePath, newFileNamePath)

  # hetScore data -----------

  oldLohPerArm <- file.path(postProcessingDir, 'reports', paste0(sampleId, '_lohPerArm.csv'))
  newHetScorePerArm <- file.path(outputDir, paste0(sampleId, '_hetScorePerArm.csv'))

  oldLohPerBin <- file.path(postProcessingDir, 'reports', paste0(sampleId, '_loh.wig.gz'))
  newHetScorePerBin <- file.path(outputDir, paste0(sampleId, '_hetScorePerBin.wig.gz'))

  copyFileToNewLocation(oldFileNamePath=oldLohPerArm, newFileNamePath=newHetScorePerArm)
  copyFileToNewLocation(oldFileNamePath=oldLohPerBin, newFileNamePath=newHetScorePerBin)




  ### refAlt count data ----
  # convert .Rdata to .Rds
  mainChromsNoY=1:23
  for (i in mainChromsNoY) {
    # i=1
    snpFull=bmdTools::loadRdata(file.path(postProcessingDir,     'loh', paste0(sampleId, '_snpVals_',i,'.Rdata')))
    countBPFull=bmdTools::loadRdata(file.path(postProcessingDir, 'loh', paste0(sampleId, '_countBP_',i,'.Rdata')))
    ichrChar=convertChromToCharacter(i)

    iRefAltCount=data.frame('chr'=ichrChar, 'pos'=snpFull, 'ref'=countBPFull$ref, 'alt'=countBPFull$alt)

    iFile=file.path(outputDir, paste0(sampleId,'_','refAltCount_', ichrChar,'.Rds'))
    loginfo('%i writing %s',i, iFile)

    saveRDS(iRefAltCount, file=iFile )
  }

  # first time this did create the data/ directory
  # saves data as .rda which I don't think I want
  my_pkg_data <- sample(1000)
  usethis::use_data(my_pkg_data)




  # without chr column
  # [m071478@mforgers3 data]$ du -hs
  # 108M    .
  # with chr column
  # [m071478@mforgers3 data]$ du -hs
  # 108M    .

  ### frequency array 1Kb, 30Kb, 100Kb ---------------
  folder=66301
  sampleId=bmdSvPipeline::getSampleId(folder)
  postProcessingDir=bmdSvPipeline::getPostProcessingDir(folder)
  outputDir=file.path(bmdTools::mainDir, 'NextGen/johnsonsh/Rprojects/BACDAC/inst/extdata')

  cnvBinnedFile <- file.path(postProcessingDir, 'cnv/cnvBinned.Rdata')
  if(file.exists(cnvBinnedFile)){
    cnvBinnedData <- bmdTools::loadRdata(cnvBinnedFile)
  }else{
    logerror('cant find cnvBinned file: %s', cnvBinnedFile)
  }
  outputWsz <- 1000
  readDepthPer1kbBin <- bmdSvPipeline:::getFreqArrayFromCnvBinned(cnvBinnedData, newWindowSize=outputWsz, maxChrom=24)
  oneKbFile=file.path(outputDir, paste0(sampleId,'_','readDepthPer1kbBin.Rds'))
  loginfo('writing %s',oneKbFile)
  saveRDS(readDepthPer1kbBin, file=oneKbFile )

  wszLinear <- 30000
  readDepthPer30kbBin <- bmdSvPipeline:::getFreqArrayFromCnvBinned(cnvBinnedData, newWindowSize=wszLinear, maxChrom=24)
  thirtyKbFile=file.path(outputDir, paste0(sampleId,'_','readDepthPer30kbBin.Rds'))
  loginfo('writing %s',thirtyKbFile)
  saveRDS(readDepthPer30kbBin, file=thirtyKbFile )

  wszPeaks <- 100000
  readDepthPer100kbBin <- bmdSvPipeline:::getFreqArrayFromCnvBinned(cnvBinnedData, newWindowSize=wszPeaks, maxChrom=22)
  hundredKbFile=file.path(outputDir, paste0(sampleId,'_','readDepthPer100kbBin.Rds'))
  loginfo('writing %s',hundredKbFile)
  saveRDS(readDepthPer100kbBin, file=hundredKbFile )


}

