if(FALSE){
  build_readme() # very best way to render README.Rmd

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

  ### transfer files from one location to the package location ----------
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
    }else{
      logwarn('file does not exist: %s',oldFileNamePath)
    }
  }


  ### segmentation data ----
  # convert <sampleId>_cnvIntervals.csv to <sampleId>_segmentation.csv
  #      file.path(bmdTools::mainDir, 'NextGen/Projects/MethodDev/MD66301/GRCh38/svar-1/cnv/TCGA-14-1402-02A_ds_cnvIntervals.csv')
  # outputDir = system.file('extdata', package = "BACDAC")
  outputDir=file.path(bmdTools::mainDir, 'NextGen/johnsonsh/Rprojects/BACDAC/inst/extdata')
  newFileNamePath=file.path(outputDir,  paste0(sampleId,'_segmentation.csv'))
  oldFileNamePath=file.path(postProcessingDir, 'cnv', paste0(sampleId,'_cnvIntervals.csv'))
  copyFileToNewLocation(oldFileNamePath, newFileNamePath)

  # hetScore data -----------
  # files made from bmdSvPipeline have different headers DO NOT USE

  devPath='/research/labs/experpath/vasm/shared/NextGen/johnsonsh/Routput/BACDAC'
  outputDir=file.path(bmdTools::mainDir, 'NextGen/johnsonsh/Rprojects/BACDAC/inst/extdata')

  oldLohPerArm <- file.path(devPath, 'reports', paste0(sampleId, '_hetScorePerArm.csv'))
  newHetScorePerArm <- file.path(outputDir, paste0(sampleId, '_hetScorePerArm.csv'))
  copyFileToNewLocation(oldFileNamePath=oldLohPerArm, newFileNamePath=newHetScorePerArm)

  oldLohPerBin <- file.path(postProcessingDir, 'reports', paste0(sampleId, '_loh.wig.gz'))
  newHetScorePerBin <- file.path(outputDir, paste0(sampleId, '_hetScorePerBin.wig.gz'))
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

numberOfBinsInEachChrom=function(){
  # what is the number of bins for each chromosome
  frq30 =readDepthPer30kbBin$readDepthArray
  bin30 =readDepthPer30kbBin$goodWindowArray

  frq100 =readDepthPer100kbBin$readDepthArray
  bin100 =readDepthPer100kbBin$goodWindowArray

  coords <- getLinearCoordinates(1:24)

  start30=binnedPosStart(coords@chromStart, binSize = 30000)
  start100=binnedPosStart(coords@chromStart, binSize = 100000)

  bin30Sum=array(dim=24)
  for(i in 1:24){
    bin30Sum[i]=   sum( bin30>=start30[i] &
                          bin30<start30[i+1])
  }

  bin100Sum=array(dim=24)
  for(i in 1:24){
    bin100Sum[i]=   sum( bin100>=start100[i] &
                           bin100<start100[i+1])
  }

  sum(bin100Sum)
  length(bin100)

  sum(bin30Sum)
  length(bin30)

  cbind(bin30Sum, bin100Sum)

}
