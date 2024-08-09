if(FALSE){
  folder=66301
  sampleId=bmdSvPipeline::getSampleId(folder)
  postProcessingDir=bmdSvPipeline::getPostProcessingDir(folder)
  outputDir='/research/labs/experpath/vasm/shared/NextGen/johnsonsh/Routput/BACDAC'
  outputDir='/research/labs/experpath/vasm/shared/NextGen/johnsonsh/Rprojects/BACDAC/inst/extdata'

  # we will be writing to this path, make sure it exists # TODO: do we need to check that the path is writable?
  if(!dir.exists(file.path(outputDir))){
    dir.create(file.path(outputDir))
  }
  if(!dir.exists(file.path(outputDir, 'reports'))){
    dir.create(file.path(outputDir, 'reports'))
  }

  # hetScorePerArmFile <- file.path(outputDir, 'reports', paste0(sampleId, '_hetScorePerArm.csv'))
  # hetScore30Kb_wigFile <- file.path(outputDir, 'reports', paste0(sampleId, '_hetScore30Kb.wig.gz'))
  # hetScoreWithReadDepthReport <- file.path(outputDir, 'reports', paste0(sampleId, '_HetScoreWithReadDepthReport.pdf'))

  ### segmentation data ----
  # convert <sampleId>_cnvIntervals.csv to <sampleId>_segmentation.csv
  # /research/labs/experpath/vasm/shared/NextGen/Projects/MethodDev/MD66301/GRCh38/svar-1/cnv/TCGA-14-1402-02A_ds_cnvIntervals.csv
  outputDir = system.file('extdata', package = "BACDAC")
  oldFileNamePath=file.path(postProcessingDir, 'cnv', paste0(sampleId,'_cnvIntervals.csv'))
  # newFileNamePath=file.path(outputDir, 'data', paste0(sampleId,'_segmentation.csv'))
  newFileNamePath=file.path(outputDir,  paste0(sampleId,'_segmentation.csv'))
  cpCmd = paste('cp', oldFileNamePath, newFileNamePath)
  if(file.exists(oldFileNamePath)){
    params <- c(oldFileNamePath, newFileNamePath)
    params <- shQuote(params)
    loginfo("Executing '%s %s'", 'cp', paste(params, collapse = " "))
    ret <- system2('cp', params)
    if(ret!=0) {
      stop(sprintf("cp '%s' returned nonzero value %d\nPlease review ", paste(params, collapse = " "), ret))
    }
  }

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




}

