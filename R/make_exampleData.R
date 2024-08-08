if(FALSE){
  folder=66301
  sampleId=bmdSvPipeline::getSampleId(folder)
  postProcessingDir=bmdSvPipeline::getPostProcessingDir(folder)
  outputDir=postProcessingDir
  snpFull=bmdTools::loadRdata(file.path(postProcessingDir,     'loh', paste0(sampleId, '_snpVals_1.Rdata')))
  countBPFull=bmdTools::loadRdata(file.path(postProcessingDir, 'loh', paste0(sampleId, '_countBP_1.Rdata')))

  # we will be writing to this path, make sure it exists # TODO: do we need to check that the path is writable?
  dir.exists(file.path(outputDir, 'reports'))
  hetScorePerArmFile <- file.path(outputDir, 'reports', paste0(sampleId, '_hetScorePerArm.csv'))
  hetScore30Kb_wigFile <- file.path(outputDir, 'reports', paste0(sampleId, '_hetScore30Kb.wig.gz'))

  hetScoreWithReadDepthReport <- file.path(outputDir, 'reports', paste0(sampleId, '_HetScoreWithReadDepthReport.pdf'))


  segmentationFile

 }

