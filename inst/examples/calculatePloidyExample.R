\dontrun{
  library(bmdSvPipeline)
  basicConfig("INFO")

  #### ENTER YOUR INFO HERE ################
  folder=58077                        # 5-digit number of the sample
  noPdf=TRUE                          # TRUE= print to screen, FALSE=print to pdf (i e. outputDir/dev/ploidy)
  outputDir = tempdir();              # output folder for pdfs etc.
  ##########################################

  ### load all the necessary inputs for calculatePloidy---------
  sampleId <- getSampleId(folder=folder)
  postProcessingDir <- getPostProcessingDir(folder=folder)
  inputs = bmdSvPipeline:::loadInputs_calculatePloidy(numfolder = folder,sampleId=sampleId, postProcessingDir = postProcessingDir)

  ### call calculatePloidy, the function to do all the ploidy work ---------
  loginfo('calculate ploidy for %s %s', sampleId, folder)
  result=calculatePloidy(postProcessingDir = postProcessingDir, sampleId=sampleId, outputDir = outputDir,
                         rgdObject = inputs$rgdObject, cnvBinnedData = inputs$cnvBinnedData, cnvIntervals=inputs$cnvIntervals, centroArray = inputs$centroArray, lohdata = inputs$lohdata,
                         noPdf=noPdf, skipExtras=TRUE )

  print(result)

  mainPeakIndex = which(result$peakInfo$rankByHeight==1)
  loginfo('Main peak is %sN',result$peakInfo[mainPeakIndex,'nCopy'])
  loginfo('tumor percentage: %s ',round(result$percentTumor) )



}
