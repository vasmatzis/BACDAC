#' check segmentation data
#'
#' helper function to check for required columns: 'chr', 'start', 'end', 'rd' and recommend column: 'cnvState'
#' * chr, start, end provide the location of the segment
#' * rd is the read depth of the segment
#' * cnvState an integer 1=loss, 2=normal, 3=gain
#'
#' @param segmentation genome segmented into predicted intervals of equal read depth
#' @returns data.frame with only those required and recommended columns, others are removed
#' @export
#'
checkSegmentation=function(segmentation){
    # check for required columns
    requiredColumns=c('chr', 'start', 'end','rd')
    missingColumnKey=which(!requiredColumns %in% names(segmentation))
    if(length(missingColumnKey)>0){
      logerror('missing required column: %s',requiredColumns[missingColumnKey])
    }

    recommendedColumns=c('cnvState')
    missingColumnKey=which(!recommendedColumns %in% names(segmentation))
    if(length(missingColumnKey)>0){
      logerror('missing recommended column: %s',recommendedColumns[missingColumnKey])
    }

    returnColumns=which(names(segmentation) %in% c(requiredColumns,recommendedColumns))

    return(as.data.frame(segmentation[,returnColumns]) )
}



# https://zenodo.org/records/13619655

#' load hsNormMat
#'
#' if file does not exist, will download from zenodo. Will check contents and return matrix if
#' successfully loaded.
#'
#' @param hsNormMatFile full path to hsNormMatFile
#' @export
loadHsNormMat=function(hsNormMatFile){
  if(!file.exists(hsNormMatFile)){
    if(!dir.exists(dirname(hsNormMatFile))){
      dir.create(dirname(hsNormMatFile))
    }
    download.file(url='https://zenodo.org/records/13619655/files/hetScoreNormMat.Rds?download=1',
                  destfile = hsNormMatFile)
  }

  # check again, should be there now.
  if(!file.exists(hsNormMatFile)){
    stop(sprintf('required file does not exist: %s:',hsNormMatFile ))
  }else{
    hsNormMat = readRDS(hsNormMatFile)
    if(nrow(hsNormMat)==101046){
      loginfo('%s loaded',basename(hsNormMatFile))
    }else{
      logwarn('expecting matrix with 101046 rows but found %i, file not loaded properly',
              nrow(hsNormMat))
    }
  }
  return(hsNormMat)
}

#' load testVals
#'
#' if file does not exist, will download from zenodo. Will check contents and return matrix if
#' successfully loaded.
#'
#' @param testValsFile full path to testValsFile
#' @export
loadTestVals=function(testValsFile){
  if(!file.exists(testValsFile)){
    if(!dir.exists(dirname(testValsFile))){
      dir.create(dirname(testValsFile))
    }
    download.file(url='https://zenodo.org/records/13619655/files/testVals.Rds?download=1',
                  destfile = testValsFile)
  }

  # check again, should be there now.
  if(!file.exists(testValsFile)){
    stop(sprintf('required file does not exist: %s:',testValsFile ))
  }else{
    testVals = readRDS(testValsFile)
    if(nrow(testVals)==997){
      loginfo('%s loaded',basename(testValsFile))
    }else{
      logwarn('expecting matrix with 997 rows but found %i, file not loaded properly',
              nrow(testVals))
    }
  }
  return(testVals)
}
