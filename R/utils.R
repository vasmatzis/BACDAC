#' load segmentation file
#' check for required columns
#'
#' @param segmentationFile full path to file to load
#'
loadSegmentationFile=function(segmentationFile){
  if(file.exists(segmentationFile)){
    segmentation <- read.csv(segmentationFile,comment.char = '#', header = TRUE)

    # check for required columns
    requiredColumns=c('chr', 'start', 'end','rd', 'cnvState')
    missingColumnKey=which(!requiredColumns %in% names(segmentation))
    if(length(missingColumnKey)>0){
      logerror('missing required column: %s',requiredColumns[missingColumnKey])
    }
    return(as.data.frame(segmentation) )


  }else{
    logerror('file not found: %s', segmentationFile)
  }
}
