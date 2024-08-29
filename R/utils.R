#' load segmentation file
#' check for required columns
#'
#' @param segmentationFile full path to file to load
#'
#' @export
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

#' chromosome segment colors
#' one color for each chromosome 1-22
#' @keywords internal
getSegmentColors=function(){
  ### one color for each chromosome 1-22
  segColors <- c(RColorBrewer::brewer.pal(9, 'Set1')[-c(6,9)] ,                     # remove the yellow and gray
                 'red', 'blue', 'cyan', 'gray45','magenta',                         # to get to 22
                 RColorBrewer::brewer.pal(8, 'Set2')[-c(8)], 'gray75',              # replace the gray with a slightly lighter color
                 'purple', "#00FF92FF"
  )

  if(FALSE){
    segColorsOld <- c(
      'black',
      RColorBrewer::brewer.pal(9, 'Set1')[-c(6)] ,                     # remove the yellow and gray
      palette()[-7],
      '#FF0000FF', '#00FF92FF', '#FFDB00FF',  '#0092FFFF',  '#4900FFFF',  '#FF00DBFF'
    )
    plot(1:22, type='n');    #abline(v=7.5);    abline(v=13.5)
    points(x=1:22, y=rep(5,22), pch=1:22, cex=2,col=segColors)
    points(x=1:22, y=rep(3,22), pch=1:22, cex=2,col=segColorsOld)
  }

  return(segColors)
}
