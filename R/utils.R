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
