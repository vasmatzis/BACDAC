
#' allele specific cnv
#'
#' 1) find CN for the segment by rounding the cnLevel to an integer
#' 2) find which closest star, within that CN, is closest to the hetScore of the segment
#' 3) the closest star is then the segment's major:minor allele
#' @param starLookUp object returned by \code{makeStarLookUpTable}
#' @inheritParams commonParameters
#' @export
#' @keywords internal
allelicCNV <- function(starLookUp, segmentData){
  # see alleleCNVtesting() for testing different strategies

  dfColumns <-  c('Chr','Start','End','Size', 'Copy_Number', 'Major_Copy_Number', 'Minor_Copy_Number')


  # segmentData=result$segmentData

  # get rid of some of those columns that are annoying me
  segmentDataOut <- segmentData[,c('chr','start','end','size','rd', 'nrd','lohScoreMedian','valid','cnLevel')]

  # round the copy number level to get an integer value
  segmentDataOut[['copy_number']] <- round(segmentDataOut$cnLevel)     # integer


  # for each segment, based on the CN and hetScore, list the major and minor allele as integers
  for(i in 1:nrow(segmentDataOut)){
    # i <- 56; segmentDataOut[i,]

    if(!is.na(segmentDataOut[i, 'valid']) &&
       segmentDataOut[i, 'valid'] == 1){ # hetScore mean and median are both not zero
      iHetScore <- segmentDataOut[i, 'lohScoreMedian']

      # find which hetScore in the starLookUp table is the closest to the segment hetScore
      iCN <- segmentDataOut[i,'copy_number']
      key <- which(starLookUp[,'cn']==iCN)
      keyInd <- which.min(abs(iHetScore-starLookUp[key,'hetScore']))
      #assign the major and minor values from the starLookUp table to the segment
      if(length(keyInd)==1){
        segmentDataOut[i,'major_copy_number'] <- starLookUp[key[keyInd],'major']
        segmentDataOut[i,'minor_copy_number'] <- starLookUp[key[keyInd],'minor']
      }
    }
  }
  # clean up columns
  # names(df)[names(df) == 'old.var.name'] <- 'new.var.name'
  names(segmentDataOut)[names(segmentDataOut) == "lohScoreMedian"] <- "heterozygosity_score"
  names(segmentDataOut)[names(segmentDataOut) == "cnLevel"] <- "copy_clonality"
  names(segmentDataOut)[names(segmentDataOut) == "nrd"] <- "normalized_read_depth"
  names(segmentDataOut)[names(segmentDataOut) == "rd"] <- "read_depth"


  zero=which(segmentDataOut$minor_copy_number==0)
  segmentDataOut[zero,]



  return(segmentData=segmentDataOut)
}

#' get LOH content metrics used to identify high-ploidy
#'
#' @param allelicSegData segment data returned from \code{calculatePloidy} and augmented with major and minor allele specific copy number
#' @return a list with the three metrics used by us (A) and other authors (B and C)
#' lohContentA_maj2_min0
#' lohContentB_maj1_min0
#' lohContentC_maj2
#' @export
#' @keywords internal
getLohContent <- function(allelicSegData){
  numTotalSegments <- nrow(allelicSegData)
  hemiDel_keys  <-  which(allelicSegData$minor==0 & allelicSegData$copy_number == 1)
  cnLOH_keys <- which(allelicSegData$minor==0 & allelicSegData$copy_number == 2)
  gainLOH_keys <- which(allelicSegData$minor==0 & (allelicSegData$copy_number ==3 | allelicSegData$copy_number ==4))
  ampLOH_keys <- which(allelicSegData$minor==0 & allelicSegData$copy_number >=5)
  allLOH_keys <- c(cnLOH_keys, gainLOH_keys,ampLOH_keys)

  sum(allelicSegData[hemiDel_keys,'size'] )/sum(allelicSegData[,'size']) # weighted by length of segment
  sum(allelicSegData[cnLOH_keys,  'size'] )/sum(allelicSegData[,'size']) # weighted by length of segment
  sum(allelicSegData[gainLOH_keys,'size'] )/sum(allelicSegData[,'size']) # weighted by length of segment
  sum(allelicSegData[ampLOH_keys, 'size'] )/sum(allelicSegData[,'size']) # weighted by length of segment
  sum(allelicSegData[allLOH_keys, 'size'] )/sum(allelicSegData[,'size']) # weighted by length of segment

  sum(allelicSegData[cnLOH_keys,'size'] )/sum(allelicSegData[,'size'])
  sum(allelicSegData[cnLOH_keys,'size'] * allelicSegData[cnLOH_keys,'copy_number']  )/sum(allelicSegData[,'size']) # also weighted by copy_number

  # cnLohRatioA = lohContentA_maj2_min0
  # cnLohRatioB = lohContentB_maj1_min0
  # mcnGtEq2Ratio = lohContentC_maj2

  #svatools methodA: minor allele = 0 and major  >=2 aka 2N+LOH
  minorEq0_cnGtEq2_keys <- which(allelicSegData$minor==0 & allelicSegData$copy_number >= 2)
  lohContentA_maj2_min0_NW        <- length(minorEq0_cnGtEq2_keys)/ numTotalSegments  # not weighted, WRONG
  lohContentA_maj2_min0_pW        <- sum(allelicSegData[minorEq0_cnGtEq2_keys,'size']* allelicSegData[minorEq0_cnGtEq2_keys,'copy_number'] )/sum(allelicSegData[,'size']) # weighted by length of segment and copy number
  lohContentA_maj2_min0           <- sum(allelicSegData[minorEq0_cnGtEq2_keys,'size'] )/sum(allelicSegData[,'size'])                                                   # weighted by length of segment

  #svatools methodB: minor allele = 0 (will include CN=1) aka all LOH including haploid LOH
  minorEq0_keys  <- which(allelicSegData$minor==0 )
  lohContentB_maj1_min0_NW <- length(minorEq0_keys)/numTotalSegments  # not weighted
  lohContentB_maj1_min0    <- sum(allelicSegData[minorEq0_keys,'size'])/sum(allelicSegData[,'size']) # weighted by length of segment


  #bielski method: major allele, copy number (MCN) >=2
  majorGtEq2_keys <- which(allelicSegData$major>=2)
  lohContentC_maj2_NW  <- length(majorGtEq2_keys)/numTotalSegments  # not weighted
  lohContentC_maj2    <- sum(allelicSegData[majorGtEq2_keys,'size'])/sum(allelicSegData[,'size']) # weighted by length of segment
  return(list(lohContentA_maj2_min0=lohContentA_maj2_min0,
              lohContentB_maj1_min0=lohContentB_maj1_min0,
              lohContentC_maj2=lohContentC_maj2)
  )
}





#' create the starLookUp table
#'
#' lists the hetScore for each possible allele fraction
#'
#' @param starCloudResult object returned by \code{plotStarsInTheClouds}
#' @param percentTumor tumor percent listed in the object returned by \code{calculatePloidy}
#' @export
#' @keywords internal
#' @returns
#' A data.frame with 5 columns: hetScore, nrd, cn, major, minor.
#' One row copy number 1 to max copy
makeStarLookUpTable <- function(starCloudResult,percentTumor){
  # percentTumor= calcPloidyResult$percentTumor

  # create the starLookUp table
  hetScore <- starCloudResult$starVals / starCloudResult$medStarVals
  nrd <- starCloudResult$plotStarRange
  cn <- round(calcCopyNumber(NRD=nrd, tau=percentTumor/100))
  starLookUp <- data.frame(hetScore,
                           nrd=round(nrd,3),
                           cn,
                           major=0,    # initialize to 0
                           minor=0 )   # initialize to 0

  # assign major and minor allele columns for each CN and hetScore
  # the LOH level will be the first row, the heterozygous level the last row, within the rows for a given copy number level.
  for(icn in 1:max(cn)){
    # icn=6
    minor <- seq(from=0, to=icn/2)
    major <-  icn-minor
    starLookUp[which(starLookUp[,'cn']==icn),'major'] <- major
    starLookUp[which(starLookUp[,'cn']==icn),'minor'] <- minor
  }
  return(starLookUp)
}

# convert the constellation y axis limits which are in NRD units and returned from the function into read depth values
convertYlimitsToRD=function(starCloudResult, wsz, expReadsIn2NPeak_1bp){
  minYplotNRD = starCloudResult$plotAxisLimits$nrdAxisLims[1]
  maxYplotNRD = starCloudResult$plotAxisLimits$nrdAxisLims[2]
  minYplotRD=  calcRD(nrd=minYplotNRD, wsz=readDepthPer30kbBin$windowSize, expReadsIn2NPeak_1bp)
  maxYplotRD=  calcRD(nrd=maxYplotNRD, wsz=readDepthPer30kbBin$windowSize, expReadsIn2NPeak_1bp)
  yAxisLimits=c(minYplotRD, maxYplotRD)
  return(yAxisLimits)
}
