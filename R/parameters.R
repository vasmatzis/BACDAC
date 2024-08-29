#' A dummy function holding a list of all common parameters so
#' the documentation can be used from one place.
#'
#' @param sampleId                   sample Identifier
#' @param alternateId                optional secondary sample identifier
#' @param segmentation               identified regions of the genome with constant read depth. Contains chromosome, start, end, expected CNV, actual CNV and other values we do not need.
#' @param segmentationBinSize        bin size used for the read depth in the segmentation data
#' @param readDepthPer30kbBin        read depth for 30 kb bins, with bin index in linear coordinates
#' @param readDepthPer100kbBin       read depth for 100 kb bins, with bin index in linear coordinates
#' @param hetScoreData               heterozygosity scores determined per 30 kb bin over a 1 Mb region
#' @param numChroms                  number of chromosomes in the reference genome to consider
#' @param dPeaksCutoff               dPeaksCutoff min grid height for a peak to be considered a digital peak
#' @param penaltyCoefForAddingGrids  penalty for adding additional peaks to the digital peak alignment
#' @param minGridHeight              minimum value that can be assigned to the gridHeights
#' @param grabDataPercentManual      portion of main peak data to grab, other peaks will be scaled based on read depth (x location), set to -1 to base off of mainPeak width
#' @param origMaxPercentCutoffManual peaks smaller than this portion of the max peak are not considered; set to -1 to use default value
#' @param pause                      pause execution until user prompts to continue, available interactively only, useful during testing
#' @param noPdf                      if present, do not create pdfs
#' @param skipExtras                 logical to turn on/off plots used for testing and debugging
#' @param minPeriodManual            user provided \code{minPeriod} for digital grid, default is -1 to indicate no user input
#' @param maxPeriodManual            user provided \code{maxPeriod} for digital grid, default is -1 to indicate no user input
#' @param forceFirstDigPeakCopyNum   value to force copy number of first digital peak, use only when ploidy calculation is wrong
#' @param minReasonableSegmentSize   initial smallest segment size to include in ploidy test segments; want to keep as large as possible to avoid 0N segments, but will decrease size if not enough segments are found
#' @param outputDir                  output directory
#' @param heterozygosityScoreThreshold peaks with a hetScore mode above this value are considered heterozygous, typically 0.98, but may vary depending on NGS library quality and preparation
#' @param hsNormMat                  hetScores from  a database of 23 Normals,  101046 x 23 rows, one row for each 30 kb segment of the genome, 1-22, X and a part of Y. Columns are values for each of the 23 Normals for each segment
#' @param segmentData                segments from segmentation that are longer than minReasonableSegmentSize and then broken up into segments no smaller than 3 Mbs
#' @param allelicSegments            segmentData from \code{calculatePloidy} then augmented in \code{plotStarsInTheClouds} with major_copy_number and minor_copy_number
#' @param peakInfo summary table of info for each peak found in \code{peaksByDensity}
#' @param n00                       Precision for creating initial grid

#' @param ideogram             The \link{ideogram} object. If not present, we attempt to use referenceGenomeDescriptor

commonParameters <- function(
 sampleId,
 alternateId,
 segmentation,
 segmentationBinSize,
 readDepthPer30kbBin,
 readDepthPer100kbBin,
 centroArray,
 hetScoreData,
 numChroms,
 dPeaksCutoff,
 penaltyCoefForAddingGrids,
 minGridHeight,
 grabDataPercentManual,
 origMaxPercentCutoffManual,
 pause,
 noPdf,
 skipExtras,
 minPeriodManual,
 maxPeriodManual,
 forceFirstDigPeakCopyNum,
 minReasonableSegmentSize,
 outputDir,
 heterozygosityScoreThreshold,
 hsNormMat
) {

}


