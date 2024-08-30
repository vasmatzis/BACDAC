#' A dummy function holding a list of all common parameters so
#' the documentation can be used from one place.
#'
#' @param sampleId                   sample Identifier
#' @param alternateId                optional secondary sample identifier
#' @param inputDir                   full path to directory with input files
#' @param outputDir                  full path to directory to output all files and intermediate objects created
#' @param segmentation               identified regions of the genome with constant read depth. Data.frame with required columns:
#'                                   chr, start, end, rd; optional: cnvState (1=loss, 2=normal, 3=gain) for color coded linear linear genome plot
#' @param segmentationBinSize        bin size used for the read depth in the segmentation data
#' @param segmentData                segments from segmentation that are longer than minReasonableSegmentSize and then broken up into segments no smaller than 3 Mbs
#' @param allelicSegments            segmentData from \code{calculatePloidy} then augmented in \code{plotStarsInTheClouds} with major_copy_number and minor_copy_number
#' @param readDepthPer30kbBin       list with two equal length arrays:
#'  * `readDepthArray` read depth for 30 kb bins, normalized for GC content and other artifacts
#'  * `goodWindowArray` linear genome position of each window (bin), masked windows have been removed
#' @param readDepthPer100kbBin       list with two equal length arrays:
#'  * `readDepthArray` read depth for 100 kb bins, normalized for GC content and other sequencing artifacts
#'  * `goodWindowArray` linear genome position of each window (bin) (masked windows have been removed)
#' @param readDepthBinSize           bp size of the bins in the read depth array, expecting 30000 or 100000
#' @param hetScoreData               heterozygosity scores determined per 30 kb bin over a 1 Mb region
#' @param dPeaksCutoff               dPeaksCutoff min grid height for a peak to be considered a digital peak
#' @param penaltyCoefForAddingGrids  penalty for adding additional peaks to the digital peak alignment
#' @param minGridHeight              minimum value that can be assigned to the gridHeights
#' @param grabDataPercentManual      portion of main peak data to grab, other peaks will be scaled based on read depth (x location),
#'                                   default is -1, and value will be base on width of the main (dominate) read depth peak
#' @param origMaxPercentCutoffManual peaks smaller than this portion of the max peak are not considered; set to -1 to use default value
#' @param noPdf                      if present (TRUE), pdf files will not be generated, instead plots are drawn on default device
#' @param minPeriodManual            user provided \code{minPeriod} for digital grid, default is -1 to indicate no user input
#' @param maxPeriodManual            user provided \code{maxPeriod} for digital grid, default is -1 to indicate no user input
#' @param forceFirstDigPeakCopyNum   value to force copy number of first digital peak, use only when ploidy calculation is wrong
#' @param minReasonableSegmentSize   initial smallest segment size to include in ploidy test segments (segmentData); want to
#'                                   keep as large as possible to avoid 0N segments, but will decrease size if not enough segments are found
#' @param heterozygosityScoreThreshold peaks with a hetScore mode above this value are considered heterozygous, typically 0.98,
#'                                   but may vary depending on NGS library quality and preparation
#' @param hsNormMat                  heterozygosity score mask, used to find genomic positions where the heterozgosity score is atypically low.
#'  * rows: hetScores for each 30 kb segment of the genome, 1-22, X and a part of Y.
#'  * columns: one column for each of the 23 normal samples
#' @param testVals                   pre-built array to assist finding each possible heterozygosity value for each copy number level
#' @param peakInfo                   summary table of info for each peak output by \code{calculatePloidy}
#' @param n00                        Precision for creating initial grid
#' @param ideogram                   object loaded with the package, contains genomic reference positions for GRCh38 cytobands
#' @param gainColor                  color for gains in the linear genome plot, default is blue
#' @param lossColor                  color for losses in the linear genome plot, default is red
#' @param allowedTumorPercent allow some tolerance to go over 100 especially for PDXs which will very near 100

commonParameters <- function(
 sampleId,
 alternateId,
 inputDir,
 outputDir,
 segmentation,
 segmentationBinSize,
 readDepthPer30kbBin,
 readDepthPer100kbBin,
 readDepthBinSize,
 hetScoreData,
 dPeaksCutoff,
 penaltyCoefForAddingGrids,
 minGridHeight,
 grabDataPercentManual,
 origMaxPercentCutoffManual,
 noPdf,
 minPeriodManual,
 maxPeriodManual,
 forceFirstDigPeakCopyNum,
 minReasonableSegmentSize,
 heterozygosityScoreThreshold,
 hsNormMat,
 testVals,
 segmentData,
 allelicSegments,
 peakInfo,
 n00,
 ideogram,
 lossColor,
 gainColor,
 allowedTumorPercent
) {

}


