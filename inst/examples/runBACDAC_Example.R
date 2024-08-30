\dontrun{
library(BACDAC)
library(logging)
basicConfig("DEBUG")
mySampleId='TCGA-14-1402-02A_ds'; myAlternateId=66301
myOutputDir <- file.path('/research/labs/experpath/vasm/shared/NextGen/johnsonsh/Routput/BACDAC', myAlternateId) # tempdir()
noPdf=TRUE                          # TRUE= print to screen, FALSE=print to pdf (i e. outputDir/dev/ploidy)

# inputDir is the path to the load package data
myInputDir <- system.file('extdata', package = "BACDAC")

# read depth data
thirtyKbFile=file.path(myInputDir, paste0(mySampleId,'_','readDepthPer30kbBin.Rds'))
readDepthPer30kbBin = readRDS(file=thirtyKbFile )
hundredKbFile=file.path(myInputDir, paste0(mySampleId,'_','readDepthPer100kbBin.Rds'))
readDepthPer100kbBin = readRDS(file=hundredKbFile )

# segmentation data
segmentationFile <- file.path(myInputDir, paste0(mySampleId, '_segmentation.csv'))
mySegmentation <- read.csv(segmentationFile, comment.char = '#') # chr, start, end, rd per
mySegmentation <- checkSegmentation(mySegmentation)

mySegmentationBinSize=30000
## load two reference files
# hsNormMat/lohMat: hetScores from 23 Normals,  101046 x 23 rows, one row for each 30kb segment of the genome, 1-22, X and a part of Y. Columns are values for each of the 23 Normals for each segment
#                   used to look for places in 23 TCGA normals where more than half dropped below the a (i.e. 0.975) cutoff.
# testVals: used to find each possible heterozygosity value for each copy number level (find the right spots for the stars)
myHsNormMat <- readRDS(
  '/research/labs/experpath/vasm/shared/NextGen/Misc/pipelineInputs/hetScoreAnalysis/hetScoreNormMat.Rds')
myTestVals <- readRDS(
  '/research/labs/experpath/vasm/shared/NextGen/Misc/pipelineInputs/hetScoreAnalysis/testVals.Rds')

runBACDAC(sampleId=mySampleId,
          alternateId=myAlternateId,
          outputDir=myOutputDir,
          inputDir=myInputDir,
          noPdf,
          readDepthPer30kbBin, readDepthPer100kbBin,
          segmentation=mySegmentation,
          segmentationBinSize=mySegmentationBinSize,
          hsNormMat=myHsNormMat,
          testVals=myTestVals)
}
