\dontrun{
library(BACDAC)
library(logging)
basicConfig("DEBUG")
mySampleId='TCGA-14-1402-02A_ds'; myAlternateId=66301
myOutputDir <- tempdir()
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

## load two reference files  ---------------
# NOTE: upload these files from Zenodo and edit to specify your path here.
# https://zenodo.org/records/13619655
hsNormMatFile <- "../referenceFiles/hetScoreNormMat.Rds"
myHsNormMat = readRDS(hsNormMatFile)
testValsFile <-  "../referenceFiles/testVals.Rds"
myTestVals = readRDS(testValsFile)

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
