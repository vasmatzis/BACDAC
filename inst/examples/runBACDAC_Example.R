\dontrun{
library(BACDAC)
library(logging)
basicConfig("DEBUG")
sampleId='TCGA-14-1402-02A_ds'; alternateId=66301
outputDir <- file.path('/research/labs/experpath/vasm/shared/NextGen/johnsonsh/Routput/BACDAC', sampleId) # tempdir()
noPdf=TRUE                          # TRUE= print to screen, FALSE=print to pdf (i e. outputDir/dev/ploidy)

# inputDir is the path to the load package data
inputDir <- system.file('extdata', package = "BACDAC")

# read depth data
thirtyKbFile=file.path(inputDir, paste0(sampleId,'_','readDepthPer30kbBin.Rds'))
readDepthPer30kbBin = readRDS(file=thirtyKbFile )
hundredKbFile=file.path(inputDir, paste0(sampleId,'_','readDepthPer100kbBin.Rds'))
readDepthPer100kbBin = readRDS(file=hundredKbFile )

# segmentation data
segmentationFile <- file.path(inputDir, paste0(sampleId, '_segmentation.csv'))
segmentation= loadSegmentationFile(segmentationFile)

## load two reference files  ---------------
# hsNormMat/lohMat: hetScores from 23 Normals,  101046 x 23 rows, one row for each 30kb segment of the genome, 1-22, X and a part of Y. Columns are values for each of the 23 Normals for each segment
#                   used to look for places in 23 TCGA normals where more than half dropped below the a (i.e. 0.975) cutoff.
# testVals: used to find each possible heterozygosity value for each copy number level (find the right spots for the stars)
hsNormMat <- bmdTools::loadRdata('/research/labs/experpath/vasm/shared/NextGen/Misc/pipelineInputs/hetScoreAnalysis/lohMat.Rdata') # aka lohMat
testVals <-  bmdTools::loadRdata(file.path('/research/labs/experpath/vasm/shared/NextGen/Misc/pipelineInputs/hetScoreAnalysis/testVals.Rdata'))

runBACDAC(sampleId, alternateId,
          outputDir,
          noPdf,
          readDepthPer30kbBin, readDepthPer100kbBin,segmentation,hsNormMat,testVals)
}
