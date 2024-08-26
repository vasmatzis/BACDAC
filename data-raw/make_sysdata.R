## code to prepare `sysdata.R` dataset goes here

# usethis::use_data(make_sysdata.R, overwrite = TRUE)

usethis::use_data_raw('make_sysdata.R')

# make ideogram for sysdata.rda
cytoBandFile       = file.path(bmdTools::mainDir, "Genome/Human/referenceFiles/cytoBand_hg38.txt")
ideogram <- bmdSvPipeline::loadIdeogram(path =cytoBandFile)
head(ideogram)
ideogram$version='GRCh38'
ideogram$path=NULL

# make rgdObject for sysdata.rda
# GRCh38 reference genome descriptor used for translating some chromosome names to numbers

rgdObject=bmdSvPipeline::exampleRgd()
rgdObject$file=NULL
rgdObject$processInformation=NULL
rgdObject$folderId=NULL
rgdObject$svaFileNamePrefix=NULL
rgdObject$commandLine=NULL
rgdObject$referenceGenome$file='/Genome/Human/GRCh38_phiX/GRCh38_full_analysis_set_plus_phiX.fna'
names(rgdObject)

usethis::use_data(rgdObject,ideogram, internal = TRUE, overwrite=TRUE)
