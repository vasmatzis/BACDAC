producing input data
================

This is not an endorsement of the following tools, but a demonstration
of how the input files may be produced. Note, we have not used these
tools ourselves and cannot attest to the quality of the results.

Step-by-step instructions.

### HMMcopy_utils

The necessary input files for HMMcopy come from HMMcopy_utils.

`git clone https://github.com/shahcompbio/hmmcopy_utils.git`  
`cd hmmcopy_util`  
`cmake .`  
`make`

There are three scripts in hmmcopy_utils/bin: readCounter, gcCounter,
mapCounter

- readCounter produces the coverage wig from bam
- gcCounter produces the gc wig from reference genome fasta
- mapCounter produces the mappability wig from a mappability reference
  file

### readCounter

set window (-w) to 30000 and limit to primary chromosomes

`hmmcopy_utils/bin/readCounter -w 30000 -c chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY`
`[fileName].bam > [fileName].30kb.wig`

### gcCounter

set window (-w) to 30000 and limit to primary chromosomes  
Use the hg38 genome.

`hmmcopy_utils/bin/gcCounter -w 30000 -c chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY`
`[referenceGenomeFasta] > gc.hg38.30kb.wig`

### mapCounter

First download a reference file, the 36bp Umap multiple mapping from
UCSC

`wget http://hgdownload.soe.ucsc.edu/gbdb/hg38/hoffmanMappability/k36.Umap.MultiTrackMappability.bw`
`hmmcopy_utils/bin/apCounter -w 30000 -c chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY`
`k36.Umap.MultiTrackMappability.bw > map.hg38.k36.wig`

These commands can also be run for 100Kb The user creates the gc/map
files once. Then creates the read depth files for each bam.

Now use those files to run HMMCopy in R to produce the necessary files
for BACDAC

### HMMCopy

``` r
library(HMMcopy)
rfile <- "[fileName].30kb.wig"
gfile <-"gc.hg38.wig"
mfile <- "map.hg38.k36.wig"
covIn <- wigsToRangedData(rfile, gfile, mfile)
covNorm <- correctReadcount(covIn)
# resort to order
covNorm <- covNorm[order(as.integer(factor(covNorm$chr,levels=paste0("chr",c(1:22,"X","Y"))))),]
# refactor the chromosome column
covNorm$chr <- factor(covNorm$chr,levels=paste0("chr",c(1:22,"X","Y")))
```

The columns of covNorm are: chr, start, end, reads, gc, map, valid,
ideal, cor.gc cor.map, copy

The “valid” column is the mask.

HMMcopy has a built in program for segmentation also, though it tends to
overcall. But primarily you want the segments. Could be better for using
on the larger of the two windows, 100Kb.

### Segmentation

``` r
default_param <- HMMsegment(covNorm, getparam = TRUE)
longseg_param <- default_param
longseg_param$e <- 0.999999999999999
longseg_param$strength <- 1e30
longseg_segments <- HMMsegment(covNorm, longseg_param, verbose = FALSE)
```

The output here has `longseg_segments$segs` which gives how it has
divided up the genome
