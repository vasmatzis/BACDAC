producing input data for BACDAC
================

This is not an endorsement of the following tools, but a demonstration
of how the input files may be produced. Note, we have not tested these
tools ourselves and cannot attest to the quality of the results.

# Basic Instructions

## ref alt counts using samtools

### Make `[inputFile]`

Assuming the bam file you have is aligned to chromosomes labelled
chr1,chr2,…,chrX,chrY  
`loh_dbSnp_20180418.tsv.gz` is our file with a list of 33629538 common
snp positions. Other such files are available from a variety of sources.

`zcat loh_dbSnp_20180418.tsv.gz | tail -n +2 | awk '{if ($1==23) {test="X"} else if ($1=24) {test="Y"} else {test=$1} print "chr"test"\t"$2}' > [inputFile]`

### Samtools

Install samtools, some sudo commands might be needed

```
wget https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2
tar -xvjf samtools-1.16.1.tar.bz2
cd samtools-1.16.1
./configure --prefix=/usr/local
make  
make install
```

### Running samtools

You can then use this as input to samtools mpileup to create a pileup
file. The early commands are telling the program what to output, base
quality, alignment quality.  
Can limit depth (-d) in this case as we are getting REF/ALT balance and
not depth itself.  
`samtools mpileup -a -A -q 10 -Q 15 -d 1000 --no-output-ins --no-output-ins --no-output-del --no-output-del --no-output-ends -l [inputFile] -o [pileupFile] [bamFile]`

Can also do by chromosome using (-r) where chromNum is `chr1`  
`samtools mpileup -a -A -q 10 -Q 15 -d 1000 --no-output-ins --no-output-ins --no-output-del --no-output-del --no-output-ends -r [chromNum] -l [inputFile] -o [pileupFileChr] [bamFile]`

### Conversion of samtools output

Convert the output into a more easily digestible form:

```
samtools mpileup -a -A -q 10 -Q 15 -d 1000 --no-output-ins --no-output-ins --no-output-del --no-output-del --no-output-ends -l [inputFile] -o [pileupFile] [bamFile]
#this can also do this by chromosome using (-r) like this where chromNum would be like 'chr1'
samtools mpileup -a -A -q 10 -Q 15 -d 1000 --no-output-ins --no-output-ins --no-output-del --no-output-del --no-output-ends -r [chromNum] -l [inputFile] -o [pileupFileChr] [bamFile]
```

`[pileupFinal]` can then be converted to the format needed

```
cat [pileupFile] | awk '{ p=$5; nA=gsub("[Aa]","",p); print $1,$2,$4,$5,nA; }' | awk '{ p=$4; nC=gsub("[Cc]","",p); print $1,$2,$3,$4,$5,nC; }' | awk '{ p=$4; nG=gsub("[Gg]","",p); print $1,$2,$3,$4,$5,$6,nG; }' | awk '{ p=$4; nT=gsub("[Tt]","",p); print $1,$2,$3,$4,$5,$6,$7,nT; }' | awk '{ p=$4; nD=gsub("[*]","",p); print $1,$2,$3,$4,$5,$6,$7,$8,nD; }' | awk 'BEGIN {OFS=",";}{ p=$4; nR=gsub("[<>]","",p); print $3,$5,$6,$7,$8,$9,nR; }' > [pileupTempFile]
echo "CHROM,POS,REF,ALT,COV,A,C,G,T,DEL,REFSKIP" >> [pileupFinal]
paste -d',' <( zcat loh_dbSnp_20180418.tsv.gz | tail -n +2 | awk 'BEGIN {OFS=","};{print $1,$2,$4,$5 }') <(cat [pileupTempFile]) >> [pileupFinal]
```

```
sw=( A:C:G:T )
echo "chr,pos,ref,alt" >> [hetInputFile]
tail -n +2 [pileupFinal] | awk -v sw="${sw[*]}" 'BEGIN {FS=","};BEGIN {OFS="\t"};{ 
n=split(sw, AR, ":");
for(i=1;i<=4;i++) { 
if($3==AR[i]) {
rOut=$( 5 + i ); 
}
if($4==AR[i]) {
aOut=$( 5 + i )
}
}
print $1,$2,rOut,aOut
}' >> [hetInputFile]
```

## binned read count data and segmentation using HHMcopy

### HMMcopy_utils

The necessary input files for HMMcopy come from HMMcopy_utils.

```
git clone https://github.com/shahcompbio/hmmcopy_utils.git
cd hmmcopy_util  
cmake .  
make  
```

There are three scripts in hmmcopy_utils/bin: readCounter, gcCounter,
mapCounter

- readCounter produces the coverage wig from bam
- gcCounter produces the gc wig from reference genome fasta
- mapCounter produces the mappability wig from a mappability reference
  file

### readCounter

set window (-w) to 30000 and limit to primary chromosomes

```
hmmcopy_utils/bin/readCounter -w 30000 -c chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY` `[fileName].bam > [fileName].30kb.wig
```

### gcCounter

set window (-w) to 30000 and limit to primary chromosomes  
Use the hg38 genome.

```
hmmcopy_utils/bin/gcCounter -w 30000 -c chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY` `[referenceGenomeFasta] > gc.hg38.30kb.wig
```

### mapCounter

First download a reference file, the 36bp Umap multiple mapping from
UCSC

```
wget http://hgdownload.soe.ucsc.edu/gbdb/hg38/hoffmanMappability/k36.Umap.MultiTrackMappability.bw
hmmcopy_utils/bin/apCounter -w 30000 -c chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY` `k36.Umap.MultiTrackMappability.bw > map.hg38.k36.wig
```

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

For segmentation data

``` r
default_param <- HMMsegment(covNorm, getparam = TRUE)
longseg_param <- default_param
longseg_param$e <- 0.999999999999999
longseg_param$strength <- 1e30
longseg_segments <- HMMsegment(covNorm, longseg_param, verbose = FALSE)
```

The output here has `longseg_segments$segs` which gives how it has
divided up the genome
