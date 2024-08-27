
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BACDAC

<!-- badges: start -->
<!-- badges: end -->

BACDAC is an R package for “Ploidy Analysis using the Binomial Allelic
Content with Discretization Algorithm and Constellation plot”

Reports tumor ploidy and purity from whole-genome sequencing data
including low-pass low-tumor NGS. Inputs include read-depth and
segmentation data, and ref/alt counts for common single nucleutide
polymorphisms (SNPs). Calculates the “Heterozygosity Score” and produces
the Constellation plot to visualize allele-specific copy-number as shown
in publication xyz.

## Installation

You can install the development version of BACDAC by cloning the
respository from the Mayo Clinic dev.azure.com (will be moved to github
eventually):

<https://dev.azure.com/mclm/GBS%20GSU/_git/bmd-bacdac>

look at the clone button and copy the SSH url.

RStudio instructions:

1)  Open RStudio: File -\> New Project -\> New repository -\> git

2)  Fill in the pop up window as follows, replace the text in red with
    the name of the package:

- Repository URL: <git@ssh.dev.azure.com>:v3/mclm/GBS%20GSU/bmd-bacdac
- Project directory name: bacdac
- Create project as subdirectory of: <path/to/your>/Rprojects/

3)  Hit ‘create project’ and smile while RStudio does all the git clone
    stuff.

4)  Build

you can do this in Rstudio, by going to the Build tab and then hitting
the “Install” tab you will be directed to install any necessary package
dependancies at this time.

Updating source code

need to pull updates and build regularly during development

# on command line:

cd <path/to/your>/Rprojects/bacdac git pull R CMD INSTALL –build .

# in Rstudio:

go to git tab, click on “git pull” go to Build tab, click on “Install”

``` r
# FILL THIS IN! HOW CAN PEOPLE INSTALL YOUR DEV PACKAGE?
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(BACDAC)
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
