#' Determine ploidy state from ploidy and 2N+LOH content
#'
#'  based on SVM line separating the two clusters n=603
#'  y = -0.411x + 1.126  # Johnson/Arnav 2023, n=603
#' @param twoNpLOH 2N+LOH content determined from \code{plotStarsInTheClouds}
#' @param ploidy ploidy value determined from \code{plotStarsInTheClouds}
#' @returns a vector the same length as the input, with the value 'near-diploid' or 'high-ploidy'
#' @examples
#' getPloidyState(twoNpLOH=.2, ploidy=2.5)
#' getPloidyState(twoNpLOH=.03, ploidy=2.6)
#'
#' @export
#'
getPloidyState <- function(twoNpLOH, ploidy){
  # sideOfLine = (1/a)*x + (1/b)*y - 1  # intercept form of line: x/a + y/b = 1 where
  # a = x intercept and b = y intercept, converted to standard form of line: ax + by + c =0
  # where a, b are coefficients.

  ploidyState <- vector(mode='character', length = length(ploidy))
  sideOfLine <-  sign((-2.38)*twoNpLOH + 2.74 -ploidy)

  hiPloidyKey = which(sideOfLine<0)
  diploidKeys = which(sideOfLine>=0)

  ploidyState[hiPloidyKey]= 'high-ploidy'
  ploidyState[diploidKeys]= 'near-diploid'
  # data.frame(twoNpLOH, ploidy, ploidyState); table(ploidyState)

  return(ploidyState)
}

#' add a separation line to a ploidy vs 2N+LOH plot
#'
#' separate near-diploid from high-ploidy
#' @param col line color
#' @param lty line type
#' @param lwd line width
addPloidySvmLineToPlot=function(col='gray', lty='solid',lwd=1.5){
  ## y = -0.411x + 1.126  # Johnson/Arnav 2023, 603
  hInt=1.126 #              when x=0
  pInt=2.740 # 1.126/.411   when y=0
  lines(x=c(pInt, 0), y=c(0, hInt),col=col,lty=lty, lwd=lwd )
}

#' plot ploidy vs 2N+LOH
#'
#' includes ploidy state separation line
#' @param twoNpLOH 2N+LOH content determined from \code{plotStarsInTheClouds}
#' @param ploidy ploidy value determined from \code{plotStarsInTheClouds}
#' @param maxPloidy max x axis value for plot, default=6
#' @examples
#' plotPloidyState(twoNpLOH=.2, ploidy=2.5)
#' plotPloidyState(twoNpLOH=.03, ploidy=2.6)
#'
#' @export
plotPloidyState=function(twoNpLOH, ploidy, maxPloidy=6){
  plot(x=ploidy,y=twoNpLOH,
       xlim=c(1,maxPloidy), ylim=c(0,1),
       xlab='ploidy', ylab='2N+LOH')
  addPloidySvmLineToPlot()
}
