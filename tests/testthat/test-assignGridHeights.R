test_that("assignGridHeights works", {
  peakInfoTest=data.frame(
    peakCol=c("green3",  "blue",    "red",     "cyan" ,   "magenta"),
    scaledGrabDataPercent=c(0.1094,0.0895, 0.0830,0.0769,0.0697),
    peakReadDepth_1bp=c(0.01569796, 0.02342289, 0.02725890,  0.03178275, 0.03866112),
    peakHeight=c(0.12749177, 0.04343453, 1.00000000, 0.03383768, 0.03103891),
    rankByHeight=c( 2, 3, 1, 4, 5),
    peakReadDepth_normX=c( 1.000000, 1.492098, 1.736461, 2.024642, 2.462811),
    rankByX=c(1,2,3,4,5 ),
    bonus=c(1.0199, 0.3475, 8.0000, 0.2707, 0.2483)
  )
  n00=100
  minGridHeight=0.2

  gridHeights <- assignGridHeights(peakInfo, n00, minGridHeight )
  nonZeroGridCoords <- which(gridHeights > 0)

  expect_equal(length(gridHeights),1000)

  expect_equal(
    c(gridHeights[nonZeroGridCoords[1:7]]),
    c(0.764925, 0.968905, 1.009701, 1.019900, 1.009701, 0.968905, 0.764925)
  )

  expect_equal(
    nonZeroGridCoords,
    c(97,   98,  99, 100, 101, 102, 103,           # peak 1
      145, 146, 147, 148, 149, 150, 151, 152, 153, # peak 2
      170, 171, 172, 173, 174, 175, 176, 177, 178, # peak 3
      198, 199, 200, 201, 202, 203, 204, 205, 206, # peak 4
      242, 243, 244, 245, 246, 247, 248, 249, 250) # peak 5
  )

  expect_equal(2 * 2, 4)
})
