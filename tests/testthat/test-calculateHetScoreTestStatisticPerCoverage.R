test_that("calculateHetScoreTestStatisticPerCoverage", {
  expect_equal(2 * 2, 4)

  # First check without trimming
  simple <- calculateHetScoreTestStatisticPerCoverage(1000, 0, 0)

  expect_equal(simple[4],
              (1 * 0 +
                 4 * 1 +
                 6 * 2 +
                 4 * 1 +
                 1 * 0) / (1 + 4 + 6 + 4 + 1)
  )

  # Check with trimming
  trimmed <- calculateHetScoreTestStatisticPerCoverage(1000, 2, 1)

  expect_equal(trimmed[1], 0)
  expect_equal(trimmed[2], 0)
  expect_equal(trimmed[3], 1) # Only single bin passes through
  expect_equal(trimmed[4],
              (0 * 0 +
                 0 * 1 +
                 6 * 2 +
                 4 * 1 +
                 0 * 0) / (0 + 0 + 6 + 4 + 0)
  )
} )
