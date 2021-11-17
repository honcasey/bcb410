library(testthat)
library(cellLineConsistency)

# Generate example data - takes a while so putting out here
intersected <- PharmacoGx::intersectPSet(c(CTRPv2, GRAY),
                                         intersectOn = c("drugs", "cell.lines"))

test_that("Check if right sensitivity measure intersection", {
  pSet = intersected
  commonSensMeasures <- intersectSensMeasures(pSet = pSet)

  commonMeas <- c("aac_recomputed", "ic50_recomputed", "HS", "E_inf", "EC50")
  expect_equal(commonSensMeasures, commonMeas)
  expect_length(commonSensMeasures, 5)
  expect_type(commonSensMeasures, "character")
})

test_that("Checking for invalid user inputs", {
  pSet = intersected

  # only one PharmacoSet (nothing to intersect)
  expect_error(intersectSensMeasures(pSet = c(CTRPv2)),
               "pSet list should have at least two PharmacoSets to
           compare correlations.")

  expect_error(intersectSensMeasures(pSet = CTRPv2),
               "pSet should be of class list, which contains a list of
         intersected PharmacoSets.")
})
