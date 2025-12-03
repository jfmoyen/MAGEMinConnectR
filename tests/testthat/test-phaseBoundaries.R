test_that("Generic lower phase boundary is found at the right place", {
  expect_equal(
    LBound_ATAC4( phasename="opx",
                  Tlow = 650,
                  Thigh = 1250,
                  maxIter = 32,
                  AmountTolerance = 0.5 / 100,
                  verbose = F,
                  printPhaseBoundaryResults = T,
                  ### arguments passed to MAGEMin() ###
                  Xoxides = Xoxides,
                  X = rock,
                  Pkbar = 4.5,
                  showMinimizationResults=F)$T_C,
    837.5,
    tolerance=1e-2
  )
})

test_that("Generic upper phase boundary is found at the right place", {
  expect_equal(
    UBound_ATAC4( phasename="pl",
                  Tlow = 650,
                  Thigh = 1250,
                  maxIter = 32,
                  AmountTolerance = 0.5 / 100,
                  verbose = F,
                  printPhaseBoundaryResults = T,
                  ### arguments passed to MAGEMin() ###
                  Xoxides = Xoxides,
                  X = rock,
                  Pkbar = 4.5,
                  showMinimizationResults=F)$T_C,
    1025,
    tolerance=1e-2
  )
})


test_that("Solidus front-end does its job", {
  expect_equal(
    Solidus_ATAC4(Tlow = 650,
                  Thigh = 1250,
                  maxIter = 32,
                  AmountTolerance = 0.5 / 100,
                  verbose = F,
                  printPhaseBoundaryResults = T,
                  ### arguments passed to MAGEMin() ###
                  Xoxides = Xoxides,
                  X = rock,
                  Pkbar = 4.5,
                  showMinimizationResults=F)$T_C,
    695.9,
    tolerance=1e-2
  )
})


test_that("Liqidus is found", {
  expect_equal(
    Liquidus_ATAC4(Tlow = 800,
                  Thigh = 1500,
                  maxIter = 32,
                  AmountTolerance = 0.5 / 100,
                  verbose = F,
                  printPhaseBoundaryResults = T,
                  ### arguments passed to MAGEMin() ###
                  Xoxides = Xoxides,
                  X = rock,
                  Pkbar = 4.5,
                  showMinimizationResults=F)$T_C,
    1024.2,
    tolerance=1e-2
  )
})
