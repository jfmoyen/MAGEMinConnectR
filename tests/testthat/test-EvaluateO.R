test_that("Fe ratios are calculated as they should", {
  expect_equal(EvaluateO(FeOt=2,FeRatio=0.4),
               0.0588,
               tolerance = 1e-3)

  expect_equal(EvaluateO(FeOt=3,FeRatio=0.1),
               0.0275,
               tolerance = 1e-3)

  expect_equal(EvaluateO(FeOt=4,FeRatio=0.6),
               0.1558,
               tolerance = 1e-3)
})
