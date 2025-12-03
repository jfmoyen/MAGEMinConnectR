test_that("Simple minimization works", {
# Note that failing this test is not necessarily an R bug. It may be a change
# in MAGEMin algorithm or a recalibration of the DB.

  expect_equal(
    getPhProp(Minimize_ATAC4(Xoxides,X = rock,
                             Pkbar = 4.5,TC = 900,
                             nameSolvus = T,
                             showMinimizationResults=F)),
    c(amp=0.2222, liq=0.6237, pl=0.0921, cpx=0.0144, opx=0.0475),
    tolerance = 1e-3,
    info = "This may also be due to a change in MAGEMin algo/DB"
  )

  expect_equal(
    getSSComp(Minimize_ATAC4(Xoxides,X = rock,
                             Pkbar = 4.5,TC = 900,
                             nameSolvus = T,
                             showMinimizationResults=F),"pl","emmol"),
    c(ab=0.2049,an=0.7947,san=0.0004),
    tolerance = 1e-3,
    info = "This may also be due to a change in MAGEMin algo/DB"
  )
})

test_that("Arguments are properly passed to MAGEMin",{

  expect_match(
    paste0(names(getPhProp(Minimize_ATAC4(Xoxides,X = rock,
                                          Pkbar = 4.5,TC = 900,
                                          nameSolvus = T,
                                          showMinimizationResults=F))),collapse=" "),
    "pl"
  )

  expect_match(
    paste0(names(getPhProp(Minimize_ATAC4(Xoxides,X = rock,
                                          Pkbar = 4.5,TC = 900,
                                          nameSolvus = F,
                                          showMinimizationResults=F))),collapse=" "),
    "fsp"
  )

})

