#######################################################
#### Demonstration of the use of MAGEMin connector
#### Solidus and liquidus
#######################################################

library(MAGEMinConnectR)

#### Find the solidus of a composition ####

MAGEMin_setup(JULIA_HOME="C:\\Program Files\\Julia-1.11.5\\bin")
# Or wherever your Julia executable is located. If the environment variable JULIA_HOME
# is properly set, it can be omitted.
# The setup function also has options for datasetr to use, units etc.
# use ?MAGEMin_setup.

# MAGEMin oxides - here for the ig dataset
Xoxides <- c("SiO2", "Al2O3", "CaO", "MgO", "FeO",
             "K2O", "Na2O", "TiO2", "O", "H2O")

# WR composition to use
rock = c(SiO2 = 57.28,
         Al2O3 = 16.5,
         CaO = 7.12,
         MgO = 5.42,
         FeO = 8.48/1.111,
         K2O = 0.76,
         Na2O = 3.52,
         TiO2 = 0.67,
         O = EvaluateO(8.48/1.111,0.3), # Recommended for basaltic andesite, Middlemost 1989
         H2O = 5.0 # julia needs float to be explicitely given as decimals (5.0 not 5)
)  #ATAC-4

## Solidus
out <- findSolidus(Tlow = 650,
                   Thigh = 1000,
                   maxIter = 32,
                   AmountTolerance = 0.5 / 100,
                   verbose = T,
                   printPhaseBoundaryResults = T,
                   ### arguments passed to MAGEMin() ###
                   Xoxides = Xoxides,
                   X = rock,
                   Pkbar = 4.5,
                   showMinimizationResults=F
            )
# 695.9 Liquid: 0.002210409
# This is really only a front-end for findLowerStability("liq",.. etc ...)

## General form:
out <- findLowerStability("opx",
                   Tlow = 800,
                   Thigh = 900,
                   maxIter = 32,
                   AmountTolerance = 0.5 / 100,
                   verbose = T,
                   printPhaseBoundaryResults = T,
                   ### arguments passed to MAGEMin() ###
                   Xoxides = Xoxides,
                   X = rock,
                   Pkbar = 4.5,
                   showMinimizationResults=F
)
# 837.5

## Upper bound
out <- findUpperStability("pl",
                          Tlow = 800,
                          Thigh = 1200,
                          maxIter = 32,
                          AmountTolerance = 0.5 / 100,
                          verbose = T,
                          printPhaseBoundaryResults = T,
                          ### arguments passed to MAGEMin() ###
                          Xoxides = Xoxides,
                          X = rock,
                          Pkbar = 4.5,
                          showMinimizationResults=F)
# 1025

## Liquidus
out <- findLiquidus(Tlow = 800,
                   Thigh = 1500,
                   maxIter = 32,
                   AmountTolerance = 0.5 / 100,
                   verbose = T,
                   printPhaseBoundaryResults = T,
                   ### arguments passed to MAGEMin() ###
                   Xoxides = Xoxides,
                   X = rock,
                   Pkbar = 4.5,
                   showMinimizationResults=F
)
# 1024.2 Liquid: 0.9989859

