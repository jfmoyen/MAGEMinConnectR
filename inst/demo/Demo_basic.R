#######################################################
#### Demonstration of the use of MAGEMin connector
#### Basic usage
#######################################################


library(MAGEMinConnectR)

#### Simple case - one compo, manually defined ####

MAGEMin_setup(JULIA_HOME="C:\\Program Files\\Julia-1.11.5\\bin")

# MAGEMin oxides - here for the ig dataset
Xoxides <- c("SiO2", "Al2O3", "CaO", "MgO", "FeO",
             "K2O", "Na2O", "TiO2", "O", "H2O")

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

# The first minimization always takes a bit longer
# because julia compiles the code the first time...
out <- MAGEMin(Xoxides,X = rock,
               Pkbar = 4.5,TC=900,
               showMinimizationResults=T)

# Phase Proportions
getPhProp(out)

# Phase composition
getSSComp(out,"liq")
getSSComp(out,"pl","emmol")

#### With a buffer ####

MAGEMin_setup(JULIA_HOME="C:\\Program Files\\Julia-1.11.5\\bin",extraMAGEMinParams = 'buffer=\"qfm\"')

# MAGEMin oxides - here for the ig dataset
Xoxides <- c("SiO2", "Al2O3", "CaO", "MgO", "FeO",
             "K2O", "Na2O", "TiO2", "O", "H2O")

rock = c(SiO2 = 57.28,
         Al2O3 = 16.5,
         CaO = 7.12,
         MgO = 5.42,
         FeO = 8.48/1.111,
         K2O = 0.76,
         Na2O = 3.52,
         TiO2 = 0.67,
         O = 10, # Excess
         H2O = 5.0
)  #ATAC-4

out <- MAGEMin(Xoxides,X = rock,
               Pkbar = 4.5,TC=1000,
               extraMAGEMinParams = 'B=2.0',
               showMinimizationResults=T)

# Phase Proportions
getPhProp(out)

# Phase composition
getSSComp(out,"liq")
getSSComp(out,"pl","emmol")

