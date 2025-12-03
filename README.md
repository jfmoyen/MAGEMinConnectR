
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MAGEMinConnectR

<!-- badges: start -->

<!-- badges: end -->

MAGEMinConnectR provides a wrapper to access MAGEMin
(<https://github.com/ComputationalThermodynamics/MAGEMin_C.jl>), a julia
package to calculate phase equlibrium. Specifically, a main function,
`MAGEMin(Xoxides,X,Pkbar,TC,...)` does single-point minimization for
given P, T and composition (X), allowing the user to wrap this into
their own R code.

Most users are probably more familiar with the MAGEMin App
(<https://github.com/ComputationalThermodynamics/MAGEMin>). There are in
fact three components: \* `MAGEMin` itself, a C code that does the heavy
loading; \* `MAGEMin_C.jl`, a julia wrapper that we exploit here; \*
`MAGEMinApp.jl`, a julia App that connects to the wrapper and plots
diagrams.

Here we are using the MAGEMin_C part, i.e. we write a wrapper to access
a wrapper…

## Installation

MAGEMinConnectR requires a working installation of MAGEMin, and a way to
reach it. Thus, the user must install (and configure):

1.  julia (<https://julialang.org>)
2.  The `MAGEMin_C` package, in julia
3.  R library `JuliaCall`, that contains the connector.

Installing `MAGEMin_C` is documented on their
[webpage](https://github.com/ComputationalThermodynamics/MAGEMin),
simply using julia’s package manager:

    julia> ]
    pkg> add MAGEMin_C

However, while you are at it, it is probably more sensible to also
install the app:

    julia> ]
    pkg> add MAGEMinApp

(which will take care of dependencies and thus also install
`MAGEMin_C` - and you get to use the App if you want!)

`JuliaCall` is installed via the usual channels,
`install.packages("JuliaCall")` or via their github page
(<https://github.com/JuliaInterop/JuliaCall>). julia can also be
installed from there.

Finally, MAGEMinConnectR is a regular R library. You can download it
from the releases sections here and install it from a local file; or use

``` r
devtools::install_github("jfmoyen/MAGEMinConnectR")
```

## Main capacities

MAGEMinConnectR can do the following things:

- Run a single-point MAGEMin calculation (and retrieve the results
  in R) - with `MAGEMin()`
- Find the stability limit of a phase - with `findLowerStability()` and
  `findUpperStability()`, and, in the case of the liquid, with
  `findSolidus()`

## Example

You will find examples file in demo/Demo_xxx.R within the package
directory. A basic example of a single point minimization would be

``` r
library(MAGEMinConnectR)
#> Le chargement a nécessité le package : JuliaCall

## Setup
# You need to point to the directory where your Julia executable is located. 
# If the environment variable JULIA_HOME is properly set, it can be omitted.
# The setup function also has options for dataset to use, units etc.
# use ?MAGEMin_setup for details.
MAGEMin_setup(JULIA_HOME="C:\\Program Files\\Julia-1.11.5\\bin") 
#> Setting up julia ... may take up to a few minutes
#> Starting Julia
#> Julia version 1.11.5 at location C:\PROGRA~1\JULIA-~2.5\bin will be used.
#> Loading setup script for JuliaCall...
#> Finish loading setup script for JuliaCall.
#> Julia succesfully started.
#> Attempting to load MAGEMin
#> Initializing MAGEMin calculation
#> MAGEMin ready!


## MAGEMin oxides - here for the ig dataset
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

## Minimization
# The first minimization always takes a bit longer
# because julia compiles the code the first time...
out <- MAGEMin(Xoxides,X = rock,
               Pkbar = 4.5,TC=900,
               showMinimizationResults=T)

## Phase Proportions
getPhProp(out)
#>        amp        liq         pl        cpx        opx 
#> 0.22219649 0.62370215 0.09213501 0.01442538 0.04754097

## Phase composition
getSSComp(out,"liq") # in the default system unit, wt%
#>        SiO2       Al2O3         CaO         MgO         FeO         K2O 
#> 0.621025461 0.143612555 0.039671654 0.020143336 0.040671571 0.011870862 
#>        Na2O        TiO2           O       Cr2O3         H2O 
#> 0.043285563 0.005705379 0.000720896 0.000000000 0.073292723
getSSComp(out,"pl","emmol") # in end-member mol fraction
#>           ab           an          san 
#> 0.2049218040 0.7947160002 0.0003621958
# Note that MAGEMin automatically "sorts" solvus minerals (fsp -> ksp, pl) unless you suppress
# this option with nameSolvus=F, which is why we have pl here.
```
