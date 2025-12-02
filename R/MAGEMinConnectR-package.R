#' @keywords internal

## usethis namespace: start
## usethis namespace: end

NULL
#' @section Installation:
#' MAGEMinConnectR requires a working installation of MAGEMin, and a way to reach it. Thus, the user must install (and configure):
#' 1. julia (https://julialang.org)
#' 2. The `MAGEMin_C` package, in julia
#' 3. R library `JuliaCall`, that contains the connector.
#'
#' Installing `MAGEMin_C` is documented on their webpage \url{https://github.com/ComputationalThermodynamics/MAGEMin}, simply using julia's package manager.
#' `JuliaCall` is installed via the usual channels, `install.packages("JuliaCall")`
#'
#' @section Basic usage:
#' You will find a longer example file in demo/Demo.R within the package directory. A basic example of a single point mnimization would be
#' ``` r
#'library(MAGEMinConnectR)
#'
#'MAGEMin_setup(JULIA_HOME="C:\\Program Files\\Julia-1.11.5\\bin")
#' # Or wherever your Julia executable is located. If the environment variable JULIA_HOME
#' # is properly set, it can be omitted.
#' # The setup function also has options for datasetr to use, units etc.
#' # use ?MAGEMin_setup.
#'
#' # MAGEMin oxides - here for the ig dataset
#' Xoxides <- c("SiO2", "Al2O3", "CaO", "MgO", "FeO",
#'             "K2O", "Na2O", "TiO2", "O", "H2O")
#'
#' # WR composition to use
#' rock = c(SiO2 = 57.28,
#'         Al2O3 = 16.5,
#'         CaO = 7.12,
#'         MgO = 5.42,
#'         FeO = 8.48/1.111,
#'         K2O = 0.76,
#'         Na2O = 3.52,
#'         TiO2 = 0.67,
#'         O = EvaluateO(8.48/1.111,0.3), # Recommended for basaltic andesite, Middlemost 1989
#'         H2O = 5.0 # julia needs float to be explicitely given as decimals (5.0 not 5)
#')  #ATAC-4
#'
#' # The first minimization always takes a bit longer
#' # because julia compiles the code the first time...
#' out <- MAGEMin(Xoxides,X = rock,
#'               Pkbar = 4.5,TC=900,
#'               showResults=T)
#'
#' # Phase Proportions
#' getPhProp(out)
#' ```
#' @section Details:
#' See the github page,  \url{https://github.com/jfmoyen/MAGEMinConnectR}
#'
"_PACKAGE"



