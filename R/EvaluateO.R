EvaluateO <- function(FeOt,FeRatio){
  #' Estimate the amount of O2 for a given FeOt and FeO/Fe2O3 mass ratio
  #' @param FeOt amount of FeOt, wt.pct.
  #' @param FeRatio FeO/Fe2O3 weight ratio (Middlemost 1989)
  #' @export

  M2 <- 72  # 56 + 16
  M3 <- 160 # 56*2 + 16*3
  RM <- M2/M3

  return(FeRatio * FeOt * ( 1 - 2* RM) / (2 * FeRatio * RM +1)   )

}
