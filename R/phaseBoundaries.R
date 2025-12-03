### Find upper/lower stability bounds of phases, using the bissection method

findSolidus <- function(Tlow = 550, Thigh = 1450,
                        maxIter = 32, AmountTolerance = 0.5 / 100,
                        verbose = F, printPhaseBoundaryResults = T,
                        ### arguments passed to MAGEMin()
                        ...) {
  #' @rdname findLowerStability
  #' @export
  out <- findLowerStability(phasename = "liq",
                            Tlow = Tlow, Thigh = Thigh,
                            maxIter = maxIter, AmountTolerance = AmountTolerance,
                            verbose = verbose, 
                            printPhaseBoundaryResults = F,
                            ...)

### Print results ####
  if (printPhaseBoundaryResults) {
    if (out$lower_liq_Boundary_found) {
      cat(
        "Solidus found at",
        round(out$T_C, 1),
        "°C after",
        out$iterations,
        "iterations.\n"
      )
    } else{
      liq <- getPhProp(out)["liq"]
      if (is.na(liq)) {
        liqText <- "no liquid"
      } else{
        liqText <- paste(round(liq * 100, 1), "wt%. liquid")
      }
      cat(
        "Solidus not found. Closest approx:",
        round(out$T_C, 1),
        "°C with",
        liqText,
        "after",
        out$iterations,
        "iterations.\n"
      )
    }
  } ### end of print block
   
  # Rename to something more explicit
  names(out)[names(out) == "lower_liq_Boundary_found"] <- "solidus_found"
  
  return(out)
}

findLowerStability <- function(phasename,
                        Tlow = 550, Thigh = 1450,
                        maxIter = 32, AmountTolerance = 0.5 / 100,
                        verbose = F, printPhaseBoundaryResults = T,
                        ### arguments passed to MAGEMin()
                        ...) {
  #' Find phase boundaries
  #' 
  #' @description
    #' Various functions that allow to find lower and upper stability bounds for phases.
    #' `findLowerStability()` and `findLowerStability()` can be used on any phase, specified with `phasename=...`.
    #' `findSolidus()` is a convenience function, strictly equivalent to `findLowerStability(phasename="liq",...)`.
    #' 
  #' @param phasename (string) the name of the phase to look for, as known from MAGEMin.
  #' @param Tlow,Thigh Initial search range (see details)
  #' @param maxIter Maximum number of iterations (see details)
  #' @param AmountTolerance The value of the phase under which the boundary is considered to be reached
  #' @param verbose Report information on every step (mostly for debuging purposes)
  #' @param printPhaseBoundaryResults Report on the outcome of the calculation
  #' @param ... Further arguments that will be passed to [MAGEMin()]
  #' @returns A list containing the MAGEMin output at the phase boundary, converted to a R list.
  #'  This is the same format as in [MAGEMin()] except that 2 fields have been added, 
  #'  `$lower_<PHASE>_Boundary_found` (a boolean) and `$iterations` (the number of iterations)
  #'  @details
    #'  These functions use the bissection method to locate the position of a phase boundary, i.e. where
    #'  its abundance becomes 0. A lower phase boundary (including solidus) is defined as the lowest temperature at
    #'  which a phase is present, an upper phase boundary as the highest. This works well for
    #'  phases with simple range (for instance liquid which is never present under the solidus and always above). For phases
    #'  that are present only in part of the range (most other minerals), this may become problematic if the initial guess is not correct.
    #'  
    #'  The run time strongly depends on the number of iterations needed. They depend marginally on the original temperature range `[Tlow,Thigh]`
    #'  (because the nature of the method is that the range is halved each time, so doubling the initial start range
    #'  only adds one iteration) but much more strongly on the resolution required (`AmountTolerance`). The max number
    #'  of iterations `maxIter` is not very likely to be reached (unless you put a very small number)
    #'  
    #'  So the user has to decide on a crude but reasonably quick approximation (high tolerance), or a slower but finer one (low tolerance). 
    #'  
    #'  Fixing the initial range is somewhat more tricky. For a simple phase (liquid) there is no reason to keep it too narrow.
    #'  For a more complex stability field, you probably need to have an idea of what you are doing.
    #'  
    #'  One of the reason a computation may run for a long time is when the system does not find a solution, and
    #'  is "stuck" at one of the bounds (upper/lower). Possible solutions would be (i) to decrease the max number of iterations
    #'  or (ii) to widen the bound - remember, this adds only one or two iterations but ensures that you will find
    #'  the solution, most likely long before reaching `maxIter`.
    #'  
  #' @export
  iteration <- 0
  phprop <- 1
  Tcurrent <- (Thigh + Tlow) / 2
  
  while (iteration < maxIter && phprop > AmountTolerance){
    
    Tcurrent <- (Thigh + Tlow) / 2
    iteration <- iteration + 1
    
    out <- MAGEMin(TC = Tcurrent, ...)
    ph <- getPhProp(out)
    phprop <- ph[phasename]
    
    if (verbose) {
      cat("\n-------------\nIteration:",
          iteration,
          "\n-------------\n")
      cat("T current:", Tcurrent, "\n")
      cat(phasename,": ", phprop, "\n",sep="")
    }
    
    if (is.na(phprop)) {
      # Phase not present, too cold !
      Tlow <- Tcurrent
      phprop <- 2 # a trick to stay in the while loop...
    } else{
      # liquid present, too hot !
      Thigh <- Tcurrent
    }
    if (verbose) {
      cat("Trange:", round(Tlow, 1), "-", round(Thigh, 1), "\n")
    }
    
    ## DEBUG
    # cat("Booleans:\n")
    # cat(iteration < maxIter,"\n")
    # cat(phprop > AmountTolerance,"\n")
    
  } # end of while loop
  
  ### Final steps
  if (verbose) {
    cat("\n================================\n")
  }
  
  boundFieldName <- paste("lower",phasename,"Boundary_found",sep="_")
    
  if (phprop < AmountTolerance) {
    out[[boundFieldName]] <- TRUE
    if (printPhaseBoundaryResults) {
      cat("Lower boundary found at",
          round(Tcurrent, 1),
          "°C after",
          iteration,
          "iterations.\n")
    }
  } else{
    out[[boundFieldName]] <- FALSE
    if (printPhaseBoundaryResults) {
      if (phprop == 2) {
        phText <- paste("no",phasename)
      } else{
        phText <- paste(round(phprop * 100, 1), "wt%.")
      }
      cat(
        "Lower boundary not found. Closest approx:",
        round(Tcurrent),
        "°C with",
        phText,
        "after",
        iteration,
        "iterations.\n"
      )
    }
    
  }
  out$iterations <- iteration
  
  invisible(out)
}


findUpperStability <- function(phasename,
                               Tlow = 550, Thigh = 1450,
                               maxIter = 32, AmountTolerance = 0.5 / 100,
                               verbose = F, printPhaseBoundaryResults = T,
                               ### arguments passed to MAGEMin()
                               ...) {
  #' @rdname findLowerStability
  #' @export
  #' 
  iteration <- 0
  phprop <- 1
  Tcurrent <- (Thigh + Tlow) / 2
  
  while (iteration < maxIter && phprop > AmountTolerance){
    
    Tcurrent <- (Thigh + Tlow) / 2
    iteration <- iteration + 1
    
    out <- MAGEMin(TC = Tcurrent, ...)
    ph <- getPhProp(out)
    phprop <- ph[phasename]
    
    if (verbose) {
      cat("\n-------------\nIteration:",
          iteration,
          "\n-------------\n")
      cat("T current:", Tcurrent, "\n")
      cat(phasename,": ", phprop, "\n",sep="")
    }
    
    if (is.na(phprop)) {
      # Phase not present, too hot !
      Thigh <- Tcurrent
      phprop <- 2 # a trick to stay in the while loop...
    } else{
      # liquid present, too cold !
      Tlow <- Tcurrent
    }
    if (verbose) {
      cat("Trange:", round(Tlow, 1), "-", round(Thigh, 1), "\n")
    }
    
    ## DEBUG
    # cat("Booleans:\n")
    # cat(iteration < maxIter,"\n")
    # cat(phprop > AmountTolerance,"\n")
    
  } # end of while loop
  
  ### Final steps
  if (verbose) {
    cat("\n================================\n")
  }
  
  boundFieldName <- paste("upper",phasename,"Boundary_found",sep="_")
  
  if (phprop < AmountTolerance) {
    out[[boundFieldName]] <- TRUE
    if (printPhaseBoundaryResults) {
      cat("Upper boundary found at",
          round(Tcurrent, 1),
          "°C after",
          iteration,
          "iterations.\n")
    }
  } else{
    out[[boundFieldName]] <- FALSE
    if (printPhaseBoundaryResults) {
      if (phprop == 2) {
        phText <- paste("no",phasename)
      } else{
        phText <- paste(round(phprop * 100, 1), "wt%.")
      }
      cat(
        "Upper boundary not found. Closest approx:",
        round(Tcurrent),
        "°C with",
        phText,
        "after",
        iteration,
        "iterations.\n"
      )
    }
    
  }
  out$iterations <- iteration
  
  invisible(out)
}
