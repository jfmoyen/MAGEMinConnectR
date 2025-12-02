###################
#
# The functions than do the real job
#
###################

MAGEMin_setup <- function(JULIA_HOME=Sys.getenv("JULIA_HOME"),
                          restart = F,
                          db = "ig",
                          sys_in = "wt",
                          extraMAGEMinParams="verbose=true"){
  #' Initiate julia and MAGEMin project. This may take up to a minute.
  #' @param JULIA_HOME: Environment variable containing the path to the
  #' julia directory. Probably something like C:\\Program Files\\Julia-1.11.5\\bin
  #' @param restart: Boolean. If true, julia instance is restarted (otherwise an existing instance is used)
  #' @param db Database to use. From https://github.com/ComputationalThermodynamics/MAGEMin_C.jl/blob/main/docs/src/MAGEMin_C/MAGEMin_C.md,
  #' it can be one of mtl -> mantle (Holland et al., 2013) - mp -> metapelite (White et al., 2014) - mb -> metabasite (Green et al., 2016) -
  #' ig -> igneous (Green et al., 2025 updated from and replacing Holland et al., 2018) - igad -> igneous alkaline dry (Weller et al., 2024) -
  #' um -> ultramafic (Evans & Frost, 2021) - sb11 -> Stixrude & Lithgow-Bertelloni (2011) - sb21 -> Stixrude & Lithgow-Bertelloni (2021) -
  #' ume -> ultramafic extended (Green et al., 2016 + Evans & Frost, 2021) - mpe -> extended metapelite (White et al., 2014 + Green et al., 2016 +
  #' Franzolin et al., 2011 + Diener et al., 2007) - mbe -> extended metabasite (Green et al., 2016 + Diener et al., 2007 + Rebay et al., 2022)
  #' Generally speaking, the thermocalc-family databases (mp, mb, ig) are better behaved, in particular they can separate
  #' solvus (e.g. afs/pl)
  #' @param sys_in System units ("wt", "mol")
  #' @param extraMAGEMinParams String. Extra parameters to be passed to MAGEMin minimization,
  #' for instance 'verbose = false, buffer= \"qfm\" ' (don't forget to protect quotation marks)
  #' If you want to use a buffer, don't forget to pass an offset in the actual MAGEMincall, as 'B=1.0' for instance (QFM+1).
  #' @details
  #' Similar to all functions in this package, the main goal is to pass commands to
  #' an underlying instance of julia. Very few checks happen R-side, so if something fails
  #' you will probably get julia's error messages.
  #'
  #' The function actully does 3 things: (i) it starts a julia instance, with
  #' JuliaCall::julia_setup(); (ii) it loads the MAGEMin library, with\code{using MAGEMin},
  #' executed in julia; (iii) it initializes a MAGEMin calculation with the required parameters.
  #'
  #' If you want to initialize a julia instance without MAGEMin, use the regular
  #' JuliaCall::julia_setup().
  #'
  #' Once a julia instance is setup, it will persist until the end of the session.
  #' All subsequent computations will happen in the same julia instance.
  #' In this case, the function will merely (re)initialize MAGEMin, overwriting
  #' previous setup, and it will be much quicker.
  #' @returns Nothing. This is purely used for side effects. However the function also triggers
  #' package variables MM$JULIA_LOADED, MM$MAGEMIN_LOADED and MM$MAGEMIN_INITIALIZED (to true).
  #' @import JuliaCall
  #' @export

  if( !the$JULIA_LOADED || restart == TRUE){
    cat("Setting up julia ... may take up to a few minutes\n")
    Sys.setenv(JULIA_HOME=JULIA_HOME)
    #Sys.setenv(JULIA_NUM_THREADS=6) # Not work but probably useless here

    cat("Starting Julia\n")
    julia_setup(force=restart)

    the$JULIA_LOADED <- TRUE
    the$MAGEMIN_LOADED <- FALSE
    the$MAGEMIN_INITIALIZED <- FALSE

    cat("Julia succesfully started.\n")
  }


  cat("Attempting to load MAGEMin\n")
  julia_command("using MAGEMin_C")
  the$MAGEMIN_LOADED <- TRUE

  cat("Initializing MAGEMin calculation\n")
  JPars <- if(!is.null(extraMAGEMinParams)){paste0(',',extraMAGEMinParams)}else{NULL}

  julia_command(paste0('db = \"',db,'\";'))
  julia_command(paste0('data = Initialize_MAGEMin(db',JPars ,');'))
  julia_command(paste0('sys_in  = \"',sys_in,'\";'))

  cat("MAGEMin ready!\n")
  the$MAGEMIN_INITIALIZED <- TRUE

  invisible(NULL)
}



MAGEMin <- function(Xoxides,X,Pkbar,TC,extraMAGEMinParams=NULL,nameSolvus=T,showResults=F){
  #' Do the actual MAGEMin calculation. This takes some time, typically 200-500 ms
  #' for one run for MAGEMin proper + a bit of overhead for the various post-processing.
  #' It is maybe clever to wrap it in a progress bar
  #' @param Xoxides Character vector. Name of the oxides to use. They should exist in the MAGEMin database.
  #' @param X Numeric vector. Composition of the system for each oxide, matching Xoxides
  #' @param Pkbar numeric. Pressure, in kbar. The function makes sure that an integer is converted to a decimal
  #' (julia makes the difference)
  #' @param TC numeric. Temperature, in degree celsius
  #' @param extraMAGEMinParams String or vector of strings. Extra parameters to be passed to MAGEMin minimization,
  #' for instance 'B=1.0' (to pass a buffer offset)
  #' @param nameSolvus boolean. Should the minerals with a solvus (e.g. pl and kfs) be renamed
  #' according to their composition? This corresponds to passing the option 'name_solvus = true'
  #' to MAGEMin (in julia).
  #' @param showResults Boolean. If true, print julia's summary of the minimlization.
  #' @returns A list containing the MAGEMin output, converted to a R list. Most fields of the MAGEMin object
  #' (technically a MAGEMin_C gmin_struct) are straightforwards to convert to R. However, the julia structure contains
  #' complex types for SS_vec (solution phases), mSS_vec (metastable solutions) and PP_vec (pure phases). SS_vec and PP_vec are
  #' also parsed and converted to (sub)lists.
  #' Note that both still contain complex types, that are NOT converted to R and remain as JuliaObjects. In SS_vec
  #' these are the end-members compositions, that are known by definition.
  #' @note A side effect of the function is also to create a variable out in julia, that can be accessed
  #' with julia_console() or julia_command() if needed.
  #' julia_command("out") will invoke julia's print method and display summary info.
  #' If the julia object needs to be preserved, the best is probably to do something like
  #' foo <- julia_eval("out") (in R)
  #' @import JuliaCall
  #' @export

  # browser()
  if(!the$MAGEMIN_INITIALIZED){
    cat("MAGEMin is not initialized! Run MAGEMin_setup()  first")
    return(NULL)
  }

  julia_command(paste0('P = ',sprintf("%1.1f",Pkbar),';' ) )
  julia_command(paste0('T = ',sprintf("%1.1f",TC),';' ) )

  JOx <- paste0('[\"',paste(Xoxides,collapse='\";\"'),'\"]')
  JX <- paste0('[',paste(X,collapse=';'),']')

  solvus <- if(nameSolvus){'name_solvus = true'}else{NULL}
  JPars <- paste(c(solvus,extraMAGEMinParams),collapse=",")
  if(JPars == ""){JPars = NULL}else{JPars <- paste0(",",JPars)}

  julia_command(paste('Xoxides  = ',JOx,";"))
  julia_command(paste('X  = ',JX,";"))

  # Actual minimization
  JCall <- paste0('out = single_point_minimization(P, T, data, X=X, Xoxides=Xoxides,sys_in=sys_in',JPars,');')
  # browser()
  julia_command(JCall)

  out <- JuliaStructToList(julia_eval("out;"))

  # Extract phase info
  if(out$SS_vec$size>0){
    out$SS_vec <- lapply( 1:as.integer(out$SS_vec$size),
                          FUN = function(ii){
                            current_phase <- julia_eval(paste0('out.SS_vec[',ii,'];'))
                            return(JuliaStructToList(current_phase))
                          })
  }else{
    out$SS_vec <- NULL
  }

  if(out$PP_vec$size>0){
    out$PP_vec <- lapply( 1:as.integer(out$PP_vec$size),
                          FUN = function(ii){
                            current_phase <- julia_eval(paste0('out.PP_vec[',ii,'];'))
                            return(JuliaStructToList(current_phase))
                          })
  }else{
    out$PP_vec <- NULL
  }
  # cat("Minimization succesful\n")

  if(showResults){
    julia_command("out")
  }

  return(out)
}
