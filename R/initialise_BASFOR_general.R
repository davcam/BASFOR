### initialise_BASFOR_general.R ##

################################################################################
### 1. MODEL LIBRARY FILE & FUNCTION FOR RUNNING THE MODEL
################################################################################

#' Run the forest model BASFOR 
#' @details This is a wrapper function to run BASFOR.
#' @param ft Integer which indicates the forest type and is set to 1 for coniferous and 2 for deciduous
#' @param p Vector that contain the model parameters. The order of the parameters in the vector is important and should be consistent with the provided dataset df_params.
#' @param w Matrix that contains the model input weather data with dimensions [number of timesteps, number of weather inputs]. The seven weather inputs are: year, day of year, radiation (MJ/day), mean air temperature (deg C), precipitation (mm/day), surface wind speed (m/s) and vapour pressure (kPa)
#' @param calf Matrix timeseries of fertilization of dimension [length of timeseries, date vector of length 3]. Where the date vector is made up of (year, day of year, fertilization in kgN /day /m^2)   
#' @param calN Matrix timeseries of N depostion of dimension [length of timeseries, date vector of length 3]. Where the date vector is made up of (year, day of year, N deposition in kgN /day /m^2)
#' @param calpT Matrix timeseries of prunings of dimension [length of timeseries, date vector of length 3]. Where the date vector is made up of (year, day of year, pruning fraction)
#' @param caltT Matrix timeseries of forest thinnings of dimension [length of timeseries, date vector of length 3]. Where the date vector is made up of (year, day of year, thinning fraction)
#' @param n Length of run in days
#' @param NOUT Number of output variables
#' @return Returns a matrix of dimension [length of run in days, number of outputs]
#' @examples
#' FORTYPE        <- as.integer(1)
#' params         <- df_params[,1]
#' year_start     <- as.integer(1900); doy_start <- as.integer(1); NDAYS <- as.integer(40543)
#' df_weather     <- weather_CONIFEROUS_1
#' matrix_weather <- weather_BASFOR( year_start, doy_start, NDAYS, df_weather )
#' init           <- initialise()
#' output         <- run_model()
#' @export
#' @useDynLib BASFOR
run_model <- function(ft      = FORTYPE,
                      p       = params,
                      w       = matrix_weather,
                      calf    = init$calendar_fert,
                      calN    = init$calendar_Ndep,
                      calpT   = init$calendar_prunT,
                      caltT   = init$calendar_thinT,
                      n       = NDAYS,
                      NOUT    = init$NOUT) {

  output <- .Fortran('BASFOR',ft,p,w,calf,calN,calpT,caltT,n,NOUT,
                     matrix(0,n,NOUT)) [[10]]
  return(output)
}

################################################################################
### 2. FUNCTION FOR READING WEATHER DATA
################################################################################
#' Create weather data matrix for BASFOR 
#' @details This function converts the weather data from a data frame into an R matrix (matrix_weather) which can be passed to BASFOR via the function run_model.
#' @param y Integer giving the start year
#' @param d Integer giving the starting day of year
#' @param n Length of run in days
#' @param df_weather dataframe containing input weather data. The columns should be year, day of year, radiation (MJ/day), mean air temperature (deg C), precipitation (mm/day), surface wind speed (m/s) and vapour pressure (kPa). See for example dataset weather_CONIFEROUS_1.
#' @return Returns a matrix (matrix_weather) of dimension [length in days, 7] which is used to input weather to BASFOR in the run_model function.
#' @examples
#' year_start     <- as.integer(1900); doy_start <- as.integer(1); NDAYS <- as.integer(40543)
#' df_weather     <- weather_CONIFEROUS_1
#' matrix_weather <- weather_BASFOR( year_start, doy_start, NDAYS, df_weather )
#' @export
weather_BASFOR <- function(y = year_start,
                           d = doy_start,
                           n = NDAYS,
                           df_weather) {
  row_start             <- 1
  while( df_weather[row_start,]$year< y ) { row_start <- row_start+1 }
  while( df_weather[row_start,]$doy < d ) { row_start <- row_start+1 }
  df_weather_sim        <- df_weather[row_start:(row_start+n-1),]
  NMAXDAYS              <- as.integer(60000)
  NWEATHER              <- as.integer(7)
  matrix_weather        <- matrix( 0., nrow=NMAXDAYS, ncol=NWEATHER )
  matrix_weather[1:n,1] <- df_weather_sim$year
  matrix_weather[1:n,2] <- df_weather_sim$doy
  matrix_weather[1:n,3] <- df_weather_sim$GR
  matrix_weather[1:n,4] <- df_weather_sim$T
  matrix_weather[1:n,5] <- df_weather_sim$RAIN   
  matrix_weather[1:n,6] <- df_weather_sim$WN
  matrix_weather[1:n,7] <- df_weather_sim$VP
  return(matrix_weather)
}
   
################################################################################
### 3. OUTPUT VARIABLES
################################################################################
################################################################################
#' Initialise variables needed to run BASFOR 
#' @details This function initilaises some matrices needed to run BASFOR
#' @return Returns a list containing calendar_fert, calendar_Ndep, calendar_prunT, calendar_thinT, outputNames, outputUnits, NOUT needed by functions run_model and plot_output to run BASFOR and plot outputs.
#' @examples
#' init <- initialise()
#' @export
initialise <- function ()
{
 calendar_fert  <- matrix( -1, nrow=100, ncol=3 )
 calendar_Ndep  <- matrix( -1, nrow=100, ncol=3 )
 calendar_prunT <- matrix( -1, nrow=100, ncol=3 )
 calendar_thinT <- matrix( -1, nrow=100, ncol=3 )

 outputNames <- c(
  "Time"           , "year"          , "doy"          ,
  "T"              , "RAIN"          ,
  "CR"             , "Cwood"         , "CL"           , "CLBS"         , "CRES"      ,
  "NL"             ,
  "WC"             , "CLITT"         , "CSOM"         , "Csoil"        , "Nsoil"     ,
  "DBH"            , "H"             , "LAI"          ,
  "Rsoil"          , "NEE_gCm2d"     , "GPP_gCm2d"    , "Reco_gCm2d"   ,
  "ET_mmd"         , "NemissionN2O"  , "NemissionNO" )
  
 outputUnits <- c(
   "y"             , "y"             , "d"            ,
   "(degC)"        , "(mm d-1)"      ,
   "(t C ha-1)"    , "(t C ha-1)"    , "(t C ha-1)"   , "(t C ha-1)"   , "(t C ha-1)",
   "(kg N m-2)"    ,
   "(m3 m-3)"      , "(kg C m-2)"    , "(kg C m-2)"   , "(t C ha-1)"   , "(t N ha-1)",
   "(m)"           , "(m)"           , "(m2 m-2)"     ,
   "(kg C m-2 d-1)", "(g C m-2 d-1)" , "(g C m-2 d-1)", "(g C m-2 d-1)",
   "(mm d-1)"      , "(kg N m-2 d-1)", "(kg N m-2 d-1)" )  
  
 NOUT <- as.integer( length(outputNames) )
return( list(calendar_fert=calendar_fert, calendar_Ndep=calendar_Ndep, calendar_prunT=calendar_prunT, calendar_thinT=calendar_thinT, outputNames=outputNames, outputUnits=outputUnits, NOUT=NOUT) )
}   

################################################################################
### 4. FUNCTIONS FOR EXPORTING THE RESULTS TO FILE (pdf with plots, txt with table)
################################################################################
#' Function for plotting BASFOR model output 
#' @details This function plots BASFOR model output
#' @param list_output a list of model outputs
#' @param vars Vector of characters strings of those output variable names which are to be plotted.
#' @examples
#' output <- run_model
#' plot_output()
#' @export
plot_output <- function(
  list_output = list(output),
  vars        = init$outputNames[-(1:3)]
) {
  outputNames = init$outputNames
  outputUnits = init$outputUnits
  leg         = paste( "Run", 1:length(list_output) )
  leg_title   = "LEGEND"
  nrow_plot   = ceiling( sqrt((length(vars)+1) * 8/11) )
  ncol_plot   = ceiling( (length(vars)+1)/nrow_plot )
  lty         = rep(1,length(list_output))
  lwd         = rep(3,length(list_output))
  par( mfrow=c(nrow_plot,ncol_plot), mar=c(2, 2, 2, 1) )
  if (!is.list(list_output)) list_output <- list(list_output) ; nlist <- length(list_output)
  col_vars <- match(vars,outputNames)                         ; nvars <- length(vars)
  for (iv in 1:nvars) {
    c       <- col_vars[iv]
    g_range <- range( sapply( 1:nlist, function(il){range(list_output[[il]][,c])} ) )
    plot( list_output[[1]][,1], list_output[[1]][,c],
          xlab="", ylab="", cex.main=1,
          main=paste(outputNames[c]," ",outputUnits[c],sep=""),
          type='l', col=1, lty=lty[1], lwd=lwd[1], ylim=g_range )
    if (nlist >= 2) {
      for (il in 2:nlist) {
      points( list_output[[il]][,1], list_output[[il]][,c],
              col=il, type='l', lty=lty[il], lwd=lwd[il] )            
      }
    }
    if ( (iv%%(nrow_plot*ncol_plot-1)==0) || (iv==nvars) ) {
      plot(1,type='n', axes=FALSE, xlab="", ylab="")
      legend("bottomright", leg, lty=lty, lwd=lwd, col=1:nlist, title = leg_title)
    }
  }
}

table_output <- function(
  list_output = list(output),
  vars        = outputNames[-(1:1)],
  file_table  = paste( "output_", format(Sys.time(),"%H_%M.txt"), sep="" ),
  leg         = paste( "Run", 1:length(list_output) )
) {
  if (!is.list(list_output)) list_output <- list(list_output) ; nlist <- length(list_output)
  col_vars <- match(vars,outputNames)                         ; nvars <- length(vars)
  table_output         <- c("day", list_output[[1]][,1:1] )
  for (il in 1:nlist) {
    table_il     <- if (nvars==1) c    (vars, list_output[[il]][,col_vars]) else
                                  rbind(vars, list_output[[il]][,col_vars])
    table_output <- cbind( table_output, table_il ) 
  }
  colnames(table_output) <- c( "",rep(leg,each=nvars) )
  write.table( table_output, file_table, sep="\t", row.names=F )
}
   
################################################################################
### 5. FUNCTIONS FOR ANALYSIS

#######################
### 5.1 Function 'SA()'
#######################
SA <- function( parname_SA = "KAC",
                pmult      = 2^(-1:1),
                vars       = outputNames[-(1:3)],
                leg_title  = parname_SA,
                nrow_plot  = ceiling( sqrt((length(vars)+1) * 8/11) ),
                ncol_plot  = ceiling( (length(vars)+1)/nrow_plot ),
                lty        = rep(1,length(pmult)),
                lwd        = rep(3,length(pmult)),
                file_init  = "initialisation/initialise_BASFOR_DECIDUOUS_DE-Hai.R",
                file_plot  = paste("SA_",parname_SA,format(Sys.time(),"_%H_%M.pdf"),sep=""),
                file_table = paste("SA_",parname_SA,format(Sys.time(),"_%H_%M.txt"),sep="")
) {
  source(file_init)
  cat( "SA initialised for:", substr(basename(file_init),1,nchar(basename(file_init))-2), "\n")
  ip_SA          <- match( parname_SA, row.names(df_params) )
  par_SA_default <- params[ip_SA]
  nmult          <- length(pmult)
  list_output    <- vector( "list", nmult )
  for (im in 1:nmult) {
    params[ip_SA]     <- par_SA_default * pmult[im]
    list_output[[im]] <- run_model(p=params)
  }
  pdf( file_plot, paper="a4r", width=11, height=8 )
  plot_output( list_output, vars=vars,
               leg=as.character(pmult*par_SA_default), leg_title=parname_SA,
               nrow_plot=nrow_plot, ncol_plot=ncol_plot, lty=lty, lwd=lwd )
  dev.off()
  table_output(list_output, vars=vars,
               file_table=file_table,
               leg=paste(parname_SA,"=",pmult*par_SA_default,sep=""))
}
