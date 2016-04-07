## ------------------------------------------------------------------------
  library(BASFOR)

## ------------------------------------------------------------------------
init           <- initialise()

## ------------------------------------------------------------------------
  FORTYPE        <- as.integer(1)

## ------------------------------------------------------------------------
  year_start     <- as.integer(1900)
  doy_start      <- as.integer(1)
  NDAYS          <- as.integer(40543)

## ------------------------------------------------------------------------
  head(weather_CONIFEROUS_1)
  df_weather     <- weather_CONIFEROUS_1
  matrix_weather <- weather_BASFOR( year_start, doy_start, NDAYS, df_weather )

## ------------------------------------------------------------------------
  row.names(df_params)
  df_params      <- df_params
  parcol         <- 1  
  params         <- df_params[,parcol]

## ------------------------------------------------------------------------
  init$calendar_Ndep [ 1, ] <- c( 1700,200,  3.0 / (365 * 10000) )
  init$calendar_Ndep [ 2, ] <- c( 1900,200,  3.0 / (365 * 10000) )
  init$calendar_Ndep [ 3, ] <- c( 1910,200,  4.7 / (365 * 10000) )
  init$calendar_Ndep [ 4, ] <- c( 1920,200,  6.4 / (365 * 10000) )
  init$calendar_Ndep [ 5, ] <- c( 1930,200,  8.9 / (365 * 10000) )
  init$calendar_Ndep [ 6, ] <- c( 1940,200, 12.4 / (365 * 10000) )
  init$calendar_Ndep [ 7, ] <- c( 1950,200, 17.2 / (365 * 10000) )
  init$calendar_Ndep [ 8, ] <- c( 1960,200, 23.9 / (365 * 10000) )
  init$calendar_Ndep [ 9, ] <- c( 1970,200, 30.9 / (365 * 10000) )
  init$calendar_Ndep [10, ] <- c( 1980,200, 35.0 / (365 * 10000) )
  init$calendar_Ndep [11, ] <- c( 1990,200, 32.2 / (365 * 10000) )
  init$calendar_Ndep [12, ] <- c( 2000,200, 25.1 / (365 * 10000) )
  init$calendar_Ndep [13, ] <- c( 2010,200, 21.2 / (365 * 10000) )
  init$calendar_Ndep [14, ] <- c( 2100,200, 21.2 / (365 * 10000) )

## ------------------------------------------------------------------------
  init$calendar_thinT[ 1, ] <- c( 1925,  1, 0.400 )
  init$calendar_thinT[ 2, ] <- c( 1935,  1, 0.250 )
  init$calendar_thinT[ 3, ] <- c( 1945,  1, 0.250 )
  init$calendar_thinT[ 4, ] <- c( 1955,  1, 0.200 )
  init$calendar_thinT[ 5, ] <- c( 1965,  1, 0.150 )
  init$calendar_thinT[ 6, ] <- c( 1975,  1, 0.400 )
  init$calendar_thinT[ 7, ] <- c( 2006,  1, 0.887 )

## ------------------------------------------------------------------------
  output      <- run_model()

## ---- fig.show='hold',fig.height=8, fig.width=8, fig.cap = "BASFOR model output for a Coniferous site"----
  outputUnits <- init$outputUnits
  outputNames <- init$outputNames
  plot_output()

