  
  
  ############################################################################
  ### BacktetDAA - RPM - CONDITIONAL STATE PROBABILITY DISTRIBUTION
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     16.04.2021
  # First version:    28.05.2020
  # --------------------------------------------------------------------------
  
  
  # This file demonstrates the procedure to compute conditional distributions
  # of the state probability
  
  
  
  require(RP)
  require(slolz)
  require(DAARC)
  wd <- "R:/Asset_Management/Research_Projects/Multi_Asset/DAA/DAARC/BacktestDAA/"
  wd_warehouse <- paste0(wd, "RPM/waRehouse/")
  source( paste0(wd, "Source/custom_functions.R") )
  
  
  
  
  # --------------------------------------------------------------------------
  # BEAR AND BULL STATE SERIES
  # --------------------------------------------------------------------------
  
  BBS <- BBSRC$new()
  BBS$setCtrl()
  BBS$setData()
  BBS$runRobust()
  states <- (1 + BBS$output$states) / 2
  states <- states[isWeekday(time(states)), ]
  BBS$plotStates()
  
  
  
  # --------------------------------------------------------------------------
  # RISK INDICATOR SERIES
  # --------------------------------------------------------------------------
  
  
  VRP <- loadSensor(sensor_name = "vrp")
  Turbulence <- loadSensor(sensor_name = "turbulence_scl2_dm")
  BBTurbulence <- loadSensor(sensor_name = "bbturbulence_base_dm")
  BCPOLZ <- loadSensor(sensor_name = "bcpolz")
  Eigen <- loadSensor(sensor_name = "eigen")
  Momentum <- loadSensor(sensor_name = "momentum")
  OFRFSI <- loadSensor(sensor_name = "ofrfsi")
 
  signals <- na.omit( cbind( vrp_signal = VRP$getSignal(),
                             vrp_raw = VRP$signal$base[ ,"vrp_raw"],
                             turbulence = Turbulence$getSignal(),
                             bbturbulence = BBTurbulence$signal$base[ ,"md_delta"],
                             bcpolz = BCPOLZ$getSignal(),
                             eigen = Eigen$getSignal(),
                             momentum = Momentum$getSignal(),
                             ofrfsi = OFRFSI$getSignal() ) )
  
  signals <- OFRFSI$getSignal()
  
  
  # Use return series as 'signals'
  signals <- getMSCIData( universe = "dm", frqncy = "d" )
  signals <- signals[isWeekday(time(signals)), ]
  

  
  # --------------------------------------------------------------------------
  # CONDITIONAL DISTRIBUTION OF STATE PROBABILITIES
  # --------------------------------------------------------------------------
  
  n_lag <- 5
  TD <- trainingData( Y_train = window(states, start(states), "2020-02-01"),
                      X_train = signals,
                      n_lag = n_lag )
  today <- "2020-01-31"
  x_eval <- TD$X_train[today, ]
  
  # Kernel on X
  # debugonce( weightsFun.kernel )
  alpha <- weightsFun.kernel( data = TD, x_eval = x_eval )
  plot(alpha)
  entropy( alpha, exponential = TRUE ) / length(alpha)
  
  # Dirichlet sampling
  m <- 10^3 + 1
  P <- rdirichlet( n = m, alpha = alpha * length(alpha) )
  P_e <- rdirichlet( n = m, alpha = alpha * entropy( alpha, exponential = TRUE ) )
  P_fev <- fevBias( x  = P, q = 100 )
  
  # Expected states as RPM weighted average of state series
  y_train <- as.numeric(TD$Y_train)
  es <- apply( P, 1, function(p) { t(p) %*% y_train })
  es_e <- apply( P_e, 1, function(p) { t(p) %*% y_train })
  es_fev <- apply( P_fev, 1, function(p) { t(p) %*% y_train })
 
  ldens <- apply( cbind(es, es_e, es_fev), 2, densFUN, from = 0, to = 1 )
  plot.ldensity( ldens, fillin = FALSE, colors = 1:length(ldens) )
  
  lapply( ldens, FUN = function(x) { entropy( x$y, exponential = TRUE) / length(x$y) } )
  
  
  
  # Out-of-sample test set
  X_eval <- window(signals, "2020-01-01", "2020-06-10")
  DS <- dirichletSampling( Y_train = states,
                           X_train = signals,
                           X_eval = X_eval,
                           weights_fun = "kernel",
                           # sclfct_lp = 10 )
                           scl_by_entropy = TRUE,
                           sclfct = NULL )
  names(DS) <- rownames(X_eval)  
  ldens <- lapply( DS, densFUN, from = 0, to = 1, n = 1000 )
  
  colors <- fBasics::divPalette(n = length(ldens), "Spectral")
  plot.ldensity( ldens, fillin = FALSE, colors = colors, cex = 0.5 )
  
  barplot( unlist(lapply(DS, mean)), col = colors )
  
  
  
  
  
  
  
  
  
  
    
  