  
  
  ############################################################################
  ### SHINY - BEAR AND BULL STATES - TESTING
  ############################################################################
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version      03.07.2020
  # First version     03.07.2020
  # --------------------------------------------------------------------------
  
  
  require(BBSolz)
  require(DAARC)
  
  wd <- "R:/Asset_Management/R/Shiny/Signal_Testing/Bear_Bull_States/"
  
  
  # BBS
  BBS <- loadSensor(sensor_name = "bbs")
  BBSFuzzy <- loadSensor(sensor_name = "bbsfuzzy")
  # X <- BBS$data$X_bm
  # X_level <- cumulated(X, "discrete")
  X_level <- sp500d
  
  # Lunde Timmermann filter
  setpar_filtering_alg( tr_bear = -30, tr_bull = 30 )
  bull <- run_filtering_alg( index = X_level * 100)
  LTF <- as.timeSeries( matrix(bull, ncol = 1,
                               dimnames = list(rownames(X_level), 
                                               "LTF")) )
  
  
  # Save environment
  global_data <- new.env()
  global_data$X <- X
  global_data$X_level <- X_level
  global_data$BBS <- BBS$signal$insample
  global_data$BBSFuzzy <- BBSFuzzy$signal$insample
  global_data$LTF <- LTF
  
  
  
  saveRDS( global_data, file = paste0(wd, "Data/global_data.rds") )
  
  
  
  
  plot( global_data$LTF )
  
  plot( X_level )
  abline( v = time(global_data$LTF)[ which(global_data$LTF == 0) ], col = 2 )
  abline( v = time(global_data$LTF)[ which(global_data$LTF == 1) ], col = 3 )
  lines( X_level )
  
  
  
  
  input <- list( minphase = 5 * 4 * 2,
                 mincycle = 5 * 4 * 7 )
  
  
  
  
  
  
  
  
  
  
  