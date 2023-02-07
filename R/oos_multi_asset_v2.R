  
  
  ############################################################################
  ### GOOD AND BAD TURBULENCE - OUT-OF-SAMPLE ANALYSIS - MULTI ASSET
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     01.04.2022
  # First version:    01.04.2022
  # --------------------------------------------------------------------------
  
  
  # --------------------------------------------------------------------------
  # REQUIRE
  # --------------------------------------------------------------------------
  
  require(RP)
  require(BSS)
  require(simolz)
  require(DAARC)
  
  wd <- "R:/Asset_Management/Research_Projects/Multi_Asset/DAA/DAARC/Sensor_Testing/Univariate/Good_And_Bad_Turbulence/"
  source( paste0(wd, "Source/source_functions.R") )
  
  
  
  # --------------------------------------------------------------------------
  # Load multi-asset data
  # --------------------------------------------------------------------------
  
  
  AC <- c("JMCXAGTR", "JMCXENTR", "JMCXIMTR", 
          "MXWOMC", "NDDLWI", 
          "NDLEEGF", "RNGL", "SBWBL", "SBWGL", "SPI", "SPTR", 
          "SX5T", "TPXDDVD", "TUKXG")
  ###
  AC <- AC[-c(1:3)]
  ###
  X_mac <- rodbcGetOLZDBReturns( assetName = AC,
                                 refCcy = "USD",
                                 frqncy = "daily",
                                 na.rm = "s" )
  X_mac <- X_mac[isWeekday(time(X_mac)), ]
  X_bm <- apply( X_mac, 1, mean) 
  
  plot( as.simTS(X_bm) )
  
  
  
  # bbtewma <- loadSensor( sensor_name = "bbtewma_base_dm" )
  # bbtewma$data$capw <- NULL
  # bbtewma$spec$width <- 1000
  # bbtewma$spec$name <- "bbtewma_base_ac"
  # bbtewma$signal <- list()
  bbtewma <- loadSensor( sensor_name = "bbtewma_base_ac" )
  bbtewma$data$X <- X_mac
  bbtewma$data$X_bm <- X_bm
  # debugonce( bbtewma$computeSignalInsample )
  bbtewma$computeSignalInsample()
  # debugonce( bbtewma$computeSignal )
  bbtewma$updateSignal()
  # bbtewma$save()
  
  
  sig_is <- bbtewma$signal$insample[ ,"delta"]
  sig_oos <- bbtewma$getSignal()
  sig <- na.omit( cbind( is = sig_is, 
                         oos = sig_oos ) )
  sig_binary <- (1 + sign(ema(sig, 1))) / 2
  test <- signalTesting.byTrading( X = bbtewma$data$X_bm,
                                   sig = sig_binary,
                                   n_lag = 1, 
                                   tc = 0.001 )
  test_nc <- signalTesting.byTrading( X = bbtewma$data$X_bm,
                                      sig = sig_binary,
                                      n_lag = 1, 
                                      tc = 0 )
  X_tmp <- na.omit(cbind(bbtewma$data$X_bm, test, test_nc))
  # X_tmp <- cbind( X_tmp, test_nc )
  plot( as.simTS(X_tmp) )
  
  stats_fields <- c("cumret", "sds", "sharpe", "maxDD")
  t( descStats( X_tmp )$stats[stats_fields, ] )
  
  
  
  
  plot( sig, plot.type = "single" )
  abline( h = 0 )
  tail( sig )
  
  
  