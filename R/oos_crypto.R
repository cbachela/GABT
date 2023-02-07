  
  
  ############################################################################
  ### GOOD AND BAD TURBULENCE - OUT-OF-SAMPLE ANALYSIS - CRYPTO
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
  # Load data
  # --------------------------------------------------------------------------
  
  ticker <- c("XRPUSD",
              "XETUSD",
              "XLCUSD",
              "XBTUSD",
              "XBNUSD",
              "XTHUSD",
              "XEOUSD",
              "XMRUSD")
  X <- rodbcGetOLZDBReturns( assetName = ticker,
                             refCcy = "Local",
                             frqncy = "daily" )
  X <- na.omit( X[ ,ticker] )
  X_bm <- apply( X, 1, mean )
  
  plot( as.simTS(cbind(X, X_bm)) )
  plot( as.simTS(cbind(X, X_bm)), logscale = FALSE )
  
  descStats( cbind(X, X_bm) )
  
  
  
  
  # --------------------------------------------------------------------------
  # Compute signal (base)
  # --------------------------------------------------------------------------
  
  
  # bbtewma <- loadSensor( sensor_name = "bbtewma_base_dm" )
  # 
  # bbtewma$data$X <- X
  # bbtewma$data$X_bm <- X_bm
  # bbtewma$data$capw <- NULL
  # bbtewma$spec$name <- "bbtewma_base_crypto"
  # bbtewma$signal <- list()
  # # debugonce( bbtewma$computeSignalInsample )
  # bbtewma$computeSignalInsample()
  # # debugonce( bbtewma$computeSignal )
  # bbtewma$updateSignal()
  # # bbtewma$save()
  
  bbtewma <- loadSensor( sensor_name = "bbtewma_base_crypto" )
  bbtewma$data$X <- X
  bbtewma$data$X_bm <- X_bm
  # debugonce( bbtewma$computeSignalInsample )
  bbtewma$computeSignalInsample()
  # debugonce( bbtewma$computeSignal )
  bbtewma$updateSignal()
  # bbtewma$save()
  
  
  sig_is <- bbtewma$signal$insample[ ,"delta"]
  sig_oos <- bbtewma$getSignal()
  sig <- na.omit( cbind( is = sig_is, oos = sig_oos ) )
  sig_binary <- (1 + sign(ema(sig, 1) + 0.05) ) / 2
  
  test <- signalTesting.byTrading( X = bbtewma$data$X_bm,
                                   sig = sig_binary,
                                   n_lag = 1, 
                                   tc = 0.001 )
  X_tmp <- na.omit( cbind( bbtewma$data$X_bm, test ) )
  # X_tmp <- cbind( X_tmp, test_nc )
  plot( as.simTS(X_tmp), logscale = FALSE )
  plot( as.simTS(X_tmp) )
  
  
  stats_fields <- c("cumret", "sds", "sharpe", "maxDD")
  t( descStats( X_tmp )$stats[stats_fields, ] )
  
  
  
  plot( sig, plot.type = "single" )
  abline( h = 0 )
  tail( sig )
  
  
  
  
  
  