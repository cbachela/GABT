    
  
  ############################################################################
  ### GOOD AND BAD TURBULENCE - LONG VOLA
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     01.11.2022
  # First version:    01.11.2022
  # --------------------------------------------------------------------------
  
  
  # --------------------------------------------------------------------------
  # REQUIRE
  # --------------------------------------------------------------------------
  
  require(RP)
  require(garcholz)
  require(DAARC)
  
  wd <- "R:/Asset_Management/Research_Projects/Multi_Asset/DAA/DAARC/Sensor_Testing/Univariate/Good_And_Bad_Turbulence/"
  source( paste0(wd, "Source/source_functions.R") )
  
  
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # LOAD SENSOR
  # --------------------------------------------------------------------------

  
  bbtewma_base_sqret <- loadSensor( sensor_name = "bbtewma_base_dm_sqret" )
  signal <- bbtewma_base_sqret$getSignal()
  
  bbtewma_scl_sqret <- loadSensor( sensor_name = "bbtewma_scl_dm_sqret" )
  signal <- bbtewma_scl_sqret$getSignal()
  
  bbtewma_scl2_sqret <- loadSensor( sensor_name = "bbtewma_scl2_dm_sqret" )
  signal <- bbtewma_scl2_sqret$getSignal()
  
  
  # ------------------------
  # Signal testing
  
  y <- bbtewma_base_sqret$data$X_bm
  # sig <- (sign( signal ) + 1) / 2
  sig <- signal
  n_lag <- 2
  tc <- 0.0015
  test <- signalTesting.byTrading( X = y,
                                   sig = sig,
                                   n_lag = n_lag,
                                   tc = tc )
  colnames(test) <- colnames(sig)
  test_nc <- signalTesting.byTrading( X = y,
                                      sig = sig,
                                      n_lag = n_lag,
                                      tc = 0 )
  colnames(test_nc) <- paste0(colnames(sig), "_nc")
  
  X_tmp <- na.omit( cbind(y, test, test_nc) )
  
  
  plot( log(cumulated(X_tmp, "discrete")), plot.type = "single" )
  lines( signal, col = "blue" )
  lines( log( cumulated(X_tmp[ ,2], "discrete" ) ), col = 2 )
  lines( log( cumulated(X_tmp[ ,3], "discrete" ) ), col = 3 )
  
  
  t(descStats(X_tmp)$stats)
  
  
  
  
  
  
  # ------------------------
  # Invest in long vola
  
  VRPObj <- loadSensor( sensor_name = "vrp" )
  VRPObj$signal$base
  vix <- VRPObj$signal$base[ ,"iv"]
  vix_ret <- returns( vix, "discrete" )
  
  
  n_lag <- 2
  sig <- (sign( signal ) + 1) / 2
  # sig <- signal
  sig_longvola <- sig * (-1) + 1
  test <- signalTesting.byTrading( X = vix_ret,
                                   sig = sig_longvola,
                                   n_lag = n_lag,
                                   tc = tc )
  colnames(test) <- colnames(sig)
  test_nc <- signalTesting.byTrading( X = vix_ret,
                                      sig = sig_longvola,
                                      n_lag = n_lag,
                                      tc = 0 )
  colnames(test_nc) <- paste0(colnames(sig), "_nc")

  X_tmp <- na.omit( cbind( bm = y,
                           vix_ret = vix_ret,
                           bt = test,
                           bt_nc = test ) )

  plot( log(cumulated(X_tmp, "discrete")), plot.type = "single" )


  sim <- apply( X_tmp[ ,c("bm", "bt_nc")], 1, mean )
  sim_all <- na.omit( cbind( bm = X_tmp[ ,"bm"],
                             bm_5050 = X_tmp[ ,"bm"] / 2,
                             sim = sim ) )
  plot( as.simTS(sim_all) )
  descStats( sim_all )
  
  
  
  # Two-asset portfolio
  
  sig_longvola <- 2 / (1 + exp(signal) ) - 1
  # sig <- (sign( signal ) + 1) / 2
  # sig_longvola <- sig * (-1) + 1
  X <- na.omit( cbind( equity = y, vola = vix_ret ) )
  w_sig <- cbind( equity = 1 - abs(sig_longvola), vola = sig_longvola )
  wmat <- na.omit( lag( w_sig, 1 ) )
  colnames(wmat) <- colnames(X)
  wmat <- wmat[time(wmat) < end(X), ]  
  # debugonce( simPortfolioR )
  sim <- simPortfolio( X = X, wghts = wmat, fc = 0.01, vc = 0,
                       language = "R" )
  
  sim_all <- na.omit( cbind( bm = X_tmp[ ,"bm"], 
                             bm_5050 = X_tmp[ ,"bm"] / 2, 
                             sim = sim ) )
  plot( as.simTS(sim_all) ) 
  descStats( sim_all )
    
  plot( apply(wmat, 1, sum) )
  plot(wmat)
  
  window(cbind(X, sig_longvola), "2020-02-20", "2020-04-30")
  
  
  