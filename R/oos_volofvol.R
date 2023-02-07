    
  
  ############################################################################
  ### GOOD AND BAD TURBULENCE - VOL OF VOL
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     18.10.2022
  # First version:    18.10.2022
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
  # INSTANTIATE SENSOR
  # --------------------------------------------------------------------------

  
  Object <- Object <- BBT$new()
  Object$setCtrl( method = "base",
                  universe = "dm" )
  Object$updateData()
  lData <- Object$data
  lData$X_bm <- lData$X_bm[isWeekday(time(lData$X_bm)), ]
  lData$X <- lData$X[isWeekday(time(lData$X)), ]
  
  # Create vola returns
  X_sq_roll <- applyRoll( Data = Object$data$X^2, 
                          Width = 21,
                          By = 1,
                          FUN = function(X) { ema(X, 0.1) } )
  # X_sq_roll <- ema( lData$X^2, 0.1 )
  X <- returns( X_sq_roll, compounding = "discrete" )
  X[abs(X) == Inf] <- 0
  plot( X[ ,1:10] )

  Object$data$X <- X
  Object$computeSignalInsample()
  
  
  plot( Object$signal$insample )
  
  
  
  
  
  # ------------------------
  # Signal testing
  
  y <- Object$data$X_bm
  signal <- Object$signal$insample[ ,"delta"] # * (-1)
  sig <- (sign( signal ) + 1) / 2
  n_lag <- 2
  test <- signalTesting.byTrading( X = y,
                                   sig = sig,
                                   n_lag = n_lag,
                                   tc = 0.004 )
  colnames(test) <- colnames(sig)
  test_nc <- signalTesting.byTrading( X = y,
                                      sig = sig,
                                      n_lag = n_lag,
                                      tc = 0 )
  colnames(test_nc) <- paste0(colnames(sig), "_nc")
  
  X_tmp <- na.omit( cbind(y, test, test_nc) )
  
  
  plot( log(cumulated(X_tmp, "discrete")), plot.type = "single" )
  abline(v = time(sig)[which(sig == 1)], col = 3)
  abline(v = time(sig)[which(sig == 0)], col = 2)
  lines(log(cumulated(X_tmp[ ,1], "discrete")), col = 1)
  lines(log(cumulated(X_tmp[ ,2], "discrete")), col = 4)
  lines(log(cumulated(X_tmp[ ,3], "discrete")), col = 5)
  
  
  descStats(X_tmp)
  
  
  
  
  
  # ------------------------
  # Invest in long vola
  
  VRPObj <- loadSensor( sensor_name = "vrp" )
  VRPObj$signal$base
  vix <- VRPObj$signal$base[ ,"iv"]
  vix_ret <- returns( vix, "discrete" )
  
    
  sig_longvola <- sig * (-1) + 1
  test <- signalTesting.byTrading( X = vix_ret,
                                   sig = sig_longvola,
                                   n_lag = n_lag,
                                   tc = 0.004 )
  colnames(test) <- colnames(sig)
  test_nc <- signalTesting.byTrading( X = vix_ret,
                                      sig = sig_longvola,
                                      n_lag = n_lag,
                                      tc = 0 )
  colnames(test_nc) <- paste0(colnames(sig), "_nc")
  
  X_tmp <- na.omit( cbind(bm = y, vix_ret, test, test_nc) )
    
  plot( log(cumulated(X_tmp, "discrete")), plot.type = "single" )
  
  
  
  sim <- apply( X_tmp[ ,c("bm", "delta_nc")], 1, mean )
  sim_all <- na.omit( cbind( X_tmp[ ,"bm"], X_tmp[ ,"bm"] / 2, sim ) )
  plot( as.simTS(sim_all) )  
      
  descStats( sim_all )
  
  
  
  
  
  
  
  