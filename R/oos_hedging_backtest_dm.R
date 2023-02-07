  
  
  ############################################################################
  ### GOOD AND BAD TURBULENCE - OUT-OF-SAMPLE HEDGIND STRATEGY - EQUITY WORLD
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     04.06.2021
  # First version:    04.06.2021
  # --------------------------------------------------------------------------
  
  
  # --------------------------------------------------------------------------
  # REQUIRE
  # --------------------------------------------------------------------------
  
  require(RP)
  require(DAARC)
  
  wd <- "R:/Asset_Management/Research_Projects/Multi_Asset/DAA/DAARC/Sensor_Testing/Univariate/Good_And_Bad_Turbulence/"
  source( paste0(wd, "Source/source_functions.R") )
  
  
  
  
  
  # --------------------------------------------------------------------------
  # SENSORS
  # --------------------------------------------------------------------------
  
  Turb <- loadSensor( sensor_name = "turbulence_scl2_dm" )
  BBT <- loadSensor( sensor_name = "bbt_base_dm" )
  BBTScl2 <- loadSensor( sensor_name = "bbt_scl2_dm" )
  BBTEWMA <- loadSensor( sensor_name = "bbtewma_base_dm" )
  BBTEWMAScl2 <- loadSensor( sensor_name = "bbtewma_scl2_dm" )
  BBTSV <- loadSensor( sensor_name = "bbtsv_base_dm" )
  BBTGarch <- loadSensor( sensor_name = "bbtgarch_base_dm" )
  BBTGarchScl2 <- loadSensor( sensor_name = "bbtgarch_scl2_dm" )
  
  
  # Turb$update()
  # BBT$update()
  # BBTScl2$update()
  # BBTEWMA$update()
  # BBTEWMAScl2$update()
  # BBTSV$update()
  # BBTGarch$update()
  
  
  X_bm <- Turb$data$X_bm
  
  
  sig_turb <- Turb$signal$scl2[ ,"md"]
  plot( sig_turb )
  
  sig_bbt <- BBT$signal$base[ ,"md_delta"]
  plot( sig_bbt )
  
  sig_bbtscl2 <- BBTScl2$signal$scl2[ ,"md_delta"]
  plot( sig_bbtscl2 )
  
  sig_bbtewma <- BBTEWMA$signal$base[ ,"signal"]
  plot( sig_bbtewma )
  
  sig_bbtewmascl2 <- BBTEWMAScl2$signal$scl2[ ,"signal"]
  plot( sig_bbtewmascl2 )
  
  sig_bbtsv <- BBTSV$signal$base[ ,"signal"]
  plot( sig_bbtsv )
  
  sig_bbtgarch <- BBTGarch$signal$base[ ,"signal"]
  plot( sig_bbtgarch )
  
  sig_bbtgarchscl2 <- BBTGarchScl2$signal$scl2[ ,"signal"]
  plot( sig_bbtgarchscl2 )
  
  
  
  
  
  # --------------------------------------------------------------------------
  # BACKTESTS
  # --------------------------------------------------------------------------
  
  delta <- ema( sig_bbtsv, alpha = 0.1 )
  
  plot( delta )
  
  sig <- 1 / (1 + exp(-delta * 1000))
  plot( sig )
  
  
  
  signals <- cbind( bbt = sig_bbt,
                    bbtscl2 = sig_bbtscl2,
                    ewma = sig_bbtewma,
                    ewmascl2 = sig_bbtewmascl2,
                    garch = sig_bbtgarch,
                    garchscl2 = sig_bbtgarchscl2,
                    sv = sig_bbtsv )
  signals <- na.omit( signals )
  
  delta <- ema( signals, alpha = 0.1 )
  delta <- ema( signals, alpha = 1 )
  sig <- 1 / (1 + exp(-delta * 1000))
  
  
  
  test <- signalTesting.byTrading( X = X_bm,
                                   sig = sig,
                                   n_lag = 1,
                                   tc = 0.004 )
  colnames(test) <- colnames(sig)
  test_nc <- signalTesting.byTrading( X = X_bm,
                                   sig = sig,
                                   n_lag = 1,
                                   tc = 0 )
  colnames(test_nc) <- paste0(colnames(sig), "_nc")
  X_tmp <- na.omit( cbind( bm = X_bm, test, test_nc ) )
  # X_tmp <- na.omit( cbind( bm = X_bm, test_nc ) )
  
  
  
  colors <- fBasics::divPalette( n = ncol(sig), "Spectral")
  colors <- c(1, colors, colors)
  plot( as.simTS(X_tmp), col = colors )
  lStats <- descStats(X_tmp)
  stats_fields <- c("cumret", "sds", "sharpe", "maxDD")
  lStats$desc
  t(lStats$stats[stats_fields, ])
  
  
  
  
  # --------------------------------------------------------------------------
  # 
  # --------------------------------------------------------------------------
  
  
  
  
  
  
  
  
  
  
  
  
  
  
    
  
  
  
  
  
  
