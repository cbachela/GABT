  
  
  ############################################################################
  ### GOOD AND BAD TURBULENCE - SENSORS
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     31.01.2021
  # First version:    10.10.2020
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
  # TURBULENCE
  # --------------------------------------------------------------------------
  
  # TurbMTEWMAScl2 <- TurbMTEWMA$copy()
  # TurbMTEWMAScl2$data <- TurbScl2$data
  # TurbMTEWMAScl2$signal <- list()
  # TurbMTEWMAScl2$computeSignalInsample()
  # TurbMTEWMAScl2$updateSignal()
  # TurbMTEWMAScl2$spec$name <- "turbulence_mtewma_scl2_dm"
  # TurbMTEWMAScl2$save()
  

  Turb <- loadSensor( sensor_name = "turbulence_base_dm", b_new = TRUE )
  # TurbScl <- loadSensor( sensor_name = "turbulence_scl_dm", b_new = TRUE )
  TurbScl2 <- loadSensor( sensor_name = "turbulence_scl2_dm", b_new = TRUE )
  TurbMTEWMA <- loadSensor( sensor_name = "turbulence_mtewma_dm", b_new = TRUE )
  TurbMTEWMAScl2 <- loadSensor( sensor_name = "turbulence_mtewma_scl2_dm", b_new = TRUE )
  
 
  Turb$update()
  Turb$computeSignalInsample()
 
  TurbScl2$update()
  TurbScl2$computeSignalInsample()
  
  TurbMTEWMA$update()
  TurbMTEWMA$computeSignalInsample()

  TurbMTEWMAScl2$update()
  TurbMTEWMAScl2$computeSignalInsample()
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # BEAR / BULL TURBULENCE
  # --------------------------------------------------------------------------
  
  # Base
  
  # BBTBase <- BBTurbulence$new()
  # BBTBase$setCtrl( method = "base", universe = "dm" )
  # BBTBase$updateData()
  
  BBTBase <- loadSensor( sensor_name = "bbturbulence_base_dm" )
  BBTBase$update()
  BBTBase$computeSignalInsample()
  BBTBase$save()
  
  
  # Using weighted return series (proportional to cap-weights)
  
  # BBTScl2 <- BBTurbulence$new()
  # BBTScl2$setCtrl( method = "scl2", universe = "dm" )
  # BBTScl2$updateData()
  
  BBTScl2 <- loadSensor( sensor_name = "bbturbulence_scl2_dm" )
  BBTScl2$update()
  BBTScl2$computeSignalInsample()
  BBTScl2$save()
  
  
  
  
  # --------------------------------------------------------------------------
  # MINIMUM TORSION EWMA BBT
  # --------------------------------------------------------------------------
  
  # Ewma
  # BBTEwma <- BBTEWMA$new()
  # BBTEwma$setCtrl( method = "base", universe = "dm" )
  # BBTEwma$updateData()
  
  BBTEwma <- loadSensor( sensor_name = "bbtewma_base_dm" )
  BBTEwma$update()
  BBTEwma$computeSignalInsample()
  BBTEwma$save()
  
  
  BBTEwma01 <- loadSensor( sensor_name = "bbtewma_base_dm" )
  BBTEwma01$signal <- list() 
  BBTEwma01$spec$ewma_alpha <- 0.01
  BBTEwma01$spec$name <- "bbtewma01_base_dm"
  BBTEwma01$update()
  BBTEwma01$computeSignalInsample()
  BBTEwma01$save()

    
  # Ewma Scl2 
  # BBTEwmaScl2 <- BBTEWMA$new()
  # BBTEwmaScl2$setCtrl( method = "scl2", universe = "dm" )
  # BBTEwmaScl2$updateData()
  
  BBTEwmaScl2 <- loadSensor( sensor_name = "bbtewma_scl2_dm" )
  BBTEwmaScl2$update()  
  BBTEwmaScl2$computeSignalInsample()
  BBTEwmaScl2$save()
  
  
  BBTEwmaScl201 <- loadSensor( sensor_name = "bbtewma_scl2_dm" )
  BBTEwmaScl201$signal <- list() 
  BBTEwmaScl201$spec$ewma_alpha <- 0.01
  BBTEwmaScl201$spec$name <- "bbtewma01_scl2_dm"
  BBTEwmaScl201$update()
  BBTEwmaScl201$computeSignalInsample()
  BBTEwmaScl201$save()
  
  
  
  
  
  
    
    
    
    
  
  
  # --------------------------------------------------------------------------
  # MINIMUM TORSION GARCH BBT
  # --------------------------------------------------------------------------

  # BBTGarch <- BBTGARCH$new()
  # BBTGarch$setCtrl( method = "base", universe = "dm" )
  # BBTGarch$updateData()
  
  BBTGarch <- loadSensor( sensor_name = "bbtgarch_base_dm" )
  # debugonce( BBTGarch$computeSignal )
  BBTGarch$update()
  BBTGarch$computeSignalInsample()
  BBTGarch$save()
  
  
  # BBTGarchScl2 <- BBTGARCH$new()
  # BBTGarchScl2$setCtrl( method = "scl2", universe = "dm" )
  # BBTGarchScl2$updateData()
  
  BBTGarchScl2 <- loadSensor( sensor_name = "bbtgarch_scl2_dm" )
  BBTGarchScl2$update()  
  BBTGarchScl2$computeSignalInsample()
  BBTGarchScl2$save()
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # MINIMUM TORSION STOCHASTIC VOLATILITY BBT
  # --------------------------------------------------------------------------
  
  # BBTsv <- BBTSV$new()
  # BBTsv$setCtrl( method = "base", universe = "dm" )
  # BBTsv$spec$width <- 504
  # BBTsv$updateData()
  
  BBTsv <- loadSensor( sensor_name = "bbtsv_base_dm" )
  # debugonce( BBTsv$computeSignal )
  BBTsv$update()
  BBTsv$computeSignalInsample()
  BBTsv$save()
  
  BBTsv$computeSignal
  
  
  plot( BBTsv$signal$base )
  
  
  # BBTsvScl2 <- BBTSV$new()
  # BBTsvScl2$setCtrl( method = "scl2", universe = "dm" )
  # BBTsvScl2$spec$width <- 504
  # BBTsvScl2$updateData()
  
  BBTsvScl2 <- loadSensor( sensor_name = "bbtsv_scl2_dm" )
  BBTsvScl2$update()  
  BBTsvScl2$computeSignalInsample()
  BBTsvScl2$save()
  
  
  plot( BBTsvScl2$signal$scl2, plot.type = "single" )
  
  
  
  
  
  
  
  
  
  
 
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # BINARY HEDGING RULES
  # --------------------------------------------------------------------------

  
  sig_base <- (pchisq( q = BBTBase$signal$base[ ,"md_all"], 
                       df = ncol(BBTBase$data$X) ) - 1) * -1
  sig_base_ema <- (pchisq( q = BBTBase$signal$base[ ,"md_all_ema"], 
                           df = ncol(BBTBase$data$X) ) - 1) * -1
  sig_ema_base <- ema(sig_base, 0.1)
  sig_base_delta <- (1 + sign(BBTBase$signal$base[ ,"signal"])) / 2
  
  
  plot(sig_base)  
  plot(sig_base_ema)
  lines(round(sig_base_ema), col = 2)
  
  
  sig_scl2 <- (pchisq( q = BBTScl2$signal$scl2[ ,"md_all"], 
                       df = ncol(BBTScl2$data$X) ) - 1) * -1
  sig_scl2_ema <- (pchisq( q = BBTScl2$signal$scl2[ ,"md_all_ema"], 
                           df = ncol(BBTScl2$data$X) ) - 1) * -1
  sig_ema_scl2 <- ema(sig_scl2, 0.1)
  sig_scl2_delta <- (1 + sign(BBTScl2$signal$scl2[ ,"signal"])) / 2
  
  
  
  
  sig_mtewma <- (pchisq( q = BBTEwma$signal$base[ ,"md_bear"], 
                       df = ncol(BBTScl2$data$X) ) - 1) * -1
  sig_mtewma2 <- sig_tmp <- BBTEwma$signal$base[ ,"md_bear"]
  sig_mtewma2[sig_tmp > ncol(BBTBase$data$X), ] <- 0
  sig_mtewma2[sig_tmp <= ncol(BBTBase$data$X), ] <- 1
  sig_mtewma_delta <- (1 + sign(BBTEwma$signal$base[ ,"signal"])) / 2
  sig_mtewma_scl2_delta <- (1 + sign(BBTEwmaScl2$signal$scl2[ ,"signal"])) / 2
  th <- 0.05
  sig_mtewma_delta_2 <- (1 + sign(BBTEwma$signal$base[ ,"signal"] + th)) / 2
  sig_mtewma_scl2_delta_2 <- (1 + sign(BBTEwmaScl2$signal$scl2[ ,"signal"] + th)) / 2
  sig_mtewma_delta_ema <- (1 + sign(ema(BBTEwma$signal$base[ ,"signal"], 0.1))) / 2
  sig_mtewma_scl2_delta_ema <- (1 + sign(ema(BBTEwmaScl2$signal$scl2[ ,"signal"], 0.1))) / 2
  
  
  
  
  
  # --------------------------------------------------------------------------
  # COMPARE RAW MAHALANOBIS DISTANCES
  # --------------------------------------------------------------------------
  
  
  

  md_base <- MDBase$signal$base[ ,"md"]
  md_base_ema <- ema( md_base, 0.1 )
  md_scl2 <- MDScl2$signal$scl2[ ,"md"]
  md_scl2_ema <- ema( md_scl2, 0.1 )
  md_base_bear <- BBTBase$signal$base[ ,"md_bear"]
  md_base_bear_ema <- ema( md_base_bear, 0.1 )
  md_scl2_bear <- BBTScl2$signal$scl2[ ,"md_bear"]
  md_scl2_bear_ema <- ema( md_scl2_bear, 0.1 )
  md_bbtewma_base_bear <- BBTEwma$signal$base[ ,"md_bear"]
  md_bbtewma_base_bear_ema <- ema( md_bbtewma_base_bear, 0.1 )
  md_bbtewma_scl2_bear <- BBTEwmaScl2$signal$scl2[ ,"md_bear"]
  md_bbtewma_scl2_bear_ema <- ema( md_bbtewma_scl2_bear, 0.1 )
  
  
  MD <- cbind( base = md_base,
               base_ema = md_base_ema,
               base_bear = md_base_bear,
               base_bear_ema = md_base_bear_ema,
               base_bbtewma_bear = md_bbtewma_base_bear,
               base_bbtewma_bear_ema = md_bbtewma_base_bear_ema,
               scl2 = md_scl2,
               scl2_ema = md_scl2_ema,
               scl2_bear = md_scl2_bear,
               scl2_bear_ema = md_scl2_bear_ema,
               scl2_bbtewma_bear = md_bbtewma_scl2_bear,
               scl2_bbtewma_bear_ema = md_bbtewma_scl2_bear_ema )
  
  plot( MD[ ,1:6] )
  plot( MD[ ,7:12] )
  
  
  
  
  # Out-of-sample k-means
  
  FUN <- function(X) 
  {
    lKM <- apply( X, 2, kmeans.adj, centers = 2 )
    states <- do.call( cbind, lapply( lKM, FUN = function(x) { timeSeries(x$cluster, rownames(X) ) } ) )
    colnames(states) <- colnames(X)
    return( states )
  }
  km_roll <- applyRoll( Data = na.omit(MD), Width = 0, Gap = 252, 
                        FUN = FUN, By = 1 )
  
  lKM <- apply( na.omit(MD), 2, kmeans.adj, centers = 4 )  
  states <- do.call( cbind, lapply( lKM, FUN = function(x) { timeSeries(x$cluster, rownames(na.omit(MD)) ) } ) )
  colnames(states) <- colnames(MD)
  
  plot( states[ ,1:6] )
  plot( states[ ,7:12] )
  
  
  
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # BACKTESTS
  # --------------------------------------------------------------------------
  
  stats_fields <- c("cumret", "sds", "sharpe", "maxDD")
  tc <- 0.004
  n_lag <- 1
  X_bm <- BBTBase$data$X_bm
  
  
  
  BBTEwma <- loadSensor( sensor_name = "bbtewma_base_dm" )
  BBTEwma05 <- loadSensor( sensor_name = "bbtewma05_base_dm" )
  BBTEwma01 <- loadSensor( sensor_name = "bbtewma01_base_dm" )
  BBTEwmaScl2 <- loadSensor( sensor_name = "bbtewma_scl2_dm" )
  BBTEwmaScl205 <- loadSensor( sensor_name = "bbtewma05_scl2_dm" )
  BBTEwmaScl201 <- loadSensor( sensor_name = "bbtewma01_scl2_dm" )
  
  # BBTEwma$update()
  # BBTEwma05$update()
  # BBTEwma01$update()
  # BBTEwmaScl2$update()
  # BBTEwmaScl205$update()
  # BBTEwmaScl201$update()
  
  
  signals <- cbind( base = BBTEwma$getSignal(),
                    base05 = BBTEwma05$getSignal(),
                    base01 = BBTEwma01$getSignal(),
                    scl2 = BBTEwmaScl2$getSignal(),
                    scl205 = BBTEwmaScl205$getSignal(),
                    scl201 = BBTEwmaScl201$getSignal() )
  # signals <- cbind( signals, avg = apply(signals, 1, mean) )
  sig <- (1 + sign(signals)) / 2
  X_bm <- BBTEwma$data$X_bm
  
  test <- signalTesting.byTrading( X = X_bm,
                                   sig = sig,
                                   n_lag = n_lag,
                                   tc = 0.004 )
  colnames(test) <- colnames(sig)
  test_nc <- signalTesting.byTrading( X = X_bm,
                                      sig = sig,
                                      n_lag = n_lag,
                                      tc = 0 )
  colnames(test_nc) <- colnames(sig)
  X_tmp <- na.omit( cbind(bm = X_bm, test) )
  X_tmp_nc <- na.omit( cbind(bm = X_bm, test_nc) )
  
  plot( as.simTS(X_tmp) )
  lStats <- descStats(X_tmp)
  t(lStats$stats[stats_fields, ])
  
  plot( as.simTS(X_tmp_nc) )
  lStats <- descStats(X_tmp_nc)
  t(lStats$stats[stats_fields, ])
  
  plot( as.simTS(tail(X_tmp_nc, 300)) )
  
  
  
  
  
  
  sig <- round(cbind(base = sig_base, 
                     base_ema = sig_base_ema,
                     ema_base = sig_ema_base,
                     base_delta = sig_base_delta))
  test <- signalTesting.byTrading( X = X_bm,
                                   sig = sig,
                                   n_lag = n_lag,
                                   tc = 0.004 )
  colnames(test) <- colnames(sig)
  test_nc <- signalTesting.byTrading( X = X_bm,
                                      sig = sig,
                                      n_lag = n_lag,
                                      tc = 0 )
  colnames(test_nc) <- colnames(sig)
  X_tmp <- na.omit( cbind(bm = X_bm, test, test_nc) )
  
  plot( as.simTS(X_tmp) )
  lStats <- descStats(X_tmp)
  t(lStats$stats[stats_fields, ])
  
  
  
  
  
  sig <- round(cbind(scl2 = sig_scl2, 
                     scl2_ema = sig_scl2_ema,
                     ema_scl2 = sig_ema_scl2,
                     scl2_delta = sig_scl2_delta))
  test <- signalTesting.byTrading( X = X_bm,
                                   sig = sig,
                                   n_lag = n_lag,
                                   tc = 0.004 )
  colnames(test) <- colnames(sig)
  test_nc <- signalTesting.byTrading( X = X_bm,
                                      sig = sig,
                                      n_lag = n_lag,
                                      tc = 0 )
  colnames(test_nc) <- colnames(sig)
  X_tmp <- na.omit( cbind(bm = X_bm, test, test_nc) )
  
  plot( as.simTS(X_tmp) )
  lStats <- descStats(X_tmp)
  t(lStats$stats[stats_fields, ])
  
  
  
  
  
  sig <- round(cbind(mtewma = sig_mtewma,
                     mtewma2 = sig_mtewma2,
                     mtewma_delta = sig_mtewma_delta,
                     mtewma_scl2_delta = sig_mtewma_scl2_delta,
                     mtewma_delta_ema = sig_mtewma_delta_ema,
                     mtewma_scl2_delta_ema = sig_mtewma_scl2_delta_ema ))
                     # mtewma_delta_2 = sig_mtewma_delta_2,
                     # mtewma_scl2_delta_2 = sig_mtewma_scl2_delta_2))
  
  test <- signalTesting.byTrading( X = X_bm,
                                   sig = sig,
                                   n_lag = n_lag,
                                   tc = 0.004 )
  colnames(test) <- colnames(sig)
  test_nc <- signalTesting.byTrading( X = X_bm,
                                      sig = sig,
                                      n_lag = n_lag,
                                      tc = 0 )
  colnames(test_nc) <- colnames(sig)
  X_tmp <- na.omit( cbind(bm = X_bm, test) )
  X_tmp <- na.omit( cbind(bm = X_bm, test_nc) )
  
  
  plot( as.simTS(X_tmp) )
  lStats <- descStats(X_tmp)
  t(lStats$stats[stats_fields, ])
  
  
  plot( as.simTS(tail(X_tmp, 300) ) )
  
  
  
  
  # Combining absolute and relative turbulence
  # Heuristic: use md for out-signal and bbt for in-signal
  # If current weight = 1 and d_t \leq \gamma --> out
  # if current weight = 0 and \Delta > 0 --> in
  
  
  
  
  

  test <- signalTesting.byTrading( X = X_bm,
                                   sig = sig,
                                   n_lag = n_lag,
                                   tc = 0.004 )
  colnames(test) <- colnames(sig)
  test_nc <- signalTesting.byTrading( X = X_bm,
                                      sig = sig,
                                      n_lag = n_lag,
                                      tc = 0 )
  colnames(test_nc) <- colnames(sig)
  X_tmp <- na.omit( cbind(bm = X_bm, test, test_nc) )
  
  plot( as.simTS(X_tmp) )
  lStats <- descStats(X_tmp)
  t(lStats$stats[stats_fields, ])
  
  
  
  
  
  