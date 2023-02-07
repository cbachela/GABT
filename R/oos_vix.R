  
  
  ############################################################################
  ### GOOD AND BAD TURBULENCE - OUT-OF-SAMPLE ANALYSIS - VIX
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
  # Load VIX
  # --------------------------------------------------------------------------
  
  # VRP <- loadSensor( sensor_name = "vrp" )
  # VRP$update()
  # vix <- VRP$signal$base[ ,"iv"]
  
  X_sptr <-  rodbcGetOLZDBReturns( assetName = "SPTR",
                                   refCcy = "USD", 
                                   frqncy = "daily", 
                                   compounding = "discrete" )
  X_sptr <- X_sptr[isWeekday(time(X_sptr)), ]
  
  vix_ret <- rodbcGetOLZDBReturns( assetName = "VIX",
                                   refCcy = "USD", 
                                   frqncy = "daily", 
                                   compounding = "continuous" )
  
  vix_ret <- vix_ret[isWeekday(time(vix_ret)), ]
  vix_base <- 21.8999999999999  # this is the level of the vix as of 01.03.1990
  vix <- cumulated( vix_ret, compounding = "continuous" ) * vix_base
  
  
 
  
  # --------------------------------------------------------------------------
  # BBT base
  # --------------------------------------------------------------------------
  
  # bbt <- loadSensor( sensor_name = "bbturbulence_base_dm" )
  # bbt$data$X_bm <- returns( vix, "discrete" )
  # bbt$data$capw <- NULL
  # bbt$spec$name <- "bbturbulence_base_vix"
  # bbt$signal <- list()
  bbt <- loadSensor( sensor_name = "bbturbulence_base_vix" )
  bbt$updateData()
  bbt$data$X_bm <- returns( vix, "discrete" )
  bbt$computeSignalInsample()
  bbt$updateSignal()
  bbt$save()
  
  plot( bbt$signal$insample )
  tail( bbt$signal$base )
  
  
  
  # --------------------------------------------------------------------------
  # BBT scl2
  # --------------------------------------------------------------------------
  
  # bbt_scl2 <- loadSensor( sensor_name = "bbturbulence_scl2_dm" )
  # bbt_scl2$data$X_bm <- returns( vix, "discrete" )
  # bbt_scl2$spec$name <- "bbturbulence_scl2_vix"
  # bbt_scl2$signal <- list()
  bbt_scl2 <- loadSensor( sensor_name = "bbturbulence_scl2_vix" )
  
  bbt_scl2$updateData()
  bbt_scl2$data$X_bm <- returns( vix, "discrete" )
  bbt_scl2$computeSignalInsample()
  bbt_scl2$updateSignal()
  bbt_scl2$save()
  
  plot( bbt_scl2$signal$insample )
  tail( bbt_scl2$signal$scl2 )
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # BBTEWMA base
  # --------------------------------------------------------------------------
  
  # bbtewma <- loadSensor( sensor_name = "bbtewma_base_dm" )
  # bbtewma$data$X_bm <- returns( vix, "discrete" )
  # bbtewma$data$capw <- NULL
  # bbtewma$spec$name <- "bbtewma_base_vix"
  # bbtewma$signal <- list()
  # # debugonce( bbtewma$computeSignalInsample )
  # bbtewma$computeSignalInsample()
  # # debugonce( bbtewma$computeSignal )
  # bbtewma$updateSignal()
  # # bbtewma$save()
  
  bbtewma <- loadSensor( sensor_name = "bbtewma_base_vix" )
  bbtewma$updateData()
  bbtewma$data$X_bm <- returns( vix, "discrete" )
  bbtewma$computeSignalInsample()
  debugonce( bbtewma$computeSignal )
  bbtewma$updateSignal()
  bbtewma$save()
  
  
  
  # --------------------------------------------------------------------------
  # BBTEWMA scl2
  # --------------------------------------------------------------------------
  
  # bbtewma_scl2 <- loadSensor( sensor_name = "bbtewma_scl2_dm" )
  # bbtewma_scl2$data$X_bm <- returns( vix, "discrete" )
  # bbtewma_scl2$spec$name <- "bbtewma_scl2_vix"
  # bbtewma_scl2$signal <- list()
  # # debugonce( bbtewma_scl2$computeSignalInsample )
  # bbtewma_scl2$computeSignalInsample()
  # # debugonce( bbtewma_scl2$computeSignal )
  # bbtewma_scl2$updateSignal()
  # # bbtewma_scl2$save()
  
  bbtewma_scl2 <- loadSensor( sensor_name = "bbtewma_scl2_vix" )
  bbtewma_scl2$updateData()
  bbtewma_scl2$data$X_bm <- returns( vix, "discrete" )
  bbtewma_scl2$computeSignalInsample()
  bbtewma_scl2$updateSignal()
  bbtewma_scl2$save()
  
  
  minphase_vec <- floor( 21 * c(0.5, 0.75, 1, 2, 3, 4, 5) )
  for ( k in seq(along = minphase_vec) ) {

    # # Bry-Boschan in R
    # bbtewma_scl2$signal <- list()
    # bbtewma_scl2$data$BBS$spec$language <- "R"
    # bbtewma_scl2$data$BBS$spec$minphase <- minphase_vec[k]
    # bbtewma_scl2$data$BBS$spec$mincycle <- minphase_vec[k] * 3
    # bbtewma_scl2$data$name <- paste0( "bbtewma_scl2_vix_bbr_phase=",
    #                                   bbtewma_scl2$data$BBS$spec$minphase, "_cycle=",
    #                                   bbtewma_scl2$data$BBS$spec$mincycle )
    # bbtewma_scl2$computeSignalInsample()
    # bbtewma_scl2$updateSignal()
    # bbtewma_scl2$save()

    # Bry-Boschan in C
    bbtewma_scl2$signal <- list()
    bbtewma_scl2$data$BBS$spec$language <- "C"
    bbtewma_scl2$data$BBS$spec$minphase <- minphase_vec[k]
    bbtewma_scl2$data$BBS$spec$mincycle <- minphase_vec[k] * 3
    bbtewma_scl2$spec$name <- paste0( "bbtewma_scl2_vix_bbc_phase=",
                                      bbtewma_scl2$data$BBS$spec$minphase, "_cycle=",
                                      bbtewma_scl2$data$BBS$spec$mincycle )
    bbtewma_scl2$computeSignalInsample()
    bbtewma_scl2$updateSignal()
    bbtewma_scl2$save()

  }
  
  
  
  
  
  # --------------------------------------------------------------------------
  # BBTGARCH base
  # --------------------------------------------------------------------------
  
  # bbtgarch <- loadSensor( sensor_name = "bbtgarch_base_dm" )
  # bbtgarch$data$X_bm <- returns( vix, "discrete" )
  # bbtgarch$data$capw <- NULL
  # bbtgarch$spec$name <- "bbtgarch_base_vix"
  # bbtgarch$signal <- list()
  bbtgarch <- loadSensor( sensor_name = "bbtgarch_base_vix" )
  bbtgarch$updateData()
  bbtgarch$data$X_bm <- returns( vix, "discrete" )
  bbtgarch$computeSignalInsample()
  bbtgarch$updateSignal()
  bbtgarch$save()
  
  plot( bbtgarch$signal$insample )
  tail( bbtgarch$signal$base )
  
  
  
  
  # --------------------------------------------------------------------------
  # Compute signal on realized vola
  # --------------------------------------------------------------------------
  
  
  bbtewma <- loadSensor( sensor_name = "bbtewma_base_dm" )
  vola_roll <- applyRoll( Data = bbtewma$data$X,
                         Width = 21,
                         FUN = function(X) { apply(X, 2, sd) } )
  # var_roll <- getCondVar( garch( bbtewma$data$X ) )
  bbtewma$data$X <- var_roll
  bbtewma$data$capw <- NULL
  bbtewma$spec$name <- "bbtewma_base_vola_roll_21d"
  bbtewma$signal <- list()
  # debugonce( bbtewma$computeSignalInsample )
  bbtewma$computeSignalInsample()
  # debugonce( bbtewma$computeSignal )
  bbtewma$updateSignal()
  bbtewma$save()

  
  
  bbtewma_scl2 <- loadSensor( sensor_name = "bbtewma_scl2_dm" )
  bbtewma_scl2$data$X <- var_roll
  bbtewma_scl2$spec$name <- "bbtewma_scl2_vola_roll_21d"
  bbtewma_scl2$signal <- list()
  bbtewma_scl2$computeSignalInsample()
  bbtewma_scl2$updateSignal()
  bbtewma_scl2$save()
  
  
  
  
 
  sig <- na.omit( cbind( bbtewma$signal$insample[ ,"delta"],
                         bbtewma$getSignal(),
                         bbtewma_scl2$signal$insample[ ,"delta"],
                         bbtewma_scl2$getSignal() ) )
  sig_binary <- (1 + sign(ema(sig, 1))) / 2
  X_bm <- bbtewma$data$X_bm
  test <- signalTesting.byTrading( X = X_bm,
                                   sig = sig_binary,
                                   n_lag = 1, 
                                   tc = 0.004 )
  test_nc <- signalTesting.byTrading( X = X_bm,
                                      sig = sig_binary,
                                      n_lag = 1, 
                                      tc = 0 )
  X_tmp <- na.omit( cbind( X_bm, test ) )
  X_tmp_nc <- na.omit( cbind( X_bm, test_nc ) )
  # plot( cumulated(X_tmp, "discrete"), plot.type = "single" )
  
  stats_fields <- c("cumret", "means", "sds", "sharpe", "maxDD")
  
  plot( as.simTS(X_tmp), logscale = TRUE )
  t( descStats( X_tmp )$stats[stats_fields, ] )
  
  plot( as.simTS(X_tmp_nc), logscale = TRUE )
  t( descStats( X_tmp_nc )$stats[stats_fields, ] )
  
  plot( bbtewma$signal$insample )
 
  
  
  
  # --------------------------------------------------------------------------
  # Compute signal on ETF
  # --------------------------------------------------------------------------
  
  tickers <- c("LVOEUR EU", "VIXY US", "VIXM US")
  
  X <- rodbcGetOLZDBReturns( assetName = tickers,
                             refCcy = "Local",
                             frqncy = "daily" )
  X <- X[isWeekday(time(X)), tickers]
  
  
  plot( cumulated(na.omit(X[ ,1]), "discrete"))
  plot( as.simTS(X) )
  
  
  descStats( na.omit(X) )  
  
  
  
  # Compare VIX to VIX-ETF's
  start_date <- start(X)
  vix[start_date, ]
  
  fc <- 0.6
  fixcost <- (1 + fc)^(1 / 252) - 1
  plot( as.simTS( na.omit(cbind(returns(vix, "discrete") - fixcost, X)) ) )
  
  
  
  
  # Run backtest
  
  Dbbtewma_vixy <- loadSensor( sensor_name = "bbtewma_scl2_dm" )

  bbtewma_vixy$data$X_bm <- na.omit(X[ ,2])
  dates <- intersect( rownames(bbtewma_vixy$data$X_bm), rownames(bbtewma_vixy$data$X) )
  bbtewma_vixy$data$X <- bbtewma_vixy$data$X[dates, ]
  bbtewma_vixy$spec$name <- "bbtewma_scl2_vixy_us"
  bbtewma_vixy$signal <- list()
  # debugonce( bbtewma_vixy$computeSignalInsample )
  bbtewma_vixy$computeSignalInsample()
  # debugonce( bbtewma_vixy$computeSignal )
  bbtewma_vixy$updateSignal()
  # bbtewma_vixy$save()
  
  
  sig_is <- bbtewma_vixy$signal$insample[ ,"delta"]
  sig_oos <- bbtewma_vixy$getSignal()
  # sig <- sig_is * (-1)
  sig <- sig_is
  sig <- na.omit( cbind( is = sig_is, oos = sig_oos ) )
  sig_binary <- (1 + sign(ema(sig, 1) - 0.3)) / 2
  test_tf <- signalTesting.byTrading( X = bbtewma_vixy$data$X_bm,
                                      sig = sig_binary,
                                      n_lag = 2, 
                                      tc = 0,
                                      ticketfee = 100,
                                      nav = 10^7,
                                      fc = 0.01 )
  X_tmp <- na.omit( cbind( bbtewma_vixy$data$X_bm, test_tf ) )
  plot( as.simTS(X_tmp) )
  # plot( as.simTS(X_tmp), logscale = FALSE )
  
  
  descStats(X_tmp)
  
  
  plot( sig )  
  abline(h = 0)  

  
  X_us <- bbtewma_vixy$data$X[ ,"US"]
  Y <- apply( na.omit( cbind(US = X_us, Vola = test_tf[ ,"signal_oos"]) ), 1, mean )
  Y2 <- apply( na.omit( cbind(US = X_us * 2, Vola = test_tf[ ,"signal_oos"]) ), 1, mean )
  
  Y_tmp <- na.omit(cbind(US = X_us, Y = Y, Y2 = Y2, US_half = X_us / 2))
  plot( as.simTS(Y_tmp) )  
  descStats(Y_tmp)
  
  
  
  
  
  
 
 
  
  