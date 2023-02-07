  
  
  ############################################################################
  ### GOOD AND BAD TURBULENCE - OUT-OF-SAMPLE - PREDICTING TURBULENCE
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     18.03.2021
  # First version:    18.03.2021
  # --------------------------------------------------------------------------
  
  
  # DESCRIPTION:
  #
  # We investigate the oos performance of predicting the Mahalanobis distance
  # with either garch or sv.
  
  
  
  
  # --------------------------------------------------------------------------
  # REQUIRE
  # --------------------------------------------------------------------------
  
  require(stochvol)
  require(RP)
  require(DAARC)
  
  wd <- "R:/Asset_Management/Research_Projects/Multi_Asset/DAA/DAARC/Sensor_Testing/Univariate/Good_And_Bad_Turbulence/"
  source( paste0(wd, "Source/daarc_functions.R") )
  
  
  
  
  
  # --------------------------------------------------------------------------
  # DATA
  # --------------------------------------------------------------------------
  
  MD <- loadSensor( sensor_name = "turbulence_base_dm" )
  X <- MD$data$X[isWeekday(time(MD$data$X)), ]
  X_bm <- MD$data$X_bm[isWeekday(time(MD$data$X_bm)), ]
  
  # Alternatively, use eqw as bm
  X_eqw <- apply( X, 1, mean )
  
  universe <- "dm"
  oos_dates <- rownames(X)[ -c(1:(252*3)) ]
  
  
  
  
  # --------------------------------------------------------------------------
  # BBQ
  # --------------------------------------------------------------------------
  
  BBS <- BBSRC$new()
  BBS$setCtrl( minphase = 21,
               mincycle = 21 * 3,
               theta = 0.1,
               logarithmic = FALSE,
               e = 0,
               k_peak = 10,
               k_trough = 10,
               l_peak = 10,
               l_trough = 10,
               language = "R",
               algo = "Bry-Boschan" )
  BBS$data <- list( X_level = cumulated(X_bm, "discrete") )
  
  minphase_vec <- round( 21 * c(0.5, 0.75, 1, 2, 3, 4, 5) )
  
  
  
  
  # --------------------------------------------------------------------------
  # IN-SAMPLE
  # --------------------------------------------------------------------------
  
  BBT <- BBTurbulence$new()
  BBT$setCtrl( universe = "dm",
               method = "base",
               verbose = TRUE )
  # BBT$spec$BBS <- BBS
  BBT$data <- list( X = X,
                    X_bm = X_bm )
  BBT$computeSignalInsample()
  
 
  
  
  
  # --------------------------------------------------------------------------
  # OUT-OF-SAMPLE
  # --------------------------------------------------------------------------
  
  BBTgarch <- BBTGARCH$new()
  BBTgarch$setCtrl( universe = "dm",
                    method = "base",
                    verbose = TRUE )
  BBTgarch$spec$BBS <- BBS
  BBTgarch$spec$garch_ctrl <- garchCtrl(steps = 21)
  BBTgarch$data <- list( X = X,
                         X_bm = X_bm )
  
  # Bry-Boschan in R
  k = 3
  BBTgarch$spec$BBS$spec$language <- "R"
  BBTgarch$spec$BBS$spec$minphase <- minphase_vec[k]
  BBTgarch$spec$BBS$spec$mincycle <- minphase_vec[k] * 3
  BBTgarch$spec$name <- paste0( "bbtgarchpred_bbr_", universe, "_phase=", 
                                BBTgarch$spec$BBS$spec$minphase, "_cycle=",
                                BBTgarch$spec$BBS$spec$mincycle, "_oos" )
  BBTgarch$computeSignalInsample()
  # debugonce( BBTgarch$computeSignalPrediction )
  BBTgarch$updateSignal()
  BBTgarch$save( path = paste0(wd, "waRehouse/Garch/") )
  
 
  
  # --------------------------------------------------------------------------
  # MEAN SQUARED ERROR
  # --------------------------------------------------------------------------
  
  idx_is <- which(names(BBTgarch$signal) == "insample")
  lSig <- lapply( BBTgarch$signal[ -idx_is ], FUN = function(x) { x[ ,"md_bull_scl"] } )
  signals <- do.call( cbind, lSig )
  signal_is <- BBT$signal$insample[ ,"bull_scl"]
  dates <- intersect( rownames(signals), rownames(signal_is) )
  signals <- signals[dates, ]
  signal_is <- signal_is[dates, ]
  
  signal_is_lag <- lapply( 1:ncol(signals), FUN = function(i) { lag(signal_is, -(i-1)) } )
  signal_is_lag <- do.call( cbind, signal_is_lag )
  
  se <- na.omit( (signals - signal_is_lag)^2 )
  mse <- apply( se, 2, mean )
  barplot( mse )
  
  plot( se, plot.type = "single" )
  
  
  
  
  # --------------------------------------------------------------------------
  # TRADING STRATEGIES
  # --------------------------------------------------------------------------
  
  tail(BBTgarch$signal$base)
  head(BBTgarch$signal$pred1)
  
  
  idx_is <- which(names(BBTgarch$signal) == "insample")
  lSig <- lapply( BBTgarch$signal[ -idx_is ], FUN = function(x) { x[ ,"signal"] } )
  signals <- do.call( cbind, lSig )
  signals <- na.omit( cbind( insample = BBTgarch$signal$insample[ ,"delta"], signals ) )
  
  tc <- 0.004
  n_lag <- 2
  penalty <- 0
  
  
  # Binary signals
  TD <- trainingData( Y_train = BBTgarch$data$X_bm,
                      X_train = (1 + sign( ema( signals, 0.1 ) )) / 2 )
  
  # debugonce(signalTesting.byTrading)
  test <- signalTesting.byTrading( X = TD$Y_train,
                                   sig = TD$X_train,
                                   n_lag = n_lag,
                                   tc = tc,
                                   penalty = penalty )
  colnames(test) <- colnames(TD$X_train)
  test_nc <- signalTesting.byTrading( X = TD$Y_train,
                                      sig = TD$X_train,
                                      n_lag = n_lag,
                                      tc = 0,
                                      penalty = penalty )
  colnames(test_nc) <- colnames(TD$X_train)
  
  X_tmp <- na.omit( cbind( bm = TD$Y_train, test ) )
  X_tmp_nc <- na.omit( cbind( bm = TD$Y_train, test_nc ) )
  
  lStats <- descStats( X_tmp )
  lStats_nc <- descStats( X_tmp_nc )
  
  stats_fields <- c("cumret", "sds", "maxDD")
  t(lStats_nc$stats[stats_fields, ])    
  t(lStats$stats[stats_fields, ]) 
  
  DD <- drawDownStats( X_tmp )
  DD
  
  
  colors <- c(1, fBasics::rainbowPalette(n = ncol(X_tmp)))
  plot( as.simTS(X_tmp[ ,1:10]), col = colors )
  plot( as.simTS( tail( X_tmp[ ,1:10], 500 ) ), col = colors )
  plot( as.simTS( tail( X_tmp, 500 ) ) )
  
  
  plot.timeSeries( log(cumulated(X_tmp, "discrete")),  col = colors )
  plot.timeSeries( log(cumulated(X_tmp_nc, "discrete")), col = colors )
  
  
  
  
  
  
  
  ####
  
  lSig <- lapply( BBTgarch$signal, FUN = function(x) { x[ ,"md_bull"] } )
  signals <- do.call( cbind, lSig )
 
  lSig <- lapply( BBTgarch$signal, FUN = function(x) { x[ ,"signal"] } )
  signals <- do.call( cbind, lSig )
  
  colors <- c(1, fBasics::divPalette(n = ncol(signals)-1, "RdYlGn"))
  plot( scale(signals, FALSE, TRUE), plot.type = "single", col = colors )
  
  
  
  tail(signals)
  
  # sig_name <- "signal"
  # tmp <- BBTgarch$signal$base[ ,sig_name]
  # for ( k in 2:length(BBTgarch$signal) ) {
  #   tmp <- cbind( tmp, lag( BBTgarch$signal[[k]][ ,sig_name], -(k-1) ) )
  # }
  # head(tmp)
  
  sig_name <- "signal"
  tmp <- do.call( cbind, lapply( BBTgarch$signal, FUN = function(x) { x[ ,sig_name] } ) )
  
  
  
  # plot( scale(tmp, FALSE, TRUE)[ ,c(1, 2)], plot.type = "single", col = colors, type = "o" )
  # 
  # head(cbind(pred1 = BBTgarch$signal$pred1[ ,sig_name],
  #            pred2 = BBTgarch$signal$pred2[ ,sig_name]))
  
  
  
  
  startdate <- "2010-02-01"
  # enddate <- "2020-03-01"
  # startdate <- "2020-03-30"
  enddate <- "2020-04-15"
  par( mfrow = c(2, 1) )
  plot( window(scale(tmp, FALSE, TRUE)[ ,c(1, 10, 20)], startdate, enddate), 
        plot.type = "single", col = colors, type = "o" )
  abline( h = 0 )
  plot( log(cumulated(window(X_bm, startdate, enddate), "discrete")), type = "o" )
  
  
  
  dev.off()
  is <- BBTgarch$signal$insample[ ,"delta"]
  colors <- c("blue", "black", fBasics::rainbowPalette(n = ncol(X_tmp)))
  plot( window(scale(cbind(is, tmp), FALSE, TRUE)[ ,c(1, 2, 3, 10, 20)], startdate, enddate), 
        plot.type = "single", col = colors, type = "o" )
  abline( h = 0 )
  
  
  
  
  
