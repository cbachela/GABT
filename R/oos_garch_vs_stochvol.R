  
  
  ############################################################################
  ### GOOD AND BAD TURBULENCE - OUT-OF-SAMPLE - GARCH VS STOCHVOL VS EWMA
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     18.03.2021
  # First version:    18.03.2021
  # --------------------------------------------------------------------------
  
  
  # DESCRIPTION:
  #
  # We investigate the oos performance of the different model parametrizations
  # during the Covid-19 crash.
  
  
  
  
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
  
  universe <- "dm"
  
  
  
  # --------------------------------------------------------------------------
  # DATES
  # --------------------------------------------------------------------------
  
  oos_dates_1 <- rownames( window(X_bm, "2011-06-01", "2012-03-31") )
  oos_dates_2 <- rownames( window(X_bm, "2015-03-01", "2016-12-30") )
  oos_dates_3 <- rownames( window(X_bm, "2018-09-01", "2019-03-30") )
  oos_dates_4 <- rownames( window(X_bm, "2020-02-01", "2020-06-30") )
  
  # oos_dates <- c(oos_dates_1[1:2], oos_dates_2[1:2], oos_dates_3[1:2], oos_dates_4[1:2])
  # oos_dates <- c(oos_dates_1, oos_dates_2, oos_dates_3, oos_dates_4)
  # oos_dates <- oos_dates_4
  oos_dates <- rownames(X_bm)[ -c(1:(252*3)) ]
  
  
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
  
  # minphase_vec <- round( 21 * c(0.5, 0.75, 1, 2, 3, 4, 5) )
  minphase_vec <- floor(21 * c(0.5, 0.75, 1, 1.5, 2))
  
  
  
  # --------------------------------------------------------------------------
  # EWMA
  # --------------------------------------------------------------------------
  
  BBTewma <- BBTEWMA$new()
  BBTewma$setCtrl( universe = "dm",
                   method = "base",
                   verbose = TRUE )
  BBTewma$spec$ewma_alpha <- 0.1
  BBTewma$spec$BBS <- BBS
  BBTewma$data <- list( X = X,
                         X_bm = X_bm )
  
  # for ( k in seq(along = minphase_vec) ) {
  #   
  #   # Bry-Boschan in R
  #   BBTewma$spec$BBS$spec$language <- "R"
  #   BBTewma$spec$BBS$spec$minphase <- minphase_vec[k]
  #   BBTewma$spec$BBS$spec$mincycle <- minphase_vec[k] * 3
  #   BBTewma$spec$name <- paste0( "bbtewma_bbr_", universe, "_phase=", 
  #                                 BBTewma$spec$BBS$spec$minphase, "_cycle=",
  #                                 BBTewma$spec$BBS$spec$mincycle, "_oos" )
  #   # debugonce( BBTewma$computeSignal )
  #   BBTewma$computeSignal( dates = oos_dates )
  #   BBTewma$save( path = paste0( wd, "waRehouse/Garch_vs_Stochvol/OOS/" ) )
  #   
  #   # Bry-Boschan in C
  #   BBTewma$spec$BBS$spec$language <- "C"
  #   BBTewma$spec$BBS$spec$minphase <- minphase_vec[k]
  #   BBTewma$spec$BBS$spec$mincycle <- minphase_vec[k] * 3
  #   BBTewma$spec$name <- paste0( "bbtewma_bbc_", universe, "_phase=", 
  #                                 BBTewma$spec$BBS$spec$minphase, "_cycle=",
  #                                 BBTewma$spec$BBS$spec$mincycle, "_oos" )
  #   BBTewma$computeSignal( dates = oos_dates )
  #   BBTewma$save( path = paste0( wd, "waRehouse/Garch_vs_Stochvol/OOS/") )
  #   
  # }
  
  # Robust Bry-Boschan in R
  BBTewma$spec$BBS$spec$language <- "R"
  BBTewma$spec$robust <- TRUE
  BBTewma$spec$minphase_vec <- minphase_vec
  BBTewma$spec$mincycle_vec <- minphase_vec * 3
  BBTewma$spec$name <- paste0( "bbtewma_bbr_robust_", universe, "_oos" )
  BBTewma$computeSignal( dates = oos_dates )
  BBTewma$save( path = paste0( wd, "waRehouse/Garch_vs_Stochvol/OOS/") )
  
  # Robust Bry-Boschan in C
  BBTewma$spec$BBS$spec$language <- "C"
  BBTewma$spec$robust <- TRUE
  BBTewma$spec$minphase_vec <- minphase_vec
  BBTewma$spec$mincycle_vec <- minphase_vec * 3
  BBTewma$spec$name <- paste0( "bbtewma_bbc_robust_", universe, "_oos" )
  BBTewma$computeSignal( dates = oos_dates )
  BBTewma$save( path = paste0( wd, "waRehouse/Garch_vs_Stochvol/OOS") )
  
  # Lunde-Timmermann
  BBTewma$spec$BBS$spec$language <- "C"
  BBTewma$spec$BBS$spec$algo <- "Lunde-Timmermann"
  BBTewma$spec$name <- paste0( "bbtewma_ltc_", universe, "_theta=", 
                                BBTewma$spec$BBS$spec$theta, "_oos" )
  BBTewma$computeSignal( dates = oos_dates )
  BBTewma$save( path = paste0( wd, "waRehouse/Garch_vs_Stochvol/OOS/") )
  
  
  
  
  # --------------------------------------------------------------------------
  # GARCH
  # --------------------------------------------------------------------------
  
  BBTgarch <- BBTGARCH$new()
  BBTgarch$setCtrl( universe = "dm",
                    method = "base",
                    verbose = TRUE )
  BBTgarch$spec$BBS <- BBS
  BBTgarch$data <- list( X = X,
                         X_bm = X_bm )
  
  # for ( k in seq(along = minphase_vec) ) {
  #   
  #   # Bry-Boschan in R
  #   BBTgarch$spec$BBS$spec$language <- "R"
  #   BBTgarch$spec$BBS$spec$minphase <- minphase_vec[k]
  #   BBTgarch$spec$BBS$spec$mincycle <- minphase_vec[k] * 3
  #   BBTgarch$spec$name <- paste0( "bbtgarch_bbr_", universe, "_phase=", 
  #                                 BBTgarch$spec$BBS$spec$minphase, "_cycle=",
  #                                 BBTgarch$spec$BBS$spec$mincycle, "_oos" )
  #   # debugonce( BBTgarch$computeSignal )
  #   BBTgarch$computeSignal( dates = oos_dates )
  #   BBTgarch$save( path = paste0( wd, "waRehouse/Garch_vs_Stochvol/OOS/" ) )
  #   
  #   # Bry-Boschan in C
  #   BBTgarch$spec$BBS$spec$language <- "C"
  #   BBTgarch$spec$BBS$spec$minphase <- minphase_vec[k]
  #   BBTgarch$spec$BBS$spec$mincycle <- minphase_vec[k] * 3
  #   BBTgarch$spec$name <- paste0( "bbtgarch_bbc_", universe, "_phase=", 
  #                                 BBTgarch$spec$BBS$spec$minphase, "_cycle=",
  #                                 BBTgarch$spec$BBS$spec$mincycle, "_oos" )
  #   BBTgarch$computeSignal( dates = oos_dates )
  #   BBTgarch$save( path = paste0( wd, "waRehouse/Garch_vs_Stochvol/OOS/") )
  #   
  # }
  
  # Robust
  # Bry-Boschan in R
  BBTgarch$spec$BBS$spec$language <- "R"
  BBTgarch$spec$robust <- TRUE
  BBTgarch$spec$minphase_vec <- minphase_vec
  BBTgarch$spec$mincycle_vec <- minphase_vec * 3
  BBTgarch$spec$name <- paste0( "bbtgarch_bbr_robust_", universe, "_oos" )
  BBTgarch$computeSignal( dates = oos_dates )
  BBTgarch$save( path = paste0( wd, "waRehouse/Garch_vs_Stochvol/OOS/") )
  
  # Bry-Boschan in C
  BBTgarch$spec$BBS$spec$language <- "C"
  BBTgarch$spec$robust <- TRUE
  BBTgarch$spec$minphase_vec <- minphase_vec
  BBTgarch$spec$mincycle_vec <- minphase_vec * 3
  BBTgarch$spec$name <- paste0( "bbtgarch_bbc_robust_", universe, "_oos" )
  BBTgarch$computeSignal( dates = oos_dates )
  BBTgarch$save( path = paste0( wd, "waRehouse/Garch_vs_Stochvol/OOS/") )
  
  
  # Lunde-Timmermann
  BBTgarch$spec$BBS$spec$language <- "C"
  BBTgarch$spec$BBS$spec$algo <- "Lunde-Timmermann"
  BBTgarch$spec$name <- paste0( "bbtgarch_ltc_", universe, "_theta=", 
                                BBTgarch$spec$BBS$spec$theta, "_oos" )
  BBTgarch$computeSignal( dates = oos_dates )
  BBTgarch$save( path = paste0( wd, "waRehouse/Garch_vs_Stochvol/OOS/") )
  
  
  
  
  # --------------------------------------------------------------------------
  # STOCHVOL
  # --------------------------------------------------------------------------
  
  BBTsv <- BBTSV$new()
  BBTsv$setCtrl( universe = "dm",
                 method = "base",
                 verbose = TRUE )
  BBTsv$spec$BBS <- BBS
  BBTsv$data <- list( X = X,
                      X_bm = X_bm )
  
  # for ( k in seq(along = minphase_vec) ) {
  #   
  #   # Bry-Boschan in R
  #   BBTsv$spec$BBS$spec$language <- "R"
  #   BBTsv$spec$BBS$spec$minphase <- minphase_vec[k]
  #   BBTsv$spec$BBS$spec$mincycle <- minphase_vec[k] * 3
  #   BBTsv$spec$name <- paste0( "bbtsv_bbr_", universe, "_phase=", 
  #                                BBTsv$spec$BBS$spec$minphase, "_cycle=",
  #                                BBTsv$spec$BBS$spec$mincycle, "_oos" )
  #   # debugonce( BBTsv$computeSignal )
  #   BBTsv$computeSignal( dates = oos_dates )
  #   BBTsv$save( path = paste0( wd, "waRehouse/Garch_vs_Stochvol/OOS/" ) )
  #   
  #   # Bry-Boschan in C
  #   BBTsv$spec$BBS$spec$language <- "C"
  #   BBTsv$spec$BBS$spec$minphase <- minphase_vec[k]
  #   BBTsv$spec$BBS$spec$mincycle <- minphase_vec[k] * 3
  #   BBTsv$spec$name <- paste0( "bbtsv_bbc_", universe, "_phase=", 
  #                                BBTsv$spec$BBS$spec$minphase, "_cycle=",
  #                                BBTsv$spec$BBS$spec$mincycle, "_oos" )
  #   BBTsv$computeSignal( dates = oos_dates )
  #   BBTsv$save( path = paste0( wd, "waRehouse/Garch_vs_Stochvol/OOS/") )
  #   
  # }
  
  # Robust
  # Bry-Boschan in R
  BBTsv$spec$BBS$spec$language <- "R"
  BBTsv$spec$robust <- TRUE
  BBTsv$spec$minphase_vec <- minphase_vec
  BBTsv$spec$mincycle_vec <- minphase_vec * 3
  BBTsv$spec$name <- paste0( "bbtsv_bbr_robust_", universe, "_oos" )
  BBTsv$computeSignal( dates = oos_dates )
  BBTsv$save( path = paste0( wd, "waRehouse/Garch_vs_Stochvol/OOS/") )
  
  # Bry-Boschan in C
  BBTsv$spec$BBS$spec$language <- "C"
  BBTsv$spec$robust <- TRUE
  BBTsv$spec$minphase_vec <- minphase_vec
  BBTsv$spec$mincycle_vec <- minphase_vec * 3
  BBTsv$spec$name <- paste0( "bbtsv_bbc_robust_", universe, "_oos" )
  BBTsv$computeSignal( dates = oos_dates )
  BBTsv$save( path = paste0( wd, "waRehouse/Garch_vs_Stochvol/OOS/") )
  
  # Lunde-Timmermann
  BBTsv$spec$BBS$spec$language <- "C"
  BBTsv$spec$BBS$spec$algo <- "Lunde-Timmermann"
  BBTsv$spec$name <- paste0( "bbtsv_ltc_", universe, "_theta=", 
                               BBTsv$spec$BBS$spec$theta, "_oos" )
  BBTsv$computeSignal( dates = oos_dates )
  BBTsv$save( path = paste0( wd, "waRehouse/Garch_vs_Stochvol/OOS/") )
  
  
  
  
  # --------------------------------------------------------------------------
  # ANALYZE
  # --------------------------------------------------------------------------
  
  analyze <- function()
  {
    
    require(stochvol)
    require(DAARC)
    wd <- "R:/Asset_Management/Research_Projects/Multi_Asset/DAA/DAARC/Sensor_Testing/Univariate/Good_And_Bad_Turbulence/"
    source( paste0(wd, "Source/daarc_functions.R") )
    
    
    filenames <- list.files(path = paste0(wd, "waRehouse/Garch_vs_Stochvol/OOS/"),
                            pattern = ".rds")
    strategy_names <- unlist( lapply(filenames, FUN = function(x) { substring(x, 1, nchar(x)-4) }) )
    lSensor <- lapply( strategy_names, FUN = function(name) { try(loadSensor(name, wd =  paste0(wd, "waRehouse/Garch_vs_Stochvol/OOS/"))) } )
    names(lSensor) <- strategy_names
    
    
    # EWMA
    names_ewma <- strategy_names[ grepl( "bbtewma_", strategy_names ) ]
    # ~~~~~~~~~~~~~~~~ -15
    lSig_ewma <- lapply( lSensor[ names_ewma[-15] ], FUN = function(x) { if ( !inherits(x, "try-error")) x$getSignal() } )
    signals_ewma <- do.call( cbind, lSig_ewma )
    colnames(signals_ewma) <- names(lSig_ewma)
  
    
    # GARCH
    names_garch <- strategy_names[ grepl( "bbtgarch_", strategy_names ) ]
    lSig_garch <- lapply( lSensor[ names_garch ], FUN = function(x) { if ( !inherits(x, "try-error")) x$getSignal() } )
    signals_garch <- do.call( cbind, lSig_garch )
  
    
    # Stochvol
    names_sv <- strategy_names[ grepl( "bbtsv_", strategy_names ) ]
    lSig_sv <- lapply( lSensor[ names_sv ], FUN = function(x) { if ( !inherits(x, "try-error")) x$getSignal() } )
    lSig_sv_prob <- lapply( lSensor[ names_sv ], FUN = function(x) { if ( !inherits(x, "try-error")) x$signal$base[ ,"prob"] } )
    names(lSig_sv_prob) <- paste0( names(lSig_sv), "_prob" )
    signals_sv_delta <- do.call( cbind, lSig_sv )
    signals_sv_prob <- do.call( cbind, lSig_sv_prob ) - 0.5   #//
    signals_sv <- cbind( signals_sv_delta, signals_sv_prob )
    
    
    Obj <- loadSensor( sensor_name = "bbtsv_bbc_robust_dm_oos",
                       b_new = FALSE,
                       wd = paste0(wd, "waRehouse/Garch_vs_Stochvol/OOS/") )
    
    edit( Obj$computeSignal )
    plot( Obj$signal$base )
    
    debugonce( Obj$computeSignal )
    Obj$updateSignal()
    
                         
    
    
    signals <- na.omit( cbind( signals_ewma, signals_garch, signals_sv ) )
    # signals <- na.omit( signals_ewma )
    # signals <- timeSeries( t( na.omit( t( signals_ewma ) ) ), time(signals_ewma) )
    
    tc <- 0.004
    n_lag <- 2
    penalty <- 0
    
    
    # Binary signals
    TD <- trainingData( Y_train = lSensor[[1]]$data$X_bm,
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
    
    
    colors <- c(1, fBasics::rainbowPalette(n = ncol(X_tmp)))
    plot( as.simTS(X_tmp), col = colors )
    plot( as.simTS( tail( X_tmp[ ,1:10], 500 ) ), col = colors )
    plot( as.simTS( tail( X_tmp, 500 ) ) )
    
    
    plot.timeSeries( log(cumulated(X_tmp, "discrete")),  col = colors )
    plot.timeSeries( log(cumulated(X_tmp_nc, "discrete")), col = colors )
    
    
    
    idx_ewma <- grepl( "bbtewma_", colnames(lStats$stats) )
    idx_garch <- grepl( "bbtgarch_", colnames(lStats$stats) )
    idx_sv <- grepl( "bbtsv_", colnames(lStats$stats) )
    
    field_name <- "cumret"
    boxplot( cbind(ewma = lStats$stats[field_name, idx_ewma],
                   garch = lStats$stats[field_name, idx_garch],
                   stochvol = lStats$stats[field_name, idx_sv]) )
    
    
    
    
    
    s1 <- lSensor[["bbtewma_bbc_dm_phase=10_cycle=30_oos"]]$signal$base
    s2 <- lSensor[["bbtewma_bbc_dm_phase=42_cycle=126_oos"]]$signal$base
    s3 <- lSensor[["bbtewma_bbc_dm_phase=84_cycle=252_oos"]]$signal$base
    
    # name <- "states"
    name <- "md_bear"
    plot( tail( cbind( s1[ ,name], s2[ ,name], s3[ ,name] ), 500 ),
          plot.type = "single" )
    abline( h = 0 )
    
    
    
    s1 <- lSensor[["bbtewma_bbc_dm_phase=10_cycle=30_oos"]]$data$BBS$output$states
    s2 <- lSensor[["bbtewma_bbc_dm_phase=42_cycle=126_oos"]]$data$BBS$output$states
    s3 <- lSensor[["bbtewma_bbc_dm_phase=84_cycle=252_oos"]]$data$BBS$output$states
    
    plot( tail(cbind(s1, s2, s3), 500) )
  
      
    
  }
  
  
  # --------------------------------------------------------------------------
  xxx <- function()
  {
    
    require(stochvol)
    require(DAARC)
    wd <- "R:/Asset_Management/Research_Projects/Multi_Asset/DAA/DAARC/Sensor_Testing/Univariate/Good_And_Bad_Turbulence/"
    source( paste0(wd, "Source/daarc_functions.R") )
    
    
    filenames <- list.files(path = paste0(wd, "waRehouse/Garch_vs_Stochvol/OOS/"),
                            pattern = ".rds")
    strategy_names <- unlist( lapply(filenames, FUN = function(x) { substring(x, 1, nchar(x)-4) }) )
    lSensor <- lapply( strategy_names, FUN = function(name) { try(loadSensor(name, wd =  paste0(wd, "waRehouse/Garch_vs_Stochvol/OOS/"))) } )
    names(lSensor) <- strategy_names
    
    Name <- "bbtsv_bbc_robust_dm_oos"
    Obj <- lSensor[[ Name ]]
    
   plot(Obj$signal$base)    
    
    
    
  }
  
  
  
  