  
  
  ############################################################################
  ### GOOD AND BAD TURBULENCE - IN SAMPLE ANALYSIS - BBQ SENSITIVITY TESTS
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     08.01.2021
  # First version:    10.10.2020
  # --------------------------------------------------------------------------
  
  
  # Description:
  #
  # We analyze the sensitivity of our relative turbulence measure, aka. pseudo
  # quadratic discrimininant analysis, to the parametrization of the BBQ algo.
  # We distinguish between the following model specifications:
  # Data weighting:
  # - base (no weighting)
  # - scl2 (proportional to cap-weighted^2)
  # Smoothing:
  # - ewma
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # REQUIRE
  # --------------------------------------------------------------------------
  
  require(stochvol)
  require(RP)
  require(DAARC)
  
  wd <- "R:/Asset_Management/Research_Projects/Multi_Asset/DAA/DAARC/Sensor_Testing/Univariate/Good_And_Bad_Turbulence/"
  source( paste0(wd, "Source/source_functions.R") )
  
  
  
  # --------------------------------------------------------------------------
  # BBQ PARAMETERS
  # --------------------------------------------------------------------------

  BBSObj <- DAARC::BBS$new()
  BBSObj$setCtrl( universe = "dm" )
  BBSObj$spec$k_peak <- BBSObj$spec$l_peak <- 
      BBSObj$spec$k_trough <- BBSObj$spec$l_trough <- 10
  BBSObj$spec$e <- 0
  BBSObj$spec$language <- "R"
  BBSObj$spec$theta <- 0.2
  
  minphase_vec <- floor( 21 * c(0.5, 0.75, 1, 1.5, 2, 3) )
  mincycle_vec <- minphase_vec * 3
  
  
  
  # --------------------------------------------------------------------------
  # INITIALIZE
  # --------------------------------------------------------------------------
  
  BBT <- loadSensor(sensor_name = "bbturbulence_base_dm" )
  BBTScl2 <- loadSensor(sensor_name = "bbturbulence_scl2_dm" )
  BBTewma <- loadSensor(sensor_name = "bbtewma_base_dm")
  BBTewmaScl2 <- loadSensor(sensor_name = "bbtewma_scl2_dm")
  
  BBT$update()
  BBTScl2$update()
  lData <- BBT$data
  wmat <- BBTScl2$data$wmat
  
  BBT0 <- BBTurbulence$new()
  BBT0$setCtrl( universe = "dm")
  BBT0$data <- lData
  
  BBTewma0 <- BBTEWMA$new()
  BBTewma0$setCtrl(universe = "dm", method = "base")
  BBTewma0$spec$ewma_alpha <- 0.1
  BBTewma0$data <- lData
  
  BBTsv <- BBTSV$new()
  BBTsv$setCtrl( universe = "dm", method = "base" )
  BBTsv$data <- lData
  
  BBTgarch <- BBTGARCH$new()
  BBTgarch$setCtrl( universe = "dm", method = "base" )
  BBTgarch$data <- lData
  
  
  
  # --------------------------------------------------------------------------
  # RUN LOOP
  # --------------------------------------------------------------------------
  
  for ( i in seq(along = minphase_vec) ) {
    
    Obj <- BBT0$copy()
    Obj$spec$method <- "base"
    Obj$spec$BBS <- BBSObj
    Obj$spec$BBS$spec$minphase <- minphase_vec[i]
    Obj$spec$BBS$spec$mincycle <- mincycle_vec[i]
    Obj$spec$name <- paste0("bbturbulence_base_dm_phase=",  
                            Obj$spec$BBS$spec$minphase, "_cycle=",
                            Obj$spec$BBS$spec$mincycle)
    Obj$computeSignalInsample()
    Obj$save( path = paste0(wd, "/waRehouse/BBQ_Sensitivity/") )
    
    Obj <- BBT0$copy()
    Obj$spec$method <- "scl2"
    Obj$spec$BBS <- BBSObj
    Obj$spec$BBS$spec$minphase <- minphase_vec[i]
    Obj$spec$BBS$spec$mincycle <- mincycle_vec[i]
    Obj$spec$name <- paste0("bbturbulence_scl2_dm_phase=",  
                            Obj$spec$BBS$spec$minphase, "_cycle=",
                            Obj$spec$BBS$spec$mincycle)
    Obj$data$wmat <- wmat
    Obj$computeSignalInsample()
    Obj$save( path = paste0(wd, "/waRehouse/BBQ_Sensitivity/") )
    
    
    # Minimum torsion ewma                 
    Obj <- BBTewma0$copy()
    Obj$spec$method <- "base"
    Obj$spec$BBS <- BBSObj
    Obj$spec$BBS$spec$minphase <- minphase_vec[i]
    Obj$spec$BBS$spec$mincycle <- mincycle_vec[i]
    Obj$spec$name <- paste0("bbtewma_base_dm_phase=",  
                            Obj$spec$BBS$spec$minphase, "_cycle=",
                            Obj$spec$BBS$spec$mincycle)
    Obj$computeSignalInsample()
    Obj$save( path = paste0(wd, "/waRehouse/BBQ_Sensitivity/") )
    
    Obj <- BBTewma0$copy()
    Obj$spec$method <- "scl2"
    Obj$spec$BBS <- BBSObj
    Obj$spec$BBS$spec$minphase <- minphase_vec[i]
    Obj$spec$BBS$spec$mincycle <- mincycle_vec[i]
    Obj$spec$name <- paste0("bbtewma_scl2_dm_phase=",  
                            Obj$spec$BBS$spec$minphase, "_cycle=",
                            Obj$spec$BBS$spec$mincycle)
    Obj$data$wmat <- wmat
    Obj$computeSignalInsample()
    Obj$save( path = paste0(wd, "/waRehouse/BBQ_Sensitivity/") )
    
    
    # Minimum torsion stochvol
    Obj <- BBTsv$copy()
    Obj$spec$method <- "base"
    Obj$spec$BBS <- BBSObj
    Obj$spec$BBS$spec$minphase <- minphase_vec[i]
    Obj$spec$BBS$spec$mincycle <- mincycle_vec[i]
    Obj$spec$name <- paste0("bbtsv_base_dm_phase=",  
                            Obj$spec$BBS$spec$minphase, "_cycle=",
                            Obj$spec$BBS$spec$mincycle)
    Obj$computeSignalInsample()
    Obj$save( path = paste0(wd, "/waRehouse/BBQ_Sensitivity/") )
    
    # Obj <- BBTsv$copy()
    # Obj$spec$method <- "scl2"
    # Obj$spec$BBS <- BBSObj
    # Obj$spec$BBS$spec$minphase <- minphase_vec[i]
    # Obj$spec$BBS$spec$mincycle <- mincycle_vec[i]
    # Obj$spec$name <- paste0("bbtsv_scl2_dm_phase=",  
    #                         Obj$spec$BBS$spec$minphase, "_cycle=",
    #                         Obj$spec$BBS$spec$mincycle)
    # Obj$data$wmat <- wmat
    # Obj$computeSignalInsample()
    # Obj$save( path = paste0(wd, "/waRehouse/BBQ_Sensitivity/") )
    
    # Minimum torsion garch
    Obj <- BBTgarch$copy()
    Obj$spec$method <- "base"
    Obj$spec$BBS <- BBSObj
    Obj$spec$BBS$spec$minphase <- minphase_vec[i]
    Obj$spec$BBS$spec$mincycle <- mincycle_vec[i]
    Obj$spec$name <- paste0("bbtgarch_base_dm_phase=",  
                            Obj$spec$BBS$spec$minphase, "_cycle=",
                            Obj$spec$BBS$spec$mincycle)
    Obj$computeSignalInsample()
    Obj$save( path = paste0(wd, "/waRehouse/BBQ_Sensitivity/") )
    
  }
  
  
  # Robust BBQ
  
  BBSF <- BBSFuzzy$new()
  BBSF$spec <- BBSObj$spec
    
  Obj <- BBT0$copy()
  Obj$spec$method <- "base"
  Obj$spec$BBS <- BBSF
  Obj$spec$BBS$spec$minphase_vec <- minphase_vec
  Obj$spec$BBS$spec$mincycle_vec <- mincycle_vec
  Obj$spec$name <- "bbturbulence_base_dm_robust"
  # debugonce( Obj$computeSignalInsample )
  Obj$computeSignalInsample()
  Obj$save( path = paste0(wd, "/waRehouse/BBQ_Sensitivity/") )
  
  Obj <- BBT0$copy()
  Obj$spec$method <- "scl2"
  Obj$spec$BBS <- BBSF
  Obj$spec$BBS$spec$minphase_vec <- minphase_vec
  Obj$spec$BBS$spec$mincycle_vec <- mincycle_vec
  Obj$spec$name <-"bbturbulence_scl2_dm_robust"
  Obj$data$wmat <- wmat
  Obj$computeSignalInsample()
  Obj$save( path = paste0(wd, "/waRehouse/BBQ_Sensitivity/") )
  
  # Minimum torsion ewma
  Obj <- BBTewma0$copy()
  Obj$spec$method <- "base"
  Obj$spec$BBS <- BBSF
  Obj$spec$BBS$spec$minphase_vec <- minphase_vec
  Obj$spec$BBS$spec$mincycle_vec <- mincycle_vec
  Obj$spec$name <- "bbtewma_base_dm_robust"
  Obj$computeSignalInsample()
  Obj$save( path = paste0(wd, "/waRehouse/BBQ_Sensitivity/") )
  
  Obj <- BBTewma0$copy()
  Obj$spec$method <- "scl2"
  Obj$spec$BBS <- BBSF
  Obj$spec$BBS$spec$minphase_vec <- minphase_vec
  Obj$spec$BBS$spec$mincycle_vec <- mincycle_vec
  Obj$spec$name <-"bbtewma_scl2_dm_robust"
  Obj$data$wmat <- wmat
  Obj$computeSignalInsample()
  Obj$save( path = paste0(wd, "/waRehouse/BBQ_Sensitivity/") )
  
  # Minimum torsion stochvol
  Obj <- BBTsv$copy()
  Obj$spec$method <- "base"
  Obj$spec$BBS <- BBSF
  Obj$spec$BBS$spec$minphase_vec <- minphase_vec
  Obj$spec$BBS$spec$mincycle_vec <- mincycle_vec
  Obj$spec$name <- "bbtsv_base_dm_robust"
  Obj$computeSignalInsample()
  Obj$save( path = paste0(wd, "/waRehouse/BBQ_Sensitivity/") )
  
  # Obj <- BBTsv$copy()
  # Obj$spec$method <- "scl2"
  # Obj$spec$BBS <- BBSF
  # Obj$spec$BBS$spec$minphase_vec <- minphase_vec
  # Obj$spec$BBS$spec$mincycle_vec <- mincycle_vec
  # Obj$spec$name <-"bbtsv_scl2_dm_robust"
  # Obj$data$wmat <- wmat
  # Obj$computeSignalInsample()
  # Obj$save( path = paste0(wd, "/waRehouse/BBQ_Sensitivity/") )
  
  # Minimum torsion garch
  Obj <- BBTgarch$copy()
  Obj$spec$method <- "base"
  Obj$spec$BBS <- BBSF
  Obj$spec$BBS$spec$minphase_vec <- minphase_vec
  Obj$spec$BBS$spec$mincycle_vec <- mincycle_vec
  Obj$spec$name <- "bbtgarch_base_dm_robust"
  Obj$computeSignalInsample()
  Obj$save( path = paste0(wd, "/waRehouse/BBQ_Sensitivity/") )
  
  
  
  # --------------------------------------------------------------------------
  # BACKTESTS
  # --------------------------------------------------------------------------

  filenames <- list.files(path = paste0(wd, "waRehouse/BBQ_Sensitivity/"),
                          pattern = ".rds")
  Names <- unlist( lapply(filenames, FUN = function(x) { substring(x, 1, nchar(x)-4) }) )
  Names <- Names[ -which(grepl("_ff10_", Names)) ]
  lSensor <- lapply( Names, FUN = function(name) { try(loadSensor(name, wd =  paste0(wd, "waRehouse/BBQ_Sensitivity/"))) } )
  lSig <- lapply( lSensor, FUN = function(x) { if ( !inherits(x, "try-error")) x$signal$insample[ ,"delta"] } )
  signals <- do.call( cbind, lSig )
  colnames(signals) <- Names
  
  tc <- 0.004
  n_lag <- 1
  penalty <- 0
  
  
  # Binary signals
  TD <- trainingData( Y_train = BBT0$data$X_bm,
                      X_train = (1 + sign(ema(signals, 0.1))) / 2 )
  
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
  # colnames(lStats$stats) <- colnames(lStats_nc$stats) <- c("bm", rep(paste0(minphase_vec, "_", mincycle_vec), 4) )
  # colnames(lStats$stats) <- colnames(lStats_nc$stats) <- c("bm", Names )
  colnames(lStats$stats) <- colnames(lStats_nc$stats) <- c("bm", rep(paste0(minphase_vec, "_", mincycle_vec), 4) )
  
  
  
  
  colors <- c(1, fBasics::rainbowPalette(n = ncol(X_tmp)-1))
  plot( as.simTS(X_tmp), col = colors )
  plot( as.simTS(X_tmp_nc), col = colors )
  
  
  stats_fields <- c("cumret", "sds", "sharpe", "maxDD")
  t(lStats$stats[stats_fields, ])
  t(lStats_nc$stats[stats_fields, ])
  
  # cbind( lStats$stats["cumret", ], lStats0$stats["cumret", ])
  # barplot( lStats$stats["cumret", ] - lStats0$stats["cumret", ])
  
  
  # Turnover
  turnover <- attr(test, "turnover")
  barplot( t(apply(aggYearly(turnover, compounding = "continuous"), 2, mean)) )
  
  
  plot( turnover[ ,ncol(turnover)] )
  plot( aggYearly(turnover[ ,c(1, ncol(turnover))], "continuous"), plot.type = "single" )
  
  
  
  plot.stats( lStats, sortby = NULL )
  
  
  n <- length(minphase_vec) + 1
  par( mfrow = c(2, 2) )
  barplot( lStats$stats["cumret", 2:(n+1)], main = "bbtewma", ylim = range(lStats$stats["cumret", ]))
  abline( h = lStats$stats["cumret", 1] )
  barplot( lStats$stats["cumret", (n+2):(2*n+1)], main = "bbtewma_scl2", ylim = range(lStats$stats["cumret", ]))
  abline( h = lStats$stats["cumret", 1] )
  barplot( lStats$stats["cumret", (2*n+2):(3*n+1)], main = "bbt", ylim = range(lStats$stats["cumret", ]))
  abline( h = lStats$stats["cumret", 1] )
  barplot( lStats$stats["cumret", (3*n+2):(4*n+1)], main = "bbt_scl2", ylim = range(lStats$stats["cumret", ]))
  abline( h = lStats$stats["cumret", 1] )
  
  
  par( mfrow = c(2, 2) )
  barplot( lStats_nc$stats["cumret", 2:(n+1)], main = "bbtewma", ylim = range(lStats_nc$stats["cumret", ]) )
  abline( h = lStats_nc$stats["cumret", 1] )
  barplot( lStats_nc$stats["cumret", (n+2):(2*n+1)], main = "bbtewma_scl2", ylim = range(lStats_nc$stats["cumret", ]))
  abline( h = lStats_nc$stats["cumret", 1] )
  barplot( lStats_nc$stats["cumret", (2*n+2):(3*n+1)], main = "bbt", ylim = range(lStats_nc$stats["cumret", ]))
  abline( h = lStats_nc$stats["cumret", 1] )
  barplot( lStats_nc$stats["cumret", (3*n+2):(4*n+1)], main = "bbt_scl2", ylim = range(lStats_nc$stats["cumret", ]))
  abline( h = lStats_nc$stats["cumret", 1] )
  
  
  dev.off()
  barplot( t(cbind(lStats$stats["cumret", ],
                  lStats_nc$stats["cumret", ])), beside = TRUE )
  
  
  
  
      
  
  
  
  
  
