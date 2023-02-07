  
  
  ############################################################################
  ### GOOD AND BAD TURBULENCE - IN SAMPLE ANALYSIS - BBQ SENSITIVITY TESTS
  ### FAMA FRENCH INDUSTRY DATASET
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     13.03.2021
  # First version:    10.10.2020
  # --------------------------------------------------------------------------
  
  
  # Description:
  #
  # We analyze the sensitivity of our relative turbulence measure, aka. pseudo
  # quadratic discrimininant analysis, to the parametrization of the BBQ algo.
 
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # REQUIRE
  # --------------------------------------------------------------------------
  
  require(stochvol)
  require(DAARC)
  wd <- "R:/Asset_Management/Research_Projects/Multi_Asset/DAA/DAARC/Sensor_Testing/Univariate/Good_And_Bad_Turbulence/"
  source( paste0(wd, "Source/daarc_functions.R") )
  
  
  
  
  
  # --------------------------------------------------------------------------
  # LOAD DATA
  # --------------------------------------------------------------------------
  
  # debugonce(get.file)
  X <- get.file( wd = paste0(wd, "Data/"), filename = "FF_Industry_10", sep = ";" )
  
  # Benchmark 
  X_bm <- apply( X, 1, mean )
  colnames(X_bm) <- "bm"
  Y <- cumulated(X_bm, "discrete")
  
  
  
  # --------------------------------------------------------------------------
  # BBQ ALGO SETTINGS
  # --------------------------------------------------------------------------
  
  BBS <- BBSRC$new()
  BBS$setCtrl( minphase = 5 * 4 * 1,
               mincycle = 5 * 4 * 3,
               e = 0,
               k_peak = 10,
               k_trough = 10,
               l_peak = 10,
               l_trough = 10,
               theta = 0.1,
               language = "C",
               logarithmic = FALSE )
  BBS$data <- list(X_level = Y)
  BBS$run()
  
  
  
  
  

  # --------------------------------------------------------------------------
  # MINIMUM TORSION SV BBT
  # --------------------------------------------------------------------------
  
  
  
  BBTsv <- BBTSV$new()
  BBTsv$setCtrl( method = "base", universe = "dm" )
  BBTsv$spec$universe <- "ff10"
  BBTsv$spec$name <- "bbtsv_ff10"
  BBTsv$spec$iso <- colnames(X)
  BBTsv$spec$BBS <- BBS
  BBTsv$data <- list( X_bm = X_bm,
                      X = X )

  # # debugonce( BBTsv$computeSignalInsample )
  # BBTsv$computeSignalInsample()
  #   
  # 
  # 
  # signal <- ema( BBTsv$signal$insample[ ,"delta"], 0.1 )
  # sig <- (1 + sign(signal)) / 2
  # test <- signalTesting.byTrading( X = X_bm,
  #                                  sig = sig,
  #                                  n_lag = 2,
  #                                  tc = 0.004 )
  # X_tmp <- na.omit( cbind( bm = X_bm, bt = test) ) 
  # plot( log(cumulated(X_tmp, "discrete")), plot.type = "single" )
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # RUN LOOP
  # --------------------------------------------------------------------------
  
  minphase_vec <- floor(21 * c(0.75, 1, 1.5, 2, 3))
  mincycle_vec <- minphase_vec * 3
  
  for ( i in seq(along = minphase_vec) ) {
    
    BBS_tmp <- BBS$copy()
    BBS_tmp$spec$language <- "R"
    BBS_tmp$spec$minphase <- minphase_vec[i]
    BBS_tmp$spec$mincycle <- mincycle_vec[i]
    Obj <- BBTsv$copy()
    Obj$spec$BBS <- BBS_tmp
    Obj$spec$name <- paste0("bbtsv_ff10_phase=",  
                            Obj$spec$BBS$spec$minphase, "_cycle=",
                            Obj$spec$BBS$spec$mincycle, "_R")
    Obj$computeSignalInsample()
    Obj$save( path = paste0(wd, "/waRehouse/BBQ_Sensitivity/") )
    
    
    BBS_tmp <- BBS$copy()
    BBS_tmp$spec$language <- "C"
    BBS_tmp$spec$minphase <- minphase_vec[i]
    BBS_tmp$spec$mincycle <- mincycle_vec[i]
    Obj <- BBTsv$copy()
    Obj$spec$BBS <- BBS_tmp
    Obj$spec$name <- paste0("bbtsv_ff10_phase=",  
                            Obj$spec$BBS$spec$minphase, "_cycle=",
                            Obj$spec$BBS$spec$mincycle, "_C")
    Obj$computeSignalInsample()
    Obj$save( path = paste0(wd, "/waRehouse/BBQ_Sensitivity/") )
    
    
  }
  
  
  # Robust BBQ
  
  BBS_tmp <- BBS$copy()
  BBS_tmp$spec$language <- "R"
  Obj <- BBTsv$copy()
  Obj$spec$robust <- TRUE
  Obj$spec$minphase_vec <- minphase_vec
  Obj$spec$mincycle_vec <- mincycle_vec
  Obj$spec$BBS <- BBS_tmp
  Obj$spec$name <-"bbtsv_ff10_robust_R"
  # debugonce( Obj$computeSignalInsample )
  Obj$computeSignalInsample()
  Obj$save( path = paste0(wd, "/waRehouse/BBQ_Sensitivity/") )
    
  
  BBS_tmp <- BBS$copy()
  BBS_tmp$spec$language <- "C"
  Obj <- BBTsv$copy()
  Obj$spec$robust <- TRUE
  Obj$spec$minphase_vec <- minphase_vec
  Obj$spec$mincycle_vec <- mincycle_vec
  Obj$spec$BBS <- BBS_tmp
  Obj$spec$name <- "bbtsv_ff10_robust_C"
  Obj$computeSignalInsample()
  Obj$save( path = paste0(wd, "/waRehouse/BBQ_Sensitivity/") )
  
  
  Obj$data$BBS$plotStates()
  
  
  # --------------------------------------------------------------------------
  # BACKTESTS
  # --------------------------------------------------------------------------

  # Parameters
  tc <- 0.004
  n_lag <- 1
  penalty <- 0
  
  # Load sensors
  filenames <- list.files(path = paste0(wd, "waRehouse/BBQ_Sensitivity/"),
                          pattern = ".rds")
  Names <- unlist( lapply(filenames, FUN = function(x) { substring(x, 1, nchar(x)-4) }) )
  Names <- Names[ grepl("_ff10_", Names) ]
  lSensor <- lapply( Names, FUN = function(name) { try(loadSensor(name, wd =  paste0(wd, "waRehouse/BBQ_Sensitivity/"))) } )
  lSig <- lapply( lSensor, FUN = function(x) { if ( !inherits(x, "try-error")) x$signal$insample[ ,"delta"] } )
  signals <- do.call( cbind, lSig )
  colnames(signals) <- Names
  
  
  
  
  # Binary signals
  TD <- trainingData( Y_train = BBTsv$data$X_bm,
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
  # colnames(lStats$stats) <- colnames(lStats_nc$stats) <- c("bm", rep(paste0(minphase_vec, "_", mincycle_vec), 4) )
  
  
  
  
  colors <- c(1, fBasics::rainbowPalette(n = ncol(X_tmp)-1))
  plot( as.simTS(X_tmp), col = colors )
  plot( as.simTS(X_tmp_nc), col = colors )
  
  
  stats_fields <- c("cumret", "means", "sds", "sharpe", "maxDD")
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
  
  
  dev.off()
  barplot( t(cbind(lStats$stats["cumret", ],
                  lStats_nc$stats["cumret", ])), beside = TRUE )
  
  
  
  
      
  
  
  
  
  
