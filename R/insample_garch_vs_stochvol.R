  
  
  ############################################################################
  ### GOOD AND BAD TURBULENCE - INSAMPLE - GARCH VS STOCHVOL VS EWMA
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     18.03.2021
  # First version:    18.03.2021
  # --------------------------------------------------------------------------
  
  
  
  
  
  
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
  # EWMA
  # --------------------------------------------------------------------------
  
  BBTewma <- BBTEWMA$new()
  BBTewma$setCtrl( universe = "dm",
                    method = "base" )
  BBTewma$spec$ewma_alpha <- 0.1
  BBTewma$spec$BBS <- BBS
  BBTewma$data <- list( X = X,
                         X_bm = X_bm )
  
  for ( k in seq(along = minphase_vec) ) {
    
    # Bry-Boschan in R
    BBTewma$spec$BBS$spec$language <- "R"
    BBTewma$spec$BBS$spec$minphase <- minphase_vec[k]
    BBTewma$spec$BBS$spec$mincycle <- minphase_vec[k] * 3
    BBTewma$spec$name <- paste0( "bbtewma_bbr_", universe, "_phase=", 
                                  BBTewma$spec$BBS$spec$minphase, "_cycle=",
                                  BBTewma$spec$BBS$spec$mincycle )
    # debugonce( BBTewma$computeSignalInsample )
    BBTewma$computeSignalInsample()
    BBTewma$save( path = paste0( wd, "waRehouse/Garch_vs_Stochvol/" ) )
    
    # Bry-Boschan in C
    BBTewma$spec$BBS$spec$language <- "C"
    BBTewma$spec$BBS$spec$minphase <- minphase_vec[k]
    BBTewma$spec$BBS$spec$mincycle <- minphase_vec[k] * 3
    BBTewma$spec$name <- paste0( "bbtewma_bbc_", universe, "_phase=", 
                                  BBTewma$spec$BBS$spec$minphase, "_cycle=",
                                  BBTewma$spec$BBS$spec$mincycle )
    BBTewma$computeSignalInsample()
    BBTewma$save( path = paste0( wd, "waRehouse/Garch_vs_Stochvol/") )
    
  }
  
  # Robust Bry-Boschan in R
  BBTewma$spec$BBS$spec$language <- "R"
  BBTewma$spec$robust <- TRUE
  BBTewma$spec$minphase_vec <- minphase_vec
  BBTewma$spec$mincycle_vec <- minphase_vec * 3
  BBTewma$spec$name <- paste0( "bbtewma_bbr_robust_", universe )
  # debugonce( BBTewma$computeSignalInsample )
  BBTewma$computeSignalInsample()
  BBTewma$save( path = paste0( wd, "waRehouse/Garch_vs_Stochvol/") )
  
  # Robust Bry-Boschan in C
  BBTewma$spec$BBS$spec$language <- "C"
  BBTewma$spec$robust <- TRUE
  BBTewma$spec$minphase_vec <- minphase_vec
  BBTewma$spec$mincycle_vec <- minphase_vec * 3
  BBTewma$spec$name <- paste0( "bbtewma_bbc_robust_", universe )
  BBTewma$computeSignalInsample()
  BBTewma$save( path = paste0( wd, "waRehouse/Garch_vs_Stochvol/") )
  
  # Lunde-Timmermann
  BBTewma$spec$BBS$spec$language <- "C"
  BBTewma$spec$BBS$spec$algo <- "Lunde-Timmermann"
  BBTewma$spec$name <- paste0( "bbtewma_ltc_", universe, "_theta=", 
                                BBTewma$spec$BBS$spec$theta )
  BBTewma$computeSignalInsample()
  BBTewma$save( path = paste0( wd, "waRehouse/Garch_vs_Stochvol/") )
  
  
  
  # --------------------------------------------------------------------------
  # GARCH
  # --------------------------------------------------------------------------
  
  BBTgarch <- BBTGARCH$new()
  BBTgarch$setCtrl( universe = "dm",
                    method = "base" )
  BBTgarch$spec$BBS <- BBS
  BBTgarch$data <- list( X = X,
                         X_bm = X_bm )
  
  for ( k in seq(along = minphase_vec) ) {
    
    # Bry-Boschan in R
    BBTgarch$spec$BBS$spec$language <- "R"
    BBTgarch$spec$BBS$spec$minphase <- minphase_vec[k]
    BBTgarch$spec$BBS$spec$mincycle <- minphase_vec[k] * 3
    BBTgarch$spec$name <- paste0( "bbtgarch_bbr_", universe, "_phase=", 
                                  BBTgarch$spec$BBS$spec$minphase, "_cycle=",
                                  BBTgarch$spec$BBS$spec$mincycle )
    # debugonce( BBTgarch$computeSignalInsample )
    BBTgarch$computeSignalInsample()
    BBTgarch$save( path = paste0( wd, "waRehouse/Garch_vs_Stochvol/" ) )
    
    # Bry-Boschan in C
    BBTgarch$spec$BBS$spec$language <- "C"
    BBTgarch$spec$BBS$spec$minphase <- minphase_vec[k]
    BBTgarch$spec$BBS$spec$mincycle <- minphase_vec[k] * 3
    BBTgarch$spec$name <- paste0( "bbtgarch_bbc_", universe, "_phase=", 
                                  BBTgarch$spec$BBS$spec$minphase, "_cycle=",
                                  BBTgarch$spec$BBS$spec$mincycle )
    BBTgarch$computeSignalInsample()
    BBTgarch$save( path = paste0( wd, "waRehouse/Garch_vs_Stochvol/") )
    
  }
  
  # Robust Bry-Boschan in R
  BBTgarch$spec$BBS$spec$language <- "R"
  BBTgarch$spec$robust <- TRUE
  BBTgarch$spec$minphase_vec <- minphase_vec
  BBTgarch$spec$mincycle_vec <- minphase_vec * 3
  BBTgarch$spec$name <- paste0( "bbtgarch_bbr_robust_", universe )
  BBTgarch$computeSignalInsample()
  BBTgarch$save( path = paste0( wd, "waRehouse/Garch_vs_Stochvol/") )
  
  # Robust Bry-Boschan in C
  BBTgarch$spec$BBS$spec$language <- "C"
  BBTgarch$spec$robust <- TRUE
  BBTgarch$spec$minphase_vec <- minphase_vec
  BBTgarch$spec$mincycle_vec <- minphase_vec * 3
  BBTgarch$spec$name <- paste0( "bbtgarch_bbc_robust_", universe )
  BBTgarch$computeSignalInsample()
  BBTgarch$save( path = paste0( wd, "waRehouse/Garch_vs_Stochvol/") )
  
  # Lunde-Timmermann
  BBTgarch$spec$BBS$spec$language <- "C"
  BBTgarch$spec$BBS$spec$algo <- "Lunde-Timmermann"
  BBTgarch$spec$name <- paste0( "bbtgarch_ltc_", universe, "_theta=", 
                                BBTgarch$spec$BBS$spec$theta )
  BBTgarch$computeSignalInsample()
  BBTgarch$save( path = paste0( wd, "waRehouse/Garch_vs_Stochvol/") )
  
  
  
  
  # --------------------------------------------------------------------------
  # STOCHVOL
  # --------------------------------------------------------------------------
  
  BBTsv <- BBTSV$new()
  BBTsv$setCtrl( universe = "dm",
                   method = "base" )
  BBTsv$spec$BBS <- BBS
  BBTsv$data <- list( X = X,
                      X_bm = X_bm )
  
  for ( k in seq(along = minphase_vec) ) {
    
    # Bry-Boschan in R
    BBTsv$spec$BBS$spec$language <- "R"
    BBTsv$spec$BBS$spec$minphase <- minphase_vec[k]
    BBTsv$spec$BBS$spec$mincycle <- minphase_vec[k] * 3
    BBTsv$spec$name <- paste0( "bbtsv_bbr_", universe, "_phase=", 
                                 BBTsv$spec$BBS$spec$minphase, "_cycle=",
                                 BBTsv$spec$BBS$spec$mincycle )
    # debugonce( BBTsv$computeSignalInsample )
    BBTsv$computeSignalInsample()
    BBTsv$save( path = paste0( wd, "waRehouse/Garch_vs_Stochvol/" ) )
    
    # Bry-Boschan in C
    BBTsv$spec$BBS$spec$language <- "C"
    BBTsv$spec$BBS$spec$minphase <- minphase_vec[k]
    BBTsv$spec$BBS$spec$mincycle <- minphase_vec[k] * 3
    BBTsv$spec$name <- paste0( "bbtsv_bbc_", universe, "_phase=", 
                                 BBTsv$spec$BBS$spec$minphase, "_cycle=",
                                 BBTsv$spec$BBS$spec$mincycle )
    BBTsv$computeSignalInsample()
    BBTsv$save( path = paste0( wd, "waRehouse/Garch_vs_Stochvol/") )
    
  }
  
  # Robust Bry-Boschan in R
  BBTsv$spec$BBS$spec$language <- "R"
  BBTsv$spec$robust <- TRUE
  BBTsv$spec$minphase_vec <- minphase_vec
  BBTsv$spec$mincycle_vec <- minphase_vec * 3
  BBTsv$spec$name <- paste0( "bbtsv_bbr_robust_", universe )
  BBTsv$computeSignalInsample()
  BBTsv$save( path = paste0( wd, "waRehouse/Garch_vs_Stochvol/") )
  
  # Robust Bry-Boschan in C
  BBTsv$spec$BBS$spec$language <- "C"
  BBTsv$spec$robust <- TRUE
  BBTsv$spec$minphase_vec <- minphase_vec
  BBTsv$spec$mincycle_vec <- minphase_vec * 3
  BBTsv$spec$name <- paste0( "bbtsv_bbc_robust_", universe )
  BBTsv$computeSignalInsample()
  BBTsv$save( path = paste0( wd, "waRehouse/Garch_vs_Stochvol/") )
  
  # Lunde-Timmermann
  BBTsv$spec$BBS$spec$language <- "C"
  BBTsv$spec$BBS$spec$algo <- "Lunde-Timmermann"
  BBTsv$spec$name <- paste0( "bbtsv_ltc_", universe, "_theta=", 
                               BBTsv$spec$BBS$spec$theta )
  BBTsv$computeSignalInsample()
  BBTsv$save( path = paste0( wd, "waRehouse/Garch_vs_Stochvol/") )
  
  
  
  
  # --------------------------------------------------------------------------
  # ANALYZE
  # --------------------------------------------------------------------------
  
  analyze <- function()
  {
    
    require(stochvol)
    require(DAARC)
    wd <- "R:/Asset_Management/Research_Projects/Multi_Asset/DAA/DAARC/Sensor_Testing/Univariate/Good_And_Bad_Turbulence/"
    source( paste0(wd, "Source/daarc_functions.R") )
    
    
    filenames <- list.files(path = paste0(wd, "waRehouse/Garch_vs_Stochvol/"),
                            pattern = ".rds")
    strategy_names <- unlist( lapply(filenames, FUN = function(x) { substring(x, 1, nchar(x)-4) }) )
    lSensor <- lapply( strategy_names, FUN = function(name) { try(loadSensor(name, wd =  paste0(wd, "waRehouse/Garch_vs_Stochvol/"))) } )
    names(lSensor) <- strategy_names
    
    
    # EWMA
    names_ewma <- strategy_names[ grepl( "bbtewma_", strategy_names ) ]
    lSig_ewma <- lapply( lSensor[ names_ewma ], FUN = function(x) { if ( !inherits(x, "try-error")) x$signal$insample[ ,"delta"] } )
    signals_ewma <- do.call( cbind, lSig_ewma )
    
    
    # GARCH
    names_garch <- strategy_names[ grepl( "bbtgarch_", strategy_names ) ]
    lSig_garch <- lapply( lSensor[ names_garch ], FUN = function(x) { if ( !inherits(x, "try-error")) x$signal$insample[ ,"delta"] } )
    signals_garch <- do.call( cbind, lSig_garch )
    
    
    # Stochvol
    names_sv <- strategy_names[ grepl( "bbtsv_", strategy_names ) ]
    lSig_sv <- lapply( lSensor[ names_sv ], FUN = function(x) { if ( !inherits(x, "try-error")) x$signal$insample[ ,"delta"] } )
    lSig_sv_prob <- lapply( lSensor[ names_sv ], FUN = function(x) { if ( !inherits(x, "try-error")) x$signal$insample[ ,"prob"] } )
    names(lSig_sv_prob) <- paste0( names(lSig_sv), "_prob" )
    signals_sv_delta <- do.call( cbind, lSig_sv )
    signals_sv_prob <- do.call( cbind, lSig_sv_prob ) - 0.5   #//
    signals_sv <- cbind( signals_sv_delta, signals_sv_prob )
    
    
    
    # signals <- cbind( signals_ewma, signals_garch, signals_sv )
    signals <- signals_garch
    
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
    plot( as.simTS(X_tmp[ ,1:10]), col = colors )
    plot( as.simTS( tail( X_tmp[ ,1:10], 500 ) ), col = colors )
    plot( as.simTS( tail( X_tmp, 500 ) ) )
    
    
    
    plot.timeSeries( log(cumulated(X_tmp, "discrete")),  col = colors )
    plot.timeSeries( log(cumulated(X_tmp_nc, "discrete")), col = colors )
    
    
    
    
    
    
    
    s1 <- lSensor[["bbtewma_ltc_dm_theta=0.1"]]$signal$insample
    s2 <- lSensor[["bbtgarch_ltc_dm_theta=0.1"]]$signal$insample
    s3 <- lSensor[["bbtsv_ltc_dm_theta=0.1"]]$signal$insample
    
    name <- "states"
    plot( cbind( s1[ ,name], s2[ ,name], s3[ ,name] ), plot.type = "single" )
    
    name <- "delta"
    plot( tail( cbind( s1[ ,name], s2[ ,name], s3[ ,name] ), 500 ), plot.type = "single" )
    abline( h = 0 )
    
    
    
    s1 <- lSensor[["bbtewma_bbc_dm_phase=21_cycle=63"]]$signal$insample
    s2 <- lSensor[["bbtewma_bbc_robust_dm"]]$signal$insample
    
    name <- "states"
    plot( cbind( s1[ ,name], s2[ ,name] ), plot.type = "single" )
      
    
  }
  
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  xxx <- function()
  {
    
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
    BBS$run()
    bbstates <- BBS$output$states
    
    
    
    # --------------------------------------------------------------------------
    # BEARS, BULLS AND MINIMUM TORSION
    # --------------------------------------------------------------------------
    
    X_train <- X
    
    # Compute regime specific stats
    
    # All
    covmat <- cov(X_train)
    
    # Bear
    p_bear <- (1 - bbstates) / 2
    stats_bear <- cov.wt(x = X_train, 
                         wt = as.numeric( p_bear / sum(p_bear) ),
                         cor = FALSE, 
                         center = TRUE, 
                         method = "ML")
    mu_bear <- stats_bear$center
    covmat_bear <- make.positive.definite( stats_bear$cov )
    
    # Bull
    p_bull <- (1 + bbstates) / 2
    stats_bull <- cov.wt(x = X_train, 
                         wt = as.numeric( p_bull / sum(p_bull) ),
                         cor = FALSE, 
                         center = TRUE, 
                         method = "ML")
    mu_bull <- stats_bull$center
    covmat_bull <- make.positive.definite( stats_bull$cov )
    
    # Minimum torsion series
    MT <-  mtca(X_train)
    S_all <- getter( MT, "sources" )
    S_bear <- timeSeries( scale(X_train, mu_bear, FALSE) %*% 
                            torsion(covmat_bear, method = "mt"), time(X_train) )
    S_bull <- timeSeries( scale(X_train, mu_bull, FALSE) %*% 
                            torsion(covmat_bull, method = "mt"), time(X_train) )
    
    
    # --------------------------------------------------------------------------
    # DIRICHLET SAMPLING IN MT SPACE
    # --------------------------------------------------------------------------
    
    n_sim <- 10^3
    h_array_bear <- array( NA, dim = c(nrow(S_all), ncol(S_all), n_sim) )
    h_array_bull <- h_array_bear
    
    for ( j in 1:ncol(S_all) ) {
      
      xsq <- S_bear[ ,j]^2
      xsq_ema <- as.numeric( ema( xsq, alpha = 0.1 ) )
      alpha <- xsq_ema / sum(xsq_ema)
      P <- rdirichlet( n = n_sim, alpha = alpha * length(alpha) )
      h_bear <- apply( P, 1, function(p) { p * xsq_ema * length(p) } )
      h_array_bear[ ,j, ] <- scale( h_bear, FALSE, rep(var(S_bear[ ,j]), n_sim) )
      
      xsq <- S_bull[ ,j]^2
      xsq_ema <- as.numeric( ema( xsq, alpha = 0.1 ) )
      alpha <- xsq_ema / sum(xsq_ema)
      P <- rdirichlet( n = n_sim, alpha = alpha * length(alpha) )
      h_bull <- apply( P, 1, function(p) { p * xsq_ema * length(p) } )
      h_array_bull[ ,j, ] <- scale( h_bull, FALSE, rep(var(S_bull[ ,j]), n_sim) )
      
    }
    
    md_bear_mat <- do.call( cbind, lapply(1:n_sim, FUN = function(k) { apply( h_array_bear[ ,,k], 1, sum ) }) )
    md_bear_mat <- timeSeries( md_bear_mat, time(S_bear) )
    md_bear_mat_scl <- scale( md_bear_mat, FALSE, TRUE )
    
    md_bull_mat <- do.call( cbind, lapply(1:n_sim, FUN = function(k) { apply( h_array_bull[ ,,k], 1, sum ) }) )
    md_bull_mat <- timeSeries( md_bull_mat, time(S_bull) )
    md_bull_mat_scl <- scale( md_bull_mat, FALSE, TRUE )
    
    md_delta_mat <- md_bear_mat_scl - md_bull_mat_scl
    prob <- apply( md_delta_mat, 1, function(x) { sum( x > 0) / n_sim } )
    
    
    plot( prob )
    plot( apply( md_delta_mat, 1, mean) ) 
    plot( apply( md_bear_mat, 1, mean) ) 
    plot( apply( md_bull_mat, 1, mean) ) 
    
    
    
    
    # Binary signals
    signals <- cbind( delta = apply( md_delta_mat, 1, mean),
                      prob = prob - 0.5 )
    TD <- trainingData( Y_train = X_bm,
                        X_train = (1 + sign( ema( signals, 1 ) )) / 2 )
    
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
    t(lStats$stats[stats_fields, ])    
    
    plot( as.simTS(X_tmp) ) 
    plot( as.simTS(X_tmp_nc) ) 
    
      
    
  
    
    
  }
  
  
  
  
  
  