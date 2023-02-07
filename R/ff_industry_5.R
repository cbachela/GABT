  
  
  ############################################################################
  ### GOOD AND BAD TURBULENCE - FAMA FRENCH INDUSTRY PORTFOLIOS
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     09.03.2021
  # First version:    09.03.2021
  # --------------------------------------------------------------------------
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # REQUIRE
  # --------------------------------------------------------------------------
  
  require(stochvol)
  require(RP)
  require(DAARC)
  
  wd <- "R:/Asset_Management/Research_Projects/Multi_Asset/DAA/DAARC/Sensor_Testing/Univariate/Good_And_Bad_Turbulence/"


  
  
  # --------------------------------------------------------------------------
  # LOAD DATA
  # --------------------------------------------------------------------------
  
  MD <- loadSensor( "turbulence" )
  X <- MD$data$X
  
  # debugonce(get.file)
  X <- get.file( wd = paste0(wd, "Data/"), filename = "FF_Industry_5", sep = ";" )
  
  head(X)
  plot(X)
  
  
  # Benchmark 
  X_bm <- apply( X, 1, mean )
  colnames(X_bm) <- "bm"
  
  
  
  descStats(cbind(X, X_bm))
  
  
 
  
  # --------------------------------------------------------------------------
  # BBQ ALGO SETTINGS
  # --------------------------------------------------------------------------
  
  BBSF <- BBSFuzzy$new()
  BBSF$setCtrl()
  BBSF$spec$universe <- "ff5"
  BBSF$spec$language <- "R"
  BBSF$spec$round <- 0
  BBSF$data <- list( X_bm = X_bm,
                     X = X_bm )
  
  BBSF$computeSignalInsample()
  head(BBSF$signal$insample)
  
  BBSF$getSignal
  BBSF$spec$round
  
  BBSF$plotStates( type = "is" )
  
  mu_bear <- meanGeo( X = BBSF$data$X_bm[which(BBSF$getSignal(method = "insample") == -1), ] )
  mu_bull <- meanGeo( X = BBSF$data$X_bm[which(BBSF$getSignal(method = "insample") == 1), ] )

  mu_bear; mu_bull
  
  
  BBSFC <- BBSF$copy()
  BBSFC$spec$language <- "C"
  BBSFC$computeSignalInsample()
  
  sig1 <- BBSF$getSignal( method = "insample" )
  sig2 <- BBSFC$getSignal( method = "insample" )
  plot( cbind( sig1, sig2) )
  
  apply( cbind(sig1, sig2), 2, sum )
  
  
  
  # --------------------------------------------------------------------------
  # TURBULENCE
  # --------------------------------------------------------------------------
  
  MDBase <- Turbulence$new()
  MDBase$setCtrl()
  MDBase$spec$b_scl <- FALSE
  MDBase$spec$width <- Inf
  MDBase$spec$universe <- "ff5"
  MDBase$data <- list(X = X, 
                      S = X,
                      X_bm = X_bm)
  
  
  # MDBase <- loadSensor( sensor_name = "turbulence_base_ff5" )
  # MDBase$update()

  MDBase$computeSignalInsample()
  # MDBase$updateSignal()
  
  signal <- MDBase$signal$insample[ ,"md"]
  signal <- MDBase$signal$base[ ,"md"]
  plot(signal)
  
  
  plot( cbind(log(cumulated(X_bm, "discrete")), 
              MDBase$signal$insample[ ,1:5]) )
  
  
  plot( MDBase$signal$insample[ ,"sig_md_ema0.157"])
  
  
  tailleft( MDBase$signal$base )
  
  
  
  
  signal <- ema( MDBase$signal$insample[ ,"sig_md"], 1)
  sig <- round(signal)
  test <- signalTesting.byTrading( X = X_bm,
                                   sig = sig,
                                   n_lag = 1,
                                   tc = 0 )
  X_tmp <- na.omit( cbind( bm = X_bm, bt = test) ) 
  plot( log(cumulated(X_tmp, "discrete")), plot.type = "single" )
  
  plot( log(cumulated(tail(X_tmp, 10000), "discrete")), plot.type = "single" )
  
  descStats(X_tmp)
  drawDownStats(X_tmp)
  
  
  
  # # Conditional density
  # DS <- dirichletSampling( Y_train = X_bm,
  #                          X_train = signal,
  #                          n_lag = 0,
  #                          weights_fun = "l1",
  #                          sclfct = NULL )
  # ldens <- lapply( DS, FUN = densFUN )
  # plot.ldensity( ldens )
  
  
  
  
  
  # --------------------------------------------------------------------------
  # BEAR / BULL TURBULENCE
  # --------------------------------------------------------------------------
  
  # Base
  BBTBase <- BBTurbulence$new()
  BBTBase$setCtrl( method = "base" )
  BBTBase$spec$universe <- "ff5"
  BBTBase$spec$iso <- colnames(X)
  BBTBase$spec$BBS <- BBSF
  BBTBase$data <- list( X_bm = X_bm,
                        X = X,
                        BBS = BBSF )
  # BBTBase <- loadSensor( sensor_name = "bbturbulence_base_ff5" )
  # BBTBase$computeSignalInsample()
  # BBTBase$update()
  
  # debugonce(BBTBase$computeSignalInsample)
  BBTBase$computeSignalInsample()
  
  
  edit(BBTBase$computeSignalBase)
  
  
  
  plot( ema(BBTBase$signal$insample[ ,"delta"], 0.1) )
  abline(h = 0)
  
  signal <- ema( BBTBase$signal$insample[ ,"delta"], 1)
  sig <- (1 + sign(signal)) / 2
  test <- signalTesting.byTrading( X = X_bm,
                                   sig = sig,
                                   n_lag = 1,
                                   tc = 0 )
  X_tmp <- na.omit( cbind( bm = X_bm, bt = test) ) 
  plot( log(cumulated(X_tmp, "discrete")), plot.type = "single" )
  
  plot( log(cumulated(tail(X_tmp, 10000), "discrete")), plot.type = "single" )
  
  descStats(X_tmp)
  drawDownStats(X_tmp)
  
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # MINIMUM TORSION EWMA BBT
  # --------------------------------------------------------------------------
  
  # Ewma
  BBTEwma <- BBTEWMA$new()
  BBTEwma$setCtrl( method = "base", universe = "dm" )
  BBTEwma$spec$universe <- "ff5"
  BBTEwma$spec$name <- "bbtewma_base_ff5"
  BBTEwma$spec$ewma_alpha <- 0.1
  BBTEwma$spec$iso <- colnames(X)
  BBTEwma$spec$name <- paste0(BBTEwma$spec$name, "_a", BBTEwma$spec$ewma_alpha)
  BBTEwma$spec$BBS <- BBSF
  BBTEwma$data <- list( X_bm = X_bm,
                        X = X )
  # BBTEwma <- loadSensor( sensor_name = "bbtewma_base_ff5_a0.1" )
  # BBTEwma$computeSignalInsample()
  # BBTEwma$update()
  
  BBTEwma$computeSignalInsample()
  edit(BBTEwma$computeSignal)

  
  
  signal <- ema( BBTEwma$signal$insample[ ,"delta"], 1)
  sig <- (1 + sign(signal)) / 2
  test <- signalTesting.byTrading( X = X_bm,
                                   sig = sig,
                                   n_lag = 1,
                                   tc = 0 )
  X_tmp <- na.omit( cbind( bm = X_bm, bt = test) ) 
  plot( log(cumulated(X_tmp, "discrete")), plot.type = "single" )
  
  plot( log(cumulated(tail(X_tmp, 10000), "discrete")), plot.type = "single" )
  
  descStats(X_tmp)
  drawDownStats(X_tmp)
  
  
  plot(signal)
  abline(h = 0)    
  
  plot( BBTEwma$data$BBS$getSignal(method = "insample") )
  
  
  
  
  
  # --------------------------------------------------------------------------
  # MINIMUM TORSION GARCH BBT
  # --------------------------------------------------------------------------
  
  BBSF$spec$language <- "C"
  BBSF$spec$minphase_vec <- 21 * c(1, 2, 3, 4)
  BBSF$spec$mincycle_vec <- BBSF$spec$minphase_vec * 3
  BBTgarch <- DAARC::BBTGARCH$new()
  BBTgarch$setCtrl( method = "base", universe = "dm" )
  BBTgarch$spec$universe <- "ff5"
  BBTgarch$spec$name <- "BBTgarch_ff5"
  BBTgarch$spec$iso <- colnames(X)
  BBTgarch$spec$BBS <- BBSF
  BBTgarch$data <- list( X_bm = X_bm,
                        X = X )
  # BBTgarch <- loadSensor( sensor_name = "bbtgarch_ff5" )
  # BBTgarch$computeSignalInsample()
  # BBTgarch$update()
  
  BBTgarch$computeSignalInsample()
  
  
  
  signal <- ema( BBTgarch$signal$insample[ ,"delta"], 0.1 )
  sig <- (1 + sign(signal)) / 2
  test <- signalTesting.byTrading( X = X_bm,
                                   sig = sig,
                                   n_lag = 1,
                                   tc = 0.004 )
  X_tmp <- na.omit( cbind( bm = X_bm, bt = test) ) 
  plot( log(cumulated(X_tmp, "discrete")), plot.type = "single" )
  
  plot( log(cumulated(tail(X_tmp, 10000), "discrete")), plot.type = "single" )
  
  descStats(X_tmp)
  drawDownStats(X_tmp)
  
  
  plot(signal)
  abline(h = 0)    
  
  
  plot( BBTgarch$signal$insample )
  
  
  
  
  # --------------------------------------------------------------------------
  # MINIMUM TORSION SV BBT
  # --------------------------------------------------------------------------
  
  BBSF$spec$language <- "C"
  BBSF$spec$minphase_vec <- 21 * c(1, 2, 3, 4)
  # BBSF$spec$minphase_vec <- 21 * c(0.5, 0.75, 1, 2)
  BBSF$spec$mincycle_vec <- BBSF$spec$minphase_vec * 3
  
  BBTsv <- BBTSV$new()
  BBTsv$setCtrl( method = "base", universe = "dm" )
  BBTsv$spec$universe <- "ff5"
  BBTsv$spec$name <- "bbtsv_ff5"
  BBTsv$spec$iso <- colnames(X)
  BBTsv$spec$BBS <- BBSF
  BBTsv$data <- list( X_bm = X_bm,
                      X = X )
  # BBTsv <- loadSensor( sensor_name = "bbtsv_base_ff5" )
  # BBTsv$computeSignalInsample()
  # BBTsv$update()
  
  # debugonce( BBTsv$computeSignalInsample )
  BBTsv$computeSignalInsample()
  BBTsv$updateSignal()
  # BBTsv$save()
  
  
  plot( BBTsv$signal$insample[ ,c("delta", "prob")] )
  
  
  
  
  signal <- ema( BBTsv$signal$insample[ ,"delta"], 0.1)
  signal <- ema( BBTsv$signal$insample[ ,"prob"] - 0.5, 1 )
  signal <- ema( BBTsv$getSignal(), 0.1 )
  sig <- (1 + sign(signal)) / 2
  test <- signalTesting.byTrading( X = X_bm,
                                   sig = sig,
                                   n_lag = 1,
                                   tc = 0.004 )
  X_tmp <- na.omit( cbind( bm = X_bm, bt = test) ) 
  test_nc <- signalTesting.byTrading( X = X_bm,
                                      sig = sig,
                                      n_lag = 1,
                                      tc = 0 )
  X_tmp_nc <- na.omit( cbind( bm = X_bm, bt_nc = test_nc) )
  X_tmp_all <- na.omit( cbind( bm = X_bm, bt = test, bt_nc = test_nc ) )
  
  plot( log(cumulated(X_tmp, "discrete")), plot.type = "single" )
  
  # plot( as.simTS(X_tmp) )
  # plot( as.simTS(X_tmp_nc) )
  plot( as.simTS(window(X_tmp_all, "1932-01-01", end(X_tmp_all))) )
  
  
  
  descStats(X_tmp_all)
  meanGeo(X_tmp_all, scalefactor = 252)
  drawDownStats(X_tmp_all)
  
  
  plot(signal)
  abline(h = 0)    

  plot( BBTsv$spec$BBS$getSignal( "insample" ) )
  
  
  
  # # Conditional density
  # DS <- dirichletSampling( Y_train = X_bm,
  #                          X_train = signal,
  #                          n_lag = 0,
  #                          weights_fun = "l1",
  #                          sclfct = NULL )
  # ldens <- lapply( DS, FUN = densFUN )
  # plot.ldensity( ldens )
  
  
  
  
  
  # --------------------------------------------------------------------------
  # RELATIVE PERFORMANCE
  # --------------------------------------------------------------------------
  
  X_delta <- (X_tmp[ ,"bt"] - X_tmp[ ,"bm"]) / (1 + X_tmp[ ,"bm"])
  colnames(X_delta) <- "Outperformance"
  
  plot( as.simTS(X_delta) )
  
  
  # Relative drawdown
  drawDownPlot( X = X_delta )
  
  DD <- drawDownStats( X = X_delta )
  DD
  
  
  # Rolling 4-year excess returns
  FUN <- function(X) 
  { 
    tmp <- apply( X, 2, function(x) { exp( sum( log( 1 + x ) ) ) - 1 } )
    tmp[2] - tmp[1]
  }
  X_delta_roll_1 <- applyRoll( Data = X_tmp, 
                               Width = 252 * 1,
                               By = 1,
                               FUN = FUN )
  plot(X_delta_roll_1)

    
  
  
  
  # --------------------------------------------------------------------------
  # GARCH vs SV
  # --------------------------------------------------------------------------
  
  BBS <- BBSRC$new()
  BBS$setCtrl( language = "C",
               minphase = 21,
               mincycle = 21 * 3 )
  BBS$data <- list( X_level = cumulated(X_bm, "discrete") )
  BBS$runRobust( )
  bbstates <- BBS$output$states
  
  
  dates <- intersect(rownames(X), rownames(bbstates))
  X_train <- X[dates, ]
  bbstates <- bbstates[dates, ]
  
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
  S_all <- getter( mtca(X_train), "sources" )
  S_bear <- timeSeries( scale(X_train, mu_bear, FALSE) %*% 
                          torsion(covmat_bear, method = "mt"), time(X_train) )
  S_bull <- timeSeries( scale(X_train, mu_bull, FALSE) %*% 
                          torsion(covmat_bull, method = "mt"), time(X_train) )
  # S_all_wght <- S_all
  # S_bear_wght <- S_bear
  # S_bull_wght <- S_bull

  # Stochastic volatility MD in MT space
  sv_all <- sv_bear <- sv_bull <- S_all * NA
  n_sim <- 1000
  sv_array_bear <- array( NA, dim = c(nrow(S_all), ncol(S_all), n_sim) )
  sv_array_bull <- sv_array_bear
  for ( j in 1:ncol(S_all) ) {
    
    # draws <- NULL
    # draws <- svsample( S_all[ ,j], draws = 1000, burnin = 100 )
    # sv_all[ ,j] <- apply( exp(draws$latent), 2, mean )
    
    draws <- NULL
    draws <- svsample( S_bear[ ,j], draws = n_sim, burnin = 100 )
    h_bear <- t(exp(draws$latent))
    sv_bear[ ,j] <- apply( h_bear, 1, mean )
    # sv_bear[ ,j] <- apply( h_bear, 1, median )
    sv_array_bear[ ,j, ] <- scale( h_bear, FALSE, rep(var(S_bear[ ,j]), n_sim) )
    
    draws <- NULL
    draws <- svsample( S_bull[ ,j], draws = n_sim, burnin = 100 )
    h_bull <- t(exp(draws$latent))
    sv_bull[ ,j] <- apply( h_bull, 1, mean )
    # sv_bull[ ,j] <- apply( h_bull, 1, median )
    sv_array_bull[ ,j, ] <- scale( h_bull, FALSE, rep(var(S_bull[ ,j]), n_sim) )
  }
  
  md_bear_mat <- do.call( cbind, lapply(1:n_sim, FUN = function(k) { apply( sv_array_bear[ ,,k], 1, sum ) }) )
  md_bear_mat <- timeSeries( md_bear_mat, time(sv_bear) )
  md_bear_mat_scl <- scale( md_bear_mat, FALSE, TRUE )

  md_bull_mat <- do.call( cbind, lapply(1:n_sim, FUN = function(k) { apply( sv_array_bull[ ,,k], 1, sum ) }) )
  md_bull_mat <- timeSeries( md_bull_mat, time(sv_bull) )
  md_bull_mat_scl <- scale( md_bull_mat, FALSE, TRUE )

  md_delta_mat <- md_bear_mat_scl - md_bull_mat_scl
  prob <- apply( md_delta_mat, 1, function(x) { sum( x > 0) / n_sim } )
  
  # md_all <- apply( scale(sv_all, FALSE, apply(S_all, 2, var)), 1, sum)
  md_bear <- apply( scale(sv_bear, FALSE, apply(S_bear, 2, var)), 1, sum)
  md_bear <- timeSeries( md_bear, time(sv_bear) )
  md_bull <- apply( scale(sv_bull, FALSE, apply(S_bull, 2, var)), 1, sum)
  md_bull <- timeSeries( md_bull, time(sv_bull) )
  
  MD <- cbind( bear = md_bear, 
               bull = md_bull )
  MD_scl <- scale( MD, FALSE, TRUE )
  md_delta <- MD_scl[ ,"bear"] - MD_scl[ ,"bull"]
  
  
  
  # GARCH MD in MT space
  cvar_bear <- getCondVar( garch( S_bear ) )
  cvar_bull <- getCondVar( garch( S_bull ) ) 
  cvar_all <- getCondVar( garch( S_all ) )
  
  md_bear_garch <- apply( scale(cvar_bear, FALSE, apply(S_bear, 2, var)), 1, sum )
  md_bear_garch <- timeSeries( md_bear_garch, time(cvar_bear) )
  md_bull_garch <- apply( scale(cvar_bull, FALSE, apply(S_bull, 2, var)), 1, sum )
  md_bull_garch <- timeSeries( md_bull_garch, time(cvar_bull) )
  md_all_garch <- apply( scale(cvar_all, FALSE, apply(S_all, 2, var)), 1, sum )
  md_all_garch <- timeSeries( md_all_garch, time(cvar_all) )
  
  MD_garch <- cbind( bear = md_bear_garch, 
                     bull = md_bull_garch,
                     all = md_all_garch )
  MD_garch_scl <- scale( MD_garch, FALSE, TRUE )
  md_garch_delta <- MD_garch_scl[ ,"bear"] - MD_garch_scl[ ,"bull"]
  
  
  
  
  plot( cbind(sv = md_delta, garch = md_garch_delta), plot.type = "single" )
  
  plot( cbind(MD_garch_scl[ ,"bear"], MD_scl[ ,"bear"]), plot.type = "single" )
  
  plot( prob )
  
  plot( x = as.numeric(md_garch_delta), y = as.numeric(X_bm) )
  
  
  # Conditional density
  DS <- dirichletSampling( Y_train = X_bm,
                           X_train = md_garch_delta,
                           n_lag = 1,
                           weights_fun = "l1",
                           sclfct = NULL )
  ldens <- lapply( DS, FUN = density )
  plot.ldensity( ldens )
  
  
  
  
  
  # GARCH
  fit <- garch(X[ ,2])
  cvol <- getCondVar(fit)

  # SV
  draws <- svsample( X[ ,2], draws = 10^3, burnin = 100 )
  h <- timeSeries( apply( exp(draws$latent), 2, mean), time(X) )
  
  draws <- svsample( X[ ,2], draws = 10^4, burnin = 1000 )
  h2 <- timeSeries( apply( exp(draws$latent), 2, mean), time(X) )
  
  plot( head(cbind(cvol, h, h2), 100), plot.type = "single" )

    
  
  
  plot( cbind(BBTsv$signal$insample[ ,"bear"], 
              BBTgarch$signal$insample[ ,"bear"]), plot.type = "single" )
  plot( ema(cbind(BBTsv$signal$insample[ ,"delta"], 
              BBTgarch$signal$insample[ ,"delta"]), 0.1),
        plot.type = "single" )
  
  
  
    
  signal <- ema( cbind(BBTsv$signal$insample[ ,"delta"],
                       BBTgarch$signal$insample[ ,"delta"]), 1 )
  # signal <- ema( cbind(md_delta, md_garch_delta), 1 )
  sig <- (1 + sign(signal)) / 2
  test <- signalTesting.byTrading( X = X_bm,
                                   sig = sig,
                                   n_lag = 0,
                                   tc = 0.00 )
  X_tmp <- na.omit( cbind( bm = X_bm, bt = test) ) 
  plot( log(cumulated(X_tmp, "discrete")), plot.type = "single" )
  
  
  
  
  
  
  
  
  