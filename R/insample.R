  

  ############################################################################
  ### GOOD AND BAD TURBULENCE - IN-SAMPLE ANALYSIS
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     22.10.2020
  # First version:    10.10.2020
  # --------------------------------------------------------------------------
  
  
  # 1) Define four turbulence measures: 
  #    md_all, md_bear, md_bull, md_delta.
  # 2) Check persistence of turbulence measures.
  # 3) Check performance conditional on turbulence measures.
  # 4) Repeat 1) - 3) for different data frequencies: 
  #    daily, weekly, monthly
  #    by i) aggregating daily signals and ii) computing signals on aggregated returns.
  # 5) Repeat 1) - 4) for different parametrizations of BBQ algo.
  # 6) Repeat 1) - 4) for fuzzy BBQ.
  
  
  # Open questions, tasks:
  # - What distance measure to use for the dirichletSampling
  # - how to condition on absolute and relative distance simultaneously
  # - what if we weight the MD by cap-weights and not cap-weights^2
  
  
  # --------------------------------------------------------------------------
  # REQUIRE
  # --------------------------------------------------------------------------
  
  require(RP)
  require(DAARC)
  
  wd <- "R:/Asset_Management/Research_Projects/Multi_Asset/DAA/DAARC/Sensor_Testing/Univariate/Good_And_Bad_Turbulence/"
  source( paste0(wd, "Source/source_functions.R") )
  
  
  
  # --------------------------------------------------------------------------
  # PARAMETERS
  # --------------------------------------------------------------------------
  
  minphase_m <- 2
  mincycle_m <- 5
  minphase_w <- minphase_m * 4
  mincycle_w <- mincycle_m * 4
  minphase_d <- minphase_m * 4 * 5
  mincycle_d <- mincycle_m * 4 * 5
  theta <- 0.15
  
  
    
  # --------------------------------------------------------------------------
  # 1) DEFINE TURBULENCE MEASURES
  # --------------------------------------------------------------------------
  
  # --------------------------------------------------------------------------
  # Construct absolute and relative Mahalanobis distances
  # --------------------------------------------------------------------------
  
  # MD <- loadSensor(sensor_name = "turbulence_scl2")
  BBT <- loadSensor(sensor_name = "bbturbulence_scl2_dm")
  BBT$signal <- list()
  BBT$spec$theta <- theta
  
  # #
  # # BBT$data$X <- BBT$data$X^2
  # BBT$data$X <- getCondVar(garch(BBT$data$X))^.5
  # #
 
  
  # Using daily returns
  BBTD <- BBT$copy()
  BBTD$spec$minphase <- minphase_d
  BBTD$spec$mincycle <- mincycle_d
  BBTD$spec$l_peak <- BBTD$spec$l_trough <- 
    BBTD$spec$k_peak <- BBTD$spec$k_trough <- minphase_d
  # debugonce(BBTD$computeSignalInsample)
  BBTD$computeSignalInsample()
  mdd <- BBTD$signal$insample
  mdd <- mdd[isWeekday(time(mdd)), ] 
  mdd_scl <- apply(mdd, 2, function(x) { x / max(x) } )
  
  
  # Using weekly returns
  BBTW <- BBT$copy()
  BBTW$data$X <- aggWeekly(BBT$data$X, day = "Tue", compounding = "discrete")
  BBTW$data$X_bm <- aggWeekly(BBT$data$X_bm, day = "Tue", compounding = "discrete")
  # BBTW$data$wmat <- applyRoll( BBT$data$wmat, 
  #                              By = 5, 
  #                              FUN = function(X) { apply(X, 2, mean) },
  #                              charvec = rownames(BBTW$data$X) )
  BBTW$data$wmat <- applyRoll( BBT$data$wmat, 
                               Width = 5, 
                               FUN = function(X) { apply(X, 2, mean) },
                               charvec = rownames(BBTW$data$X) )
  BBTW$spec$minphase <- minphase_w
  BBTW$spec$mincycle <- mincycle_w
  BBTW$spec$l_peak <- BBTW$spec$l_trough <- 
    BBTW$spec$k_peak <- BBTW$spec$k_trough <- minphase_w
  BBTW$computeSignalInsample()
  mdw <- BBTW$signal$insample
  mdw_scl <-  apply(mdw, 2, function(x) { x / max(x) } )
  
  
  # Using monthly returns
  BBTM <- BBT$copy()
  BBTM$data$X <- aggMonthly(BBTM$data$X, compounding = "discrete")
  BBTM$data$X_bm <- aggMonthly(BBTM$data$X_bm, compounding = "discrete")
  # BBTM$data$wmat <- applyRoll( BBT$data$wmat, 
  #                              By = 21, 
  #                              FUN = function(X) { apply(X, 2, mean) },
  #                              charvec = rownames(BBTM$data$X) )
  BBTM$data$wmat <- applyRoll( BBT$data$wmat, 
                               Width = 21, 
                               FUN = function(X) { apply(X, 2, mean) },
                               charvec = rownames(BBTM$data$X) )
  BBTM$spec$minphase <- minphase_m
  BBTM$spec$mincycle <- mincycle_m
  BBTM$spec$l_peak <- BBTM$spec$l_trough <- 
      BBTM$spec$k_peak <- BBTM$spec$k_trough <- minphase_m
  # debugonce(BBTM$computeSignalInsample)
  BBTM$computeSignalInsample()
  mdm <- BBTM$signal$insample
  mdm_scl <-  apply(mdm, 2, function(x) { x / max(x) } )
  

  # --------------------------------------------------------------------------
  # Rolling returns
  # --------------------------------------------------------------------------

  X_roll <- applyRoll( Data = BBT$data$X,
                       Width = 21,
                       By = 1,
                       FUN = function(X) { exp( apply(log(1 + X), 2, sum)) - 1 } )
  BBTMroll <- BBT$copy()
  BBTMroll$data$X <- X_roll
  BBTMroll$computeSignalInsample()
  mdm_roll <- BBTMroll$signal$insample

  BBT$computeSignalInsample()
  plot( cbind(BBT$signal$insample[ ,"delta"], 
              mdm_roll[ ,"delta"]), plot.type = "single" )
  
  
  s1 <- mdm_roll[ ,"delta"]
  s2 <- BBT$signal$insample[ ,"delta"]
  sig <- s1[ ,rep(1, 2)] * 0
  sig[s1 > 0, 1] <- 1
  sig[s2 > 0, 2] <- 1
  test <- signalTesting.byTrading( X = BBT$data$X_bm,
                                   sig = sig,
                                   n_lag = 1,
                                   tc = 0.004 )
  X_tmp <- na.omit( cbind(BBT$data$X_bm, test) )
  plot( as.simTS(X_tmp) )
  
  
  
  
    
  # --------------------------------------------------------------------------
  # Compare signal on aggregated return to aggregated signals
  # --------------------------------------------------------------------------
  
  # Aggregate daily distances to weekly
  mddma5 <- filter( mdd, rep(1/5, 5), side = 1 )
  colnames(mddma5) <- paste0( colnames(mdd), "ma5" )

  mddma21 <- filter( mdd, rep(1/21, 21), side = 1 )
  colnames(mddma21) <- paste0( colnames(mdd), "ma21" )
  
  
  # Monthly
  Name <- "all"
  plot( cbind(mdm[ ,Name], mddma21[rownames(mdm), paste0(Name, "ma21")]),
        plot.type = "single" )
  
  
  # Weekly
  Name <- "delta"
  plot( cbind(mdw[ ,Name], mddma5[rownames(mdw), paste0(Name, "ma5")]),
        plot.type = "single" )

  
  
  
  # --------------------------------------------------------------------------
  # Add garch implied volatility and squared returns as vola proxy
  # --------------------------------------------------------------------------
  
  # daily
  md_d <- cbind( mdd, 
                 x_sq = BBTD$data$X_bm[rownames(mdd), ]^2,
                 garch = getCondVar(garch(BBTD$data$X_bm[rownames(mdd), ])) )
  # weekly
  md_w <- cbind( mdw, 
                 bear_scl_ma5 = mddma5[rownames(mdw), "bear_sclma5"],
                 bull_scl_ma5 = mddma5[rownames(mdw), "bull_sclma5"],
                 delta_ma5 = mddma5[rownames(mdw), "deltama5"],
                 x_sq = BBTW$data$X_bm[rownames(mdw), ]^2,
                 garch = getCondVar(garch(BBTW$data$X_bm[rownames(mdw), ])) )
  # monthly
  md_m <- cbind( mdm, 
                 bear_scl_ma21 = mddma21[rownames(mdm), "bear_sclma21"],
                 bull_scl_ma21 = mddma21[rownames(mdm), "bull_sclma21"],
                 delta_ma21 = mddma21[rownames(mdm), "deltama21"],
                 x_sq = BBTM$data$X_bm[rownames(mdm), ]^2,
                 garch = getCondVar(garch(BBTM$data$X_bm[rownames(mdm), ])) )
  
  # Scale
  md_d_scl <-  apply(md_d, 2, function(x) { x / max(x) } )
  md_w_scl <-  apply(md_w, 2, function(x) { x / max(x) } )
  md_m_scl <-  apply(md_m, 2, function(x) { x / max(x) } )

  
  

  
  plot(md_d[ ,c("bear", "bear_scl")], plot.type = "single")
  
  
  
  # Omit certain signals for further analysis
  strategies <- c("states", "bear_scl", "bull_scl", "all", "delta", "x_sq", "garch")
  strategies_w <- c("states", "bear_scl", "bull_scl", "all", "delta", "x_sq", "garch",
                    "bear_scl_ma5", "bull_scl_ma5", "delta_ma5")
  strategies_m <- c("states", "bear_scl", "bull_scl", "all", "delta", "x_sq", "garch",
                    "bear_scl_ma21", "bull_scl_ma21", "delta_ma21")
  
  mdd <- md_d[ ,strategies]
  mdw <- md_w[ ,strategies_w]
  mdm <- md_m[ ,strategies_m]
  mdd_scl <- md_d_scl[ ,strategies]
  mdw_scl <- md_w_scl[ ,strategies_w]
  mdm_scl <- md_m_scl[ ,strategies_m]
  
  
  
  # --------------------------------------------------------------------------
  # 2) SIGNAL PERSISTENCE
  # --------------------------------------------------------------------------
  
  lacf_d <- apply( mdd, 2, acf )
  lacf_w <- apply( mdw, 2, acf )
  lacf_m <- apply( mdm, 2, acf )
  
  plot( lacf_d[[4]] )
  
  
  descStats( X = mdd, spec = descStatsSpec(what = "portmanteau") )
  descStats( X = mdw, spec = descStatsSpec(what = "portmanteau") )
  descStats( X = mdm, spec = descStatsSpec(what = "portmanteau") )
  
  
  # Compare with distance series computed on synthetic iidN data
  X <- log(1 + BBT$data$X)
  X_iidN <- rmvnorm( n = nrow(X), 
                     mean = apply(X, 2, mean), 
                     sigma = cov(X),
                     method = "eigen" )
  X_iidN <- timeSeries(X_iidN, time(X))
  md_iidN <- mahalanobisTS(X_iidN)
  
  # Compare with distance series computed on synthetic garch data
  fit <- garch(X)
  # debugonce(rmvgarch)
  X_garchN <- rmvgarch( fit, 
                        startDate = start(X),
                        nSample = 1000,   # BUG! other values won't work.
                        M1 = apply(X, 2, mean),
                        M2 = cov(X) )
  md_garchN <- mahalanobisTS(X_garchN)
  
  # md_synth <- cbind(iidN = md_iidN,
  #                   garchN = md_garchN)
  # plot( md_synth )
  # acf(md_iidN)
  plot( md_garchN )
  
  descStats( X = md_iidN, spec = descStatsSpec(what = "portmanteau") )
  descStats( X = md_garchN, spec = descStatsSpec(what = "portmanteau") )
  
  
 
  md_tmp <- mdm
  
  plot( x = md_tmp[ ,"delta"], y = md_tmp[ ,"all"] )
  
  
  color <- rev( fBasics::divPalette(n = nrow(md_tmp), "RdYlGn") )
  colors <- color[ order(as.numeric(md_tmp[ ,"delta"])) ]
  plot( x = as.numeric(md_tmp[ ,"bear"]), 
        y = as.numeric(md_tmp[ ,"bull"]), 
        xlab = "bear", ylab = "bull",
        xlim = range(md_tmp[ ,c("bear", "bull")]),
        ylim = range(md_tmp[ ,c("bear", "bull")]),
        col = colors )   # colors dont show correctly
  
  cor( md_tmp[ ,4:5] )
  cor( md_tmp[ ,2:3] )
  
  
  
  
  plot( mdm[ ,"bear_scl"] / mdm[ ,"bull_scl"] )
  plot( mdm[ ,"bull_scl"] / mdm[ ,"bear_scl"] )
  plot( BBTM$signal$insample[ ,"bear"] / BBTM$signal$insample[ ,"bull"] )
  
  
  KCDLR <- kcdens( Y_train = TDM$Y_train,
                   X_train = mdm[ ,"bear_scl"] / mdm[ ,"bull_scl"],
                   # X_eval = X_eval,
                   n_lag = 0,
                   kcd_fun = "liracine",
                   bw_method = "normal-reference" )
  
  matplot( x = KCDLR$y_eval,
           y = t(KCDLR$cdens),
           type = "l" )
  
  
  
  
  # --------------------------------------------------------------------------
  # 3) CONDITIONAL PERFORMANCE
  # --------------------------------------------------------------------------
  
  # --------------------------------------------------------------------------
  # Posterior means
  # --------------------------------------------------------------------------
  
  TDD <- trainingData( Y_train = BBTD$data$X_bm,
                       X_train = mdd_scl )
  TDW <- trainingData( Y_train = BBTW$data$X_bm,
                       X_train = mdw_scl )
  TDM <- trainingData( Y_train = BBTM$data$X_bm,
                       X_train = mdm_scl )
  lDSD <- list()
  for ( j in 1:ncol(TDD$X_train) ) {
    
    # debugonce(dirichletSampling)
    X_eval <- seq(from = min(TDD$X_train[ ,j]), to = max(TDD$X_train[ ,j]), length.out = 5)
    X_eval2 <- quantile( TDD$X_train[ ,j], seq(from = 0, to = 1, length.out = 5) )
    cbind(X_eval, X_eval2)
    
    lDSD[[j]] <- dirichletSampling( Y_train = TDD$Y_train,
                                    X_train = TDD$X_train[ ,j],
                                    # X_eval = X_eval,
                                    weights_fun = "kernel",
                                    correct_bias = FALSE,
                                    sclfct = NULL )
    
  }
  names(lDSD) <- colnames(TDD$X_train)
  
  
  lDSW <- list()
  lDSM <- list()
  for ( j in 1:ncol(TDW$X_train) ) {
    
    lDSW[[j]] <- dirichletSampling( Y_train = TDW$Y_train,
                                    X_train = TDW$X_train[ ,j],
                                    # X_eval = X_eval,
                                    weights_fun = "kernel",
                                    # centers = 5,
                                    correct_bias = FALSE,
                                    sclfct = NULL )
    lDSM[[j]] <- dirichletSampling( Y_train = TDM$Y_train,
                                    X_train = TDM$X_train[ ,j],
                                    weights_fun = "kernel",
                                    correct_bias = FALSE,
                                    sclfct = NULL )
  }
  names(lDSW) <- colnames(TDW$X_train)
  names(lDSM) <- colnames(TDM$X_train)
  

  lDS <- lDSW
  
  # from <- quantile( unlist(lDS), 0 )
  # to <- quantile( unlist(lDS), 1 )
  from <- min( unlist(lDS) ) * 1.3
  to <- max( unlist(lDS) ) * 1.3
  n <- 999
  lldens <- lapply( lDS, FUN = function(x) { 
                lapply(x, density, from = from, to = to, n = n) } )
  for ( i in 1:length(lldens) ) {
    plot.ldensity( lldens[[i]], main = paste0("Posterior means conditional on ", names(lDS)[[i]]) )
  }
  
  
  
  
  
   
  
  lDS <- lDSD
  lMu <- lapply( lDS, FUN = function(x) { unlist(lapply(x, FUN = mean)) } )
  Mu <- do.call( cbind, lMu )
  
  colors <- rev(fBasics::divPalette(n = nrow(Mu), "RdYlGn"))
  colors <- rev(fBasics::divPalette(n = nrow(Mu), "Spectral"))
  barplot( Mu, beside = TRUE, col = colors )
  names(lMu)
  
  
  
  # Performance after signal
  # debugonce(pas)
  lPAS <- pas( X = TDD$Y_train,
               sig = TDD$X_train[ ,"delta"], 
               n_lags = 21 )
  
  str( lPAS[[1]] )
   
  
  i <- 1
  plot( x = as.numeric(lPAS[[i]]$X_train),
        y = as.numeric(lPAS[[i]]$Y_train) )
  abline(h = 0, col = "grey")
  
  
  
  
  
  
  plot( x = TDD$X_train[ ,"x_sq"], y = TDD$Y_train  )
  
  
  KCDH <- kcdens( Y_train = TDM$Y_train,
                  X_train = TDM$X_train[ ,5],
                  kcd_fun = "hyndman" )
  
  plot(KCDH$fit)
  KCD$x_eval
  
  
  plot( TDM$X_train[ ,4:5] )
  
  
  a <- seq(from = 0, to = 1, length.out = 10)
  b <- seq(from = -1, to = 1, length.out = 10)
  expand.grid(a, b)
  
  KCDLR <- kcdens( Y_train = TDM$Y_train,
                   X_train = TDM$X_train[ ,2:3],
                   # X_eval = X_eval,
                   n_lag = 1,
                   kcd_fun = "liracine",
                   bw_method = "normal-reference" )
  
  matplot( x = KCDLR$y_eval,
           y = t(KCDLR$cdens),
           type = "l",
           col = colors )
  
  mu_kcdlr <- apply(KCDLR$cdens, 1, function(p) { p %*% KCDLR$y_eval } )
  names(mu_kcdlr) <- rownames(KCDLR$x_eval)
  barplot( mu_kcdlr )
  
  
  
  cbind( KCDLR$x_eval, mu_kcdlr )
  
  
  
  # --------------------------------------------------------------------------
  # Regimes
  # --------------------------------------------------------------------------
  
  # Calm, Turbulent, Good Turbulent, Bad Turbulent
  
  
  mdd_regimes <- mdd[ ,rep(1, 4)] * 0
  colnames(mdd_regimes) <- c("calm", "turbulent", "bad_turbulent", "good_turbulent")
  mdd_regimes[mdd[ ,"all"] < quantile(mdd[ ,"all"], 0.9), 1] <- 1
  mdd_regimes[mdd[ ,"all"] > quantile(mdd[ ,"all"], 0.9), 2] <- 1
  mdd_regimes[mdd[ ,"delta"] < quantile(mdd[ ,"delta"], 0.1), 3] <- 1
  mdd_regimes[mdd[ ,"delta"] > quantile(mdd[ ,"delta"], 0.9), 4] <- 1
  
  plot(mdd_regimes)
  
  
  mdw_regimes <- mdw[ ,rep(1, 4)] * 0
  colnames(mdw_regimes) <- c("calm", "turbulent", "bad_turbulent", "good_turbulent")
  mdw_regimes[mdw[ ,"all"] < quantile(mdw[ ,"all"], 0.9), 1] <- 1
  mdw_regimes[mdw[ ,"all"] > quantile(mdw[ ,"all"], 0.9), 2] <- 1
  mdw_regimes[mdw[ ,"delta"] < quantile(mdw[ ,"delta"], 0.1), 3] <- 1
  mdw_regimes[mdw[ ,"delta"] > quantile(mdw[ ,"delta"], 0.9), 4] <- 1
  
  
  plot(mdw_regimes)
  
  
  
  mdm_regimes <- mdm[ ,rep(1, 4)] * 0
  colnames(mdm_regimes) <- c("calm", "turbulent", "bad_turbulent", "good_turbulent")
  mdm_regimes[mdm[ ,"all"] < quantile(mdm[ ,"all"], 0.9), 1] <- 1
  mdm_regimes[mdm[ ,"all"] > quantile(mdm[ ,"all"], 0.9), 2] <- 1
  mdm_regimes[mdm[ ,"delta"] < quantile(mdm[ ,"delta"], 0.1), 3] <- 1
  mdm_regimes[mdm[ ,"delta"] > quantile(mdm[ ,"delta"], 0.9), 4] <- 1
  
  plot(mdm_regimes)
  
  
  
  
  
  possible_regimes <- rbind( c(1, 0, 0, 0), 
                             c(1, 0, 1, 0), 
                             c(1, 0, 0, 1),
                             c(0, 1, 0, 0),
                             c(0, 1, 1, 0),
                             c(0, 1, 0, 1) )
  
  # Create binary timeseries for each regime 
  md_regimes <- mdw_regimes
  regimes <- t( apply( md_regimes, 1,
                       function(s) { apply(possible_regimes, 1, 
                                           function(r) { as.numeric(all(s == r)) } ) } ) )
  regimes <- timeSeries(regimes, time(md_regimes) )
  # colnames(regimes) <- paste0("regime_", 1:ncol(regimes))
  colnames(regimes) <- c("calm", "calm-bad", "calm-good",
                         "turbulent", "turbulent-bad", "turbulent-good")
  
  # Create categorical timeseries containing regimes 1 to 8
  market_regimes <- apply( regimes, 1, function(x) { which(x == 1) } )
  
  plot( regimes )
  plot( market_regimes )
  
  
  Y_train <- lag(TDW$Y_train, 0)  # ~~~~~~~~~~~~~~~~ 
  Y_train[is.na(Y_train)] <- 0
  y <- as.numeric(Y_train)
  lmu_regimes <- list()
  lcdens_regimes <- list()
  for ( j in 1:ncol(regimes) ) {
    
    P <- rdirichlet( n = 10^3, 
                     alpha = regimes[ ,j] / sum(regimes[ ,j]) )
    lcdens_regimes[[j]] <- P %*% y
    P <- rdirichlet( n = 10^3, 
                     alpha = regimes[ ,j] / sum(regimes[ ,j]) * nrow(regimes) )
    lmu_regimes[[j]] <- P %*% y
  }
  
  ldens <- lapply( lcdens_regimes, density )
  names(ldens) <- colnames(regimes)
  colors <- 1:ncol(regimes)
  plot.ldensity( ldens, colors = colors, main = "Conditional Posterior Predictive Distributions" )
  
  
  ldens <- lapply( lmu_regimes, density )
  names(ldens) <- colnames(regimes)
  colors <- 1:ncol(regimes)
  plot.ldensity( ldens, colors = colors, main = "Posterior Conditional Means" )
  
  
  
  DS <- dirichletSampling( Y_train = TDW$Y_train,
                           X_train = market_regimes,
                           correct_bias = FALSE,
                           sclfct = NULL )
  ldens <- lapply( DS, density )
  plot.ldensity( ldens )
  
  
  KCD <- kcdens( Y_train = TDW$Y_train,
                 X_train = market_regimes,
                 kcd_fun = "liracine" )
  
  matplot( x = KCD$y_eval,
           y = t(KCD$cdens),
           type = "l" )             
  
  
  
  
  # --------------------------------------------------------------------------
  # 
  # --------------------------------------------------------------------------
  
  
  
  
  
  
  
  
  
  
  
 