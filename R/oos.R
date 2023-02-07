  
  
  ############################################################################
  ### GOOD AND BAD TURBULENCE - OUT-OF-SAMPLE ANALYSIS
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     16.10.2020
  # First version:    10.10.2020
  # --------------------------------------------------------------------------
  
  
  
  
  
  
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
  BBT <- loadSensor(sensor_name = "bbturbulence")
  BBT$signal <- list()
  BBT$spec$theta <- theta
  
  
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
  
  
  # # Aggregate daily distances to weekly
  # mddma5 <- filter( mdd, rep(1/5, 5), side = 1 )
  # colnames(mddma5) <- c("mddma5")
  # 
  # mddma21 <- filter( mdd, rep(1/21, 21), side = 1 )
  # colnames(mddma21) <- c("mddma21")
  
  

  
  
  
  # --------------------------------------------------------------------------
  # Add garch implied volatility and squared returns as vola proxy
  # --------------------------------------------------------------------------
  
  # daily
  md_d <- cbind( mdd, 
                 x_sq = BBTD$data$X_bm[rownames(mdd), ]^2,
                 garch = getCondVar(garch(BBTD$data$X_bm[rownames(mdd), ])) )
  # weekly
  md_w <- cbind( mdw, 
                 # mddma5,
                 x_sq = BBTW$data$X_bm[rownames(mdw), ]^2,
                 garch = getCondVar(garch(BBTW$data$X_bm[rownames(mdw), ])) )
  # monthly
  md_m <- cbind( mdm, 
                 # mddma21,
                 x_sq = BBTM$data$X_bm[rownames(mdm), ]^2,
                 garch = getCondVar(garch(BBTM$data$X_bm[rownames(mdm), ])) )
  
  # Scale
  md_d_scl <-  apply(md_d, 2, function(x) { x / max(x) } )
  md_w_scl <-  apply(md_w, 2, function(x) { x / max(x) } )
  md_m_scl <-  apply(md_m, 2, function(x) { x / max(x) } )
  
  
  
  plot( md_d_scl )
  plot( md_w_scl )
  plot( md_m_scl )
  
  plot(md_d[ ,c("bear", "bear_scl")], plot.type = "single")
  
  
  
  # Omit certain signals for further analysis
  strategies <- c("states", "bear_scl", "bull_scl", "all", "delta", "x_sq", "garch")
  # strategies_w <- c("states", "bear_scl", "bull_scl", "all", "delta", "x_sq", "garch")
  
  mdd <- md_d[ ,strategies]
  mdw <- md_w[ ,strategies]
  mdm <- md_m[ ,strategies]
  mdd_scl <- md_d_scl[ ,strategies]
  mdw_scl <- md_w_scl[ ,strategies]
  mdm_scl <- md_m_scl[ ,strategies]
  
  
  
 
  
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
  lDSW <- list()
  lDSM <- list()
  
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
    lDSW[[j]] <- dirichletSampling( Y_train = TDW$Y_train,
                                    X_train = TDW$X_train[ ,j],
                                    # X_eval = X_eval,
                                    weights_fun = "kernel",
                                    # centers = 5,
                                    correct_bias = FALSE,
                                    sclfct = 1 )
    lDSM[[j]] <- dirichletSampling( Y_train = TDM$Y_train,
                                    X_train = TDM$X_train[ ,j],
                                    weights_fun = "kernel",
                                    correct_bias = FALSE,
                                    sclfct = NULL )
  }
  
  
  lDS <- lDSM
  names(lDS) <- colnames(TDD$X_train)
  
  
  # from <- quantile( unlist(lDS), 0 )
  # to <- quantile( unlist(lDS), 1 )
  from <- min( unlist(lDS) ) * 1.3
  to <- max( unlist(lDS) ) * 1.3
  n <- 999
  lldens <- lapply( lDS, FUN = function(x) { lapply(x, density, from = from, to = to, n = n) } )
  plot.ldensity( lldens[[1]], main = "Posterior means conditional on market regime" )
  plot.ldensity( lldens[[2]], main = "Posterior means conditional on bear Turbulence" )
  plot.ldensity( lldens[[3]], main = "Posterior means conditional on bull Turbulence" )
  plot.ldensity( lldens[[4]], main = "Posterior means conditional on absolute turbulence" )
  plot.ldensity( lldens[[5]], main = "Posterior means conditional on relative turbulence", fillin = FALSE )
  plot.ldensity( lldens[[5]], main = "Posterior means conditional on relative turbulence")
  plot.ldensity( lldens[[6]], main = "Posterior means conditional on X^2" )
  plot.ldensity( lldens[[7]], main = "Posterior means conditional on garch vola" )
  
  
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
  # CLASSIFICATION
  # --------------------------------------------------------------------------

  
  
  dates <- rownames(BBT$data$X)
  dates_train <- head(dates, floor(length(dates) * 0.5))
  dates_test <- dates[ !(dates %in% dates_train) ]
  
  
  # Daily
  
  # Construct signal on training data
  Obj <- BBT$copy()
  X_w <- aggWeekly(BBT$data$X, 
                   day = "Tue",
                   compounding = "discrete")
  Obj$data$X <- BBT$data$X[dates_train, ]
  Obj$data$X_bm <- BBT$data$X_bm[dates_train, ]
  Obj$data$wmat <- BBT$data$wmat[dates_train, ] 
  Obj$spec$minphase <- minphase_d
  Obj$spec$mincycle <- mincycle_d
  Obj$spec$l_peak <- Obj$spec$l_trough <- 
      Obj$spec$k_peak <- Obj$spec$k_trough <- minphase_d
  Obj$computeSignalInsample()
  mdd <- Obj$signal$insample
  mdd_scl <-  apply(mdd, 2, function(x) { x / max(x) } )
  
  # Test data
  X_test <- BBT$data$X[dates_test, ]
  mu <- Obj$signal$insample_stats$mu
  sigma <- Obj$signal$insample_stats$sigma
  mu_bear <- Obj$signal$insample_stats$mu_bear
  sigma_bear <- Obj$signal$insample_stats$sigma_bear
  mu_bull <- Obj$signal$insample_stats$mu_bull
  sigma_bull <- Obj$signal$insample_stats$sigma_bull
  wmat <- BBT$data$wmat
  charvec <- intersect( dates_test, rownames(wmat) )
  FUN <- function(today)
  {
    x <- as.numeric(X_test[today, ])
    wghts <- as.numeric(wmat[today, ])
    md_bear <- weightedMahalanobis(x = x,
                                   center = mu_bear, 
                                   covmat = sigma_bear, 
                                   wghts = wghts, 
                                   scl = TRUE)
    md_bull <- weightedMahalanobis(x = x,
                                   center = mu_bull, 
                                   covmat = sigma_bull, 
                                   wghts = wghts, 
                                   scl = TRUE)
    md_all <- weightedMahalanobis(x = x,
                                  center = mu,
                                  covmat = sigma,
                                  wghts = wghts, 
                                  scl = TRUE)
    md <- c(md_bear, md_bull, md_all)
    return( md )
  }
  lMD <- lapply( charvec, FUN = FUN )
  MD <- timeSeries( do.call( rbind, lMD ), charvec )
  colnames(MD) <- c("bear", "bull", "all")
  MD_scl <- scale( MD, FALSE, TRUE )
  md_delta <- MD_scl[ ,"bear"] - MD_scl[ ,"bull"]
  MD_scl <- cbind(MD_scl, delta = md_delta)
  
  plot( MD_scl )
  
  
  
  states_train <- Obj$signal$insample[ ,"states"]
  states_test <- (1 + sign(ema(md_delta, 0.1))) / 2
  states_test_true <- (1 + sign(BBTD$signal$insample[rownames(states_test), "states"])) / 2
  
  ROC <- roc( y_hat = states_test,
              y_true = states_test_true, 
              th = 1 )
  ROC    
  
  
  plot( cbind(states_test, states_test_true) )
  
  
  test <- signalTesting.byTrading( X = BBTD$data$X_bm,
                                   sig = states_test,
                                   n_lag = 1,
                                   tc = 0.004 )
  X_tmp <- cbind(BBTD$data$X_bm[rownames(test), ], test = test )  
  plot( as.simTS(X_tmp) )  
  
  descStats(X_tmp)
  
  
  
  
  # Weekly

  # Construct signal on training data
  Obj <- BBT$copy()
  X_w <- aggWeekly(BBT$data$X, 
                   day = "Tue",
                   compounding = "discrete")
  Obj$data$X <- X_w[intersect(dates_train, rownames(X_w)), ]
  Obj$data$X_bm <- aggWeekly(BBT$data$X_bm[dates_train, ], 
                             day = "Tue", 
                             compounding = "discrete")
  Obj$data$wmat <- applyRoll( BBT$data$wmat, 
                              Width = 5,
                              FUN = function(X) { apply(X, 2, mean) },
                              charvec = rownames(X_w) )
  Obj$spec$minphase <- minphase_w
  Obj$spec$mincycle <- mincycle_w
  Obj$spec$l_peak <- Obj$spec$l_trough <- 
    Obj$spec$k_peak <- Obj$spec$k_trough <- minphase_w
  Obj$computeSignalInsample()
  mdw <- Obj$signal$insample
  mdw_scl <-  apply(mdw, 2, function(x) { x / max(x) } )

  
  # Test data
  X_test <- X_w[intersect(dates_test, rownames(X_w)), ]
  mu <- Obj$signal$insample_stats$mu
  sigma <- Obj$signal$insample_stats$sigma
  mu_bear <- Obj$signal$insample_stats$mu_bear
  sigma_bear <- Obj$signal$insample_stats$sigma_bear
  mu_bull <- Obj$signal$insample_stats$mu_bull
  sigma_bull <- Obj$signal$insample_stats$sigma_bull
  wmat <- Obj$data$wmat
  charvec <- intersect( dates_test, rownames(Obj$data$wmat) )
  FUN <- function(today)
  {
    x <- as.numeric(X_test[today, ])
    wghts <- as.numeric(wmat[today, ])
    md_bear <- weightedMahalanobis(x = x,
                                   center = mu_bear, 
                                   covmat = sigma_bear, 
                                   wghts = wghts, 
                                   scl = TRUE)
    md_bull <- weightedMahalanobis(x = x,
                                   center = mu_bull, 
                                   covmat = sigma_bull, 
                                   wghts = wghts, 
                                   scl = TRUE)
    md_all <- weightedMahalanobis(x = x,
                                  center = mu,
                                  covmat = sigma,
                                  wghts = wghts, 
                                  scl = TRUE)
    md <- c(md_bear, md_bull, md_all)
    return( md )
  }
  lMD <- lapply( charvec, FUN = FUN )
  MD <- timeSeries( do.call( rbind, lMD ), charvec )
  colnames(MD) <- c("bear", "bull", "all")
  MD_scl <- scale( MD, FALSE, TRUE )
  md_delta <- MD_scl[ ,"bear"] - MD_scl[ ,"bull"]
  
  MD_scl <- cbind(MD_scl, delta = md_delta)
  
  plot( MD_scl )
  
  
  
  states_train <- Obj$signal$insample[ ,"states"]
  states_test <- (1 + sign(ema(md_delta, 0.1))) / 2
  states_test_true <- (1 + sign(BBTW$signal$insample[rownames(states_test),"states"])) / 2
  
  ROC <- roc( y_hat = states_test,
              y_true = states_test_true, 
              th = 1 )
  ROC    
  
  
  plot( cbind(states_test, states_test_true) )
  
  
  test <- signalTesting.byTrading( X = BBTW$data$X_bm,
                                   sig = states_test,
                                   n_lag = 1,
                                   tc = 0.004 )
  X_tmp <- cbind(BBTW$data$X_bm[rownames(test), ], test )  
  plot( as.simTS(X_tmp) )  
  
  
  
  
  
  
  
  
