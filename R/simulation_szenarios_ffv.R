  
  
  ############################################################################
  ### BBTurbulence - FULLY FLEXIBLE VIEWS SZENARIOS SIMULATION 
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     17.12.2020
  # First version:    10.10.2020
  # --------------------------------------------------------------------------
  
  # Generate synthetic data as re-weighted real data.
  
  
  
  # --------------------------------------------------------------------------
  # REQUIRE
  # --------------------------------------------------------------------------
  
  require(stochvol)
  require(rugarch)
  require(DAARC)
  
  wd <- "R:/Asset_Management/Research_Projects/Multi_Asset/DAA/DAARC/Sensor_Testing/Univariate/Good_And_Bad_Turbulence/"
  source( paste0(wd, "Source/source_functions.R") )
  
  
  
  # --------------------------------------------------------------------------
  # 
  # --------------------------------------------------------------------------
  
  # BBT <- loadSensor( sensor_name = "bbtsv_scl_dm",
  #                    wd = "R:/Asset_Management/Research_Projects/Multi_Asset/DAA/Data/Sensors/" )
  BBT <- loadSensor( sensor_name = "bbtewma_base_dm",
                     wd = "R:/Asset_Management/Research_Projects/Multi_Asset/DAA/Data/Sensors/" )
  X <- BBT$data$X[isWeekday(time(BBT$data$X)), ]

  
  # Empirical probabilities
  p_emp <- rep(1/nrow(X), nrow(X))
  
  
  # Generate persistant signal
  set.seed(1234)
  P <- X * NA
  P[1, ] <- runif(ncol(X), 0, 1)
  for ( i in 2:nrow(X) ) {
    P[i, ] <- P[i-1, ] * 0.99 + rnorm(ncol(X), 0, 1)
  }
  P <- ema(P, 0.01) * 100
  P <- apply( P, 2, function(p) { (p + abs(min(p))) / sum(p + abs(min(p))) } ) 
  
  
  plot(P[ ,1])
  
  
 
  
  # Generate persistant signal based on garch volas of underlyings
  variance_model = list(model = "sGARCH",
                        garchOrder = c(1, 1),
                        submodel = NULL, 
                        external.regressors = NULL, 
                        variance.targeting = FALSE)
  mean_model = list(armaOrder = c(0, 0), 
                    include.mean = TRUE, 
                    archm = FALSE,
                    archpow = 1, 
                    arfima = FALSE, 
                    external.regressors = NULL, 
                    archex = FALSE)
  distribution_model = "norm"
  seed <- 1234
  start_pars <- list()
  fixed_pars <- list(mu = 0.0001,
                     omega = 0.0001,
                     alpha1 = 0.05,
                     beta1 = 0.94)
  spec <- ugarchspec( variance.model = variance_model,
                      mean.model = mean_model,
                      distribution.model = distribution_model,
                      start.pars = start_pars,
                      fixed.pars = fixed_pars )

  P <- X * NA
  for ( j in 1:ncol(X) ) {
    seed <- round(runif(1, 0, 10^3))
    fixed_pars$mu <- rnorm(1, 0, 0.001)
    path <- ugarchpath( spec = spec,
                        n.sim = nrow(X),
                        n.start = 1,
                        m.sim = 1,
                        rseed = seed )
    sim <- as.timeSeries( as.numeric(path@path$seriesSim),
                          time(X) )
    cvol <- getCondVar(garch(sim))
    P[ ,j] <- cvol / sum(cvol)
  }
  
  plot( P[ ,1:10] )
  
  
  
  
  # Generate alternative universes
  lZ <- list()
  lBBT <- list()
  for ( j in 1:ncol(P) ) {
    
    Z <- X * P[ ,rep(j, ncol(X))] * nrow(P) # * c(-1, 1)[round(runif(1, 1, 2))]
    Z_bm <- apply(Z, 1, mean)
    lZ[[j]] <- Z
  }
  
  sim_tmp <- do.call( cbind, lapply(lZ, FUN = function(x) { apply(x, 1, mean) } )  )
  
  plotSimTS( as.simTS(sim_tmp) )
  plotSimTS( as.simTS(apply(sim_tmp, 1, mean)) )
  
  W <- rdirichlet( n = 1, alpha = rep(1, ncol(sim_tmp)) )
  sim <- timeSeries( sim_tmp %*% t(W), time(sim_tmp) )
  plotSimTS( as.simTS(sim) )
  
  
  plotSimTS( as.simTS(lZ[[1]]) )
  
  
  
  
  for ( j in 1:ncol(P) ) {
    
    BBTSYNT <- BBTEWMA$new()
    # BBTSYNT <- BBTSV$new()
    # BBTSYNT <- BBTGARCH$new()
    BBTSYNT$setCtrl( method = "base", universe = "dm" )
    BBTSYNT$data <- list( X = lZ[[j]],
                          X_bm = apply(lZ[[j]], 1, mean),
                          BBS = BBTSYNT$spec$BBS )
    # debugonce( BBTSYNT$computeSignalInsample )
    BBTSYNT$computeSignalInsample()
    lBBT[[j]] <- BBTSYNT
  }
  
  
  lX_tmp <- list()
  for ( j in 1:length(lBBT) ) {
    
    Z_bm <- apply(lZ[[j]], 1, mean)
    sig <- (1 + sign( ema( lBBT[[j]]$signal$insample[ ,"delta"], 0.1 ) ) ) / 2
    test <- signalTesting.byTrading( X = Z_bm,
                                     sig = sig,
                                     n_lag = 1,
                                     tc = 0.001,
                                     penalty = 0.1 )
    colnames(test) <- colnames(sig)
    test_nc <- signalTesting.byTrading( X = Z_bm,
                                        sig = sig,
                                        n_lag = 1,
                                        tc = 0,
                                        penalty = 0.1 )
    colnames(test_nc) <- colnames(sig)
    X_tmp <- na.omit(cbind(Z_bm, test, test_nc))
    lX_tmp[[j]] <- X_tmp
    
    # plot( as.simTS(X_tmp), col = colors, main = "Insample - Before and after costs" )
    
  }
  
  
  stats_fields <- c("cumret", "means", "sds", "sharpe", "maxDD")
  for ( j in 1:length(lX_tmp) ) {
    plotSimTS( as.simTS(lX_tmp[[j]]), main = "Insample - Before and after costs" )
    t( descStats( lX_tmp[[j]] )$stats[stats_fields, ] )
    # drawDownStats(X_tmp)
  }
  
  
  
  lStats <- lapply( lX_tmp, FUN = descStats )
  lapply( lStats, FUN = function(x) { x$stats[stats_fields, ] } )
  
  
  sim <- do.call( cbind, lapply(lX_tmp, FUN = function(x) { x[ ,1] } )  )
  plot( as.simTS(sim) )
  
  
  
  
  
  # --------------------------------------------------------------------------
  # SAVE
  # --------------------------------------------------------------------------
  
  env <- new.env()
  env$lZ <- lZ
  env$lBBT <- lBBT
  
  saveRDS( env, file = paste0(wd, "waRehouse/tmp.rds") )
  
  
  
  
  
  
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # WIP
  # --------------------------------------------------------------------------
  
  BBT <- loadSensor( sensor_name = "bbturbulence" )
  X <- BBT$data$X[isWeekday(time(BBT$data$X)), ]
  X_bm <- BBT$data$X_bm[isWeekday(time(BBT$data$X_bm)), ]
  
  # Empirical probabilities
  p_emp <- rep(1/nrow(X), nrow(X))

 
  # Generate autocorrelated random probability vectors
  set.seed(1234)
  AR <- 0.99
  P <- X * NA
  P[1, ] <- runif(ncol(X), 0, 1)
  for ( i in 2:nrow(X) ) {
    P[i, ] <- P[i-1, ] * AR + rnorm(ncol(X), 0, 1)
  }
  P <- ema(P, 0.01) * 100
  P <- apply( P, 2, function(p) { (p + abs(min(p))) / sum(p + abs(min(p))) } ) 
  
  plot( P[ ,1] )
  
  
  
  # Generate synthetic returns
  lZ <- list()
  Z <- X * 0
  for ( j in 1:ncol(P) ) {
    plusminus <- sample( x = c(-1, 1), size = nrow(P), replace = TRUE, prob = c(0.3, 0.7))
    lZ[[j]] <- X * P[ ,rep(j, ncol(X))] * nrow(P) * plusminus
    Z <- Z + lZ[[j]]
  }
  Z <- Z / ncol(P)
  
  
    
  plot( as.simTS(X) )
  plot( as.simTS(Z) )
  plot( as.simTS(lZ[[1]]) )
  plot( as.simTS(lZ[[10]]) )
  
  
  lCov <- lapply( lZ, FUN = cov )
  lMu <- lapply( lZ, FUN = meanGeo )
  
  plot( as.CovMatrix(cov(X)) )
  plot( as.CovMatrix(lCov[[8]]) )

  barplot( unlist( lapply(lMu, FUN = function(x) { x[1] } ) ) )
  abline( h = meanGeo(X[ ,1]) )  
  
  
  
  
  # Backtest
  
  lBBT <- list()
  for ( j in 1:ncol(P) ) {
    
    BBTSYNT <- BBTEWMA$new()
    # BBTSYNT <- BBTSV$new()
    # BBTSYNT <- BBTGARCH$new()
    BBTSYNT$setCtrl( method = "base", universe = "dm" )
    BBTSYNT$data <- list( X = lZ[[j]],
                          X_bm = apply(lZ[[j]], 1, mean) )
    BBTSYNT$computeSignalInsample()
    lBBT[[j]] <- BBTSYNT
  }
  
  
  lX_bt <- list()
  for ( j in 1:length(lBBT) ) {
    
    
    Z_bm <- lBBT[[j]]$data$X_bm
    sig <- (1 + sign( ema( lBBT[[j]]$signal$insample[ ,"delta"], 0.1 ) ) ) / 2
    test <- signalTesting.byTrading( X = Z_bm,
                                     sig = sig,
                                     n_lag = 1,
                                     tc = 0.004,
                                     penalty = 0 )
    colnames(test) <- colnames(sig)
    test_nc <- signalTesting.byTrading( X = Z_bm,
                                        sig = sig,
                                        n_lag = 1,
                                        tc = 0,
                                        penalty = 0 )
    colnames(test_nc) <- colnames(sig)
    X_bt <- na.omit( cbind( Z_bm, test, test_nc ) )
    lX_bt[[j]] <- X_bt
    
    # plot( as.simTS(X_bt), col = colors, main = "Insample - Before and after costs" )
    
  }
  
  
  # for ( j in 1:length(lX_tmp) ) {
  #   plot( as.simTS(lX_tmp[[j]]), col = colors, main = "Insample - Before and after costs" )
  #   t( descStats( lX_tmp[[j]] )$stats[stats_fields, ] )
  #   # drawDownStats(X_tmp)
  # }
  
  
  stats_fields <- c("cumret", "means", "sds", "sharpe", "maxDD")
  lStats <- lapply( lX_bt, FUN = descStats )
  lapply( lStats, FUN = function(x) { x$stats[stats_fields, ] } )
  
  mu <- do.call( cbind, lapply( lStats, FUN = function(x) { x$stats["means", ] } ))
  sds <- do.call( cbind, lapply( lStats, FUN = function(x) { x$stats["sds", ] } ))
  sharpe <- do.call( cbind, lapply( lStats, FUN = function(x) { x$stats["sharpe", ] } ))
  maxdd <- do.call( cbind, lapply( lStats, FUN = function(x) { x$stats["maxDD", ] } ))
  
  
  barplot( mu[2:3, ] - mu[1, ], beside = TRUE, col = 1:2 )
  barplot( sds[2:3, ] - sds[1, ], beside = TRUE, col = 1:2 )
  barplot( sharpe[2:3, ] - sharpe[1, ], beside = TRUE, col = 1:2 )
  barplot( maxdd[2:3, ] - maxdd[1, ], beside = TRUE, col = 1:2 )
  barplot( maxdd, beside = TRUE, col = 1:3 )
  
  
  plot( as.simTS(lX_bt[[1]]) )
  
  
  sim <- do.call( cbind, lapply(lX_bt, FUN = function(x) { x[ ,1] } )  )
  plot( as.simTS(sim) )
  
  plot( as.simTS(lX_bt[[1]]) )
  
  
  
  