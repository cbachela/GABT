  
  
  ############################################################################
  ### BBT ON ROBECO VOLA PERCENTILE DATA
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     02.06.2021
  # First version:    02.06.2021
  # --------------------------------------------------------------------------
    
  
  
  
  # --------------------------------------------------------------------------
  # REQUIRE
  # --------------------------------------------------------------------------
  
  require(DAARC)
  
  wd <- "R:/Asset_Management/Research_Projects/Multi_Asset/DAA/DAARC/Sensor_Testing/Univariate/Good_And_Bad_Turbulence/"
  source( paste0(wd, "Source/source_functions.R") )
  
  
  # --------------------------------------------------------------------------
  # Data
  # --------------------------------------------------------------------------
  
  X_all <- get.file( wd = paste0(wd, "Data/"), 
                     filename = "Vola_Percentile_Portfolios_Robeco" )
  X <- X_all[ ,-which(colnames(X_all) == "Rm.Rf")]
  
  # Remove great depression
  # X_all <- X_all[rownames(X_all) > "1939-12-31", ]
  ##
  
  X_bm <- X_all[ ,"Rm.Rf"]
  X_eqw <- apply( X, 1, mean )
  
  colors <- c(fBasics::divPalette(n = ncol(X), "RdYlGn"), 1)
  plot( as.simTS(X_all), col = colors )
  lStats <- descStats( X_all )
  lStats$stats
  
  barplot( lStats$stats["sharpe", ])
  
  
  
  # Load Fama-French industry portfolios
  wd_ff <- "R:/Asset_Management/Research_Projects/Multi_Asset/DAA/DAARC/Sensor_Testing/Univariate/Good_And_Bad_Turbulence/"
  X_ff <- get.file( wd = paste0(wd_ff, "Data/"), 
                    filename = "FF_Industry_5", sep = ";" )
  X_ff_eqw <- apply( X_ff, 1, mean )
  
  # Aggregate to monthly
  X_ff <- aggMonthly( X = X_ff, compounding = "discrete" )
  X_ff_eqw <- aggMonthly( X = X_ff_eqw, compounding = "discrete" )
  
  
  dates <- intersect( rownames(X_all), rownames(X_ff) )
  X_ff <- X_ff[dates, ]
  X_ff_eqw <- X_ff_eqw[dates, ]
  X <- X[dates, ]
  X_bm <- X_bm[dates, ]
  X_eqw <- X_eqw[dates, ]
  
  
  
  
  plot( as.simTS(na.omit(cbind( X_eqw, X_ff_eqw))), plot.type = "single" )
  
  
  
  
  # --------------------------------------------------------------------------
  # Sensor
  # --------------------------------------------------------------------------
  
  BBSSensor <- BBS$new()
  BBSSensor$setCtrl()
  BBSSensor$spec$minphase = 3
  BBSSensor$spec$mincycle = 8
  BBSSensor$spec$theta = 0.1
  BBSSensor$spec$e = 0
  BBSSensor$spec$language = "R"
  BBSSensor$spec$logarithmic = FALSE
  BBSSensor$spec$k_peak = 2
  BBSSensor$spec$k_trough = 2
  BBSSensor$spec$l_peak = 2
  BBSSensor$spec$l_trough = 2
  
  # BBT <- BBTurbulence$new()
  BBT <- BBTEWMA$new()
  BBT$setCtrl( method = "base", universe = "dm" )
  BBT$spec$universe <- "robeco"
  BBT$spec$name <- "BBTulence_base_robeco"
  BBT$spec$BBS <- BBSSensor
  # BBT$spec$iso <- colnames(X)
  # BBT$data$X <- X
  # BBT$data$X_bm <- X_bm
  BBT$spec$iso <- colnames(X_ff)
  BBT$data$X <- X_ff
  BBT$data$X_bm <- X_ff_eqw
  
  # debugonce( BBT$computeSignalInsample )
  BBT$computeSignalInsample()
  
  
  plot( BBT$signal$insample )
  plot( BBT$signal$insample[ ,"delta"] )
  abline( h = 0 )
  
  plot( BBT$signal$insample[ ,c("bear_scl", "bull_scl")], plot.type = "single" )
  plot( BBT$signal$insample[ ,c("bear", "bull")], plot.type = "single" )
  
  delta <- BBT$signal$insample[ ,c("bear")] - BBT$signal$insample[ ,c("bull")]
  plot( delta )
  abline( h = 0 )
  
  
  
  
  delta_ema <- ema( BBT$signal$insample[ ,"delta"], alpha = 0.1 )
  # delta_ema <- ema( delta, alpha = 1 )
  sig <- (1 + sign(delta_ema)) / 2
  
  plot( delta_ema )
  abline( h = 0)
  plot( sig )
  
  Y <- X_eqw
  test <- signalTesting.byTrading( X = Y,
                                   sig = sig, 
                                   tc = 0.004,
                                   n_lag = 1 )
  X_tmp <- na.omit( cbind( Y, bt = test ) )
  
  plot( as.simTS(X_tmp) )
  descStats(X_tmp)
  
  
  
  # --------------------------------------------------------------------------
  # Sensor - daily data
  # --------------------------------------------------------------------------
  
  BBSSensor <- BBS$new()
  BBSSensor$setCtrl()
  
  # BBT <- BBTurbulence$new()
  BBT <- BBTEWMA$new()
  BBT$setCtrl( method = "base", universe = "dm" )
  BBT$spec$universe <- "robeco"
  BBT$spec$name <- "BBTulence_base_robeco"
  BBT$spec$BBS <- BBSSensor
  # BBT$spec$iso <- colnames(X)
  # BBT$data$X <- X
  # BBT$data$X_bm <- X_bm
  BBT$spec$iso <- colnames(X_ff)
  BBT$data$X <- X_ff
  BBT$data$X_bm <- X_ff_eqw
  
  # debugonce( BBT$computeSignalInsample )
  BBT$computeSignalInsample()
  
  delta <- BBT$signal$insample[ ,"delta"]
  plot( delta )
  abline( h = 0 )
  
  
  # delta_m <- aggMonthly( delta, compounding = "continuous" )
  delta_ema <- ema( delta, alpha = 1 )
  sig <- (1 + sign(delta_ema)) / 2
  
  plot( delta_ema )
  abline( h = 0)
  plot( sig )
  
  Y <- X_ff_eqw
  test <- signalTesting.byTrading( X = Y,
                                   sig = sig, 
                                   tc = 0.004,
                                   n_lag = 1 )
  X_tmp <- na.omit( cbind( Y, bt = test ) )
  
  plot( as.simTS(X_tmp) )
  descStats(X_tmp)
  
  
  
  
  #########
  
  
  X_eval <- seq(from = min(delta), to = max(delta), length.out = 10)
  W <- weightsFun( data = list(X_train = delta), 
                   x_eval = X_eval,
                   method = "l1" )
  
  dim(W)
  head(W)
  
  
  Mu <- t(X) %*% W
  Mu  
  barplot( Mu[1, ] )
  barplot( Mu[10, ] )
  
  
  
  DS <- dirichletSampling( Y_train = X_ff_eqw,
                           X_train = delta,
                           sclfct = 1,
                           weights_fun = "kernel" )
  ldens <- lapply( DS, densFUN )
  plot.ldensity( ldens )
  
  
  
  
  #####
  
  
  KM <- DAA::kmeans.adj( x = delta, centers = 20 )
  KM$cluster[which(KM$cluster > 15), ] <- 15
  KM$cluster[which(KM$cluster < 6), ] <- 6
  KM$cluster <- KM$cluster - 5
  
  wmat <- X[rownames(KM$cluster), ] * 0
  for ( today in rownames(KM$cluster) ) {
    wmat[today, KM$cluster[today, ]] <- 1
  }
  
  head(wmat)
  
  plot( KM$cluster )
  
  
  
  sim_bt <- simPortfolio( X = X, wghts = wmat, fc = 0, vc = 0 )  
  
  sim <- na.omit( cbind( X_bm, sim_bt, X ) )
  
  descStats(sim)
  
  colors <- c(1, "magenta", fBasics::divPalette(n = ncol(X), "RdYlGn"))
  plot( as.simTS(sim), col = colors )
  
  weightsBarPlot( wmat )
  plot( wmat, plot.type = "multiple" )
  
  
  
  
  
  