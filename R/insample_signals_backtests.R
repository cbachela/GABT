  
  
  ############################################################################
  ### BACKTEST - PPP ON HEURISTIC BBTURBULENCE STRATEGIES
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     16.10.2020
  # First version:    16.10.2020
  # --------------------------------------------------------------------------

  # DESCRIPTION:
  #
  # We run a series of univariate binary backtests, i.e., either be 
  # fully invested in an index or not at all, using the signalTesting.byTrading
  # function. The binary state variable is defined by the sign of the
  # smoothed BBTurbulence signal, where the smoothing is done with an 
  # exponentially weighted moving average. We loop over a grid of half-live 
  # parameters generating a backtest for each parameter setup and feed the 
  # backtest simulations to an out-of-sample PPP optimization for automated
  # parameter tuning.
  # Note: 
  # The heuristic strategies are excellent before and terrible after costs.
  
  
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
  # COMPUTE BBTurbulence SIGNAL
  # --------------------------------------------------------------------------
  
  BBT <- loadSensor(sensor_name = "bbturbulence")
  BBTD <- BBT$copy()
  BBTD$spec$minphase <- minphase_d
  BBTD$spec$mincycle <- mincycle_d
  BBTD$spec$l_peak <- BBTD$spec$l_trough <- 
    BBTD$spec$k_peak <- BBTD$spec$k_trough <- minphase_d
  # debugonce(BBTD$computeSignalInsample)
  BBTD$computeSignalInsample()
  mdd <- BBTD$signal$insample
  mdd <- mdd[isWeekday(time(mdd)), ]
  
  # mdd <- BBT$signal$base
  # mdd <- mdd[isWeekday(time(mdd)), ]
  
  
  # --------------------------------------------------------------------------
  # HEURISTIC STRATEGIES
  # --------------------------------------------------------------------------
  
  X_bm <- BBT$data$X_bm[isWeekday(time(BBT$data$X_bm)), ]
  # TD <- trainingData( Y_train = aggWeekly(X_bm, day = "Tue", "discrete"),
  #                     X_train = mdw[ ,"delta"] )
  TD <- trainingData( Y_train = X_bm,
                      X_train = mdd[ ,"delta"] )
  
  # Loop over grid of half-live parameter
  # tau_vec <- round( seq(from = 0.001, to = 0.2, length.out = 50), 3 )
  tau_vec <- round( seq(from = 0.01, to = 1, length.out = 50)^2, 5 )
  n_lag <- 1
  tc <- 0.004
  
  
  # FUN <- function(x) { (1 + sign(ema(signal, alpha = x))) / 2 } 
  FUN <- function(x) { ema( TD$X_train, alpha = x ) }
  lSignal_ema <- lapply( tau_vec, FUN = FUN )
  sigs1_ema <- do.call( cbind, lSignal_ema )
  th <- 0
  sigs1 <- apply( sigs1_ema, 2, function(x) { (1 + sign(x - th)) / 2 } )
  
  sims1 <- signalTesting.byTrading( X = TD$Y_train,
                                    sig = sigs1,
                                    n_lag = n_lag,
                                    tc = tc )
  colnames(sims1) <- paste0("tau=", tau_vec)
  sims1_nc <- signalTesting.byTrading( X = TD$Y_train,
                                       sig = sigs1,
                                       n_lag = n_lag,
                                       tc = 0 )
  colnames(sims1_nc) <- paste0("tau=", tau_vec)
  stats1 <- descStats( sims1 )
  stats1_nc <- descStats( sims1_nc )
  
  
  # Alternative 2, 
  # use sigmoidal transformation of signal to get a weight between zero and one.
  FUN <- function(x) { 1 / (1 + exp(-ema(TD$X_train, alpha = x) * 100)) }
  lSig <- lapply( tau_vec, FUN = FUN )
  sigs2 <- do.call( cbind, lSig )
  sims2 <- signalTesting.byTrading( X = TD$Y_train,
                                    sig = sigs2,
                                    n_lag = n_lag,
                                    tc = tc )
  colnames(sims2) <- paste0("tau=", tau_vec)
  sims2_nc <- signalTesting.byTrading( X = TD$Y_train,
                                       sig = sigs2,
                                       n_lag = n_lag,
                                       tc = 0 )
  colnames(sims2_nc) <- paste0("tau=", tau_vec)
  stats2 <- descStats( sims2 )
  stats2_nc <- descStats( sims2_nc )
  
  
  sim <- cbind( sims1, sims2)
  sig <- cbind( sigs1, sigs2 )
  
  
  
  colors <- fBasics::divPalette(n = ncol(sims1), "RdYlGn")
  plot( log(cumulated(sims1, "discrete")), plot.type = "single", col = colors )
  lines( log(cumulated(TD$Y_train, "discrete")) )
  
  plot( log(cumulated(sims1_nc, "discrete")), plot.type = "single", col = colors )
  lines( log(cumulated(TD$Y_train, "discrete")) )
  
  
  
  
  
  
  
  # Alternative 3, 
  # apply smoothing on bear/bull Mahalanobis distance 
  # before taking the difference
  FUN <- function(x) { ema(lMDBB$MD[ ,"bear_scl"], x) - ema(lMDBB$MD[ ,"bull_scl"], x) }
  lSig <- lapply( tau_vec, FUN = FUN )
  sigs3 <- do.call( cbind, lSig )
  sigs3 <- apply( sigs3, 2, function(x) { (1 + sign(x)) / 2 } )
  sims3 <- signalTesting.byTrading( X = y,
                                    sig = sigs3,
                                    n_lag = n_lag,
                                    tc = tc )
  colnames(sims3) <- paste0("tau=", tau_vec)
  
  stats3 <- descStats(sims3)
  
  
  
  # Note:
  # Alternatives 1 and 3 are equivalent. 
  # Alternative 2 works less good than 1, with and without costs.
  
  
  
  # Alternative 4:
  # smooth the data before computing the distances
  
  FUN <- function(x) { ema(signal_roll, alpha = x) } 
  lSignal_ema <- lapply( tau_vec, FUN = FUN )
  sigs4_ema <- do.call( cbind, lSignal_ema )
  th <- 0
  sigs4 <- apply( sigs4_ema, 2, function(x) { (1 + sign(x - th)) / 2 } )
  
  sims4 <- signalTesting.byTrading( X = y,
                                    sig = sigs4,
                                    n_lag = n_lag,
                                    tc = tc )
  colnames(sims4) <- paste0("tau=", tau_vec)
  sims4_nc <- signalTesting.byTrading( X = y,
                                       sig = sigs4,
                                       n_lag = n_lag,
                                       tc = 0 )
  colnames(sims4_nc) <- paste0("tau=", tau_vec)
  
  stats4 <- descStats( sims4 )
  stats4_nc <- descStats( sims4_nc )
  
  
  
  
  # Alternative 5:
  # smooth weights
  
  # th <- 0
  # sig5 <- (1 * sign(signal - th)) / 2
  sig5 <- 1 / (1 + exp(-signal * 100))
  plot(sig5)
  FUN <- function(x) { ema(sig5, alpha = x) } 
  lSignal <- lapply( tau_vec, FUN = FUN )
  sigs5 <- do.call( cbind, lSignal )
  sims5 <- signalTesting.byTrading( X = y,
                                    sig = sigs5,
                                    n_lag = n_lag,
                                    tc = tc )
  colnames(sims5) <- paste0("tau=", tau_vec)
  sims5_nc <- signalTesting.byTrading( X = y,
                                       sig = sigs5,
                                       n_lag = n_lag,
                                       tc = 0 )
  colnames(sims5_nc) <- paste0("tau=", tau_vec)
  
  stats5 <- descStats( sims5 )
  stats5_nc <- descStats( sims5_nc )
  
  
  colors <- fBasics::divPalette(n = ncol(sims5), "RdYlGn")
  plot( log(cumulated(sims5, "discrete")), plot.type = "single", col = colors )
  lines( log(cumulated(y, "discrete")) )
  
  
  
  
  
  # Analyze
  
  range(stats1$stats["cumret", ])
  range(stats2$stats["cumret", ])
  range(stats3$stats["cumret", ])
  range(stats4$stats["cumret", ])
  range(stats5$stats["cumret", ])
  
  range(stats1_nc$stats["cumret", ])
  range(stats2_nc$stats["cumret", ])
  range(stats4_nc$stats["cumret", ])
  range(stats4_nc$stats["cumret", ])
  range(stats5_nc$stats["cumret", ])
  
  
  
  
  # Alternative 6
  # Run BBQ on cumulated signal
  
  signal <- Object$signal$base[ ,"md_delta"]  # out-of-sample signal
  # sig_level <- cumulated(signal / 100, "discrete")
  sig_level <- timeSeries( cumsum(signal), time(signal) )
  
  # signal <- signal_roll
  # sig_level <- log(cumulated(signal / 100, "discrete"))
  # 
  # BCP <- loadSensor( sensor_name = "bcpolz" )
  # # signal <- BCP$getSignal() - 1
  # signal <- BCP$signal$base$mumat[ ,1]
  # sig_level <- log(cumulated(signal / 100, "discrete"))
  
  
  theta <- quantile(signal / 100, 0.9)
  theta <- quantile(returns(sig_level, "discrete"), 0.9)
  
  # In-sample
  BBS <- bbs( X = sig_level,
              minphase = 5,
              mincycle = 21 * 3,
              theta = theta,
              logarithmic = FALSE,
              e = 0 )
  
  # theta <- quantile(signal / 100, 0.9) * 100 * 10
  theta <- quantile(returns(sig_level, "discrete"), 0.9) * 10
  setpar_filtering_alg(tr_bull = theta,
                       tr_bear = theta)
  bull_lt <- run_filtering_alg(index =  sig_level )
  BBS <- timeSeries( matrix(sign(as.numeric(bull_lt)-0.5), ncol = 1), time(sig_level) ) 
  
  plot(BBS)
  
  
  sig_bbs_is <- (1 + BBS) / 2
  test <- signalTesting.byTrading( X = y,
                                   # sig = sig_bbs_is,
                                   sig = BBS,
                                   n_lag = n_lag,
                                   tc = tc )
  X_tmp <- na.omit(cbind(bm = y, test = test))
  descStats(X_tmp)
  plot(as.simTS(X_tmp))
  
  
  
  # Out-of-sample
  
  dates <- rownames(signal)[ -c(1:252) ]
  sig_bbs_oos <- signal[dates, ] * NA
  
  tic <- Sys.time()
  for ( today in dates ) {
    
    Y <- sig_level[rownames(sig_level) <= today, ]
    BBS_tmp <- bbs( X = Y,
                    minphase = 5,
                    mincycle = 21 * 3,
                    # theta = theta,
                    logarithmic = FALSE,
                    e = 0 )
    sig_bbs_oos[today, ] <- tail(BBS_tmp, 1)
    # theta <- 2
    # setpar_filtering_alg(tr_bull = theta,
    #                      tr_bear = theta)
    # bull_lt <- run_filtering_alg(index =  Y )
    # sig_bbs_oos[today, ] <- ifelse( isTRUE(tail(bull_lt, 1)), 1, -1 )
    
  }
  (toc <- Sys.time() - tic)
  
  test <- signalTesting.byTrading( X = y,
                                   sig = (1 + sig_bbs_oos) / 2,
                                   # sig = sig_bbs_oos,
                                   n_lag = n_lag,
                                   tc = tc )
  X_tmp <- na.omit(cbind(bm = y, test = test))
  descStats(X_tmp)
  plot(as.simTS(X_tmp))
  
  
  plot( sig_bbs_oos )
  
  
  
  
  
  par( mfrow = c(2, 1) )
  
  plot( sig_level )
  abline( v = time(BBS)[ which(BBS == -1) ], col = 2 )
  abline( v = time(BBS)[ which(BBS == 1) ], col = 3 )
  lines( sig_level, lwd = 2, col = 4 )
  
  plot( log(X_level) )
  abline( v = time(BBS)[ which(BBS == -1) ], col = 2 )
  abline( v = time(BBS)[ which(BBS == 1) ], col = 3 )
  lines( log(X_level), lwd = 2, col = 4 )
  
  dev.off()
  
  
  plot( sig_level )
  abline( v = time(sig_bbs_oos)[ which(sig_bbs_oos == -1) ], col = 2 )
  abline( v = time(sig_bbs_oos)[ which(sig_bbs_oos == 1) ], col = 3 )
  lines( sig_level, lwd = 2, col = 4 )
  
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Run momPortfolio backtest on simulations
  # --------------------------------------------------------------------------
  
  require(AAA)
  require(SIDRC)
  wd <- "R:/Asset_Management/Research_Projects/Multi_Asset/SID/BacktestSIDRC/"
  source("R:/Asset_Management/Research_Projects/Multi_Asset/SID/BacktestSID/Source/BacktestSID_returns.R" )
  source("R:/Asset_Management/Research_Projects/Multi_Asset/SID/BacktestSID/Source/BacktestSID_portfolio.R" )
  source("R:/Asset_Management/Research_Projects/Multi_Asset/SID/BacktestSID/Source/BacktestSID_zzz.R" )
  source( paste0(wd, "Source/BacktestSIDRC_zzz.R") )
  
  
  
  Constraints <- constraints( selection = colnames(sim) )
  addConstraint(Constraints) <- budgetConstraint()
  addConstraint(Constraints) <- boxConstraint()
  
  
  BT2 <- BacktestSIDRC$new()
  # debugonce(aggWeekly)
  # BT$data$X <- aggWeekly( sim, compounding = "discrete", day = "Wed" )
  BT2$data$X_est <- sim
  BT2$data$X_sim <- sim
  rebdates <- rownames(BT2$data$X_est)[ seq(from = 252, to = nrow(BT2$data$X_est), by = 1) ]
  BT2$setCtrl( width = 252,
               alpha = 0.001,
               tau = 26 * 5,
               riskaversion = 1,
               Constraints = Constraints,
               NAMES_ADAPTIVE_EST = colnames(sim),
               # wd_warehouse = paste0(wd, "Momentum/MOFM/waRehouse/"),
               rebdates = rebdates )
  BT2$runSID( method = "mom" )
  
  
  # Extract weights of single asset investment strategy
  
  theta_mat <- BT2$output$weights
  # strat_array <- list2array( lWmat_strat )
  strat_mat <- sig
  
  wmat <- theta_mat[ ,1] * NA
  for ( today in rownames(theta_mat) ) {
    
    coeff <- as.numeric( theta_mat[today, ] )
    wmat_tmp <- strat_mat[today, ]
    wmat[today, ] <- coeff %*% t(wmat_tmp)
  }
  
  
  plot( wmat )
  weightsBarPlot(tail(wmat, 200))
  weightsBarPlot(tail(theta_mat, 200))
  
  boxplot( as.data.frame(theta_mat) )
  
  
  sim_mom <- signalTesting.byTrading( X = y,
                                      sig = wmat,
                                      n_lag = 1,
                                      tc = 0.004 )
  sim_mom_nc <- signalTesting.byTrading( X = y,
                                         sig = wmat,
                                         n_lag = 1,
                                         tc = 0 )
  
  
  sim_all <- na.omit( cbind( bm = X_bm,
                             mom = sim_mom, 
                             mom_nc = sim_mom_nc, 
                             sim ) )
  
  
  plot( as.simTS(sim_all) )
  plot( as.simTS(sim_all[ ,1:20]) )
  
  plot( log(cumulated(sim_all, "discrete")), plot.type = "single" )
  
  
  
  # Statistics
  stats_fields <- c("cumret", "sds", "sharpe", "maxDD")
  lStats <- descStats( sim_all )
  stats <- t( lStats$stats[stats_fields, ] )
  
  stats
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # PPP
  # --------------------------------------------------------------------------
  
  Names <- c("Cash", "BM")
  strategy_names <- colnames(sim)
  strat_mat <- sig
  lWmat <- list()
  for ( i in 1:ncol(strat_mat) ) {
    tmp <- cbind(1 - strat_mat[ ,i], strat_mat[ ,i])
    colnames(tmp) <- Names
    lWmat[[i]] <- tmp
  }
  names(lWmat) <- strategy_names
  w_array <- list2array( L = lWmat )
  Constraints <- constraints(selection = Names)
  addConstraint(Constraints) <- boxConstraint()
  addConstraint(Constraints) <- budgetConstraint()
  
  BT3 <- BacktestSIDRC$new()
  X_tmp <- cbind(X_bm * 0, X_bm)
  colnames(X_tmp) <- Names
  BT3$data$X_est <- X_tmp
  BT3$data$X_sim <- X_tmp
  BT3$data$w_array <- w_array
  dates <- intersect(rownames(BT3$data$X_est), rownames(sig))
  rebdates <- dates[ seq(from = 252, to = length(dates), by = 1) ]
  BT3$setCtrl( width = 26 * 5,
               alpha = 0.05,
               riskaversion = 1,
               lambda_1 = 1,
               lambda_2 = 0,
               lambda_3 = 1,
               strategy_names = strategy_names,
               n_lag = 0,
               tc_vec = setNames( c(0, 0.04), Names ),
               # tc_vec = setNames( c(0, 0.004), Names ),
               hide_cash = FALSE,
               substract_wbar = TRUE,
               NAMES_ADAPTIVE_EST = Names,
               Constraints = Constraints,
               # wd_warehouse = paste0(wd, "Momentum/MOFM/waRehouse/"),
               rebdates = rebdates )
  
  # debugonce(BT3$pppPortfolio)
  BT3$runSID( method = "ppp" )
  
  
  plot( BT3$output$weights )
  weightsBarPlot( tail(BT3$output$weights, 100) )
  colors <- fBasics::divPalette(n = ncol(BT3$output$theta), "RdYlGn")
  weightsBarPlot( tail(BT3$output$theta, 100), col = colors )
  
  
  plot( tail(sigs1, 100), plot.type = "single", col = colors )
  
  
  
  
  bt_ppp <- signalTesting.byTrading(X = X_bm,
                                    sig =  BT3$output$weights[ ,"BM"],
                                    n_lag = 1,
                                    tc = 0.004)
  X_tmp <- na.omit(cbind(X_bm, bt_ppp))
  
  plot( log(cumulated(X_tmp, "discrete")), plot.type = "single" )
  # abline(v = time(sig)[which(sig == 1)], col = 3)
  # abline(v = time(sig)[which(sig == 0)], col = 2)
  lines(log(cumulated(X_tmp[ ,1], "discrete")), col = 1)
  lines(log(cumulated(X_tmp[ ,2], "discrete")), col = 4)
  
  descStats(X_tmp)
  
  
  sim_all <- na.omit( cbind(bm = X_bm, 
                            ppp = test,
                            mom = sim_mom, 
                            mom_nc = sim_mom_nc) )
  plot( as.simTS(sim_all) )
  
  
  
  
  
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # SAVE
  # --------------------------------------------------------------------------
  
  env <- new.env()
  
  env$sim_all <- sim_all
  env$BT_MOM <- BT2
  env$BT_PPP <- BT3
  
  wd_save <- "R:/Asset_Management/Research_Projects/Multi_Asset/DAA/DAARC/..."
  saveRDS( env, file = paste0(wd_save, "....rds") )
  
  
  
  
  
  
  
  
  
  
