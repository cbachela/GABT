  
  
  ############################################################################
  ### GOOD AND BAD TURBULENCE - OUT-OF-SAMPLE ANALYSIS - MULTI ASSET
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     25.10.2020
  # First version:    25.10.2020
  # --------------------------------------------------------------------------
  
  
  # --------------------------------------------------------------------------
  # REQUIRE
  # --------------------------------------------------------------------------
  
  require(MAE)
  require(SIDRC)
  require(DAARC)
  
  wd <- "R:/Asset_Management/Research_Projects/Multi_Asset/DAA/DAARC/Sensor_Testing/Univariate/Good_And_Bad_Turbulence/"
  source( paste0(wd, "Source/source_functions.R") )
  source( "R:/Asset_Management/Research_Projects/Multi_Asset/SID/BacktestSID/Source/BacktestSID_returns.R" )
  source( "R:/Asset_Management/Research_Projects/Multi_Asset/SID/BacktestSID/Source/BacktestSID_portfolio.R" )
  source( "R:/Asset_Management/Research_Projects/Multi_Asset/SID/BacktestSID/Source/BacktestSID_zzz.R" )
  source( "R:/Asset_Management/Research_Projects/Multi_Asset/SID/BacktestSIDRC/Source/BacktestSIDRC_zzz.R" )
  
  
  
  
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
  
  
  # Weekly BBS indicator in-sample
  BBS <- loadSensor("bbs", b_new = TRUE)
  BBSW <- BBS$copy()
  BBSW$setCtrl()
  BBSW$spec$minphase <- minphase_w
  BBSW$spec$mincycle <- mincycle_w
  BBSW$spec$ccy <- "Local"
  BBSW$spec$width <- 52 * 3
  BBSW$spec$name <- "bbs_w"
  BBSW$data <- list()
  BBSW$updateData()
  BBSW$data$X <- aggWeekly( X = BBSW$data$X,
                            day = "Tue",
                            compounding = "discrete" )
  BBSW$data$X_bm <- aggWeekly( X = BBSW$data$X_bm,
                               day = "Tue",
                               compounding = "discrete" )
  BBSW$signal <- list()
  BBSW$updateSignal()
  states_is_w <- (1 + BBSW$signal$insample) / 2
  
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # MULTI-ASSET DATASET
  # --------------------------------------------------------------------------
  
  MAE <- mae( b_update = FALSE )
  X_est <- MAE$X_est
  
  
  # --------------------------------------------------------------------------
  # 1) DEFINE TURBULENCE MEASURES
  # --------------------------------------------------------------------------
  
  # --------------------------------------------------------------------------
  # Construct absolute and relative Mahalanobis distances
  # --------------------------------------------------------------------------
  
  # BBT <- loadSensor( sensor_name = "bbturbulence_scl2_dm", 
  #                    b_new = FALSE )
  # 
  # # Using weekly returns
  # BBTW <- BBT$copy()
  # BBTW$spec$name <- paste0(BBT$spec$name, "_w")
  # BBTW$spec$width <- floor(BBT$spec$width / 5)
  # BBTW$signal <- list()
  # BBTW$data$X <- aggWeekly(BBT$data$X, day = "Tue", compounding = "discrete")
  # BBTW$data$X_bm <- aggWeekly(BBT$data$X_bm, day = "Tue", compounding = "discrete")
  # BBTW$data$wmat <- applyRoll( BBT$data$wmat, 
  #                              Width = 5, 
  #                              FUN = function(X) { apply(X, 2, mean) },
  #                              charvec = rownames(BBTW$data$X) )
  # BBTW$spec$minphase <- minphase_w
  # BBTW$spec$mincycle <- mincycle_w
  # BBTW$spec$l_peak <- BBTW$spec$l_trough <- 
  #   BBTW$spec$k_peak <- BBTW$spec$k_trough <- minphase_w
  # # debugonce(BBTW$update)
  # BBTW$updateSignal()
  # mdw <- BBTW$signal$scl2
  # mdw <- mdw[ ,c("md_bear", "md_bull", "md_all", "md_delta")]
  # mdw_scl <-  apply(mdw, 2, function(x) { x / max(x) } )
  
  BBTW <- loadSensor( sensor_name = "bbturbulence_scl2_dm_w", 
                      b_new = FALSE )
  BBTW$updateSignal()
  mdw <- BBTW$signal$scl2[ ,c("md_all", "md_delta")]
  
  
  # --------------------------------------------------------------------------
  # 2) CONDITIONAL PERFORMANCE
  # --------------------------------------------------------------------------
  
  # --------------------------------------------------------------------------
  # Posterior means
  # --------------------------------------------------------------------------
  
  TD <- trainingData( Y_train = X_est[ ,-1],
                      X_train = mdw )
  lDS <- list()
  
  for ( j in 1:ncol(TD$Y_train) ) {
    tmp <- list()
    for ( k in 1:ncol(TD$X_train) ) {
      tmp[[k]] <- dirichletSampling( Y_train = TD$Y_train[ ,j],
                                     X_train = TD$X_train[ ,k],
                                     n_lag = 0,
                                     # X_eval = X_eval,
                                     weights_fun = "kernel",
                                     # centers = 5,
                                     correct_bias = FALSE,
                                     sclfct = 1 )
    }
    names(tmp) <- colnames(TD$X_train)
    lDS[[j]] <- tmp
  }
  names(lDS) <- colnames(TD$Y_train)
  
  
  
  # from <- quantile( unlist(lDS), 0 )
  # to <- quantile( unlist(lDS), 1 )
  from <- min( unlist(lDS) ) * 1.3
  to <- max( unlist(lDS) ) * 1.3
  n <- 999
  lldens <- lapply( lDS[["Equity_DM"]], FUN = function(x) { lapply(x, density, from = from, to = to, n = n) } )
  for ( i in 1:length(lldens) ) {
    plot.ldensity( lldens[[i]], 
                   main = paste0("Posterior means conditional on ", 
                                 names(lldens)[i]) )
  }
 
  
  
  
  lMu <- list()
  for ( j in 1:length(lDS) ) {
    tmp <- lapply( lDS[[j]], FUN = function(x) { unlist(lapply(x, FUN = mean)) } )
    lMu[[j]] <- do.call( cbind, tmp )
  }
  
  colors <- rev(fBasics::divPalette(n = nrow(lMu[[1]]), "RdYlGn"))
  colors <- rev(fBasics::divPalette(n = nrow(lMu[[1]]), "Spectral"))
  j <- 4; barplot( lMu[[j]], beside = TRUE, col = colors )
  names(lMu)
  
  
  
  
  
  
  KCDH <- kcdens( Y_train = TD$Y_train[ ,"Equity_DM"],
                  X_train = TD$X_train[ ,"md_delta"],
                  kcd_fun = "hyndman",
                  n_lag = 1 )
  plot(KCDH$fit)

  

  KCDLR <- kcdens( Y_train = TD$Y_train[ ,"Equity_DM"],
                   X_train = TD$X_train[ ,"md_delta"],
                   # X_eval = X_eval,
                   n_lag = 0,
                   kcd_fun = "liracine",
                   bw_method = "normal-reference" )
  
  matplot( x = KCDLR$y_eval,
           y = t(KCDLR$cdens),
           type = "l")
  
  mu_kcdlr <- apply(KCDLR$cdens, 1, function(p) { p %*% KCDLR$y_eval } )
  names(mu_kcdlr) <- rownames(KCDLR$x_eval)
  barplot( mu_kcdlr )
  
  
  cbind( KCDLR$x_eval, mu_kcdlr )
  
  
  
  # Logistic regression
  
  TD <- trainingData( Y_train = states_is_w,
                      # X_train = mdw_scl[ ,c("md_bear", "md_bull")],
                      X_train = mdw_scl[ ,c("md_all", "md_delta")],
                      n_lag = 0 )
  
  reg <- regression( Y_train = TD$Y_train,
                     X_train = TD$X_train,
                     type = "logit" )
  summary(reg$reg)
  
  plot( reg$z, reg$y )
  plot( reg$y )
  
  
  # ROC Curves
  ROC <- roc( y_hat = reg$y, 
              y_true = TD$Y_train, 
              th = 0.5 )
  ROC
  
  th_vec <- round( seq(from = 0, to = 1, length.out = 10), 2 )
  lROC <- lapply(th_vec, FUN = function(x) { roc(th = x, y_hat = reg$y, y_true = TD$Y_train)})
  ROC <- do.call( rbind, lROC )
  rownames(ROC) <- th_vec
  ROC
  
  ROC <- pROC::roc(response = as.numeric(TD$Y_train),
                   # predictor = as.numeric(reg$y) * 0 + rnorm(nrow(reg$y)),
                   predictor = as.numeric(reg$y),
                   auc = TRUE)
  plot(ROC)
  ROC$auc
  
  
  
  
  # OLS regression
  
  TD <- trainingData( Y_train = BBS$data$X_bm,
                      # X_train = mdw_scl[ ,c("md_bear", "md_bull")],
                      X_train = mdw_scl[ ,c("md_all", "md_delta")],
                      n_lag = 1 )
  
  reg <- regression( Y_train = TD$Y_train,
                     X_train = TD$X_train,
                     type = "ols" )
  summary(reg$reg)
  
  
  
  # Elastic net regression
  
  # Differently smoothened covariates
  tau_vec <- round( seq(from = 0.99, to = 0.01, length.out = 50), 2 )
  tau_vec <- c(0.9999, 0.999, tau_vec)
  halflive <- -log(2) / log(tau_vec)
  FUN <- function(x) { ema( mdw[ ,"md_delta"], alpha = x ) }
  lSignal_ema <- lapply( tau_vec, FUN = FUN )
  sig_ema <- do.call( cbind, lSignal_ema )
  colnames(sig_ema) <- paste0("tau", tau_vec)
  
  TD <- trainingData( X_train = sig_ema,
                      Y_train = states_is_w,
                      # BBTW$data$X_bm,
                      n_lag = 0 )
  
  reg <- regression( Y_train = TD$Y_train,
                     X_train = TD$X_train,
                     type = "elnet",
                     elnet_alpha = 1 )  # elnet_alpha = 1 --> LASSO, i.e., L1 penalty
  
  summary(reg$reg)
  
  coef(reg$reg)
  plot( reg$coeffmat )
  
  plot( halflive )
  plot( x = halflive, y = as.numeric(reg$coeffmat)[-1] )
  

  
  idx <- which( as.numeric(reg$coeffmat)[-1] > 0 )
  reg <- regression( Y_train = TD$Y_train,
                     X_train = TD$X_train[ ,idx],
                     type = "logit" )
  summary(reg$reg)
  
  plot( reg$z, reg$y )
  
  plot( reg$y )
  
  plot( sig_ema[ ,idx], plot.type = "single" )
  plot( apply(sig_ema[ ,idx], 1, mean) )
  abline(h = 0)
  
  
  
  
  # --------------------------------------------------------------------------
  # 3) BACKTEST
  # --------------------------------------------------------------------------

  
  # --------------------------------------------------------------------------
  # INSTANTIATE SIDRC CLASS, UPDATE DATA AND SIGNALS
  # --------------------------------------------------------------------------
  
  SID <- SIDRC$new()
  SID$setCtrl()
  SID$setConstraints()
  # debugonce(SID$appendData)
  SID$updateData()  # ~~~~~~~~~~~~~~~~ make append
  lData <- SID$data
  SID$data <- lData
  # SID$updateSensors()
  
  
  # --------------------------------------------------------------------------
  # INSTANTIATE BACKTEST CLASS, UPDATE DATA AND SIGNALS
  # --------------------------------------------------------------------------
  
  BT0 <- BacktestSIDRC$new()
  signals <- TD$X_train[ ,c("md_bear", "md_bull", "md_all", "md_delta")]
  BT0$data <- list( X_est = SID$data$X_est,
                    X_sim = SID$data$X_sim,
                    signals = signals )
  
  
  # --------------------------------------------------------------------------
  # LOOP OVER SENSORS
  # --------------------------------------------------------------------------
  
  rebdates <- intersect( rownames(BT0$data$X_est), rownames(BT0$data$signals) )
  rebdates <- rebdates[ -c(1:103) ]
  BT0$setCtrl( width = Inf,
               alpha = 0.05,
               riskaversion = 1,
               n_lag = 0,
               NAMES_ADAPTIVE_EST = getConstraints(SID$constraints, "selection"),
               wd_warehouse = paste0(wd, "waRehouse/"),
               rebdates = rebdates )
  
  # sensor_names <- BT$spec$SIGNAL_NAMES
  sensor_names <- colnames(signals)
  b_relentropy <- c(TRUE, FALSE)
  
  for ( k in seq(along = sensor_names) ) {
    for ( j in seq(along = b_relentropy) ) {
      
      if ( isTRUE( b_relentropy[j] ) ) {
        name <- paste0( "pgck_minklic(", sensor_names[k], ")")
      } else {
        name <- paste0( "pgck(", sensor_names[k], ")")
      }
      
      BT_old <- try( loadBacktest( wd = BT$spec$wd_warehouse,
                                   name = name ) )
      if ( !inherits(BT_old, "try-error") ) {
        BT <- BT_old$copy()
        BT_old$data <- BT0$data
        BT_old$spec$rebdates <- BT0$spec$rebdates
      } else {
        BT <- BT0$copy()
        BT$output <- list()
        BT$spec$signal_names <- sensor_names[k]
        BT$spec$b_relentropy <- b_relentropy[j]
        BT$spec$name <- name
      }
      BT$runSID( method = "pgck" )
      
      BT$save( wd = BT$spec$wd_warehouse,
               name = name,
               without_data = FALSE )
      
    }
  }
  
  
  p <- 1
  FUN <- function(x)
  {
    y <- x * 0
    y[x > 0] <- x[x > 0]
    ans <- y^p / sum(y^p) 
    return(ans)
  }
  tmp <- t( apply( BT$output$mu_post, 1, FUN ) )
  wmat <- timeSeries( tmp, rownames(BT$output$mu_post) )
  sim <- simPortfolio( X = BT$data$X_sim,
                                         wghts = wmat )
  plot( sim )
  
  plot(wmat, plot.type = "single")
  
  plot(BT$output$mu_post, plot.type = "single")
  
  plot(BT$output$weights, plot.type = "single")
  
  
  
  
  # --------------------------------------------------------------------------
  # ANALYZE
  # --------------------------------------------------------------------------
  
  analyze <- function()
  {
    
    fc <- 0.01
    vc <- 0.004
    wd_warehouse <- paste0(wd, "waRehouse/")
    file_names <- list.files( path = wd_warehouse, 
                              pattern = ".rds" )
    portfolio_names <- substr(file_names, 1, nchar(file_names) - 4)
    
    lBT <- lapply( portfolio_names, 
                   loadBacktest, 
                   wd = wd_warehouse, 
                   b_new = TRUE )
    names(lBT) <- portfolio_names
    lSim <- lapply( lBT, FUN = function(x) { try( x$simulate(fc = fc, vc = vc) ) } )
    lSim_novc <- lapply( lBT, FUN = function(x) { try( x$simulate(fc = fc, vc = 0) ) } )
    idx <- unlist( lapply(lSim, FUN = function(x) { !inherits(x, "try-error") } ) )
    sim_strat <- do.call( cbindsimTS, lSim[ idx ] )
    sim_strat_novc <- do.call( cbindsimTS, lSim_novc[ idx ] )
    Names <- portfolio_names[ idx ]
    colnames(sim_strat) <- colnames(sim_strat_novc) <- Names
    sim_strat_w <- aggWeekly(sim_strat, 
                             compounding = "discrete", 
                             day = "Tue")
    sim_strat_novc_w <- aggWeekly(sim_strat_novc, 
                                  compounding = "discrete", 
                                  day = "Tue")
    
    
    stats_fields <- c("cumret", "sds", "sharpe", "maxDD")
    t(descStats(na.omit(sim_strat))$stats[stats_fields, ])
    t(descStats(na.omit(sim_strat_novc))$stats[stats_fields, ])
    
    
    plot( sim_strat )
    plot( sim_strat_novc )
    
    
    weightsBarPlot( tail(lBT[[1]]$output$weights, 300) )
    names(lBT)
    
    head(tail(lBT[[1]]$output$weights, 35))
    
    
    
    # Benchmark
    
    
  }
  
  
  
  
  