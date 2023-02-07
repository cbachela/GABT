  
  
  ############################################################################
  ### BacktetDAA - RPM1
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     29.06.2020
  # First version:    28.05.2020
  # --------------------------------------------------------------------------
  
  
  # This file runs DAA backtests using the RPM methodology.
  
  
  
  require(RP)
  require(slolz)
  require(DAARC)
  wd <- "R:/Asset_Management/Research_Projects/Multi_Asset/DAA/DAARC/BacktestDAA/"
  wd_warehouse <- paste0(wd, "RPM/waRehouse/")
  source( paste0(wd, "Source/custom_functions.R") )
  
  

  
  # --------------------------------------------------------------------------
  # BEARS AND BULLS
  # --------------------------------------------------------------------------
  
  Object <- loadSensor( sensor_name = "bbs" )
  Object$update()
  
  X_bm <- Object$data$X_bm
  X_level <- log(cumulated(X_bm, "discrete"))
  BBS_is <- Object$signal$insample
  states_is <- (1 + BBS_is) / 2
  
  # Compare in- and out-of-sample states
  BBS <- Object$getSignal()
  par( mfrow = c(2, 1) )
  plot(X_level)
  abline(v = time(BBS)[which(BBS == 1)], col = "green")
  abline(v = time(BBS)[which(BBS == 0)], col = "red")
  lines(X_level)
  plot(X_level)
  abline(v = time(BBS_is)[which(BBS_is == 1)], col = "green")
  abline(v = time(BBS_is)[which(BBS_is == -1)], col = "red")
  lines(X_level)
  
  dev.off()
  
  
  
  # ---
  
  mincycle <- 5 * 4 * 9 
  minphase <- 5 * 4 * 3
  k_peak <- k_trough <- l_peak <- l_trough <- minphase * 1
  e <- 0
  
  BBS <- bbs(X = X_level,
             mincycle = mincycle, 
             minphase = minphase,
             k.peak = k_peak,
             l.peak = l_peak,
             k.trough = k_trough,
             l.trough = l_trough,
             e = e)
  states_is <- (1 + BBS) / 2
  
  
  par( mfrow = c(2, 1) )
  
  plot(X_level)
  abline(v = time(BBS)[which(BBS == 1)], col = "green")
  abline(v = time(BBS)[which(BBS == -1)], col = "red")
  lines(X_level)
  
  plot(X_level)
  abline(v = time(BBS_is)[which(BBS_is == 1)], col = "green")
  abline(v = time(BBS_is)[which(BBS_is == -1)], col = "red")
  lines(X_level)
  
  dev.off()
  
  
  
  # Omit last phase
  states <- states_is
  states <- states[ rownames(states) <= "2020-02-12", ]

  
 
  
  
  
  
  # --------------------------------------------------------------------------
  # SIGNALS
  # --------------------------------------------------------------------------
  
  VRP <- loadSensor(sensor_name = "vrp")
  Turbulence <- loadSensor(sensor_name = "turbulence")
  BBTurbulence <- loadSensor(sensor_name = "bbturbulence")
  BCPOLZ <- loadSensor(sensor_name = "bcpolz")
  Eigen <- loadSensor(sensor_name = "eigen")
  Momentum <- loadSensor(sensor_name = "momentum")
  OFRFSI <- loadSensor(sensor_name = "ofrfsi")
  
  
  md_delta <- BBTurbulence$signal$base[ ,"md_delta"]
  md_delta_ema <- ema(md_delta, alpha = 0.02)
  md_delta_sigmoid <- 1 / (1 + exp(-md_delta * 100))
  
  plot( md_delta_sigmoid )
  
  
  signals <- na.omit( cbind( vrp_signal = VRP$getSignal(),
                             vrp_raw = VRP$signal$base[ ,"vrp_raw"],
                             turbulence = Turbulence$getSignal(),
                             bbturbulence = BBTurbulence$signal$base[ ,"md_delta"],
                             # bbturbulence_sigmoid = md_delta_sigmoid,
                             # bbturbulence_ema_sigmoid = 1 / (1 + exp(-md_delta_ema * 100)),
                             bcpolz = BCPOLZ$getSignal(),
                             eigen = Eigen$getSignal(),
                             momentum = Momentum$getSignal(),
                             ofrfsi = OFRFSI$getSignal() ) )
  
  
  # signals <- OFRFSI$getSignal()
  # colnames(signals) <- "ofrfsi"
  
  
  
  # Define training dates
  training_dates <- intersect( rownames(signals), rownames(states) )
  states <- states[training_dates, ]
  BBS <- BBS[training_dates, ]

  
  
  
  # --------------------------------------------------------------------------
  # BACKTESTS
  # --------------------------------------------------------------------------
  
  source( paste0(wd, "Source/custom_functions.R") )
  
  BacktestDAAWIP <- setRefClass( Class = "BacktestDAAWIP", 
                                 contains = "BacktestDAA",
                                 methods = list() )
  
  BacktestDAAWIP$methods( rpm = BacktestDAA.rpm_mvu )
  BacktestDAAWIP$methods( rpm_disc = BacktestDAA.rpm_mvu_disc )
  BacktestDAAWIP$methods( runDAA = BacktestDAA.runDAA )
  
  
  # Instantiate Backtest class
  BT <- BacktestDAAWIP$new()
  BT$setCtrl()
  BT$updateData()
  
  
  # Define rebalancing dates
  X <- BT$data$X[ ,"BM"]
  rebdates <- rownames(X)[ seq(from = 252, to = nrow(X), by = 1) ]
  
  # signal_names <- colnames(signals)
  signal_names <- c("bbturbulence.md_delta", 
                    "bbturbulence_sigmoid.md_delta",
                    "bbturbulence_ema_sigmoid")
  
  for ( signal_name in signal_names ) {
      
    BT <- BacktestDAAWIP$new( )
    BT$spec$rebdates <- intersect( rebdates, rownames(signals) )
    # FUN <- function(X) { mom.ewma( X = X, tau = 21 ) }
    # sig_roll <- applyRoll( Data = signals[ ,signal_name], 
    #                        Width = 252, 
    #                        By = 1, 
    #                        FUN = FUN )
    # FUN = function(X) { as.numeric(X[nrow(X), ]) - as.numeric(X[1, ]) }
    # sig_roll_delta <- applyRoll( Data = sig_roll,
    #                              Width = 21,
    #                              By = 1,
    #                              FUN = FUN )
    # tmp <- abs(sig_roll[rownames(sig_roll_delta), ]) * sig_roll_delta
    # plot( tmp )
    # 
    # BT$data <- list( X = X,
    #                  signal = sig_roll,
    #                  states = states,
    #                  TD = trainingData( Y_train = BBS,
    #                                     X_train = sig_roll ))
    BT$data <- list( X = X,
                     signal = signals[ ,signal_name],
                     states = states,
                     TD = trainingData( Y_train = BBS,
                                        X_train = signals[ ,signal_name] ))
                                       

    BT$spec$mvu_lambda <- 1
    BT$runDAA( rebdates = BT$spec$rebdates, method = "rpm" )
    # name <- paste0("rpm_mvu1_", signal_name, "_ewma21")
    name <- paste0("rpm_mvu1_", signal_name)
    BT$save( wd = wd_warehouse,
             name = name,
             without_data = FALSE )
    
    BT$output <- list()
    BT$spec$mvu_lambda <- 1
    BT$runDAA( rebdates = BT$spec$rebdates, method = "rpm_disc" )
    # name <- paste0("rpm_mvu1_disc_", signal_name, "_ewma21")
    name <- paste0("rpm_mvu1_disc_", signal_name)
    BT$save( wd = wd_warehouse, 
             name = name,
             without_data = FALSE )
    
    
  }
  
  
  plot( BT$output$weights )
  
  
  
  
  DS <- dirichletSampling( Y_train = BT$data$TD$Y_train, # BT$data$X[ ,"BM"],
                           X_train = tmp,
                           sclfct = 1 )
  ldens <- lapply( DS, densFUN )
  plot.ldensity( ldens )
  
  
  
  
  
  #####
  # EWMA
  sig <- BT$data$signal
  sig_ema <- ema( sig, alpha = 0.01 )
  FUN <- function(X) { mom.ewma(X = X, tau = 21 ) }
  sig_ewma <- applyRoll( Data = sig, Width = 252, FUN = FUN, By = 1 )
  
  sig_all <- cbind(ema = sig_ema, ewma = sig_ewma )
  plot(sig_all)
  plot(sig_all, plot.type = "single" )
  
  
  
  
  
  BT1 <- loadBacktest( wd = wd_warehouse,
                       name = "rpm_mvu1_disc_ofrfsi" )
  BT2 <- loadBacktest( wd = wd_warehouse,
                       name = "rpm_mvu2_disc_ofrfsi" )
  
  wghts <- cbind(BT1$output$weights, BT2$outpout$weights)
  plot( wghts )
  
  
  
  
  
  
  
  
  plot( BT$output$weights )
  
  
  
  BT0 <- BacktestBase$new()
  BT1 <- loadBacktest( wd = wd_warehouse, 
                       name = "rpm_mvu_disc_ofrfsi", 
                       b_new = TRUE )
  BT <- BT0
  BT$data <- BT1$data
  BT$spec <- BT1$spec
  BT$output <- BT1$output
  BT$spec$verbose <- TRUE
  
  
  
  tail( BT$output$weights, 100 )
  plot( cbind(BT$output$weights, BT$output$weights2), plot.type = "single" )
  
  plot( BT$output$lambda )
  plot( apply(BT$output$exp_state, 1, mean), plot.type = "single" )
  colors <- fBasics::divPalette(n = ncol(BT$output$exp_util), "RdYlGn")
  plot( BT$output$exp_util, plot.type = "single", col = colors )
  plot( BT$output$exp_regret, plot.type = "single", col = colors )
  
  

  
  
  
  
  sig <- BT$output$weights
  sig <- (1 + sign(apply(BT$output$exp_state, 1, mean))) / 2
  sig <- (1 + sign(BT$output$exp_util[ ,2])) / 2
  sig <- (1 - sign(BT$data$signal)) / 2

  y <- BT$data$X[ ,"BM"]
  test <- signalTesting.byTrading(X = y,
                                  sig = sig,
                                  n_lag = 1,
                                  tc = 0)
  X_tmp <- na.omit(cbind(bm = y, test))
  
  plot( log(cumulated(X_tmp, "discrete")), plot.type = "single" )

  abline(v = time(sig)[which(sig == 1)], col = 3)
  abline(v = time(sig)[which(sig == 0)], col = 2)
  lines(log(cumulated(X_tmp[ ,1], "discrete")), col = 1)
  lines(log(cumulated(X_tmp[ ,2], "discrete")), col = 4)
  
  
  t(descStats(X_tmp)$stats[c("cumret", "sds", "sharpe", "maxDD"), ])
  
  
  drawDownPlot( X_tmp )
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Augment univariate signal to multivariate set by considering also past
  # signal values. Boil back down to univariate signal using Mahalanobis dist.
  # --------------------------------------------------------------------------
  
  require(AAA)
  
  signal <- signals[ ,signal_name]
  
  n_days <- 21
  tmp <- lapply( 0:n_days, FUN = function(x) { lag(signal, x) } )
  X <- na.omit( do.call( cbind, tmp ) )
  today <- "2020-06-17"
  today <- "2020-03-10"
  today <- "2008-10-10"
  X_eval <- X[today, ]
  # X_eval <- apply(X, 2, mean)
  covmat <- cov( X )
  md <- mahalanobis.ts( X = X, center = X_eval, covmat = covmat )
  eps <- density(md)$bw
  p <- exp( -1 / (2 * eps) * md )
  d <- abs(md - as.numeric(md[today, ]))
  eps_d <- density(d)$bw
  pd <- exp( -1 / (2 * eps_d) * d )
  
  plot( md )
  plot( p )
  plot( d )
  plot( pd ) 
  barplot(X_eval)
    
  
  # debugonce(importanceSampling)
  IS <- importanceSampling(Y_train = X_bm,
                           X_train = md)
  plot.is( IS )

  
  DS <- dirichletSampling(Y_train = X_bm,
                          X_train = md,
                          sclfct = nrow(md)) 
  ldens <- lapply( DS, density )
  plot.ldensity( ldens )
  
  barplot( unlist( lapply( DS, mean ) ) )
  
  
  DS <- dirichletSampling(Y_train = X_bm,
                          X_train = md,
                          sclfct = 1 ) 
  ldens <- lapply( DS, density )
  plot.ldensity( ldens )
  
  
  
  
  #################
  
  n_days <- 5
  tmp <- lapply( 0:n_days, FUN = function(x) { lag(signal, x) } )
  X <- na.omit( do.call( cbind, tmp ) )
  dates <- c("2020-06-17", "2020-03-10", "2008-10-10")
  # dates <- tail(rownames(X), 80)
  X_eval <- X[dates, ]
  
  DS <- dirichletSampling( Y_train = X_bm,
                           X_train = X,
                           X_eval = X_eval, 
                           weights_fun = "kernel",
                           sclfct = 1 )
  ldens <- lapply( DS, densFUN )
  plot.ldensity( ldens )
  
  barplot( unlist( lapply( DS, mean) ) ) 
    
  
  # VS.
  X_eval <- apply(X, 2, mean)
  covmat <- cov( X )
  md <- mahalanobis.ts( X = X, center = X_eval, covmat = covmat )
  
  DS <- dirichletSampling( Y_train = X_bm,
                           X_train = md,
                           X_eval = as.numeric(md[dates, ]),
                           weights_fun = "l1",
                           sclfct = 1 )
  ldens <- lapply( DS, densFUN )
  plot.ldensity( ldens )
  
  barplot( unlist( lapply( DS, mean) ) ) 
  
  
  
  

    