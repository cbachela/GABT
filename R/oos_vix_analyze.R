  
  
  ############################################################################
  ### GOOD AND BAD TURBULENCE - OUT-OF-SAMPLE ANALYSIS - VIX
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     01.04.2022
  # First version:    01.04.2022
  # --------------------------------------------------------------------------
  
  
  # --------------------------------------------------------------------------
  # REQUIRE
  # --------------------------------------------------------------------------
  
  require(RP)
  require(BSS)
  require(simolz)
  require(DAARC)
  require(visolz)
  
  wd <- "R:/Asset_Management/Research_Projects/Multi_Asset/DAA/DAARC/Sensor_Testing/Univariate/Good_And_Bad_Turbulence/"
  source( paste0(wd, "Source/source_functions.R") )
  
  
  
  # --------------------------------------------------------------------------
  # Load VIX
  # --------------------------------------------------------------------------
  
  # VRP <- loadSensor( sensor_name = "vrp" )
  # VRP$update()
  # vix <- VRP$signal$base[ ,"iv"]
  
  X_sptr <-  rodbcGetOLZDBReturns( assetName = "SPTR",
                                   refCcy = "USD", 
                                   frqncy = "daily", 
                                   compounding = "discrete" )
  X_sptr <- X_sptr[isWeekday(time(X_sptr)), ]
  
  vix_ret <- rodbcGetOLZDBReturns( assetName = "VIX",
                                   refCcy = "USD", 
                                   frqncy = "daily", 
                                   compounding = "continuous" )
  
  vix_ret <- vix_ret[isWeekday(time(vix_ret)), ]
  vix_base <- 21.8999999999999  # this is the level of the vix as of 01.03.1990
  vix <- cumulated( vix_ret, compounding = "continuous" ) * vix_base
  
  
 
  
  
  # --------------------------------------------------------------------------
  # BBT base
  # --------------------------------------------------------------------------
  
  bbt <- loadSensor( sensor_name = "bbturbulence_base_vix" )
  
  
  sig_is <- bbt$signal$insample[ ,"delta"]
  sig_oos <- bbt$getSignal()
  sig <- na.omit( cbind( is = sig_is, 
                         oos = bbt$signal$base[ ,"md_delta"],
                         oos_ema = sig_oos ) )
  ###
  # sig <- scale(sig, TRUE, FALSE)
  ###
  # sig <- sig * (-1)
  sig_binary <- (1 + sign(ema(sig, 1))) / 2
  apply( sig_binary, 2, mean )
  
  
  # plot( sig )
  # plot( bbt$signal$base[ ,c("md_bear_scl", "md_bull_scl", "md_all_scl", "md_delta")] )
  
  X_bm <- bbt$data$X_bm
  test <- signalTesting.byTrading( X = X_bm,
                                   sig = sig_binary,
                                   n_lag = 1, 
                                   tc = 0.004 )
  test_nc <- signalTesting.byTrading( X = X_bm,
                                      sig = sig_binary,
                                      n_lag = 1, 
                                      tc = 0 )
  X_tmp <- na.omit( cbind( X_bm, test ) )
  X_tmp_nc <- na.omit( cbind( X_bm, test_nc ) )
  # plot( cumulated(X_tmp, "discrete"), plot.type = "single" )
  
  stats_fields <- c("cumret", "means", "sds", "sharpe", "maxDD")
  
  plot( as.simTS(X_tmp), logscale = TRUE )
  t( descStats( X_tmp )$stats[stats_fields, ] )
  
  plot( as.simTS(X_tmp_nc), logscale = TRUE )
  t( descStats( X_tmp_nc )$stats[stats_fields, ] )
  

  
  BBSObj <- bbt$data$BBS$copy()
  BBSObj$data$X_level <- cumulated( bbt$data$X_bm, "discrete" )
  BBSObj$data$X <- bbt$data$X_bm
  BBSObj$runRobust()
  BBSObj$output$states <- sign( sig[ ,"oos"] )
  # debugonce( BBSObj$plotStates )
  BBSObj$plotStates( logarithmic = FALSE )
  # debugonce( BBSObj$phaseStats )
  BBSObj$phaseStats( n_lag = 0 )
  BBSObj$phaseStats( n_lag = 1 )
  do.call( cbind, lapply( 0:3, FUN = function(i) { BBSObj$phaseStats( n_lag = i )[[1]] } ) )
  do.call( cbind, lapply( 0:3, FUN = function(i) { BBSObj$phaseStats( n_lag = i )[[2]] } ) )
  
  
  
  # --------------------------------------------------------------------------
  # BBT scl2
  # --------------------------------------------------------------------------
  
  bbt_scl2 <- loadSensor( sensor_name = "bbturbulence_scl2_vix" )
  
  sig_is <- bbt_scl2$signal$insample[ ,"delta"]
  sig_oos <- bbt_scl2$getSignal()
  sig <- na.omit( cbind( is = sig_is, 
                         oos = bbt_scl2$signal$scl2[ ,"md_delta"],
                         oos_ema = sig_oos ) )
  ###
  # sig <- scale(sig, TRUE, FALSE)
  ###
  # sig <- sig * (-1)
  sig_binary <- (1 + sign(ema(sig, 1))) / 2
  apply( sig_binary, 2, mean )
  
  
  X_bm <- bbt_scl2$data$X_bm
  X_bm <- X_sptr
  test <- signalTesting.byTrading( X = X_bm,
                                   sig = sig_binary,
                                   n_lag = 1, 
                                   tc = 0.004 )
  test_nc <- signalTesting.byTrading( X =X_bm,
                                      sig = sig_binary,
                                      n_lag = 1, 
                                      tc = 0 )
  X_tmp <- na.omit( cbind( X_bm, test ) )
  X_tmp_nc <- na.omit( cbind( X_bm, test_nc ) )
  # plot( cumulated(X_tmp, "discrete"), plot.type = "single" )
  
  stats_fields <- c("cumret", "means", "sds", "sharpe", "maxDD")
  
  plot( as.simTS(X_tmp), logscale = TRUE )
  t( descStats( X_tmp )$stats[stats_fields, ] )
  
  plot( as.simTS(X_tmp_nc), logscale = TRUE )
  t( descStats( X_tmp_nc )$stats[stats_fields, ] )
  
  
  
  
  BBSObj <- bbt_scl2$data$BBS$copy()
  BBSObj$data$X_level <- cumulated( bbt_scl2$data$X_bm, "discrete" )
  BBSObj$data$X <- bbt_scl2$data$X_bm
  BBSObj$runRobust()
  BBSObj$output$states <- sign( sig[ ,"oos"] )
  # debugonce( BBSObj$plotStates )
  BBSObj$plotStates( logarithmic = FALSE )
  # debugonce( BBSObj$phaseStats )
  BBSObj$phaseStats( n_lag = 0 )
  BBSObj$phaseStats( n_lag = 1 )
  do.call( cbind, lapply( 0:3, FUN = function(i) { BBSObj$phaseStats( n_lag = i )[[1]] } ) )
  do.call( cbind, lapply( 0:3, FUN = function(i) { BBSObj$phaseStats( n_lag = i )[[2]] } ) )
  
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # BBTEWMA base
  # --------------------------------------------------------------------------
  
  bbtewma <- loadSensor( sensor_name = "bbtewma_base_vix" )

  
  sig_is <- bbtewma$signal$insample[ ,"delta"]
  sig_oos <- bbtewma$getSignal()
  sig <- na.omit( cbind( is = sig_is, 
                         bbs = bbtewma$signal$base[ ,"states"],
                         oos = sig_oos ) )
  ###
  sig <- scale(sig, TRUE, FALSE)
  ###
  # sig <- sig * (-1)
  sig_binary <- (1 + sign(ema(sig, 1))) / 2
  apply( sig_binary, 2, mean )

  
  X_bm <- bbtewma$data$X_bm
  # X_bm <- X_sptr
  test <- signalTesting.byTrading( X = X_bm,
                                   sig = sig_binary,
                                   n_lag = 1, 
                                   tc = 0.004 )
  test_nc <- signalTesting.byTrading( X =X_bm,
                                      sig = sig_binary,
                                      n_lag = 1, 
                                      tc = 0 )
  X_tmp <- na.omit( cbind( X_bm, test ) )
  X_tmp_nc <- na.omit( cbind( X_bm, test_nc ) )
  # plot( cumulated(X_tmp, "discrete"), plot.type = "single" )
  
  stats_fields <- c("cumret", "means", "sds", "sharpe", "maxDD")
  
  plot( as.simTS(X_tmp), logscale = TRUE )
  t( descStats( X_tmp )$stats[stats_fields, ] )
  
  plot( as.simTS(X_tmp_nc), logscale = TRUE )
  t( descStats( X_tmp_nc )$stats[stats_fields, ] )
  
  
  
  BBSObj <- bbtewma$data$BBS$copy()
  BBSObj$data$X_level <- cumulated( bbtewma$data$X_bm, "discrete" )
  BBSObj$data$X <- bbtewma$data$X_bm
  BBSObj$runRobust()
  BBSObj$output$states <- sign( sig[ ,"oos"] )
  # BBSObj$output$states <- sign( bbtewma$signal$base[ ,"states"] )
  # BBSObj$output$states <- sign( bbtewma$signal$insample[ ,"states"])
  # debugonce( BBSObj$plotStates )
  # BBSObj$plotStates( logarithmic = FALSE )
  # debugonce( BBSObj$phaseStats )
  BBSObj$phaseStats( n_lag = 0 )
  BBSObj$phaseStats( n_lag = 1 )
  do.call( cbind, lapply( 0:3, FUN = function(i) { BBSObj$phaseStats( n_lag = i )[[1]] } ) )
  do.call( cbind, lapply( 0:3, FUN = function(i) { BBSObj$phaseStats( n_lag = i )[[2]] } ) )
  
  
  
  # --------------------------------------------------------------------------
  # BBTEWMA scl2
  # --------------------------------------------------------------------------

  bbtewma_scl2 <- loadSensor( sensor_name = "bbtewma_scl2_vix" )
  # bbtewma_scl2 <- loadSensor( sensor_name = "bbtewma_scl2_vix_bbc_phase=10_cycle=30" )
  # bbtewma_scl2 <- loadSensor( sensor_name = "bbtewma_scl2_vix_bbc_phase=15_cycle=45" )
  
  # Signals
  sig_is <- bbtewma_scl2$signal$insample[ ,"delta"]
  sig_oos <- bbtewma_scl2$getSignal()
  sig <- na.omit( cbind( is = sig_is, 
                         bbs = bbtewma_scl2$signal$scl2[ ,"states"],
                         oos = sig_oos ) )
  ###
  sig <- scale( sig, TRUE, FALSE )
  ###
  # sig <- sig * (-1)  # ~~~~~~~~~~~~~~~~~~~~~~
  sig_binary <- (1 + sign(ema(sig, 1) - 0)) / 2
  apply( sig_binary, 2, mean )
  
  # Backtest
  X_bm <- bbtewma_scl2$data$X_bm
  # X_bm <- X_sptr
  test <- signalTesting.byTrading( X = X_bm,
                                   sig = sig_binary,
                                   n_lag = 1, 
                                   tc = 0.004 )
  test_nc <- signalTesting.byTrading( X =X_bm,
                                      sig = sig_binary,
                                      n_lag = 1, 
                                      tc = 0 )
  X_tmp <- na.omit( cbind( X_bm, test ) )
  X_tmp_nc <- na.omit( cbind( X_bm, test_nc ) )
  # plot( cumulated(X_tmp, "discrete"), plot.type = "single" )
  
  stats_fields <- c("cumret", "means", "meansAr", "sds", "sharpe", "maxDD")
  
  plot( as.simTS(X_tmp), logscale = TRUE )
  plot( as.simTS(head(tail(X_tmp, 700), 200)), logscale = TRUE )
  t( descStats( X_tmp )$stats[stats_fields, ] )
  
  plot( as.simTS(X_tmp_nc), logscale = TRUE )
  t( descStats( X_tmp_nc )$stats[stats_fields, ] )
  
  
  tmp <- window( na.omit(cbind(sig_binary, cumulated(X_bm, "discrete"))), "2020-02-10", "2020-03-31" )
  tmp <- window( na.omit(cbind(sig_binary, cumulated(X_bm, "discrete"))), "2008-02-10", "2008-11-30" )
  plot(tmp, type = "o")
  
  
  ldens <- apply( X_tmp, 2, density )
  slolz:::plot.ldensity( ldens, fillin = FALSE, xlim = c(-0.4, 0.4), ylim = c(0, 30) )
  abline( v = 0, col = "grey" )
  
  
  
  
  BBSObj <- bbtewma_scl2$data$BBS$copy()
  BBSObj$data$X_level <- cumulated( bbtewma_scl2$data$X_bm, "discrete" )
  BBSObj$data$X <- bbtewma_scl2$data$X_bm
  BBSObj$runRobust()
  # BBSObj$output$states <- sign( sig[ ,"oos"] )
  # BBSObj$output$states <- sign( sig[ ,"bbs"] )
  BBSObj$output$states <- sign( sig[ ,"is"] )
  # BBSObj$plotStates( logarithmic = FALSE )
  BBSObj$phaseStats( n_lag = 0 )
  BBSObj$phaseStats( n_lag = 1 )
  do.call( cbind, lapply( 0:3, FUN = function(i) { BBSObj$phaseStats( n_lag = i )[[1]] } ) )
  do.call( cbind, lapply( 0:3, FUN = function(i) { BBSObj$phaseStats( n_lag = i )[[2]] } ) )
  
  
  
  
  n_lag <- 1
  annualize <- TRUE
  states <- BBSObj$output$states
  if ( is.null(states)  ) {
    stop("No states in output list")
  }
  X <- BBSObj$data$X
  if ( is.null(states)  ) {
    stop("No dataset 'X' in data list")
  }
  dates <- intersect( rownames(states), rownames(X) )
  X <- X[dates, ]
  states <- states[dates, ]
  
  dates_up <- dates[ intersect(1:length(dates), which(states == 1) + n_lag) ]
  dates_down <- dates[ intersect(1:length(dates), which(states == -1) + n_lag) ]
  X1 <- X[dates_up, ]
  X2 <- X[dates_down, ]
  
  if ( isTRUE( annualize) ) {
    frqncy <- xts::periodicity(X)
    scalefactor <- switch(frqncy$scale,
                          daily   = 252,
                          weekly  = 52,
                          monthly = 12, 
                          yearly  = 1)
  } else {
    scalefactor <- 1
  }
  
  lBoot1 <- apply( X1, 2, FUN = boot, statistic = meanGeoBoot, R = 10^3 )
  ldens1 <- lapply( lBoot1, FUN = function(x) { density( x$t * scalefactor ) } ) 
  lBoot2 <- apply( X2, 2, FUN = boot, statistic = meanGeoBoot, R = 10^3 )
  ldens2 <- lapply( lBoot2, FUN = function(x) { density( x$t * scalefactor ) } )
  
  ans <- list( "states==1" = ldens1,
               "states==-1" = ldens2 )
  
  slolz:::plot.ldensity( lapply( ans, FUN = function(x) { x[[1]] } ) )
  
  
  
  
  
  
  
  
  # Ticketfee
  # debugonce( signalTesting.byTrading )
  test_tf <- signalTesting.byTrading( X = X_bm,
                                      sig = sig_binary,
                                      n_lag = 1, 
                                      tc = 0,
                                      ticketfee = 100,
                                      nav = 10^7,
                                      fc = 0 )
  # X_tmp <- na.omit( cbind( X_bm, test, test_tf ) )
  X_tmp <- na.omit( cbind( X_bm, test_tf ) )
  plot( as.simTS(X_tmp) )
  
  descStats(X_tmp)
  
  
  
  TD <- trainingData( Y_train = bbtewma_scl2$data$X_bm,
                      X_train = bbtewma_scl2$getSignal(), 
                      n_lag = 1 )
  
  plot( x = TD$X_train, y = TD$Y_train )
  cor( TD$DF )
  
  
  
  X_train <- lag( bbtewma_scl2$data$X_bm, -1)
  colnames(X_train) <- colnames(bbtewma_scl2$data$X_bm)
  Y_train <- (1 + sign(ema(bbtewma_scl2$getSignal() * (-1), 1) - 0)) / 2
  reg <- regression( Y_train = Y_train,
                     X_train = X_train,
                     type = "logit" )
  summary(reg$reg)
  plot( x = reg$z, y = reg$y )
  
  
  
  
  
  
  X <- bbtewma_scl2$data$X_bm
  descStats(X, descStatsSpec(compounding = "discrete"))  
  
  TD <- trainingData( Y_train = bbtewma_scl2$data$X_bm,
                      X_train = sig[ ,"oos"], 
                      n_lag = 1 )
  
  # y_hat <- (1 + sign(ema(TD$X_train - 0, 1))) / 2
  # y_true <- (1 + sign(TD$Y_train)) / 2
  y_hat <- TD$X_train
  y_true <- TD$Y_train
  roc( y_hat = y_hat, y_true = y_true, th = 0.5 )
  roc( y_hat = sig, y_true = bbtewma_scl2$data$X_bm, th = 0.1 )
  
  rocFUN <- function(x) { roc( y_hat = y_hat, y_true = y_true, th = x ) }
  lROC <- lapply( seq(from = -1, to = 1, by = 0.01), FUN = rocFUN )
  ROC <- do.call( rbind, lROC )
  plot( ROC[ ,1], ROC[ ,2] )
  
  
  
  require(pROC)
  pROC_obj <- roc( TD$Y_train, 
                   TD$X_train,
                   smoothed = TRUE,
                   # arguments for ci
                   ci=TRUE, ci.alpha=0.9, stratified=FALSE,
                   # arguments for plot
                   plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                   print.auc=TRUE, show.thres=TRUE )
  
  
  sens.ci <- ci.se(pROC_obj)
  plot(sens.ci, type="shape", col="lightblue")
  ## Warning in plot.ci.se(sens.ci, type = "shape", col = "lightblue"): Low
  ## definition shape.
  plot(sens.ci, type="bars")
  
  
  
  
  
  lPAS <- pas( X = bbtewma_scl2$data$X_bm,
               sig = bbtewma_scl2$getSignal() * (-1),
               n_lags = 0:5,
               compounding = "discrete" )
  
  plot( x = lPAS[[1]]$X_train, y = lPAS[[1]]$Y_train )
  points( x = as.numeric(lPAS[[6]]$X_train), y = as.numeric(lPAS[[6]]$Y_train), col = 2 )
  
  unlist( lapply( 1:length(lPAS), FUN = function(i) { cor( lPAS[[i]]$DF )[1, 2] } ) )
  
  
  
  
  # No-trade-zone
  
  signal <- sig_oos * 0
  # signal[sig > 0.1, ] <- 1
  # signal[sig < 0.1, ] <- -1
  signal[sig > 0, ] <- 1
  mean(signal)
  
  test_tf <- signalTesting.byTrading( X = bbtewma_scl2$data$X_bm,
                                      sig = signal,
                                      n_lag = 1, 
                                      tc = 0,
                                      ticketfee = 100,
                                      nav = 10^7,
                                      fc = 0.1 )
  X_tmp <- na.omit( cbind( bbtewma_scl2$data$X_bm, test_tf ) )
  plot( as.simTS(X_tmp) )
  
  
  
  
  # --------------------------------------------------------------------------
  # BBTGARCH base
  # --------------------------------------------------------------------------
  
  bbtgarch <- loadSensor( sensor_name = "bbtgarch_base_vix" )
  
  
  sig_is <- bbtgarch$signal$insample[ ,"delta"]
  sig_oos <- bbtgarch$getSignal()
  sig <- na.omit( cbind( is = sig_is, 
                         bbs = bbtgarch$signal$base[ ,"states"],
                         oos = sig_oos ) )
  ###
  # sig <- scale(sig, TRUE, FALSE)
  ###
  # sig <- sig * (-1)
  sig_binary <- (1 + sign(ema(sig, 1))) / 2
  apply( sig_binary, 2, mean )
  
  
  X_bm <- bbtgarch$data$X_bm
  test <- signalTesting.byTrading( X = X_bm,
                                   sig = sig_binary,
                                   n_lag = 1, 
                                   tc = 0.004 )
  test_nc <- signalTesting.byTrading( X =X_bm,
                                      sig = sig_binary,
                                      n_lag = 1, 
                                      tc = 0 )
  X_tmp <- na.omit( cbind( X_bm, test ) )
  X_tmp_nc <- na.omit( cbind( X_bm, test_nc ) )
  # plot( cumulated(X_tmp, "discrete"), plot.type = "single" )
  
  stats_fields <- c("cumret", "means", "sds", "sharpe", "maxDD")
  
  plot( as.simTS(X_tmp), logscale = TRUE )
  t( descStats( X_tmp )$stats[stats_fields, ] )
  
  plot( as.simTS(X_tmp_nc), logscale = TRUE )
  t( descStats( X_tmp_nc )$stats[stats_fields, ] )
  
  
  
  
  BBSObj <- bbtgarch$data$BBS$copy()
  BBSObj$data$X_level <- cumulated( bbtgarch$data$X_bm, "discrete" )
  BBSObj$data$X <- bbtgarch$data$X_bm
  BBSObj$runRobust()
  BBSObj$output$states <- sign( sig[ ,"oos"] )
  # BBSObj$output$states <- sign( sig[ ,"bbs"] )
  # BBSObj$output$states <- sign( sig[ ,"is"] )
  # BBSObj$plotStates( logarithmic = FALSE )
  BBSObj$phaseStats( n_lag = 0 )
  BBSObj$phaseStats( n_lag = 1 )
  do.call( cbind, lapply( 0:3, FUN = function(i) { BBSObj$phaseStats( n_lag = i )[[1]] } ) )
  do.call( cbind, lapply( 0:3, FUN = function(i) { BBSObj$phaseStats( n_lag = i )[[2]] } ) )
  
  
  
  n_lag <- 2
  annualize <- TRUE
  states <- BBSObj$output$states
  if ( is.null(states)  ) {
    stop("No states in output list")
  }
  X <- BBSObj$data$X
  if ( is.null(states)  ) {
    stop("No dataset 'X' in data list")
  }
  dates <- intersect( rownames(states), rownames(X) )
  X <- X[dates, ]
  states <- states[dates, ]
  
  dates_up <- dates[ intersect(1:length(dates), which(states == 1) + n_lag) ]
  dates_down <- dates[ intersect(1:length(dates), which(states == -1) + n_lag) ]
  X1 <- X[dates_up, ]
  X2 <- X[dates_down, ]
  
  if ( isTRUE( annualize) ) {
    frqncy <- xts::periodicity(X)
    scalefactor <- switch(frqncy$scale,
                          daily   = 252,
                          weekly  = 52,
                          monthly = 12, 
                          yearly  = 1)
  } else {
    scalefactor <- 1
  }
  
  lBoot1 <- apply( X1, 2, FUN = boot, statistic = meanGeoBoot, R = 10^3 )
  ldens1 <- lapply( lBoot1, FUN = function(x) { density( x$t * scalefactor ) } ) 
  lBoot2 <- apply( X2, 2, FUN = boot, statistic = meanGeoBoot, R = 10^3 )
  ldens2 <- lapply( lBoot2, FUN = function(x) { density( x$t * scalefactor ) } )
  
  ans <- list( "states==1" = ldens1,
               "states==-1" = ldens2 )
  
  slolz:::plot.ldensity( lapply( ans, FUN = function(x) { x[[1]] } ) )
  
  
  
  
  
  
  
  
  
  
  
  
  
  ###############################
  
  bbtewma <- loadSensor( sensor_name = "bbtewma_base_vix" )
  bbtewma_scl2 <- loadSensor( sensor_name = "bbtewma_scl2_vix" )
  tmp <- na.omit( cbind( base = bbtewma$signal$base[ ,"states"],
                         scl2 = bbtewma_scl2$signal$scl2[ ,"states"] ) )
  plot(tmp[ ,1] - tmp[ ,2] )
  
  
  bbtewma$computeSignalInsample()
  bbtewma_scl2$computeSignalInsample()
  
  s1 <- bbtewma$signal$insample
  s2 <- bbtewma_scl2$signal$insample
  tmp <- cbind( s1[ ,"delta"], s2[ ,"delta"] )
  
  head(s1)
  head(s2)
  plot( tmp )
  plot( tmp, plot.type = "single" )
  
  
  
  X_bm <- bbtewma$data$X_bm
  # sig <- tmp * (-1)
  sig <- tmp
  sig <- scale( sig, TRUE, FALSE )
  sig_binary <- (1 + sign(ema(sig, 1))) / 2
  apply( sig_binary, 2, mean )
  test <- signalTesting.byTrading( X = X_bm,
                                   sig = sig_binary,
                                   # sig = (1 + bbtewma$signal$insample[ ,"states"]) / 2,
                                   n_lag = 1, 
                                   tc = 0.004 )
  X_tmp <- na.omit( cbind( X_bm, test ) )
  plot( as.simTS(X_tmp), logscale = TRUE )
  stats_fields <- c("cumret", "means", "sds", "sharpe", "maxDD")
  t( descStats( X_tmp )$stats[stats_fields, ] )
  
  
  
  BBSObj <- bbtewma$data$BBS$copy()
  BBSObj$data$X_level <- cumulated( bbtewma$data$X_bm, "discrete" )
  BBSObj$data$X <- bbtewma$data$X_bm
  BBSObj$runRobust()
  BBSObj$output$states <- sig_binary[ ,2] * 2 - 1
  BBSObj$phaseStats( n_lag = 0 )
  BBSObj$phaseStats( n_lag = 1 )
  
  
  
  tmp <- cbind( is = bbtewma_scl2$signal$insample[ ,"states"],
                sig = sig_binary[ ,2] * 2 - 1 )
  plot( tmp )
  
  plot( tail(tmp, 300) )
  
 
  
  
  
  # --------------------------------------------------------------------------
  # Compute signal on squared returns
  # --------------------------------------------------------------------------
  
  
  bbtewma <- loadSensor( sensor_name = "bbtewma_base_dm" )

  X_sptr_sq <- ema( X_sptr^2, 0.01 )[-c(1:10), ]
  bbtewma$data$X_bm <- returns( X_sptr_sq, "discrete" )
  bbtewma$data$capw <- NULL
  bbtewma$spec$name <- "bbtewma_base_sptrsq"
  bbtewma$signal <- list()
  # debugonce( bbtewma$computeSignalInsample )
  bbtewma$computeSignalInsample()
  # debugonce( bbtewma$computeSignal )
  bbtewma$updateSignal()
  # bbtewma$save()

  
  BBSObj <- bbtewma$data$BBS$copy()
  BBSObj$data$X_level <- cumulated( bbtewma$data$X_bm, "discrete" )
  BBSObj$data$X <- bbtewma$data$X_bm
  BBSObj$runRobust()
  # debugonce( BBSObj$plotStates )
  BBSObj$plotStates( logarithmic = FALSE )
  BBSObj$phaseStats()
  
  
  
  
  sig_is <- bbtewma$signal$insample[ ,"delta"]
  sig_oos <- bbtewma$getSignal()
  sig <- na.omit( cbind( is = sig_is, 
                         bbs = bbtewma$signal$base[ ,"states"],
                         oos = sig_oos ) )
  sig <- sig * (-1)
  sig_binary <- (1 + sign(ema(sig, 1))) / 2
  apply( sig_binary, 2, mean )
  
  
  X_bm <- bbtewma$data$X_bm
  # X_bm <- X_sptr
  test <- signalTesting.byTrading( X = X_bm,
                                   sig = sig_binary,
                                   n_lag = 1, 
                                   tc = 0.004 )
  test_nc <- signalTesting.byTrading( X = X_bm,
                                      sig = sig_binary,
                                      n_lag = 1, 
                                      tc = 0 )
  X_tmp <- na.omit( cbind( X_bm, test ) )
  # X_tmp <- na.omit( cbind( X_bm, test_nc ) )
  plot( as.simTS(X_tmp) )
  plot( cumulated(X_tmp, "discrete"), plot.type = "single" )
  
  stats_fields <- c("cumret", "means", "sds", "sharpe", "maxDD")
  t( descStats( X_tmp )$stats[stats_fields, ] )
  
  
  
  
  
  
  
  
 
    
  
  
  
  
 
  
  
  
  # --------------------------------------------------------------------------
  # Compute signal on ETF
  # --------------------------------------------------------------------------
  
  tickers <- c("LVOEUR EU", "VIXY US", "VIXM US")
  
  X <- rodbcGetOLZDBReturns( assetName = tickers,
                             refCcy = "Local",
                             frqncy = "daily" )
  X <- X[isWeekday(time(X)), tickers]
  
  
  plot( cumulated(na.omit(X[ ,1]), "discrete"))
  plot( as.simTS(X) )
  
  
  descStats( na.omit(X) )  
  
  
  
  # Compare VIX to VIX-ETF's
  start_date <- start(X)
  vix[start_date, ]
  
  fc <- 0.6
  fixcost <- (1 + fc)^(1 / 252) - 1
  plot( as.simTS( na.omit(cbind(returns(vix, "discrete") - fixcost, X)) ) )
  
  
  
  
  # Run backtest
  
  Dbbtewma_vixy <- loadSensor( sensor_name = "bbtewma_scl2_dm" )

  bbtewma_vixy$data$X_bm <- na.omit(X[ ,2])
  dates <- intersect( rownames(bbtewma_vixy$data$X_bm), rownames(bbtewma_vixy$data$X) )
  bbtewma_vixy$data$X <- bbtewma_vixy$data$X[dates, ]
  bbtewma_vixy$spec$name <- "bbtewma_scl2_vixy_us"
  bbtewma_vixy$signal <- list()
  # debugonce( bbtewma_vixy$computeSignalInsample )
  bbtewma_vixy$computeSignalInsample()
  # debugonce( bbtewma_vixy$computeSignal )
  bbtewma_vixy$updateSignal()
  # bbtewma_vixy$save()
  
  
  sig_is <- bbtewma_vixy$signal$insample[ ,"delta"]
  sig_oos <- bbtewma_vixy$getSignal()
  # sig <- sig_is * (-1)
  sig <- sig_is
  sig <- na.omit( cbind( is = sig_is, oos = sig_oos ) )
  sig_binary <- (1 + sign(ema(sig, 1) - 0.3)) / 2
  test_tf <- signalTesting.byTrading( X = bbtewma_vixy$data$X_bm,
                                      sig = sig_binary,
                                      n_lag = 2, 
                                      tc = 0,
                                      ticketfee = 100,
                                      nav = 10^7,
                                      fc = 0.01 )
  X_tmp <- na.omit( cbind( bbtewma_vixy$data$X_bm, test_tf ) )
  plot( as.simTS(X_tmp) )
  # plot( as.simTS(X_tmp), logscale = FALSE )
  
  
  descStats(X_tmp)
  
  
  plot( sig )  
  abline(h = 0)  

  
  X_us <- bbtewma_vixy$data$X[ ,"US"]
  Y <- apply( na.omit( cbind(US = X_us, Vola = test_tf[ ,"signal_oos"]) ), 1, mean )
  Y2 <- apply( na.omit( cbind(US = X_us * 2, Vola = test_tf[ ,"signal_oos"]) ), 1, mean )
  
  Y_tmp <- na.omit(cbind(US = X_us, Y = Y, Y2 = Y2, US_half = X_us / 2))
  plot( as.simTS(Y_tmp) )  
  descStats(Y_tmp)
  
  
  
  
  
  
  
    
      