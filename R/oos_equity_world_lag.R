  
  
  ############################################################################
  ### GOOD AND BAD TURBULENCE - OUT-OF-SAMPLE ANALYSIS 
  ### EQUITY WORLD - INCLUDING LAGGED VARIATES
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     15.11.2020
  # First version:    10.10.2020
  # --------------------------------------------------------------------------
  
  
  
  
  
  # --------------------------------------------------------------------------
  # REQUIRE
  # --------------------------------------------------------------------------
  
  require(marima)
  require(RP)
  require(DAARC)
  
  wd <- "R:/Asset_Management/Research_Projects/Multi_Asset/DAA/DAARC/Sensor_Testing/Univariate/Good_And_Bad_Turbulence/"
  source( paste0(wd, "Source/source_functions.R") )
  
  
  
  

  # --------------------------------------------------------------------------
  # ARMA MODEL FOR BEAR AND BULL MD
  # --------------------------------------------------------------------------
  
  BBT0 <- loadSensor( sensor_name = "bbturbulence_scl2_dm",
                      b_new = TRUE )
  
  
  signal <- BBT0$getSignal()
  md <- BBT0$signal$scl2[ ,c("md_bear_scl", "md_bull_scl")]
  plot(md)    
  
  acf(md)
  
  mod <- define.model( kvar = 2, ar = 1, ma = 1 )
  fit <- marima( DATA = md, 
                 ar.pattern = mod$ar.pattern, 
                 ma.pattern = mod$ma.pattern )
  fit
  ls(fit)  
  
  fitted <- timeSeries( t(fit$fitted), time(md) )
  plot(fitted)  
  
  md_delta <- fitted[ ,1] - fitted[ ,2]
  plot(md_delta*10)
  lines(BBT0$signal$scl2[ ,"md_delta"], col = 2)
  
  
  
  mod <- define.model( kvar = 1, ar = 1:21, ma = 1:21 )
  fit <- marima( DATA = BBT0$signal$scl2[ ,"md_delta"], 
                 ar.pattern = mod$ar.pattern, 
                 ma.pattern = mod$ma.pattern )
  fit
  ls(fit)  
  
  fitted <- timeSeries( t(fit$fitted), time(md) )
  plot(fitted*10)  
  lines(BBT0$signal$scl2[ ,"md_delta"], col = 2)
  
  
  # EWMA
  md_bear_ema <- ema( BBT0$signal$scl2[ ,"md_bear_scl"], 0.05 )
  md_bull_ema <- ema( BBT0$signal$scl2[ ,"md_bull_scl"], 0.05 )
  md_delta_ema <- md_bear_ema - md_bull_ema
  
  
      
  
  sig <- (1 + sign(md_delta)) / 2
  sig <- (1 + sign(BBT0$signal$scl2[ ,"md_delta"])) / 2
  sig <- (1 + sign(fitted)) / 2
  sig <- (1 + sign(scale(fitted, TRUE, FALSE) + mean(BBT0$signal$scl2[ ,"md_delta"]))) / 2
  sig <- (1 + sign(md_delta_ema)) / 2
  test <- signalTesting.byTrading( X = BBT0$data$X_bm,
                                   sig = sig, 
                                   n_lag = 1,
                                   tc = 0.004 )
  colnames(test) <- colnames(sig)
  test_nc <- signalTesting.byTrading( X = BBT0$data$X_bm,
                                      sig = sig, 
                                      n_lag = 1,
                                      tc = 0 )
  colnames(test_nc) <- paste0(colnames(sig), "_nc")
  X_tmp <- na.omit( cbind( bm = BBT0$data$X_bm, test, test_nc) )
  plot( as.simTS(X_tmp) )
  descStats(X_tmp)
  
  
  
  
  
  # --------------------------------------------------------------------------
  # LAGGED VARIABLES 
  # --------------------------------------------------------------------------
  
  BBT0 <- loadSensor( sensor_name = "bbturbulence_scl2_dm",
                      b_new = TRUE )
  
  
  signal <- BBT0$getSignal()
  signal <- BBT0$signal$scl2[ ,"md_delta"]
  plot( signal )  

  sig <- (1 + sign(signal)) / 2
  test <- signalTesting.byTrading( X = BBT0$data$X_bm,
                                   sig = sig, 
                                   n_lag = 1,
                                   tc = 0.004 )
  colnames(test) <- colnames(sig)
  test_nc <- signalTesting.byTrading( X = BBT0$data$X_bm,
                                      sig = sig, 
                                      n_lag = 1,
                                      tc = 0 )
  colnames(test_nc) <- paste0(colnames(sig), "_nc")
  X_tmp <- na.omit( cbind( bm = BBT0$data$X_bm, test, test_nc) )
  plot( as.simTS(X_tmp) )
  descStats(X_tmp)
  
  
  
  # Adding lagged covariates narrows the gap between backtests
  # with and without transaction costs because trades become sparse.
  
  # Accuracy
  states_is <- BBT0$signal$insample[ ,"states"]
  
  
  # Our goal is not to maximize accuracy, but to identify
  # regime changes 
  
  
  
  
  
  
  
  
  
  
  
  
  today <- "2012-03-05"
  BBT0$signal$scl2[today, ]
  
  barplot(BBT0$data$X[today, ])
  
  start_date <- "2011-01-01"
  end_date <- "2012-12-31"
  plot( window(BBT0$signal$scl2[ ,"md_delta"], start_date, end_date) )  
  lines( window(BBT0$signal$scl2[ ,"signal"], start_date, end_date),col = 2 )  
  abline( h = 0, col = "grey" )
  abline( v = as.timeDate(today), col = "grey" )
  
  
  
  signal <- BBT0$signal$scl2[ ,"md_delta"]
  plot( signal )
  cvol <- getCondVar( garch(signal) )^0.5
  plot(cvol)    

  
  
    
  
  BBT <- loadSensor( sensor_name = "bbturbulence_scl2_dm_lag0to21",
                      b_new = TRUE )
  plot( BBT$signal$scl2[ ,-c(1:4)] )  
  BBT$computeSignalInsample()
  signal <- BBT$signal$insample[ ,"delta"]
  signal <- BBT$signal$scl2[ ,"md_delta"]
    
  
  
  
  
  ##################
  
  BBT0 <- loadSensor( sensor_name = "bbturbulence_scl2_dm",
                      b_new = TRUE )
  BBT0$computeSignalInsample()
  
  states <- (1 + BBT0$signal$insample[ ,"states"]) / 2  
  plot(density(states)) 
  
  
  sclfct <- NULL
  DS_bear <- dirichletSampling( Y_train = states, 
                                X_train = BBT0$signal$insample[ ,"bear"],
                                weights_fun = "l1",
                                sclfct = sclfct,
                                scl_by_entropy = TRUE )
  DS_bull <- dirichletSampling( Y_train = states, 
                                X_train = BBT0$signal$insample[ ,"bull"],
                                weights_fun = "l1",
                                sclfct = sclfct,
                                scl_by_entropy = TRUE )
  DS_all <- dirichletSampling( Y_train = states, 
                               X_train = BBT0$signal$insample[ ,"all"],
                               weights_fun = "l1",
                               sclfct = sclfct,
                               scl_by_entropy = TRUE )
  
  ldens_bear <- lapply( DS_bear, density )
  ldens_bull <- lapply( DS_bull, density )
  ldens_all <- lapply( DS_all, density )
  
  plot.ldensity( ldens_bear )  
  plot.ldensity( ldens_bull )  
  plot.ldensity( ldens_all )  
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # BBT OBJECT WITH LAG 
  # --------------------------------------------------------------------------
  
  require(DAARC)
  
  BBT <- loadSensor( sensor_name = "bbturbulence_scl2_dm_lag0to21" )
  
 
  
  
  signals <- na.omit( cbind( sig1 = BBT$signal$scl2[ ,"signal"],
                             sig2 = BBT$signal$base[ ,"signal"] ) )
  plot( signals )  
  
  plot(signals, plot.type = "single")
  lines( ema(signals[ ,1], 0.05), col = 3 )
  
  
  
  # Conditional densities
  # Y_train <- BBT$data$X_bm
  Y_train <- BBT$signal$insample[ ,"states"]
  sclfct <- NULL
  scl_by_entropy <- TRUE
  weights_fun <- "l1"
  lags <- c(0, 1, 2, 3)
  lDS1 <- lDS2 <- list()
  for ( i in seq(along = lags) ) {
    DS1 <- dirichletSampling( Y_train = Y_train,
                              X_train = signals[ ,1],
                              n_lag = lags[i],
                              weights_fun = weights_fun,
                              sclfct = sclfct,
                              scl_by_entropy = scl_by_entropy )
    DS2 <- dirichletSampling( Y_train = Y_train,
                              X_train = signals[ ,2],
                              n_lag = lags[i],
                              weights_fun = weights_fun,
                              sclfct = sclfct,
                              scl_by_entropy = scl_by_entropy )
    lDS1[[i]] <- DS1
    lDS2[[i]] <- DS2
  }
  
  ldens1 <- lapply( lDS1, FUN = function(x) { lapply( x, density) } ) 
  ldens2 <- lapply( lDS2, FUN = function(x) { lapply( x, density) } ) 
  
  
  par( mfrow = c(length(ldens1), 1) )
  for ( i in 1:length(ldens1) ) {
    plot.ldensity( ldens1[[i]] )
  }
  
  par( mfrow = c(length(ldens2), 1) )
  for ( i in 1:length(ldens2) ) {
    plot.ldensity( ldens2[[i]] )
  }
  
  dev.off()
  
  
  lmu1 <- lapply( lDS1, FUN = function(x) { lapply(x, mean) } )
  mu1 <- do.call( cbind, lapply(lmu1, unlist) )
  colnames(mu1) <- lags
  lmu2 <- lapply( lDS2, FUN = function(x) { lapply(x, mean) } )
  mu2 <- do.call( cbind, lapply(lmu2, unlist) )
  colnames(mu2) <- lags
  mu <- cbind( mu1, mu2 )
  
  
  barplot( t(mu), beside = TRUE, col = 1:ncol(mu) )  
  barplot( mu, beside = TRUE, col = fBasics::divPalette(n = nrow(mu), "RdYlGn") )  
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # COMPARE BBT SENSORS WITH DIFFERENT LAGS IN THE DATA
  # --------------------------------------------------------------------------
  
  require(DAARC)
  
  BBT <- loadSensor( sensor_name = "bbturbulence_base_dm" )
  BBT5 <- loadSensor( sensor_name = "bbturbulence_base_dm_lag0to5" )
  BBT10 <- loadSensor( sensor_name = "bbturbulence_base_dm_lag0to10" )
  BBT15 <- loadSensor( sensor_name = "bbturbulence_base_dm_lag0to15" )
  BBT21 <- loadSensor( sensor_name = "bbturbulence_base_dm_lag0to21" )
  BBT30 <- loadSensor( sensor_name = "bbturbulence_base_dm_lag0to30" )
  
  BBTscl2 <- loadSensor( sensor_name = "bbturbulence_scl2_dm" )
  BBT5scl2 <- loadSensor( sensor_name = "bbturbulence_scl2_dm_lag0to5" )
  BBT10scl2 <- loadSensor( sensor_name = "bbturbulence_scl2_dm_lag0to10" )
  BBT15scl2 <- loadSensor( sensor_name = "bbturbulence_scl2_dm_lag0to15" )
  BBT21scl2 <- loadSensor( sensor_name = "bbturbulence_scl2_dm_lag0to21" )
  BBT30scl2 <- loadSensor( sensor_name = "bbturbulence_scl2_dm_lag0to21" )
  
  BBT$computeSignalInsample()
  BBT5$computeSignalInsample()
  BBT10$computeSignalInsample()
  BBT15$computeSignalInsample()
  BBT21$computeSignalInsample()
  BBT30$computeSignalInsample()
  BBTscl2$computeSignalInsample()
  
  
  colnames(BBT5scl2$signal$insample)
  
  
  # In-sample
  sig_is <- na.omit( cbind( bbt = (1 + sign(BBT$signal$insample[ ,"delta"])) / 2,
                            bbtscl2 = (1 + sign(BBTscl2$signal$insample[ ,"delta"])) / 2,
                            bbt5 = (1 + sign(BBT5$signal$insample[ ,"delta"])) / 2,
                            bbt5scl2 = (1 + sign(BBT5scl2$signal$insample[ ,"delta"])) / 2,
                            bbt10 = (1 + sign(BBT10$signal$insample[ ,"delta"])) / 2,
                            bbt10scl2 = (1 + sign(BBT10scl2$signal$insample[ ,"delta"])) / 2,
                            bbt15 = (1 + sign(BBT15$signal$insample[ ,"delta"])) / 2,
                            bbt15scl2 = (1 + sign(BBT15scl2$signal$insample[ ,"delta"])) / 2,
                            bbt21 = (1 + sign(BBT21$signal$insample[ ,"delta"])) / 2,
                            bbt21scl2 = (1 + sign(BBT21scl2$signal$insample[ ,"delta"])) / 2,
                            bbt30 = (1 + sign(BBT30$signal$insample[ ,"delta"])) / 2,
                            bbt30scl2 = (1 + sign(BBT30scl2$signal$insample[ ,"delta"])) / 2 ) )
  # Out-of-sample
  sig_oos <- na.omit( cbind( bbt = (1 + sign(BBT$signal$base[ ,"md_delta"])) / 2,
                             bbtscl2 = (1 + sign(BBTscl2$signal$scl2[ ,"md_delta"])) / 2,
                             bbt5 = (1 + sign(BBT5$signal$base[ ,"md_delta"])) / 2,
                             bbt5scl2 = (1 + sign(BBT5scl2$signal$scl2[ ,"md_delta"])) / 2,
                             bbt10 = (1 + sign(BBT10$signal$base[ ,"md_delta"])) / 2,
                             bbt10scl2 = (1 + sign(BBT10scl2$signal$scl2[ ,"md_delta"])) / 2,
                             bbt15 = (1 + sign(BBT15$signal$base[ ,"md_delta"])) / 2,
                             bbt15scl2 = (1 + sign(BBT15scl2$signal$scl2[ ,"md_delta"])) / 2,
                             bbt21 = (1 + sign(BBT21$signal$base[ ,"md_delta"])) / 2,
                             bbt21scl2 = (1 + sign(BBT21scl2$signal$scl2[ ,"md_delta"])) / 2,
                             bbt30 = (1 + sign(BBT30$signal$base[ ,"md_delta"])) / 2,
                             bbt30scl2 = (1 + sign(BBT30scl2$signal$scl2[ ,"md_delta"])) / 2 ) )
  
  sig <- sig_is
  sig <- sig_oos
  
  X_bm <- BBT5$data$X_bm
  test <- signalTesting.byTrading( X = X_bm,
                                   sig = sig,
                                   n_lag = 1,
                                   tc = 0.004 )
  colnames(test) <- colnames(sig)
  test_nc <- signalTesting.byTrading( X = X_bm,
                                      sig = sig,
                                      n_lag = 1,
                                      tc = 0 )
  colnames(test_nc) <- colnames(sig)
  X_tmp <- na.omit(cbind(X_bm, test))
  X_tmp_nc <- na.omit(cbind(X_bm, test_nc))
  
  colors <- c(1, fBasics::rainbowPalette(ncol(test)))
  stats_fields <- c("cumret", "means", "sds", "sharpe", "maxDD")
  
  plot( as.simTS(X_tmp), col = colors )
  t( descStats( X_tmp )$stats[stats_fields, ] )
  
  plot( as.simTS(X_tmp_nc), col = colors )
  t( descStats( X_tmp_nc )$stats[stats_fields, ] )
  
  
  
  
  
  