  
  
  ############################################################################
  ### GOOD AND BAD TURBULENCE - OUT-OF-SAMPLE ANALYSIS - EQUITY WORLD
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     01.11.2020
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
  
  # ...
  
  
  # --------------------------------------------------------------------------
  # 1) DEFINE TURBULENCE MEASURES
  # --------------------------------------------------------------------------
  
  # --------------------------------------------------------------------------
  # Construct absolute and relative Mahalanobis distances
  # --------------------------------------------------------------------------
  
  BBT <- loadSensor( sensor_name = "bbturbulence_scl2_dm", 
                     b_new = TRUE )
  signals <- BBT$signal$scl2[ ,c("md_all", "md_delta", "signal")]
  signals <- cbind( signals, 
                    ema(signals, 0.05),
                    ema(signals, 0.0178) )
  
  
  head( tail(signals, 180), 20 )
  
  
  dd <- drawDowns( X = BBT$data$X_bm[isWeekday(time(BBT$data$X_bm)), ] )
  
  plot( cbind(dd, ema(BBT$signal$scl2[ ,c("md_all", "md_delta")], 0.02)),
        plot.type = "multiple" )

  
  
  ###
  BBT$spec$minphase <- 1
  BBT$spec$mincycle <- 2
  BBT$spec$k_peak <- BBT$spec$k_trough <- 
    BBT$spec$l_peak <- BBT$spec$l_trough <- BBT$spec$minphase
  # debugonce( BBT$computeSignalInsample )
  BBT$computeSignalInsample()
  states <- BBT$signal$insample[ ,"states"]
  BBT$signal$insample_stats
  
  
  plot( log(cumulated(BBT$data$X_bm, "discrete")) )
  abline( v = time(states)[ which(states == 1) ], col = 3 )
  abline( v = time(states)[ which(states == -1) ], col = 2 )
  lines( log(cumulated(BBT$data$X_bm, "discrete")) )
  
  
  
 
  
  
  
  BBT$spec$minphase
  BBT$spec$mincycle
  
  
  portmanteau(X = BBT$signal$insample[ ,"all"])
  
  
  
  # --------------------------------------------------------------------------
  # HEURISTIC BACKTEST RULE
  # --------------------------------------------------------------------------
  
  sig_tmp <- BBT$signal$scl[ ,c("md_delta", "signal")]
  sig_tmp <- ema(sig_tmp, 0.1)
  
  sig_tmp <- BBT$signal$insample[ ,"delta"]
  
  sig <- sig_tmp * 0 + 1
  sig[ sig_tmp < 0 ] <- 0
  test <- signalTesting.byTrading( X = BBT$data$X_bm,
                                   sig = sig, 
                                   n_lag = 1,
                                   tc = 0.00 )
  colnames(test) <- colnames(sig)
  X_tmp <- na.omit( cbind( bm = BBT$data$X_bm, test) )
  
  plot( as.simTS(X_tmp) )
  plot( sig_tmp )
  
  
  # Signal as difference in two moving averages
  sig_fast <- ema( BBT$signal$scl2[ ,"md_delta"], 0.05)
  sig_slow <- ema( BBT$signal$scl2[ ,"md_delta"], 0.01)
  signal <- sig_fast - sig_slow
  
  plot(cbind(sig_fast, sig_slow, signal), plot.type = "single")
  abline( h = 0 )
  
  sig_tmp <- cbind(fast = sig_fast, 
                   slow = sig_slow, 
                   signal = signal)
  sig <- sig_tmp * 0 + 1
  sig[ sig_tmp < 0 ] <- 0
  test <- signalTesting.byTrading( X = BBT$data$X_bm,
                                   sig = sig, 
                                   n_lag = 1,
                                   tc = 0.004)
  colnames(test) <- colnames(sig)
  X_tmp <- na.omit( cbind( bm = BBT$data$X_bm, test) )
  
  plot( as.simTS(X_tmp) )
  
  
  signal <- BBT$getSignal()
  plot(signal)
  abline( v = time(states)[ which(states == 1) ], col = 3 )
  abline( v = time(states)[ which(states == -1) ], col = 2 )
  abline( h = 0 )
  lines(signal)
  
  
  sig_fast <- log(ema(cumulated(signal, "discrete"), 1))
  sig_slow <- log(ema(cumulated(signal, "discrete"), 0.05))
  sig_delta <- sig_fast - sig_slow
  plot( cbind(sig_fast, sig_slow, sig_delta * 10), plot.type = "single")
  abline(h = 0)
  
  sig <- (1 + sign(sig_delta)) / 2
  test <- signalTesting.byTrading( X = BBT$data$X_bm,
                                   sig = sig, 
                                   n_lag = 1,
                                   tc = 0.004 )
  colnames(test) <- colnames(sig)
  X_tmp <- na.omit( cbind( bm = BBT$data$X_bm, test) )
  plot( as.simTS(X_tmp) )
  
  
  
  # --------------------------------------------------------------------------
  # LOGISTIC REGRESSION WITH LAGGED RELATIVE TURBULENCE VALUES
  # --------------------------------------------------------------------------
  
  BBT <- loadSensor( sensor_name = "bbturbulence_scl2_dm", 
                     b_new = TRUE )
  BBT$computeSignalInsample()
  
  signal <- BBT$signal$scl2[ ,"md_delta"]
  signals <- na.omit( lag( signal, 0:100 ) )
  colnames(signals) <- paste0("lag", 0:100 )
  
  reg <- regression( Y_train = (1 + BBT$signal$insample[ ,"states"]) / 2,
                     X_train = signals,
                     type = "elnet" )
  summary(reg$reg)
  
  plot( reg$z, reg$y )
  
  plot( reg$y )
  
  
  
  sig <- cbind(reg$y, round(reg$y))
  test <- signalTesting.byTrading( X = BBT$data$X_bm,
                                   sig = sig, 
                                   n_lag = 1,
                                   tc = 0.004 )
  colnames(test) <- colnames(sig)
  X_tmp <- na.omit( cbind( bm = BBT$data$X_bm, test) )
  plot( as.simTS(X_tmp) )
  
  
  
  # --------------------------------------------------------------------------
  # REINFORCEMENT
  # --------------------------------------------------------------------------
  
  BBT <- loadSensor( sensor_name = "bbturbulence_scl2_dm", 
                     b_new = TRUE )
  # BBT$computeSignalInsample()
  states_is <- BBT$signal$insample[ ,"states"]
  signal <- BBT$getSignal()
  
  
  dates <- intersect( rownames(signal), rownames(states_is) )
  error <- sign(signal[dates, ]) - states_is[dates, ] != 0
  error_adj <- error
  error_adj_ema <- ema(error, 0.1)
  states_adj <- states_is
  
  for ( today in dates[-c(1:10)] ) {
    
    error_adj_ema <- ema( error_adj[rownames(error_adj) <= today, ], 0.1 ) 
    if ( error_adj_ema[today, ] > 0.5 ) {
      states_adj[today, ] <- states_is[today, ] * (-1)
      error_adj[today, ] <- sign(signal[today, ]) - states_adj[today, ] != 0
    }
      
  }
  
  
  plot( error_adj )
  plot( cbind(states_is, sign(signal), states_adj) )
  
  
  sig <- cbind(states_is, sign(signal), states_adj)
  test <- signalTesting.byTrading( X = BBT$data$X_bm,
                                   sig = sig, 
                                   n_lag = 1,
                                   tc = 0.004 )
  colnames(test) <- colnames(sig)
  X_tmp <- na.omit( cbind( bm = BBT$data$X_bm, test) )
  plot( as.simTS(X_tmp) )
  
  
  
  
  
  # --------------------------------------------------------------------------
  # DEEP QDA (I.E., LAGGED DISTANCES AS INPUT)
  # --------------------------------------------------------------------------
  
  BBT <- loadSensor( sensor_name = "bbturbulence_scl2_dm", 
                     b_new = TRUE )
  # BBT$computeSignalInsample()
  signal <- BBT$signal$scl2[ ,"md_delta"]
  # signal <- BBT$data$X_bm
  # signal <- BBT$data$X
  signals <- na.omit( lag( signal, 0:21) )
  
  BBT2 <- BBT$copy()
  BBT2$spec$method <- "base"
  BBT2$spec$width <- 252 * 1     # ??
  BBT2$data$X <- signals
  BBT2$data$wmat <- NULL
  # BBT2$computeSignalInsample()
  # BBT2$updateSignal()
  
  plot( BBT2$signal$insample )
  plot( BBT2$signal$scl2[ ,1:10] )
  plot( BBT2$getSignal() )
  
  
  
  
  # Compare
  BBT1 <- loadSensor( "bbturbulence_scl2_dm_lag0to21", b_new = FALSE )
  BBT2 <- loadSensor( "bbturbulence_scl2_dm_lag0to21_onsignal", b_new = FALSE )
  
  sig1 <- (1 + sign(BBT1$getSignal())) / 2
  sig2 <- (1 + sign(BBT2$getSignal())) / 2
  sig <- cbind(sig1, sig2)
  test <- signalTesting.byTrading( X = BBT$data$X_bm,
                                   sig = sig, 
                                   n_lag = 1,
                                   tc = 0.004 )
  colnames(test) <- colnames(sig)
  test_nc <- signalTesting.byTrading( X = BBT$data$X_bm,
                                   sig = sig, 
                                   n_lag = 1,
                                   tc = 0 )
  colnames(test_nc) <- paste0(colnames(sig), "_nc")
  X_tmp <- na.omit( cbind( bm = BBT$data$X_bm, test, test_nc) )
  plot( as.simTS(X_tmp) )
  descStats(X_tmp)
  
  
  
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # MAP MD BEAR AND MD BULL TO PROBABILITIES
  # --------------------------------------------------------------------------
  
  Turb <- loadSensor( "turbulence" )
  Turb$computeSignalBase
  
  BBT <- BBT0$copy()
  BBT$data$wmat <- NULL
  BBT$spec$b_scl <- FALSE
  BBT$spec$method <- "base"
  BBT$computeSignalInsample()
  # MD <- BBT0$signal$scl2[ ,c("md_bear", "md_bull")]
  MD <- BBT$signal$insample[ ,c("bear", "bull", "all")]
  # MD <- scale(MD, FALSE, TRUE)
  P <- apply( MD, 2, function(x) { pchisq(x, df = ncol(BBT0$data$X)) } )
  # P <- apply( MD, 2, function(x) { pchisq(x, df = 1) } )
  sig <- (P - 1) * (-1)
  plot( sig ) 
  
  plot( sig[ ,1] - sig[ ,2] )
  
  plot( ema(sig, 0.05 ), plot.type = "single" )
  
  
  
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
  t(descStats(X_tmp)$stats[c("cumret", "maxDD"), ])
  
  
  
      
  # --------------------------------------------------------------------------
  # 3) CONDITIONAL PERFORMANCE
  # --------------------------------------------------------------------------
  
  # --------------------------------------------------------------------------
  # Posterior means
  # --------------------------------------------------------------------------
  
  start_date <- "2005-01-01"
  TD <- trainingData( Y_train = BBT$data$X_bm,
                      X_train = signals[rownames(signals) >= start_date, ] )
  lDS_kernel <- list()
  lDS_l1 <- list()
  lDS_knn <- list()
  lDS_qnn <- list()
  lDS_cmeans <- list()
  sclfct <- NULL
  
  for ( j in 1:ncol(TD$X_train) ) {
    
    lDS_kernel[[j]] <- dirichletSampling( Y_train = TD$Y_train,
                                          X_train = TD$X_train[ ,j],
                                          weights_fun = "kernel",
                                          correct_bias = FALSE,
                                          sclfct = sclfct )
    lDS_l1[[j]] <- dirichletSampling( Y_train = TD$Y_train,
                                      X_train = TD$X_train[ ,j],
                                      weights_fun = "l1",
                                      correct_bias = FALSE,
                                      sclfct = sclfct )
    lDS_knn[[j]] <- dirichletSampling( Y_train = TD$Y_train,
                                       X_train = TD$X_train[ ,j],
                                       weights_fun = "knn",
                                       k = 100,
                                       correct_bias = FALSE,
                                       sclfct = sclfct )
    lDS_qnn[[j]] <- dirichletSampling( Y_train = TD$Y_train,
                                       X_train = TD$X_train[ ,j],
                                       weights_fun = "qnn",
                                       alpha = 0.05,
                                       correct_bias = FALSE,
                                       sclfct = sclfct )
    lDS_cmeans[[j]] <- dirichletSampling( Y_train = TD$Y_train,
                                          X_train = TD$X_train[ ,j],
                                          weights_fun = "cmeans",
                                          centers = 2,
                                          correct_bias = FALSE,
                                          sclfct = sclfct )
    
  }
  
  lDS <- lDS_l1
  names(lDS) <- colnames(TD$X_train)
  
  
  # from <- quantile( unlist(lDS), 0 )
  # to <- quantile( unlist(lDS), 1 )
  from <- min( unlist(lDS) ) * 1.3
  to <- max( unlist(lDS) ) * 1.3
  n <- 999
  ldens <- lapply( lDS, FUN = function(x) { lapply(x, density, from = from, to = to, n = n) } )
  for ( i in 1:length(ldens) ) {
    plot.ldensity( ldens[[i]], 
                   main = paste0("Posterior means conditional on ", 
                                 names(ldens)[i]) )
  }
  
  
  
  
  lDS <- lDS_l1
  lMu <- lapply( lDS, FUN = function(x) { unlist(lapply(x, FUN = mean)) } )
  Mu <- do.call( cbind, lMu )
  
  colors <- rev(fBasics::divPalette(n = nrow(Mu), "Spectral"))
  barplot( Mu, beside = TRUE, col = colors )
  names(lMu)
  
  
  DF <- as.data.frame( cbind(md = do.call( cbind, lDS[[1]] ),
                             delta = do.call( cbind, lDS[[2]] )) )
  boxplot( DF, col = colors )
  
  par( mfrow = c(1, 2) )
  boxplot( DF[ ,1:length(lDS[[1]])], col = colors,
           ylim = range(DF), main = "Absolute Turbulence" )
  abline( h = 0, col = "lightgrey" )
  boxplot( DF[ ,-c(1:length(lDS[[2]]))], col = colors,
           ylim = range(DF), main = "Relative Turbulence" )
  abline( h = 0, col = "lightgrey" )
  
  dev.off()
  
  
  
  
  plot( x = as.numeric(ema(TD$X_train[ ,2], 0.1)),
        y = as.numeric(TD$Y_train) )
  
  
  
  
  
  
  
  # Performance after signal
  
  # debugonce( pas )
  start_date <- "2005-01-01"
  TD <- trainingData( Y_train = BBT$data$X_bm[rownames(BBT$data$X_bm) > start_date, ],
                      X_train = ema(signals, 1) )
  # sigAggFUN <- parse( text = "function(X) { sum(X) }" )
  lPAS <- pas( X = TD$Y_train,
               sig = TD$X_train[ ,"md_delta"],
               # sigAggFUN = sigAggFUN,
               sigAggFUN = sum,
               n_lags = c(0, 5, 10) )
  
  str( lPAS[[1]] )
  names(lPAS[[1]])
  
  
  i <- 2
  plot( x = as.numeric(lPAS[[i]]$X_train),
        y = as.numeric(lPAS[[i]]$Y_train) )
  abline(h = 0, col = "grey")
  cor(lPAS[[i]]$DF)
  
  reg <- regression( Y_train = lPAS[[i]]$Y_train,
                     X_train = lPAS[[i]]$X_train,
                     type = "ols" )
  summary(reg$reg)
  
  mu <- unlist( lapply( lPAS, FUN = function(x) { descStats(x$Y_train)$stats["means", ] } ) )
  barplot( mu )
  
  
  
  
  # --------------------------------------------------------------------------
  # 
  # --------------------------------------------------------------------------
  
  plot( signals )
  
  
  sig <- signals[ ,2*1-1] * signals[ ,2*1]
  plot( sig )
  TD <- trainingData( Y_train = BBT$data$X_bm,
                      X_train = sig )
  DS <- dirichletSampling( Y_train = TD$Y_train,
                           # X_train = TD$X_train,
                           X_train = signals[ ,2],
                           weights_fun = "l1",
                           correct_bias = FALSE,
                           sclfct = 1 )

  ldens <- lapply( DS, density )
  plot.ldensity( ldens, fillin = FALSE )
  
  
  
  

  sig <- ema( BBT$signal$scl2[ ,c("md_delta")], alpha = 0.05 )
  plot(sig)
  head( tail(sig, 180), 20 )
  

  CM <- cluster.cmeans( sig, centers = 2 )
  memb_scl <- apply( CM$membership, 2, function(x) { x / sum(x) } )
  CM$cm$centers
  plot( CM$cluster )
  plot( CM$membership )
  
  DS <- dirichletSampling( Y_train = TD$Y_train,
                           X_train = sig,
                           weights_mat = memb_scl,
                           correct_bias = FALSE,
                           sclfct = 1 )
  
  ldens <- lapply( DS, density )
  plot.ldensity( ldens, fillin = FALSE )
  
  
  # Or:
  DS <- dirichletSampling( Y_train = TD$Y_train,
                           X_train = sig,
                           X_eval = CM$cm$centers,
                           # weights_fun = "cmeans",
                           # centers = 2,
                           weights_fun = "kernel",
                           correct_bias = FALSE,
                           sclfct = 1 )
  
  ldens <- lapply( DS, density )
  plot.ldensity( ldens, fillin = FALSE )
  
  
  
  
  
  # --------------------------------------------------------------------------
  # CLASSIFICATION
  # --------------------------------------------------------------------------

  BBS <- bbs( X = cumulated(BBT$data$X_bm, "discrete"),
              minphase = BBT$spec$minphase,
              mincycle = BBT$spec$mincycle,
              theta = BBT$spec$theta,
              language = "C" )
  states <- (1 + BBS) / 2
  
  
  start_date <- "2001-01-01"
  alpha <- 0.0178
  TD <- trainingData( Y_train = states[rownames(states) >= start_date, ],
                      X_train = ema( BBT$signal$scl2[ ,"md_delta"], alpha ) )
  KM <- cluster.kmeans(x = TD$X_train, centers = 2)
  DS <- dirichletSampling( Y_train = TD$Y_train,
                           X_train = TD$X_train,
                           X_eval = sort(KM$km$centers),
                           weights_fun = "kernel",
                           correct_bias = FALSE,
                           sclfct = 1 )
  
  ldens <- lapply( DS, density )
  plot.ldensity( ldens, fillin = FALSE )
  
  
  
  alpha <- 0.02
  TD <- trainingData( Y_train = states[rownames(states) >= start_date, ],
                      X_train = ema( BBT$signal$scl2[ ,"md_delta"], alpha ) )
  X_eval <- tail(TD$X_train, 1)
  w1 <- weightsFun.kernel( data = TD, x_eval = X_eval )
  w2 <- (1 - w1 / max(w1)) / sum( (1 - w1 / max(w1)) )
  p1 <- t( rdirichlet( n = 1, alpha = w1 * length(w1) ) )
  p2 <- t( rdirichlet( n = 1, alpha = w2 * length(w2) ) )
  
  W <- cbind(w1, w2, p1, p2)
  as.numeric( t(TD$Y_train) %*% W )
  
  
  plot( as.timeSeries(W) )
  
  
  TD <- trainingData( Y_train = states[rownames(states) >= start_date, ],
                      X_train = BBT$signal$scl2[ ,"md_delta"] )
  # tau_vec <- round( seq(from = 0.1, to = 1, length.out = 100)^2, 5 )
  # tau_vec
  # FUN <- function(x) { ema( TD$X_train, alpha = x ) }
  # lSignal <- lapply( tau_vec, FUN = FUN )
  # signals <- do.call( cbind, lSignal_ema )
  X_eval <- tail(TD$X_train, 1)
  w1 <- weightsFun.kernel( data = TD, x_eval = X_eval )
  w2 <- (1 - w1 / max(w1)) / sum( (1 - w1 / max(w1)) )
  p1 <- t( rdirichlet( n = 1, alpha = w1 * length(w1) ) )
  p2 <- t( rdirichlet( n = 1, alpha = w2 * length(w2) ) )
  
  W <- cbind(w1, w2, p1, p2)
  as.numeric( t(TD$Y_train) %*% W )
  
  
  X_eval
  range(TD$X_train)
  tail(TD$X_train)
  plot( tail(TD$X_train, 200) )
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # HIERARCHICHAL QDA
  # --------------------------------------------------------------------------
  
  BBT$computeSignalInsample()
  
  # Define prior according to logistic regression
  reg <- regression( Y_train = states,
                     X_train = ema(BBT$signal$insample[ ,c("all")], 
                                   0.02),
                     type = "logit" )
  summary(reg$reg)
  
  plot( cbind(reg$y, states), plot.type = "single" )
  plot( x = reg$z, y = reg$y )
  mean(reg$y); mean(states)
  
  
  plot( log( reg$y / (1 - reg$y) ) )
  abline(h = 0 )
  
  plot( log( (1 - reg$y) / reg$y ) )
  abline(h = 0 )
  
  plot( BBT$signal$scl2[ ,"md_delta"])
  lines( log( reg$y / (1 - reg$y) ), col = 2 )
  lines( log( (1 - reg$y) / reg$y ), col = 3 )
  
  
  
  X <- BBT$data$X
  BBS <- BBT$signal$insample[ ,"states"]
  qda_posterior <- BBS * NA
  qda_class <- BBS * NA
  for ( i in 1:nrow(X) ) {
    p <- reg$y[rownames(X)[i], ]
    QDA <- qda( x = X,
                grouping = BBS,
                prior = c(1 - p, p) )
    pred <- predict( QDA, X[i, ] )
    qda_class[i, ] <- as.numeric( pred$class )
    qda_posterior[i, ] <- as.numeric( pred$posterior[2] )
  }
  
  plot( qda_posterior )
  plot( qda_class )
  
  
  
  
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # OPTIMIZE OVER PARAMETER GRID
  # --------------------------------------------------------------------------
  
  # How to smooth the signal such that the conditional probabilities
  # for either class are as distinctive as possible. What conditioning 
  # values to use, i.e., what threshold to use?
  
  # Grid, threshold and smoothing parameter
  
  n_lag <- 1
  tc <- 0.004
  tau_vec <- round( seq(from = 0.1, to = 1, length.out = 10)^2, 5 )
  tau_vec
  FUN <- function(x) { ema( TD$X_train, alpha = x ) }
  lSignal <- lapply( tau_vec, FUN = FUN )
  signals <- do.call( cbind, lSignal_ema )
  
  
  # --------------------------------------------------------------------------
  objFun <- function( x )
  {
    sig <- ema( TD$X_train, alpha = x[1] )
    signal <- sig * 0
    # signal[ sig < x[2] ] <- 1
    signal[ sig > 0 ] <- 1
    ret <- signalTesting.byTrading( X = TD$Y_train,
                                    sig = signal,
                                    n_lag = 1,
                                    tc = 0.004 )
    # X_tmp <- na.omit( cbind(TD$Y_train, ret) )
    # plot( as.simTS(X_tmp) )
    dd <- drawDowns(ret)
    cum <- exp( sum( log( 1 + ret ) ) ) - 1
    ans <- (cum + min(dd)) * (-1)
    return( ans )
  }
  
  
  
  TD <- trainingData( Y_train = BBT$data$X_bm,
                      X_train = BBT$signal$scl2[ ,"md_all"] )
  # control <- Rsolnp:::.solnpctrl(list())
  opt <- solnp( pars = c(0.1, quantile(TD$X_train, 0.9)),
                fun = objFun,
                LB = c(0, 0) )
  opt$values
  objFun( x = opt$values )
  
  
  opt <- donlp2( par = c(0.9, 100),
                 fn = objFun,
                 par.lower = c(0.001, min(TD$X_train)),
                 par.upper = c(1, max(TD$X_train)) )
  
  opt$fx
  opt$par
  objFun( x = opt$par )
  
  
  
  
  
  
  
  TD <- trainingData( Y_train = BBT$data$X_bm,
                      X_train = BBT$signal$scl2[ ,"md_delta"] ) 
  control <- Rsolnp:::.solnpctrl(list())
  
  # debugonce(solnp)
  opt <- solnp( pars = c(1, 0),
                fun = objFun,
                control = control,
                LB = c(0.001, min(TD$X_train)),
                UB = c(1, max(TD$X_train)) )
  opt$values
  
  objFun( x = opt$values )
  
  
  opt <- optimize( f = objFun, interval = c(0, 1), maximum = FALSE )
  opt
  objFun( x = opt$minimum )
  
  
  
  
  
  # --------------------------------------------------------------------------
  objFun2 <- function( x )
  {
    sig <- ema( TD$X_train, alpha = x[1] )
    signal <- sig * 0
    # signal[ sig < x[2] ] <- 1
    signal[ sig > 0 ] <- 1
    ret <- signalTesting.byTrading( X = TD$Y_train,
                                    sig = signal,
                                    n_lag = 1,
                                    tc = 0.004 )
    # X_tmp <- na.omit( cbind(TD$Y_train, ret) )
    # plot( as.simTS(X_tmp) )
    cum <- exp( sum( log( 1 + ret ) ) ) - 1
    ans <- cum * (-1)
    return( ans )
  }
  # --------------------------------------------------------------------------
  ineqFun <- function(x)
  {
    sig <- ema( TD$X_train, alpha = x[1] )
    signal <- sig * 0
    # signal[ sig < x[2] ] <- 1
    signal[ sig > 0 ] <- 1
    ret <- signalTesting.byTrading( X = TD$Y_train,
                                    sig = signal,
                                    n_lag = 1,
                                    tc = 0.004 )
    dd <- drawDowns(ret)
    ans <- min(dd)
    return( ans )
  }
  
  TD <- trainingData( Y_train = BBT$data$X_bm,
                      X_train = BBT$signal$scl2[ ,"md_delta"] ) 
  control <- Rsolnp:::.solnpctrl(list())
  
  # debugonce(solnp)
  opt <- solnp( pars = 0.018,
                fun = objFun2,
                ineqfun = ineqFun,
                ineqLB = -0.4,
                ineqUB = 1,
                control = control,
                LB = 0.001,
                UB = 9.999 )
  opt$values
  objFun( x = opt$values )
  
  
  
  # USING DEOptim
  
  # Penalty
  penaltyFun <- function( max_dd )
  {
    penalty <- 0
    if ( max_dd < -0.4 ) {
      penalty <- 1e05
    }
    return( penalty )
  }
  
  # Objective function
  fn <- function( x ) 
  {
    sig <- ema( TD$X_train, alpha = x[1] )
    signal <- sig * 0
    signal[ sig > x[2] ] <- 1
    # signal[ sig < x[2] ] <- 1
    ret <- signalTesting.byTrading( X = TD$Y_train,
                                    sig = signal,
                                    n_lag = 1,
                                    tc = 0.004 )
    # X_tmp <- na.omit( cbind(TD$Y_train, ret) )
    # plot( as.simTS(X_tmp) )
    obj <- exp( sum( log( 1 + ret ) ) ) - 1
    dd <- drawDowns(ret)
    # penalty <- penaltyFun( max_dd = min(dd) )
    # ans <- -obj + penalty
    dd_stats <- na.omit( drawDownStats(X = ret, to = 5)[ ,"Depth"],
                         method = "z" )
    penalty <- sum( dd_stats * (5:1 / sum(5:1)) * 5 )
    ans <- -obj - penalty
    return( ans )
  }
  
  TD <- trainingData( Y_train = BBT$data$X_bm,
                      X_train = BBT$signal$scl2[ ,"md_delta"] )
                      # X_train = BBT$signal$scl2[ ,"md_all"] ) 
  initialpop <- expand.grid( list(a = seq(from = 0.05, to = 0.99, length.out = 5),
                                  # b = seq(from = quantile(TD$X_train, 0.1), 
                                  #         to = quantile(TD$X_train, 0.9), 
                                  #         length.out = 5)) )
                                  b = seq(from = -0.01,
                                          to = 0.01,
                                          length.out = 5)) )
                                  # b = seq(from = quantile(TD$X_train, 0.7), 
                                  #         to = quantile(TD$X_train, 0.99), 
                                  #         length.out = 5)) )
  DEcntrl <- DEoptim::DEoptim.control(NP = nrow(initialpop),
                                      initialpop = as.matrix(initialpop),
                                      itermax = 200)
  lower <- c(0.001, min(TD$X_train))
  upper <- c(0.999, max(TD$X_train))
  opt <- DEoptim( fn = fn, 
                  lower = lower, 
                  upper = upper, 
                  control = DEcntrl )
  
  opt$optim$bestmem
  opt$optim$bestval
  fn( x = opt$optim$bestmem )
  
  
  
  
  # Using absolute and relative distance
  
  # Penalty
  penaltyFun <- function( max_dd )
  {
    penalty <- 0
    if ( max_dd < -0.2 ) {
      penalty <- 1e05
    }
    return( penalty )
  }
  
  # Objective function
  fn <- function( x ) 
  {
    sig1 <- ema( TD$X_train[ ,1], alpha = x[1] )
    sig2 <- ema( TD$X_train[ ,2], alpha = x[2] )
    signal <- sig1 * 0 + 1
    idx <- intersect( which(sig1 > x[3]), which(sig2 < x[4]) )
    signal[idx, ] <- 0
    ret <- signalTesting.byTrading( X = TD$Y_train,
                                    sig = signal,
                                    n_lag = 1,
                                    tc = 0.004 )
    # X_tmp <- na.omit( cbind(TD$Y_train, ret) )
    # plot( as.simTS(X_tmp) )
    obj <- exp( sum( log( 1 + ret ) ) ) - 1
    dd <- drawDowns(ret)
    # penalty <- penaltyFun( max_dd = min(dd) )
    # ans <- -obj + penalty
    dd_stats <- na.omit( drawDownStats(X = ret, to = 10)[ ,"Depth"],
                         method = "z" )
    n <- sum( dd_stats < -0.1 )
    penalty <- sum( dd_stats[1:n] * (n:1 / sum(n:1)) * n )
    ans <- -obj - penalty
    return( ans )
  }
  
  TD <- trainingData( Y_train = BBT$data$X_bm,
                      X_train = BBT$signal$scl2[ ,c("md_all", "md_delta")] )
  n <- 3
  initialpop <- expand.grid( list( x1 = seq(from = 0.01, to = 0.99, length.out = n),
                                   x2 = seq(from = 0.01, to = 0.99, length.out = n),
                                   x3 = seq(from = 5, to = 10, length.out = n),
                                   x4 = seq(from = -0.01, to = 0.01, length.out = n) ) )
  DEcntrl <- DEoptim::DEoptim.control(NP = nrow(initialpop),
                                      initialpop = as.matrix(initialpop),
                                      itermax = 200)
  lower <- c(0.001, 0.001, min(TD$X_train[ ,1]), min(TD$X_train[ ,2]))
  upper <- c(0.999, 0.999, max(TD$X_train[ ,1]), max(TD$X_train[ ,2]))
  opt2 <- DEoptim( fn = fn, 
                   lower = lower, 
                   upper = upper, 
                   control = DEcntrl )
  
  opt2$optim$bestmem
  opt2$optim$bestval
  fn( x = opt2$optim$bestmem )
  
  
  apply( initialpop, 1, FUN = fn )
  
  
  
  
  ###############
  
  fn <- function( x ) 
  {
    sig_ema <- ema( TD$X_train, x[1] )
    sig_ema_sigmoid <- 1 / (1 + exp(-(sig_ema - x[2]) * x[3]))
    p_sig <- sig_ema_sigmoid / sum(sig_ema_sigmoid)
    # plot(cbind(TD$Y_train, p_sig), plot.type = "single")
    ans <- relentropy( p = TD$Y_train, 
                       q = p_sig, 
                       exponential = TRUE )
    return( ans )
  }
  
  states <- BBT$signal$insample[ ,"states"]
  TD <- trainingData( Y_train = states / sum(states),
                      X_train = BBT$signal$scl2[ ,"md_delta"] )
  n <- 4
  initialpop <- expand.grid( list( x1 = seq(from = 0.01, to = 0.99, length.out = n),
                                   x2 = seq(from = -0.01, to = 0.1, length.out = 2),
                                   x3 = seq(from = 10^2, to = 10^3, length.out = n) ) )
  initialpop <- rbind( initialpop, c(0.018, 0, 500) )
  DEcntrl <- DEoptim::DEoptim.control(NP = nrow(initialpop),
                                      initialpop = as.matrix(initialpop),
                                      itermax = 200)
  lower <- c(0.001, -1, 1)
  upper <- c(0.999, 1, 10^4)
  opt3 <- DEoptim( fn = fn, 
                   lower = lower, 
                   upper = upper, 
                   control = DEcntrl )
  
  opt3$optim$bestmem
  opt3$optim$bestval
  fn( x = opt3$optim$bestmem )
    
  tmp <- apply( initialpop, 1, fn )
  plot(tmp)  
  

  
  sig <- p_sig / max(p_sig)
  # sig <- TD$Y_train / max(TD$Y_train)
  test <- signalTesting.byTrading( X = BBT$data$X_bm,
                                   sig = sig,
                                   n_lag = 1,
                                   tc = 0.00 )
  X_tmp <- na.omit( cbind( BBT$data$X_bm, test ) )
  plot( as.simTS(X_tmp) )
  
  
  
  
  
  
  
  ####################
  
  
  plot( BBT$signal$scl2 )
  
    
  
  
  
  
  