  


  # --------------------------------------------------------------------------
  # Old
  # --------------------------------------------------------------------------
  
  X_bm_w <- aggWeekly( X = X_bm, day = "Tue", compounding = "discrete" )
  signal_w <- aggWeekly( X = signal, day = "Tue", compounding = "continuous" )  
  
  DS <- dirichletSampling( Y_train = log(1 + X_bm_w),
                           X_train = signal_w,
                           n_lag = 0,
                           weights_fun = "kernel",
                           # weights_fun = "knn", "qnn", "cmeans",
                           # weights_fun = "qnn",
                           # alpha = 0.1,
                           # weights_fun = "cmeans",
                           # centers = 2,
                           correct_bias = FALSE,
                           # apply_roll_y = apply_roll_y,
                           sclfct = 1 )
  
  ldens <- lapply( DS, densFUN )
  plot.ldensity( ldens )
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Monthly data
  # --------------------------------------------------------------------------
  
  X_bm_m <- aggMonthly( X = X_bm, compounding = "discrete" )
  signal_m <- aggMonthly( X = signal, compounding = "continuous" )  
  
  DS <- dirichletSampling( Y_train = log(1 + X_bm_m),
                           X_train = signal_m,
                           n_lag = 1,
                           weights_fun = "kernel",
                           # weights_fun = "knn", "qnn", "cmeans",
                           # weights_fun = "qnn",
                           # alpha = 0.1,
                           # weights_fun = "cmeans",
                           # centers = 2,
                           correct_bias = FALSE,
                           # apply_roll_y = apply_roll_y,
                           sclfct = NULL )
  
  ldens <- lapply( DS, densFUN )
  plot.ldensity( ldens )
  
  
  
  # --------------------------------------------------------------------------
  # Bear Bull as dependent variable
  # --------------------------------------------------------------------------
  
  BBSObj <- loadSensor( sensor_name = "bbs", b_new = TRUE )
  BBS <- BBSObj$signal$base
  
  # ...
  
  
  
  # --------------------------------------------------------------------------
  # Fuzzy Bear Bull classification to generate the signal
  # --------------------------------------------------------------------------
  
  # Generate MD_bear and MD_bull series based on a 
  # fuzzy bull membership index and it's complement.
  
  BBT <- loadSensor( sensor_name = "bbturbulence" )
  BBT$computeSignalBase
  
  # In-sample
  
  BBSF <- loadSensor( sensor_name = "bbsfuzzy" )
  I <- (1 + BBSF$signal$insample) / 2
  IC <- 1 - I
  
  BBFT <- loadSensor( sensor_name = "bbfturbulence" )
  I <- (1 + BBFT$signal$insample) / 2
  IC <- 1 - I
  
  BBFT$signal$`phase=5cycle=21`
  
  
  bbft <- do.call( cbind, lapply( BBFT$signal, function(x) { x[ ,"signal"] } ) )
  plot( ema(apply(bbft, 1, mean), 0.1) )
  plot( ema(apply(bbft, 1, sum), 0.1) )
  plot( ema(bbft, 0.1), plot.type = "single" )
  
  
  
  
  
  # edit(BBT$computeSignalBase)
  data <- BBT$data
  md_all <- data$X_bm * NA
  md_bear <- md_bull <- md_delta <- md_all
  
  X_train <- data$X
  
  # Compute mahalanobis distance on sensor signals 
  # during bear markets
  p_bear <- IC / sum(IC)
  TD <- trainingData(Y_train = p_bear, 
                     X_train = X_train)
  mu_bear <- t(TD$X_train) %*% (TD$Y_train / sum(TD$Y_train))
  scnd_mom = t(TD$X_train) %*% (TD$X_train * ((TD$Y_train / sum(TD$Y_train)) %*% matrix(1, 1, ncol(TD$X_train))) )
  scnd_mom = (scnd_mom + t(scnd_mom)) / 2
  sigma_bear = scnd_mom - mu_bear %*% t(mu_bear)
  
  
  # Compute mahalanobis distance on sensor signals 
  # during bull markets
  p_bull <- I / sum(I)
  TD <- trainingData(Y_train = p_bull, 
                     X_train = X_train)
  mu_bull <- t(TD$X_train) %*% (TD$Y_train / sum(TD$Y_train))
  scnd_mom = t(TD$X_train) %*% (TD$X_train * ((TD$Y_train / sum(TD$Y_train)) %*% matrix(1, 1, ncol(TD$X_train))) )
  scnd_mom = (scnd_mom + t(scnd_mom)) / 2
  sigma_bull = scnd_mom - mu_bull %*% t(mu_bull)
  
  
  wghts <- NULL
  # if ( !is.null(data$wmat) ) {
  #   wghts <- as.numeric(data$wmat[today, ])
  # }
  md_bear_tmp <- weightedMahalanobis(x = X_train,
                                     center = mu_bear,
                                     covmat = sigma_bear,
                                     wghts = wghts,
                                     scl = TRUE)
  md_bull_tmp <- weightedMahalanobis(x = X_train,
                                     center = mu_bull,
                                     covmat = sigma_bull,
                                     wghts = wghts,
                                     scl = TRUE)
  md_all_tmp <- weightedMahalanobis(x = X_train,
                                    # center = mu, 
                                    # covmat = covmat, 
                                    wghts = wghts,
                                    scl = TRUE)
  md_bear_tmp <- as.timeSeries(md_bear_tmp)
  md_bull_tmp <- as.timeSeries(md_bull_tmp)
  md_all_tmp <- as.timeSeries(md_all_tmp)
  
  # Scale and take difference between bear and bull distances
  MD <- cbind(bear = md_bear_tmp,
              bull = md_bull_tmp)
  MD_scl <- scale(MD, FALSE, TRUE)
  md_delta_tmp <- MD_scl[ ,"bear"] - MD_scl[ ,"bull"]
  
  plot(MD_scl, plot.type = "single" )
  plot(md_delta_tmp)
  
  
  
  
  DS <- dirichletSampling( Y_train = log(1 + X_bm),
                           X_train = md_delta_tmp,
                           n_lag = 0,
                           weights_fun = "kernel",
                           correct_bias = FALSE,
                           sclfct = 1 )
  
  ldens <- lapply( DS, densFUN )
  plot.ldensity( ldens )
  
  
  
  
  # --------------------------------------------------------------------------
  # Fuzzy Bear Bull classification to generate the signal
  # --------------------------------------------------------------------------
  
  # Generate resampled MD_bear and MD_bull series based on a 
  # fuzzy bull membership index and it's complement.
  
  n_sim <- 999
  
  # Bear
  TD <- trainingData(Y_train = p_bear, 
                     X_train = X_train)
  # sclfct <- 1
  # sclfct <- length(TD$Y_train)
  sclfct <- entropy( TD$Y_train, exponential = TRUE )
  P <- rdirichlet( n = n_sim, alpha = TD$Y_train  * sclfct )
  
  mu_bear_fuzzy <- P %*% TD$X_train
  lSigma_bear_fuzzy <- list()
  lMD_bear_fuzzy <- list()
  for ( i in 1:nrow(P) ) {
    mu <- mu_bear_fuzzy[i, ]
    p <- as.numeric(P[i, ])
    scnd_mom = t(TD$X_train) %*% (TD$X_train * p %*% matrix(1, 1, ncol(TD$X_train)))
    scnd_mom = (scnd_mom + t(scnd_mom)) / 2
    lSigma_bear_fuzzy[[i]] = scnd_mom - mu %*% t(mu)
    sigma <- make.positive.definite(lSigma_bear_fuzzy[[i]])
    lMD_bear_fuzzy[[i]] <- weightedMahalanobis(x = TD$X_train,
                                               center = mu,
                                               covmat = sigma,
                                               wghts = wghts,
                                               scl = TRUE)
  }
  md_bear_fuzzy <- timeSeries( do.call( cbind, lMD_bear_fuzzy ),
                               time(TD$X_train) )
  md_bear_fuzzy <- scale(md_bear_fuzzy, FALSE, TRUE)
  
  # Bull
  TD <- trainingData(Y_train = p_bull, 
                     X_train = X_train)
  # sclfct <- 1
  # sclfct <- length(TD$Y_train)
  sclfct <- entropy( TD$Y_train, exponential = TRUE )
  P <- rdirichlet( n = n_sim, alpha = TD$Y_train  * sclfct )
  
  mu_bull_fuzzy <- P %*% TD$X_train
  lSigma_bull_fuzzy <- list()
  lMD_bull_fuzzy <- list()
  for ( i in 1:nrow(P) ) {
    mu <- mu_bull_fuzzy[i, ]
    p <- as.numeric(P[i, ])
    scnd_mom = t(TD$X_train) %*% (TD$X_train * p %*% matrix(1, 1, ncol(TD$X_train)))
    scnd_mom = (scnd_mom + t(scnd_mom)) / 2
    lSigma_bull_fuzzy[[i]] = scnd_mom - mu %*% t(mu)
    sigma <- make.positive.definite(lSigma_bull_fuzzy[[i]])
    lMD_bull_fuzzy[[i]] <- weightedMahalanobis(x = TD$X_train,
                                               center = mu,
                                               covmat = sigma,
                                               wghts = wghts,
                                               scl = TRUE)
  }
  md_bull_fuzzy <- timeSeries( do.call( cbind, lMD_bull_fuzzy ),
                               time(TD$X_train) )
  md_bull_fuzzy <- scale(md_bull_fuzzy, FALSE, TRUE)
  
  
  # ldens <- lapply( lMD_bull_fuzzy, FUN = densFUN )
  # plot.ldensity( ldens )    
  # plot( md_bull_fuzzy, plot.type = "single" )      
  
  
  
  dens_bear <- densFUN( md_bear_fuzzy )
  dens_bull <- densFUN( md_bull_fuzzy )
  # dens_bear <- densFUN( apply(md_bear_fuzzy, 1, mean) )
  # dens_bull <- densFUN( apply(md_bull_fuzzy, 1, mean) )
  
  ldens <- list(bear = dens_bear, 
                bull = dens_bull)
  plot.ldensity( ldens )
  
  
  plot( cbind(apply(md_bear_fuzzy, 1, mean),
              apply(md_bull_fuzzy, 1, mean)), 
        plot.type = "single" )
  
  
  md_delta_fuzzy <- md_bear_fuzzy - md_bull_fuzzy
  md_delta_fuzzy_ema <- ema( md_delta_fuzzy, alpha = 0.1 )
  
  plot( apply(md_delta_fuzzy_ema, 1, mean), plot.type = "single" )
  
  
  
  ldens <- lapply( 1:nrow(md_delta_fuzzy), 
                   FUN = function(idx) { densFUN(md_delta_fuzzy[idx, ]) } )  
  names(ldens) <- rownames(md_delta_fuzzy)
  plot.ldensity( ldens )
  
  
  
  
  DS <- dirichletSampling( Y_train = log(1 + X_bm),
                           X_train = apply(md_delta_fuzzy_ema, 1, mean),
                           n_lag = 0,
                           weights_fun = "kernel",
                           correct_bias = FALSE,
                           sclfct = 1 )
  
  ldens <- lapply( DS, densFUN )
  plot.ldensity( ldens )
  
  
  
  
  
  
  
  
