  
  
  ############################################################################
  ### BBTurbulence - SZENARIOS SIMULATION 
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     16.10.2020
  # First version:    10.10.2020
  # --------------------------------------------------------------------------

  
  # --------------------------------------------------------------------------
  # REQUIRE
  # --------------------------------------------------------------------------
  
  require(rugarch)
  require(DAARC)
  
  wd <- "R:/Asset_Management/Research_Projects/Multi_Asset/DAA/DAARC/Sensor_Testing/Univariate/Good_And_Bad_Turbulence/"
  source( paste0(wd, "Source/source_functions.R") )
  
  
  
  # --------------------------------------------------------------------------
  # IN-SAMPLE BEAR BULL MAHALANOBIS ANALYSIS
  # --------------------------------------------------------------------------
  
  BBT <- loadSensor( sensor_name = "bbturbulence" )
  X_bm <- BBT$data$X_bm
  Y <- log( 1 + BBT$data$X[isWeekday(time(BBT$data$X)), ] )
  
  # Bears and Bulls on eqw series
  BBS <- bbs(X = cumulated(X_bm, "discrete"), 
             mincycle = mincycle_d, 
             minphase = minphase_d, 
             k.peak = minphase_d, 
             l.peak = minphase_d, 
             k.trough = minphase_d, 
             l.trough = minphase_d, 
             logarithmic = FALSE, 
             theta = theta, 
             e = 0)
  
  # Bear 
  p_bear <- (1 - BBS) / 2
  TD <- trainingData(Y_train = p_bear, 
                     X_train = Y)
  mu_bear <- t(TD$X_train) %*% (TD$Y_train / sum(TD$Y_train))
  scnd_mom = t(TD$X_train) %*% (TD$X_train * ((TD$Y_train/sum(TD$Y_train)) %*% 
                                                matrix(1, 1, ncol(TD$X_train))))
  scnd_mom = (scnd_mom + t(scnd_mom))/2
  sigma_bear = scnd_mom - mu_bear %*% t(mu_bear)
  md_bear <- mahalanobisTS(X = Y, 
                           center = as.numeric(mu_bear), 
                           covmat = sigma_bear,
                           scl = FALSE)
  
  # Bull
  p_bull <- (1 + BBS) / 2
  TD <- trainingData(Y_train = p_bull, 
                     X_train = Y)
  mu_bull <- t(TD$X_train) %*% (TD$Y_train / sum(TD$Y_train))
  scnd_mom = t(TD$X_train) %*% (TD$X_train * ((TD$Y_train/sum(TD$Y_train)) %*% 
                                                matrix(1, 1, ncol(TD$X_train))))
  scnd_mom = (scnd_mom + t(scnd_mom))/2
  sigma_bull = scnd_mom - mu_bull %*% t(mu_bull)
  md_bull <- mahalanobisTS(X = Y, 
                           center = as.numeric(mu_bull), 
                           covmat = sigma_bull,
                           scl = FALSE)
  # Full sample
  md_all <- mahalanobisTS(Y)
  
  MD <- cbind( bear = md_bear,
               bull = md_bull,
               all = md_all )
  MD_scl <- scale( MD, FALSE, TRUE )
  md_delta <- MD_scl[ ,"bear"] - MD_scl[ ,"bull"]
  
  ans <- cbind( MD_scl, delta = md_delta )
  plot( ans )
  
  
  
  
  # --------------------------------------------------------------------------
  # PARAMETERS
  # --------------------------------------------------------------------------
  
  # BBQ parameters
  minphase_m <- 2
  mincycle_m <- 5
  minphase_w <- minphase_m * 4
  mincycle_w <- mincycle_m * 4
  minphase_d <- minphase_m * 4 * 5
  mincycle_d <- mincycle_m * 4 * 5
  theta <- 0.15
  
  # Szenario parameters
  M1 <- apply(Y, 2, mean)
  #
  M1 <- abs(M1) * 2 # ~~~~~~~~~~~~~~ 
  #
  M2 <- cov(Y)
  n_sample = 10^3 * 10
  asTS = TRUE
  start_date = NULL
  frqncy = "day"
  seed = 1234
  variance_model = list(model = "gjrGARCH",
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
  distribution_model = "sstd"
  # start_pars = list(alpha1 = 0.05,
  #                   beta1 = 0.99,
  #                   gamma1 = 0.02)
  start_pars = list()
  # fixed_pars = list(alpha1 = 0.05,
  #                   beta1 = 0.93,
  #                   gamma1 = 0.02)
  fixed_pars = list()
  
  
  
  
  # --------------------------------------------------------------------------
  # SZENARIO SIMULATION
  # --------------------------------------------------------------------------
  
  BBT <- loadSensor( sensor_name = "bbturbulence" )
  X <- BBT$data$X
  X_bm <- BBT$data$X_bm
  
  X_eqw  <- timeSeries( X %*% rep(1/ncol(X), ncol(X)), time(X) )
  spec <- ugarchspec( variance.model = variance_model,
                      mean.model = mean_model,
                      distribution.model = distribution_model,
                      start.pars = start_pars,
                      fixed.pars = fixed_pars )
  fit <- ugarchfit( data = log(1 + X_eqw),
                    spec = spec )
  sim <- ugarchsim( fit, 
                    n.sim = nrow(X_eqw), 
                    n.start = 1, 
                    m.sim = 1, 
                    startMethod = "sample",
                    rseed = 40 )  
  sim <- asTimeSeries( sim@simulation$seriesSim, frqncy = "d" )
  sim <- exp(sim) - 1
  
  plot( as.simTS(sim) )  
  cvol <- getCondVar(garch(sim))
  
  
  # alpha <- cvol / sum(cvol)
  # alpha <- abs(sim) / sum(abs(sim))
  alpha <- weightsFun.l1( data = list(X_train = abs(sim)), x_eval = quantile(abs(sim), 0.75) )
  samples <- rdirichlet(n = 20, alpha = alpha * length(alpha) )  
  P <- as.timeSeries( t(samples), time(sim) )
  P <- ema( P, 1 )
  
  # Z_bm <- X_bm[rownames(X), ]
  Z_bm <- sim
  x <- as.numeric(Z_bm)
  tmp <- apply(P, 2, function(p) { x * p * length(p) } )
  X_synt <- as.timeSeries( tmp, time(Z_bm) )  
  
  # p <- as.numeric(P[1, ])
  # tmp <- apply( X, 2, function(x) { x * p * length(p) } )
  # X_synt <- as.timeSeries( tmp, time(X_bm) )  
  
  
  
  plot(alpha) 
  plot(cvol)      
  plot(P[ ,1:10])
  plot( as.simTS(X_synt) )
  plot( X_synt[ ,1:10] )
  plot( apply(P, 1, mean) )  
  
  cor(X_synt)
  
  
  apply( X_synt^2, 2, portmanteau )
  fit <- garch(X_synt)
  
  
  
  
  # --------------------------------------------------------------------------
  # BEAR BULL SZENARIO ANALYSIS
  # --------------------------------------------------------------------------
  
  BBTSYNT <- BBTurbulence$new()
  # BBTSYNT <- BBTEWMA$new()
  # BBTSYNT <- BBTSV$new()
  # BBTSYNT <- BBTGARCH$new()
  BBTSYNT$setCtrl( method = "base", universe = "dm" )
  BBTSYNT$spec$iso <- colnames(X_synt)
  BBTSYNT$data <- list( X = X_synt,
                        X_bm = apply(X_synt, 1, mean) )
  BBTSYNT$computeSignalInsample()
 
  plot(BBTSYNT$signal$insample)
  
  
  P <- apply( ema(BBTSYNT$signal$insample[ ,c("bear", "bull", "all")], 0.1), 2, 
              function(x) { pchisq(x, df = ncol(X)) } )
  sig <- (P - 1) * (-1)
  
  # Signal testing
  Z_bm <- BBTSYNT$data$X_bm
  sig <- (1 + sign( ema( BBTSYNT$signal$insample[ ,"delta"], 0.1 ) ) ) / 2

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
  X_tmp <- na.omit(cbind(bm = Z_bm, bt = test, bt_nc = test_nc))
  
  plot( as.simTS(X_tmp) )  
  
  
  
  # Regressions
  
  reg_logit <- regression( Y_train = (1 + BBTSYNT$signal$insample[ ,"states"]) / 2,
                           X_train = ema( na.omit(BBTSYNT$signal$insample[ ,c("delta", "bull")]), 1),
                           type = "logit" )
  summary( reg_logit$reg )
  
  plot( x = reg_logit$z, y = reg_logit$y )
  
  plot( reg_logit$y )
  
 
  
  DS <- dirichletSampling( Y_train = Z_bm,
                           X_train = BBTSYNT$signal$insample[ ,"bull"],
                           weights_fun = "kernel",
                           correct_bias = FALSE,
                           sclfct = NULL )
  ldens <- lapply( DS, FUN = densFUN )
  plot.ldensity( ldens )
  
  
  
  KCD <- kcdens( Y_train = Z_bm,
                 X_train = BBTSYNT$signal$insample[ ,"delta"],
                 kcd_fun = "hyndman" )
  plot( KCD$fit )
  
  
  
  

  
  
    
  
  
  