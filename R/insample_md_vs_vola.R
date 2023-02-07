  
  
  ############################################################################
  ### GOOD AND BAD TURBULENCE - IN-SAMPLE - MAHALANOBIS DISTANCE VS. VOLATILITY
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     22.10.2020
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
  # DATA
  # --------------------------------------------------------------------------
  
  TurbSensor <- loadSensor( "turbulence" )
  TurbSensorScl <- loadSensor( "turbulence_scl2" )
  VRPSensor <- loadSensor( "vrp" )
  BBSSensor <- loadSensor( "bbs" )
  BBTSVSensor <- loadSensor( "bbtsv_base_dm" )
  BBTSensor <- loadSensor( "bbturbulence_base_dm" )
  
  
  
  # Index
  X_bm <- TurbSensor$data$X_bm
  
  # Implied volatility (VIX)
  vix <- VRPSensor$signal$base[ ,"iv"]
  
  # Realized volatility
  rv <- VRPSensor$signal$base[ ,"rv"]

  # # Garch
  # cvol <- getCondVar( garch(X_bm) )
  
  # BBT
  bbt <- BBTSensor$signal$base[ ,"md_delta"]
  bbt_ema <- BBTSensor$signal$base[ ,"signal"]
  bbtsv <- BBTSVSensor$signal$base[ ,"signal"]
  

  
  signals <- cbind( md = TurbSensor$signal$base[ ,"md"],
                    md_scl = TurbSensorScl$signal$scl2[ ,"md"],
                    vix = vix,
                    rv = rv,
                    bbtsv = bbtsv,
                    bbt = bbt,
                    bbt_ema)
  # garch = cvol )
  colnames(signals) <- c("md", "md_scl", "vix", "rv", "bbtsv", "bbt", "bbt_ema")
  
  
  plot( tail( na.omit(signals[ ,c("bbtsv", "bbt", "bbt_ema")]), 100 ), plot.type = "single", type = "o" )
  
  
  
  # --------------------------------------------------------------------------
  # DESCRIPTIVE ANALYSIS
  # --------------------------------------------------------------------------
  
  TD <- trainingData( Y_train = X_bm,
                     X_train = signals )
  
  
  lDS <- list()
  for ( j in 1:ncol(TD$X_train) ) {
    
    lDS[[j]] <- dirichletSampling( Y_train = TD$Y_train,
                                   X_train = TD$X_train[ ,j],
                                   n_lag = 0,
                                   weights_fun = "l1",
                                   correct_bias = FALSE,
                                   sclfct = NULL )
                                   # scl_by_entropy = TRUE )
    
  }
  names(lDS) <- colnames(TD$X_train)
  

  from <- min( unlist(lDS) ) * 1.2
  to <- max( unlist(lDS) ) * 1.2
  n <- 999
  lldens <- lapply( lDS, FUN = function(x) { 
              lapply(x, density, from = from, to = to, n = n) } )
  # par(mfrow = c(length(lldens), 1))
  for ( i in 1:length(lldens) ) {
    plot.ldensity( lldens[[i]], 
                   main = paste0("Posterior means conditional on ", 
                                 names(lDS)[[i]]) )
  }
  
  
  
  
  # --------------------------------------------------------------------------
  # PREDICTIVE ANALYSIS
  # --------------------------------------------------------------------------
  
  BBSSensor$computeSignalInsample()
  states <- (1 + BBSSensor$signal$insample) / 2
  TD <- trainingData( Y_train = states,
                      X_train = signals )
  n_lag <- 5
  lDS_pred <- list()
  for ( j in 1:ncol(TD$X_train) ) {
    
    lDS_pred[[j]] <- dirichletSampling( Y_train = TD$Y_train,
                                        X_train = TD$X_train[ ,j],
                                        n_lag = n_lag,
                                        weights_fun = "l1",
                                        # weights_fun = "qnn",
                                        # alpha = 0.1,
                                        correct_bias = FALSE,
                                        sclfct = NULL ) 
                                        # scl_by_entropy = TRUE )
    
  }
  names(lDS_pred) <- colnames(TD$X_train)
  
  
  # from <- min( unlist(lDS_pred) ) * 1.3
  # to <- max( unlist(lDS_pred) ) * 1.3
  # n <- 999
  # lldens <- lapply( lDS_pred, FUN = function(x) { 
  #             lapply(x, density, from = from, to = to, n = n) } )
  lldens <- lapply( lDS_pred, FUN = function(x) { 
                    lapply(x, density) } )
  for ( i in 1:length(lldens) ) {
    plot.ldensity( lldens[[i]], 
                   fillin = FALSE,
                   main = paste0("Posterior means conditional on ", 
                                 names(lDS_pred)[[i]]) )
  }
  
  
  
  
  # --------------------------------------------------------------------------
  # PREDICTIVE ANALYSIS - SCATTERPLOT
  # --------------------------------------------------------------------------  
  
  TD <- trainingData( Y_train = log(1 + X_bm),
                      X_train = signals,
                      n_lag = 1,
                      apply_roll_y = list( fun = sum, width = 10, by = 1 ) )
  j <- 2 
  plot( x = as.numeric(TD$X_train[ ,j]), y = as.numeric(TD$Y_train), main = colnames(TD$X_train)[j] ) 
  
  
  # debugonce( pas )
  n_lags <- c(0, 5, 10, 21)
  PAS <- pas( X = log(1 + X_bm),
              sig = signals,
              n_lags = n_lags )
              # sigAggFUN = parse(text = "mean") )
 
  
  j <- "bbt"
  par( mfrow = c(length(n_lags), 1) )
  corvec <- numeric(length(n_lags))
  for ( k in seq(along = n_lags) ) {
    plot( x = as.numeric(PAS[[k]]$X_train[ ,j]),
          y = as.numeric(PAS[[k]]$Y_train), 
          main = colnames(PAS[[k]]$X_train)[j] )
    abline(h = 0)
    corvec[k] <- cor( as.numeric(PAS[[k]]$X_train[ ,j]), as.numeric(PAS[[k]]$Y_train) )
  }
  corvec
  
  
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # PREDICTIVE ANALYSIS - LOGISTIC REGRESSION
  # --------------------------------------------------------------------------  
  
  
  n_lags <- c(0, 5, 10, 21)
  states <- (1 + BBSSensor$signal$insample) / 2
  j <- "vix"
  lReg <- list()
  for ( k in seq(along = n_lags) ) {
    lReg[[k]] <- regression( Y_train = states,  
                             X_train = signals,
                             n_lag = n_lags[k],
                             type = "logit" )
  }
  
  lapply( lReg, FUN = function(x) { summary(x$reg) } )
  
  
  dev.off()  
  plot( x = lReg[[1]]$z, y = lReg[[1]]$y )
  
  
  
  
  
  
  
  
  
  
