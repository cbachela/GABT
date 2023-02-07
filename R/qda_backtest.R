  
  
  ############################################################################
  ### BBTurbulence - QDA
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     29.10.2020
  # First version:    10.10.2020
  # --------------------------------------------------------------------------
  
  # Idea: 
  # Compute relative turbulence as the difference of (scaled) good and bad
  # turbulence. Smooth the difference by augmenting current difference 
  # values with lagged ones and computing the Mahalanobis distance.
  # Obtain a distribution of state probabilities by sampling training lables

  
  
  # --------------------------------------------------------------------------
  # REQUIRE
  # --------------------------------------------------------------------------
  
  require(MASS)
  require(RP)
  require(DAARC)
  
  # wd <- "H:/Papers Cyril/PhD_Papers/Good_And_Bad_Turbulence/R/"
  wd <- "R:/Asset_Management/Research_Projects/Multi_Asset/DAA/DAARC/Sensor_Testing/Univariate/Good_And_Bad_Turbulence/"
  source( paste0(wd, "Source/source_functions.R") )
  
  
  BBSObj <- loadSensor( sensor_name = "bbs", b_new = TRUE )
  BBTObj <- loadSensor( sensor_name = "bbturbulence", b_new = TRUE )  
  X <- BBTObj$data$X[isWeekday(time(BBTObj$data$X)), ]
  X_bm <- BBTObj$data$X_bm[isWeekday(time(BBTObj$data$X_bm)), ]

  BBTObj$computeSignalInsample()
  BBS <- BBTObj$signal$insample[ ,"states"]
  
  plot( cumulated(X_bm, "discrete") )
  abline( v = time(BBS)[which(BBS == -1)], col = "red" )
  abline( v = time(BBS)[which(BBS == 1)], col = "green" )
  lines( cumulated(X_bm, "discrete") )
  
  
  
  
  # --------------------------------------------------------------------------
  # Augment indicator with lagged values
  # --------------------------------------------------------------------------
  
  BBT <- BBTObj$copy()
  md_delta <- BBT$signal$scl2[ ,"md_delta"]
  
  plot(md_delta)  
    
  tmp <- lapply( 0:10, function(n) { lag(md_delta, n) } )
  sig <- na.omit( do.call( cbind, tmp ) )
  md_augm <- mahalanobisTS( X = sig )
  
  plot(md_augm)

  
  tau_vec <- round( seq(from = 0.01, to = 1, length.out = 50)^2, 5 )
  FUN <- function(x) { ema( md_delta, alpha = x ) }
  lSignal_ema <- lapply( tau_vec, FUN = FUN )
  sigs1_ema <- do.call( cbind, lSignal_ema )
  
  plot(sigs1_ema[ ,1:10])
  
  head( tail(sigs1_ema[ ,1:10], 180), 20 )
  
  
  
  ###############
  
  minphase_w <- 1 * 4
  mincycle_w <- 5 * 4
  BBTW <- BBTObj$copy()
  BBTW$setCtrl()
  BBTW$spec$minphase <- minphase_w
  BBTW$spec$mincycle <- mincycle_w
  BBTW$spec$k_peak <- BBTW$spec$k_trough <- BBTW$spec$l_peak <- 
      BBTW$spec$k_trough <- minphase_w
  BBTW$data <- list()
  BBTW$signal <- list()
  FUN <- function(X) { exp( apply(log(1 + X), 2, mean) ) - 1 }
  Width <- 21
  X_roll <- applyRoll( BBTObj$data$X,
                       Width = Width,
                       By = 1,
                       FUN = FUN )
  X_bm_roll <- applyRoll( BBTObj$data$X_bm,
                          Width = Width,
                          By = 1,
                          FUN = FUN )
  wmat_roll <- applyRoll( BBTObj$data$wmat, 
                          Width = Width, 
                          FUN = function(X) { apply(X, 2, mean) },
                          By = 1 )
  
  dates <- rownames(BBTObj$data$X)[ -c(1:(252*3)) ]
  ans <- BBTObj$signal$insample[dates, ] * NA
  
  for ( today in dates ) {
    
    dates_tmp <- rownames(X_roll)[ rownames(X_roll) <= today ]
    idx <- rev( dates_tmp[ seq(from = length(dates_tmp), to = 1, by = -Width) ] )
    
    Obj <- BBTW$copy()
    Obj$data$X <- X_roll[idx, ]
    Obj$data$X_bm <- X_bm_roll[idx, ]
    Obj$data$wmat <- wmat_roll[idx, ]
    # debugonce( Obj$computeSignalInsample )
    Obj$computeSignalInsample()
    ans[today, ] <- tail(Obj$signal$insample, 1)
    
  }  
  
  
  
  plot( ans[ ,"delta"])
  plot( ans[ ,"states"])
  plot( ans[ ,"all"])
  
  insample_roll <- applyRoll( BBTObj$signal$insample,
                              Width = 21,
                              By = 1,
                              FUN = function(X) { apply(X, 2, mean) } )
  
  Name <- "delta"
  tmp <- na.omit( cbind( insample = BBTObj$signal$insample[ ,Name], 
                         backtest = ans[ ,Name],
                         insample_roll = insample_roll[ ,Name] ) )
  plot( tmp, plot.type = "single" )
  
  
  
  
  # Signal testing
  
  sig <- tmp * 0
  sig[ tmp > 0 ] <- 1
  test <- signalTesting.byTrading( X = X_bm,
                                   sig = sig,
                                   n_lag = 1,
                                   tc = 0.00)
  X_tmp <- na.omit( cbind(X_bm, test) )
  
  plot( cumulated(X_tmp, "discrete"), plot.type = "single" )
  legend("topleft", colnames(X_tmp), lwd = 2, col = 1:ncol(X_tmp), text.col = 1:ncol(X_tmp), bty = "n")
  
  
  
  
  
  
  
  
  
  
  
  
