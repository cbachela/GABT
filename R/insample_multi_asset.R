
  
  ############################################################################
  ### GOOD AND BAD TURBULENCE - IN-SAMPLE ANALYSIS - MULTI ASSET DATASET
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     17.10.2020
  # First version:    17.10.2020
  # --------------------------------------------------------------------------



  require(RP)
  require(DAARC)
  
  wd <- "R:/Asset_Management/Research_Projects/Multi_Asset/DAA/DAARC/Sensor_Testing/Univariate/Good_And_Bad_Turbulence/"
  source( paste0(wd, "Source/source_functions.R") )
  
  
  require(SID)
  SID <- sid( b_update = FALSE )
  X <- SID$X_est


  
  
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
  # Construct absolute and relative Mahalanobis distances
  # --------------------------------------------------------------------------
  
  BBT <- loadSensor(sensor_name = "bbturbulence")
  BBT$signal <- list()
  BBT$spec$theta <- theta
  
  
  # Using daily returns
  X_d <- SID$X_sim[ ,colnames(SID$X_est)]
  BBTD <- BBT$copy()
  BBTD$data$X <- X_d
  BBTD$data$X_bm <- X_d[ ,"Equity_DM_hedged"]
  BBTD$data$wmat <- NULL
  BBTD$spec$minphase <- minphase_d
  BBTD$spec$mincycle <- mincycle_d
  BBTD$spec$l_peak <- BBTD$spec$l_trough <- 
    BBTD$spec$k_peak <- BBTD$spec$k_trough <- minphase_d
  # debugonce( BBTD$computeSignalInsample )
  BBTD$computeSignalInsample()
  mdd <- BBTD$signal$insample
  mdd <- mdd[isWeekday(time(mdd)), ] 
  mdd_scl <- apply(mdd, 2, function(x) { x / max(x) } )
  
  plot(mdd)  
  
  
  # Using weekly returns
  X_w <- SID$X_est
  BBTW <- BBT$copy()
  BBTW$data$X <- X_w
  BBTW$data$X_bm <- X_w[ ,"Equity_DM_hedged"]
  BBTW$data$wmat <- NULL
  BBTW$spec$minphase <- minphase_w
  BBTW$spec$mincycle <- mincycle_w
  BBTW$spec$l_peak <- BBTW$spec$l_trough <- 
    BBTW$spec$k_peak <- BBTW$spec$k_trough <- minphase_w
  # debugonce( BBTW$computeSignalInsample )
  BBTW$computeSignalInsample()
  mdw <- BBTW$signal$insample
  mdw <- mdw[isWeekday(time(mdw)), ] 
  mdw_scl <- apply(mdw, 2, function(x) { x / max(x) } )

  plot(mdw)  
  
  
  Names <- c("bear_scl", "bull_scl", "all", "delta")
  plot( mdw[ ,Names], plot.type = "single" )
  
  plot( ema(mdw[ ,"delta"], 0.1) )
  abline( h = 0 )
  
  plot( ema(mdd[ ,"delta"], 0.1) )
  abline( h = 0 )
  
  
  plot( ema(mdw[ ,c("bear_scl", "bull_scl")], 0.1), plot.type = "single" )

  
  
  
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
  
  
  
  
  
  
 
  
  
  
  
  
  
  
  
  