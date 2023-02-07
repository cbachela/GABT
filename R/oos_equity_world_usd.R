  
  
  ############################################################################
  ### GOOD AND BAD TURBULENCE - DM - USD
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     20.02.2022
  # First version:    20.02.2022
  # --------------------------------------------------------------------------
  
  
  # --------------------------------------------------------------------------
  # REQUIRE
  # --------------------------------------------------------------------------
  
  require(RP)
  require(BSS)
  require(simolz)
  require(DAARC)
  
  wd <- "R:/Asset_Management/Research_Projects/Multi_Asset/DAA/DAARC/Sensor_Testing/Univariate/Good_And_Bad_Turbulence/"
  source( paste0(wd, "Source/source_functions.R") )

  
  
  
 
  # --------------------------------------------------------------------------
  # LOCAL CCY
  # --------------------------------------------------------------------------
  
  bbtewma_base <- loadSensor( sensor_name = "bbtewma_base_dm" )
  bbtewma_scl <- loadSensor( sensor_name = "bbtewma_scl_dm" )
  bbtewma_scl2 <- loadSensor( sensor_name = "bbtewma_scl2_dm" )
 
  # debugonce( bbtewma_base$computeSignal )
  bbtewma_base$computeSignalInsample()
  bbtewma_scl$computeSignalInsample()
  bbtewma_scl2$computeSignalInsample()
  bbtewma_base$update()
  bbtewma_scl$update()
  bbtewma_scl2$update()
  
  
  # BBTEWMA base with demean
  Obj <- DAARC::BBTEWMA$new()
  Obj$setCtrl( universe = "dm",
               method = "base" )
  Obj$spec$name <- "bbtewma_base_dm_qis"
  Obj$spec$demean <- TRUE
  Obj$data <- list( X = bbtewma_base$data$X,
                    BBS = bbtewma_base$data$BBS,
                    X_bm = bbtewma_base$data$X_bm,
                    capw = bbtewma_base$data$capw )
  Obj$computeSignalInsample()
  Obj$updateSignal()
  
  
  sig_base <- bbtewma_base$signal$insample[ ,"delta"]
  sig_base_demean <- Obj$signal$insample[ ,"delta"]
  sig_base <- bbtewma_base$getSignal()
  sig_base_demean <- Obj$getSignal()
  sig <- cbind( base = sig_base,
                base_demean = sig_base_demean )
  
  
  test <- signalTesting.byTrading( X = bbtewma_base$data$X_bm,
                                   sig = (1 + sign(ema(sig, 1))) / 2,
                                   n_lag = 1, 
                                   tc = 0 )
  X_tmp <- na.omit(cbind(bbtewma_base$data$X_bm, test))
  plot( as.simTS(X_tmp) )
  
  stats_fields <- c("cumret", "sds", "sharpe", "maxDD")
  t( descStats( X_tmp )$stats[stats_fields, ] )
  
  
  
  sig_base <- bbtewma_base$getSignal()
  sig_scl <- bbtewma_scl$getSignal()
  sig_scl2 <- bbtewma_scl2$getSignal()
  sig <- cbind( base = sig_base,
                scl = sig_scl,
                scl2 = sig_scl2 )
  sig <- cbind( sig, avg = apply( sig, 1, mean) )
  
  plot( tail(sig, 100), plot.type = "single", type = "o" )
  # plot( tail(sig, 10000), plot.type = "single", type = "l" )
  abline(h = 0)
  tail( sig, 100 )
  tail( bbtewma_base$data$X_bm )
  bbtewma_base$spec$ccy
  
  
  plot( bbtewma_base$signal$base )
  plot( bbtewma_scl2$signal$scl2 )
  
  
  
  sig_tmp <- cbind( sig )
  # sig_tmp <- window(sig_tmp, "2015-01-01", end(sig_tmp))
  test <- signalTesting.byTrading( X = bbtewma_base$data$X_bm,
                                   sig = (1 + sign(ema(sig_tmp, 1))) / 2,
                                   n_lag = 1, 
                                   tc = 0 )
  X_tmp <- na.omit(cbind(bbtewma_base$data$X_bm, test))
  plot( as.simTS(X_tmp) )
  
  stats_fields <- c("cumret", "sds", "sharpe", "maxDD")
  t( descStats( X_tmp )$stats[stats_fields, ] )
  
  
  
  a <- (1 + sign(ema(sig_tmp[ ,3], 1))) / 2
  sum(abs(na.omit(diff(a))))
  nrow(a)
  
  
  
  
  X <- bbtewma_base$data$X
  tail(X)
  
  plot( as.simTS(tail(X, 100)) )
  
  barplot( meanGeo(tail(X, 100)) )
  
  
  
  
  
  # --------------------------------------------------------------------------
  # USD
  # --------------------------------------------------------------------------
  
  # bbtewma_base <- loadSensor( sensor_name = "bbtewma_base_dm" )
  # bbtewma_scl <- loadSensor( sensor_name = "bbtewma_scl_dm" )
  # bbtewma_scl2 <- loadSensor( sensor_name = "bbtewma_scl2_dm" )
  # bbtewma_base$computeSignalInsample()
  # bbtewma_scl$computeSignalInsample()
  # bbtewma_scl2$computeSignalInsample()
  # bbtewma_base_usd <- bbtewma_base$copy()
  # bbtewma_scl_usd <- bbtewma_scl$copy()
  # bbtewma_scl2_usd <- bbtewma_scl2$copy()
  # bbtewma_base_usd$spec$ccy <- "USD"
  # bbtewma_scl_usd$spec$ccy <- "USD"
  # bbtewma_scl2_usd$spec$ccy <- "USD"
  # bbtewma_base_usd$spec$name <- paste0(bbtewma_base$spec$name, "_usd")
  # bbtewma_scl_usd$spec$name <- paste0(bbtewma_scl$spec$name, "_usd")
  # bbtewma_scl2_usd$spec$name <- paste0(bbtewma_scl2$spec$name, "_usd")
  # bbtewma_base_usd$signal <- list()
  # bbtewma_scl_usd$signal <- list()
  # bbtewma_scl2_usd$signal <- list()
  # bbtewma_base_usd$data <- list( BBS = bbtewma_base$data$BBS )
  # bbtewma_scl_usd$data <- list( BBS = bbtewma_base$data$BBS )
  # bbtewma_scl2_usd$data <- list( BBS = bbtewma_base$data$BBS )
  # # debugonce(bbtewma_scl_usd$computeData)
  # bbtewma_base_usd$updateData()
  # bbtewma_scl_usd$updateData()
  # bbtewma_scl2_usd$updateData()
  
  
  bbtewma_base_usd <- loadSensor( sensor_name = "bbtewma_base_dm_usd" )
  bbtewma_scl_usd <- loadSensor( sensor_name = "bbtewma_scl_dm_usd" )
  bbtewma_scl2_usd <- loadSensor( sensor_name = "bbtewma_scl2_dm_usd" )
  
  # debugonce( bbtewma_base_usd$computeSignalInsample )
  bbtewma_base_usd$computeSignalInsample()
  bbtewma_scl_usd$computeSignalInsample()
  bbtewma_scl2_usd$computeSignalInsample()
  bbtewma_base_usd$update()
  bbtewma_scl_usd$update()
  bbtewma_scl2_usd$update()
  
  
  
  sig_base <- bbtewma_base$getSignal()
  sig_scl <- bbtewma_scl$getSignal()
  sig_scl2 <- bbtewma_scl2$getSignal()
  sig <- cbind( base = sig_base,
                scl = sig_scl,
                scl2 = sig_scl2 )
  sig <- cbind( sig, avg = apply( sig, 1, mean) )
  
  sig_base_usd <- bbtewma_base_usd$getSignal()
  sig_scl_usd <- bbtewma_scl_usd$getSignal()
  sig_scl2_usd <- bbtewma_scl2_usd$getSignal()
  sig_usd <- cbind( base = sig_base_usd,
                scl = sig_scl_usd,
                scl2 = sig_scl2_usd )
  sig_usd <- cbind( sig_usd, avg = apply( sig_usd, 1, mean) )
  
  
  plot( tail(sig, 100), plot.type = "single", type = "o" )
  abline(h = 0)
  plot( tail(sig_usd, 100), plot.type = "single", type = "o" )
  abline(h = 0)
  tail( cbind(sig, sig_usd), 10 )

  
  tmp <- cbind(sig_base, sig_base_usd)
  tail(tmp)
  
  tmp <- cbind(bbtewma_base$signal$insample[ ,"delta"],
               bbtewma_base_usd$signal$insample[ ,"delta"])
  tmp
  
  tmp <- cbind( bbtewma_base$data$X_bm, 
                bbtewma_base_usd$data$X_bm)
  tmp
  plot( tail(tmp, 100), plot.type = "single" )
  abline(h = 0)
  
  headleft(bbtewma_base$data$X)
  headleft(bbtewma_base_usd$data$X)
  
  
  
  stats_fields <- c("cumret", "means", "sds", "sharpe", "maxDD")
  
  sig_tmp <- (1 + sign(ema(sig, 1))) / 2
  # sig_tmp <- window(sig_tmp, "2015-01-01", end(sig_tmp))
  # debugonce( signalTesting.byTrading )
  test <- signalTesting.byTrading( X = bbtewma_base$data$X_bm,
                                   sig = sig_tmp,
                                   n_lag = 1, 
                                   tc = 0 )
  X_tmp <- na.omit(cbind(bbtewma_base$data$X_bm, test))
  plot( as.simTS(X_tmp) )
  t( descStats( X_tmp )$stats[stats_fields, ] )
  
  
  sig_tmp <- (1 + sign(ema(sig_usd, 1))) / 2
  # sig_tmp <- window(sig_tmp, "2015-01-01", end(sig_tmp))
  # debugonce( signalTesting.byTrading )
  test <- signalTesting.byTrading( X = bbtewma_base_usd$data$X_bm,
                                   sig = sig_tmp,
                                   n_lag = 1, 
                                   tc = 0.001,
                                   penalty = 0.5 )
  X_tmp_usd <- na.omit(cbind(bbtewma_base_usd$data$X_bm, test))
  plot( as.simTS(X_tmp_usd) )
  t( descStats( X_tmp_usd )$stats[stats_fields, ] )
  
  
  
  # k-means
  KM <- DAA::kmeans.adj( x = sig_usd[ ,2], centers = 2 )
  sig_tmp <- KM$cluster / max(KM$cluster)
  test <- signalTesting.byTrading( X = bbtewma_base_usd$data$X_bm,
                                   sig = sig_tmp,
                                   n_lag = 1, 
                                   tc = 0,
                                   penalty = 0 )
  X_tmp_usd <- na.omit(cbind(bbtewma_base_usd$data$X_bm, test))
  plot( as.simTS(X_tmp_usd) )
  t( descStats( X_tmp_usd )$stats[stats_fields, ] )
  
  
  
  # quintile
  tmp <- sig_usd[ ,2]
  # quant <- quantile( tmp, c(0.4, 0.9) )
  quant <- c(-0.01, quantile( tmp, 0.99) )
  quant
  sig_tmp <- tmp * 0
  sig_tmp[ tmp > quant[1], ] <- 2
  sig_tmp[ tmp > quant[2], ] <- 5
  test <- signalTesting.byTrading( X = bbtewma_base_usd$data$X_bm,
                                   sig = sig_tmp,
                                   n_lag = 1, 
                                   tc = 0,
                                   penalty = -1 )
  X_tmp_usd <- na.omit(cbind(bbtewma_base_usd$data$X_bm, test))
  plot( as.simTS(X_tmp_usd) )
  t( descStats( X_tmp_usd )$stats[stats_fields, ] )
  
  
  plot(density(tmp))
  abline(v = quant)
  
  plot( sig_tmp )

    
  

