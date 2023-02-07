  
  ############################################################################
  ### BBTurbulence - IN-SAMPLE FACTOR ANALYSIS
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     10.10.2020
  # First version:    10.10.2020
  # --------------------------------------------------------------------------
  
  
  # --------------------------------------------------------------------------
  # REQUIRE
  # --------------------------------------------------------------------------
  
  require(RP)
  require(DAARC)
  
  wd <- "R:/Asset_Management/Research_Projects/Multi_Asset/DAA/DAARC/Sensor_Testing/Univariate/Good_And_Bad_Turbulence/"
  source( paste0(wd, "Source/source_functions.R") )
  
  
  BBSObj <- loadSensor( sensor_name = "bbs", b_new = TRUE )
  BBTObj <- loadSensor( sensor_name = "bbturbulence", b_new = TRUE )  
  X <- BBTObj$data$X
  X_bm <- BBTObj$data$X_bm
  states_is <- BBSObj$signal$insample
  md_is <- BBTObj$signal$insample
  
  
  minphase <- 5
  mincycle <- 21
  BBS <- bbs(X = cumulated(X_bm, "discrete"), 
             mincycle = mincycle, 
             minphase = minphase, 
             k.peak = minphase, 
             l.peak = minphase, 
             k.trough = minphase, 
             l.trough = minphase, 
             logarithmic = FALSE, 
             theta = Inf, 
             e = 0)
  
  plot( cumulated(X_bm, "discrete") )
  abline( v = time(BBS)[which(BBS == -1)], col = "red" )
  abline( v = time(BBS)[which(BBS == 1)], col = "green" )
  lines( cumulated(X_bm, "discrete") )
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Factor strategies
  # --------------------------------------------------------------------------
  
  X_fact <- factorIndices( universe = "europe" )
  # idx <- which(colnames(X_fact) == "capw")
  # X_delta <- (X_fact[ ,-idx] - X_fact[ ,rep(idx, ncol(X_fact) - 1)]) / 
  #   (1 + X_fact[ ,rep(idx, ncol(X_fact) - 1)])
  # 
  # plot( as.simTS(X_delta) )
  # descStats(X_delta)
  X_bm <- na.omit(X_fact[ ,"capw"])
  X_fact <- na.omit(X_fact[ ,-which(colnames(X_fact) == "capw")])
  X_bm <- X_bm[isWeekday(time(X_bm)), ]
  X_fact <- X_fact[isWeekday(time(X_fact)), ]
  
  
  BBTFact <- BBTurbulence$new()
  # BBTFact <- BBTGARCH$new()
  # BBTFact <- BBTEWMA$new()
  # BBTFact <- BBTSV$new()
  BBTFact$setCtrl( method = "base", universe = "dm" )
  BBTFact$spec$iso <- colnames(X_fact)
  BBTFact$data <- list(X = X_fact,
                       X_bm = X_bm)
  # debugonce( BBTFact$computeSignalInsample )
  BBTFact$computeSignalInsample()
  
  plot( na.omit(BBTFact$signal$insample) )
  
  
  # Signal testing
  
  sig <- (1 + sign(ema(na.omit(BBTFact$signal$insample[ ,"delta"]), 0.1))) / 2
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
  
  plot( as.simTS(X_tmp), col = colors, main = "Insample - After costs" )
  t( descStats( X_tmp )$stats[stats_fields, ] )
  drawDownStats(X_tmp)
  
  plot( as.simTS(X_tmp_nc), col = colors, main = "Insample - Before costs" )
  t( descStats( X_tmp_nc )$stats[stats_fields, ] )
  drawDownStats(X_tmp_nc)
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Sector indices
  # --------------------------------------------------------------------------
  
  env <- readRDS( file = "R:/Asset_Management/Research_Projects/Assignments/Sector_Performances/sector_performances_msci_classification.rds" )
  X_sect <- rodbcGetOLZDBReturns( assetName = c("NDDUWI", env$ticker_cyclicals, env$ticker_defensive),
                             frqncy = "daily",
                             refCcy = "USD" )
  X_bm <- X_sect[ ,"NDDUWI"]
  X_sect <- na.omit(X_sect[ ,-which(colnames(X_sect) == "NDDUWI")])
  X_bm <- X_bm[isWeekday(time(X_bm)), ]
  X_sect <- X_sect[isWeekday(time(X_sect)), ]
  
  
  BBTFact <- BBTurbulence$new()
  # BBTFact <- BBTGARCH$new()
  # BBTFact <- BBTEWMA$new()
  # BBTFact <- BBTSV$new()
  BBTFact$setCtrl( method = "base", universe = "dm" )
  BBTFact$spec$iso <- colnames(X_sect)
  BBTFact$data <- list(X = X_sect,
                       X_bm = X_bm)
  # debugonce( BBTFact$computeSignalInsample )
  BBTFact$computeSignalInsample()
  
  plot( BBTFact$signal$insample )
  
  
  # Signal testing
  
  alpha <- 0.1
  sig <- (1 + sign(ema(na.omit(BBTFact$signal$insample[ ,"delta"]), alpha))) / 2
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
  
  plot( as.simTS(X_tmp), col = colors, main = "Insample - After costs" )
  t( descStats( X_tmp )$stats[stats_fields, ] )
  drawDownStats(X_tmp)
  
  plot( as.simTS(X_tmp_nc), col = colors, main = "Insample - Before costs" )
  t( descStats( X_tmp_nc )$stats[stats_fields, ] )
  drawDownStats(X_tmp_nc)
  
  
  
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Daily
  # --------------------------------------------------------------------------
  
  signal <- BBTObj$signal$base[ ,"md_delta"]
  
  
  
  lDS <- list()
  for ( j in 1:ncol(X_delta) ) {
    
    lDS[[j]] <- dirichletSampling( Y_train = Y_train_w[ ,j],
                                   X_train = X_train_w,
                                   sclfct = NULL,
                                   # sclfct = 1,
                                   correct_bias = FALSE )
  }
  
  from <- quantile( unlist(lDS), 0.001 )
  to <- quantile( unlist(lDS), 0.999 )
  n <- 500
  lldens <- lapply( lDS, FUN = function(x) { lapply(x, density, from = from, to = to, n = n) } )
  
  
  par( mfrow = c(2, 2) )
  plot.ldensity( lldens[[1]], fillin = FALSE, main = colnames(X_delta)[1] )
  plot.ldensity( lldens[[2]], fillin = FALSE, main = colnames(X_delta)[2] )
  plot.ldensity( lldens[[3]], fillin = FALSE, main = colnames(X_delta)[3] )
  plot.ldensity( lldens[[4]], fillin = FALSE, main = colnames(X_delta)[4] )
  
  dev.off()
  
  
  
  # --------------------------------------------------------------------------
  # Weekly
  # --------------------------------------------------------------------------
  
  X_fact_w <- aggWeekly(X_fact, day = "Tue", compounding = "discrete")
  Y_train_w <- aggWeekly(X_delta, day = "Tue", compounding = "discrete")
  # idx <- which(colnames(X_fact) == "capw")
  # Y_train_w2 <- (X_fact_w[ ,-idx] - X_fact_w[ ,rep(idx, ncol(X_fact_w) - 1)]) / 
  #   (1 + X_fact_w[ ,rep(idx, ncol(X_fact_w) - 1)])
  # range(Y_train_w - Y_train_w2)  # same same
  
  signal <- BBTObj$signal$base[ ,"md_delta"]
  X_train_w <- filter( signal, rep(1/5, 5), side = 1 )
  
  lDS <- list()
  for ( j in 1:ncol(Y_train_w) ) {
    
    lDS[[j]] <- dirichletSampling( Y_train = Y_train_w[ ,j],
                                   X_train = X_train_w,
                                   # sclfct = NULL,
                                   sclfct = 1,
                                   correct_bias = FALSE )
  }
  
  from <- quantile( unlist(lDS), 0.001 )
  to <- quantile( unlist(lDS), 0.999 )
  n <- 500
  lldens <- lapply( lDS, FUN = function(x) { lapply(x, density, from = from, to = to, n = n) } )
  
  
  par( mfrow = c(2, 2) )
  plot.ldensity( lldens[[1]], fillin = FALSE, main = colnames(X_delta)[1] )
  plot.ldensity( lldens[[2]], fillin = FALSE, main = colnames(X_delta)[2] )
  plot.ldensity( lldens[[3]], fillin = FALSE, main = colnames(X_delta)[3] )
  plot.ldensity( lldens[[4]], fillin = FALSE, main = colnames(X_delta)[4] )
  
  dev.off()
  
  
  # --------------------------------------------------------------------------
  # Monthly
  # --------------------------------------------------------------------------
  
  Y_train_m <- aggMonthly(X_delta, compounding = "discrete")
  signal <- BBTObj$signal$base[ ,"md_delta"]
  X_train_m <- na.omit( filter( signal, rep(1/21, 21), side = 1 ) )
  
  lDS <- list()
  for ( j in 1:ncol(Y_train_m) ) {
    # debugonce(dirichletSampling)
    lDS[[j]] <- dirichletSampling( Y_train = Y_train_m[ ,j],
                                   X_train = X_train_m,
                                   sclfct = NULL,
                                   # sclfct = 1,
                                   n_lag = 1,
                                   correct_bias = FALSE )
  }
  
  from <- quantile( unlist(lDS), 0.001 )
  to <- quantile( unlist(lDS), 0.999 )
  n <- 500
  lldens <- lapply( lDS, FUN = function(x) { lapply(x, density, from = from, to = to, n = n) } )
  
  
  par( mfrow = c(2, 2) )
  plot.ldensity( lldens[[1]], fillin = FALSE, main = colnames(X_delta)[1] )
  plot.ldensity( lldens[[2]], fillin = FALSE, main = colnames(X_delta)[2] )
  plot.ldensity( lldens[[3]], fillin = FALSE, main = colnames(X_delta)[3] )
  plot.ldensity( lldens[[4]], fillin = FALSE, main = colnames(X_delta)[4] )
  
  dev.off()
  
  
  
  