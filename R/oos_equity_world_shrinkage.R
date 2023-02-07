  
  
  ############################################################################
  ### GOOD AND BAD TURBULENCE - DM - NONLINEAR SHRINKAGE
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

  
  
  
  VRPObj <- loadSensor( sensor_name = "vrp" )
  VRPObj$update()
  tail( VRPObj$signal$base, 10 )
  
  
  bbt_base <- loadSensor( sensor_name = "bbturbulence_base_dm" )
  bbt_scl <- loadSensor( sensor_name = "bbturbulence_scl_dm" )
  bbt_scl2 <- loadSensor( sensor_name = "bbturbulence_scl2_dm" )
  
  # debugonce( bbtewma_base$computeSignal )
  bbt_base$update()
  bbt_scl$update()
  bbt_scl2$update()
  
  sig_base <- bbt_base$signal$base[ ,"md_delta"]
  sig_scl2 <- bbt_scl2$signal$scl2[ ,"md_delta"]
  sig <- cbind( base = sig_base,
                scl2 = sig_scl2 )
  sig <- cbind( sig, avg = apply( sig, 1, mean) )
  
  plot( tail(sig, 100), plot.type = "single", type = "o" )
  abline(h = 0)
  
 
  # --------------------------------------------------------------------------
  # CAPITALIZATION WEIGHTED BENCHMARK
  # --------------------------------------------------------------------------
  
  bbtewma_base <- loadSensor( sensor_name = "bbtewma_base_dm" )
  bbtewma_scl <- loadSensor( sensor_name = "bbtewma_scl_dm" )
  bbtewma_scl2 <- loadSensor( sensor_name = "bbtewma_scl2_dm" )
 
  # debugonce( bbtewma_base$computeSignal )
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
  
  
  # In-sample
  
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
  
  
  
  # Out-of-sample
  
  sig_base <- bbtewma_base$getSignal()
  sig_scl <- bbtewma_scl$getSignal()
  sig_scl2 <- bbtewma_scl2$getSignal()
  sig <- cbind( base = sig_base,
                # scl = sig_scl,
                scl2 = sig_scl2 )
  sig <- cbind( sig, avg = apply( sig, 1, mean) )
  
  plot( tail(sig, 100), plot.type = "single", type = "o" )
  abline(h = 0)
  tail( sig, 100 )
  tail( bbtewma_base$data$X_bm )
  bbtewma_base$spec$ccy
  
  plot( as.simTS(tail(bbtewma_base$data$X_bm, 200)) )
  plot( bbtewma_base$signal$base )
  plot( bbtewma_scl2$signal$scl2 )
  
  
  
  sig_tmp <- cbind( sig )
  # sig_tmp <- window(sig_tmp, "2015-01-01", end(sig_tmp))
  test <- signalTesting.byTrading( X = bbtewma_base$data$X_bm,
                                   sig = (1 + sign(ema(sig_tmp, 1))) / 2,
                                   n_lag = 1, 
                                   tc = 0.002 )
  X_tmp <- na.omit(cbind(bbtewma_base$data$X_bm, test))
  plot( as.simTS(X_tmp) )
  
  stats_fields <- c("cumret", "sds", "sharpe", "maxDD")
  t( descStats( X_tmp )$stats[stats_fields, ] )
  
  
  
  a <- (1 + sign(ema(sig_tmp[ ,3], 1))) / 2
  sum(abs(na.omit(diff(a))))
  nrow(a)
  
  
  
  
  # Dirichlet sampling
  # sig_base <- bbtewma_base$getSignal()
  # sig_scl <- bbtewma_scl$getSignal()
  # sig_scl2 <- bbtewma_scl2$getSignal()
  # sig <- cbind( base = sig_base,
  #               scl = sig_scl,
  #               scl2 = sig_scl2 )
  # sig <- cbind( sig, avg = apply( sig, 1, mean) )
  DS <- dirichletSampling( Y_train = bbtewma_base$data$X_bm,
                           X_train = sig[ ,1],
                           # sclfct = 1,
                           sclfct = NULL,
                           n_lag = 0,
                           weights_fun = "l1" )
  ldens <- lapply( DS, densFUN )
  plot.ldensity( ldens )
  barplot( unlist( lapply( DS, mean ) ) )
  
  
  
  
  # --------------------------------------------------------------------------
  # CAPITALIZATION WEIGHTED BENCHMARK - NONLINEAR SHRINKAGE
  # --------------------------------------------------------------------------
  
  # Obj <- DAARC::BBTEWMA$new()
  # Obj$setCtrl( universe = "dm", 
  #              method = "base" )
  # Obj$spec$name <- "bbtewma_base_dm_qis"
  # Obj$spec$cov_spec <- covCtrl( method = "qis" )
  # Obj$data <- list( X = bbtewma_base$data$X,
  #                   BBS = bbtewma_base$data$BBS,
  #                   X_bm = bbtewma_base$data$X_bm,
  #                   capw = bbtewma_base$data$capw )
  # Obj$computeSignalInsample()
  # # debugonce( Obj$computeSignal )
  # Obj$updateSignal()
  # Obj$save()
  # 
  # 
  # Obj_scl <- Obj$copy()
  # Obj_scl$spec$method <- "scl"
  # Obj$spec$cov_spec <- covCtrl( method = "qis" )
  # Obj_scl$spec$name <- "bbtewma_scl_dm_qis"
  # Obj_scl$signal <- list()
  # Obj_scl$computeSignalInsample()
  # Obj_scl$updateSignal()
  # Obj_scl$save()
  # 
  # 
  # Obj_scl2 <- Obj$copy()
  # Obj_scl2$spec$method <- "scl2"
  # Obj$spec$cov_spec <- covCtrl( method = "qis" )
  # Obj_scl2$spec$name <- "bbtewma_scl2_dm_qis"
  # Obj_scl2$signal <- list()
  # Obj_scl2$computeSignalInsample()
  # Obj_scl2$updateSignal()
  # Obj_scl2$save()
  
  
  
  bbtewma_base_qis <- loadSensor( sensor_name = "bbtewma_base_dm_qis" )
  bbtewma_scl_qis <- loadSensor( sensor_name = "bbtewma_scl_dm_qis" )
  bbtewma_scl2_qis <- loadSensor( sensor_name = "bbtewma_scl2_dm_qis" )
  
  
  
  sig_pearson <- cbind( bbtewma_base$signal$insample[ ,"delta"],
                        bbtewma_scl$signal$insample[ ,"delta"],
                        bbtewma_scl2$signal$insample[ ,"delta"] )
  sig_qis <- cbind( bbtewma_base_qis$signal$insample[ ,"delta"],
                    bbtewma_scl_qis$signal$insample[ ,"delta"],
                    bbtewma_scl2_qis$signal$insample[ ,"delta"] )
  
  sig_pearson <- cbind( bbtewma_base$getSignal(),
                        bbtewma_scl$getSignal(),
                        bbtewma_scl2$getSignal() )
  sig_qis <- cbind( bbtewma_base_qis$getSignal(),
                    bbtewma_scl_qis$getSignal(),
                    bbtewma_scl2_qis$getSignal() )
  
  plot( sig_tmp )
  plot(sig_tmp, plot.type = "single" )
  
  ldens <- apply( sig_tmp, 2, densFUN )
  plot.ldensity( ldens, fillin = FALSE )
    
  
  
  sig_tmp <- cbind( sig_pearson, sig_qis )
  
  test <- signalTesting.byTrading( X = Obj$data$X_bm,
                                   sig = (1 + sign(ema(sig_tmp, 1))) / 2,
                                   n_lag = 1, 
                                   tc = 0 )
  X_tmp <- na.omit(cbind(Obj$data$X_bm, test))
  plot( as.simTS(X_tmp) )
  
  stats_fields <- c("cumret", "sds", "sharpe", "maxDD")
  t( descStats( X_tmp )$stats[stats_fields, ] )
  
  
  
  
  
  # --------------------------------------------------------------------------
  # CAPITALIZATION WEIGHTED BENCHMARK - EWMA COVMAT
  # --------------------------------------------------------------------------
  
  # Obj <- DAARC::BBTEWMA$new()
  # Obj$setCtrl( universe = "dm", 
  #              method = "base" )
  # Obj$spec$name <- "bbtewma_base_dm_ewma"
  # Obj$spec$cov_spec <- covCtrl( method = "ewma", ellipsis = list(tau = 250) )
  # Obj$data <- list( X = bbtewma_base$data$X,
  #                   BBS = bbtewma_base$data$BBS,
  #                   X_bm = bbtewma_base$data$X_bm,
  #                   capw = bbtewma_base$data$capw )
  # Obj$computeSignalInsample()
  # # debugonce( Obj$computeSignal )
  # Obj$updateSignal()
  # Obj$save()  # 
  # 
  # Obj_scl <- Obj$copy()
  # Obj_scl$spec$method <- "scl"
  # Obj$spec$cov_spec <- covCtrl( method = "ewma", ellipsis = list(tau = 250) )
  # Obj_scl$spec$name <- "bbtewma_scl_dm_ewma"
  # Obj_scl$signal <- list()
  # Obj_scl$computeSignalInsample()
  # Obj_scl$updateSignal()
  # Obj_scl$save()  # 
  # 
  # Obj_scl2 <- Obj$copy()
  # Obj_scl2$spec$method <- "scl2"
  # Obj$spec$cov_spec <- covCtrl( method = "ewma", ellipsis = list(tau = 250) )
  # Obj_scl2$spec$name <- "bbtewma_scl2_dm_ewma"
  # Obj_scl2$signal <- list()
  # Obj_scl2$computeSignalInsample()
  # Obj_scl2$updateSignal()
  # Obj_scl2$save()
  
  
  
  
  bbtewma_base_ewma <- loadSensor( sensor_name = "bbtewma_base_dm_ewma" )
  bbtewma_scl_ewma <- loadSensor( sensor_name = "bbtewma_scl_dm_ewma" )
  bbtewma_scl2_ewma <- loadSensor( sensor_name = "bbtewma_scl2_dm_ewma" )
  
  
  sig_pearson <- cbind( bbtewma_base$signal$insample[ ,"delta"],
                        bbtewma_scl$signal$insample[ ,"delta"],
                        bbtewma_scl2$signal$insample[ ,"delta"] )
  sig_qis <- cbind( bbtewma_base_qis$signal$insample[ ,"delta"],
                    bbtewma_scl_qis$signal$insample[ ,"delta"],
                    bbtewma_scl2_qis$signal$insample[ ,"delta"] )
  sig_ewma <- cbind( bbtewma_base_ewma$signal$insample[ ,"delta"],
                    bbtewma_scl_ewma$signal$insample[ ,"delta"],
                    bbtewma_scl2_ewma$signal$insample[ ,"delta"] )
  
  sig_pearson <- cbind( bbtewma_base$getSignal(),
                        bbtewma_scl$getSignal(),
                        bbtewma_scl2$getSignal() )
  sig_qis <- cbind( bbtewma_base_qis$getSignal(),
                    bbtewma_scl_qis$getSignal(),
                    bbtewma_scl2_qis$getSignal() )
  sig_ewma <- cbind( bbtewma_base_ewma$getSignal(),
                    bbtewma_scl_ewma$getSignal(),
                    bbtewma_scl2_ewma$getSignal() )
  
  colnames(sig_pearson) <- paste0(c("base", "scl", "scl2"))
  colnames(sig_qis) <- paste0(c("base", "scl", "scl2"), "_qis")
  colnames(sig_ewma) <- paste0(c("base", "scl", "scl2"), "_ewma")
  
  
  plot( sig_tmp )
  plot(sig_tmp, plot.type = "single" )
  
  ldens <- apply( sig_tmp, 2, densFUN )
  plot.ldensity( ldens, fillin = FALSE )
  
  
  
  
  sig_tmp <- cbind( sig_pearson, sig_qis, sig_ewma )
  
  test <- signalTesting.byTrading( X = Obj$data$X_bm,
                                   sig = (1 + sign(ema(sig_tmp, 1))) / 2,
                                   n_lag = 1, 
                                   tc = 0 )
  X_tmp <- na.omit(cbind(Obj$data$X_bm, test))
  plot( as.simTS(X_tmp) )
  
  stats_fields <- c("cumret", "sds", "sharpe", "maxDD")
  t( descStats( X_tmp )$stats[stats_fields, ] )
  
  
  
  
  
