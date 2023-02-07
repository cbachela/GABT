  
  
  ############################################################################
  ### GOOD AND BAD TURBULENCE - RANDOM PORTFOLIOS - DM
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
  
  source( "H:/R/RP/R/class_Polytope.R")
  source( "H:/R/RP/R/class_Simplex.R")
  
  

  
  
  
  
  
  # --------------------------------------------------------------------------
  # CAPITALIZATION WEIGHTED BENCHMARK
  # --------------------------------------------------------------------------
  
  # Object <- loadSensor( sensor_name = "bbtewma_base_dm" )
  # 
  # Obj <- DAARC::BBTEWMA$new()
  # Obj$setCtrl( universe = "dm",
  #              method = "base" )
  # Obj$spec$name <- "bbtewma_base_dm"
  # Obj$data <- list( X = Object$data$X,
  #                   BBS = Object$data$BBS,
  #                   X_bm = Object$data$X_bm,
  #                   capw = Object$data$capw )
  # Obj$computeSignalInsample()
  # Obj$updateSignal()
  # Obj$save()
  # 
  # Obj_scl <- Obj$copy()
  # Obj_scl$spec$method <- "scl"
  # Obj_scl$spec$name <- "bbtewma_scl_dm"
  # Obj_scl$signal <- list()
  # Obj_scl$computeSignalInsample()
  # Obj_scl$updateSignal()
  # Obj_scl$save()
  # 
  # 
  # Obj_scl2 <- Obj$copy()
  # Obj_scl2$spec$method <- "scl2"
  # Obj_scl2$spec$name <- "bbtewma_scl2_dm"
  # Obj_scl2$signal <- list()
  # Obj_scl2$computeSignalInsample()
  # Obj_scl2$updateSignal()
  # Obj_scl2$save()
  
  
  bbtewma_base <- loadSensor( sensor_name = "bbtewma_base_dm" )
  bbtewma_scl <- loadSensor( sensor_name = "bbtewma_scl_dm" )
  bbtewma_scl2 <- loadSensor( sensor_name = "bbtewma_scl2_dm" )
  
  bbtewma_base$update()
  bbtewma_scl$update()
  bbtewma_scl2$update()
  
  
  sig_tmp <- cbind( bbtewma_base$signal$insample[ ,"delta"],
                    bbtewma_scl$signal$insample[ ,"delta"],
                    bbtewma_scl2$signal$insample[ ,"delta"] )
  
  sig_tmp <- cbind( bbtewma_base$getSignal(),
                    bbtewma_scl$getSignal(),
                    bbtewma_scl2$getSignal() )
  
  plot( sig_tmp )
  plot(sig_tmp, plot.type = "single" )
  
  ldens <- apply( sig_tmp, 2, densFUN )
  plot.ldensity( ldens, fillin = FALSE )
    
  
  
  test <- signalTesting.byTrading( X = bbtewma_base$data$X_bm,
                                   sig = (1 + sign(ema(sig_tmp, 1))) / 2,
                                   n_lag = 1, 
                                   tc = 0 )
  X_tmp <- na.omit(cbind(bbtewma_base$data$X_bm, test))
  plot( as.simTS(X_tmp) )
  
  stats_fields <- c("cumret", "sds", "sharpe", "maxDD")
  descStats( X_tmp )$stats[stats_fields, ]
  
  
  
  
  
  # --------------------------------------------------------------------------
  # EQUALLY WEIGHTED BENCHMARK
  # --------------------------------------------------------------------------
  
  # Object <- loadSensor( sensor_name = "bbtewma_base_dm" )
  # 
  # capw <- Object$data$capw
  # eqw <- capw * 0 + 1/ncol(capw)
  # bm_eqw <- simPortfolio( X = Object$data$X,
  #                         wghts = eqw,
  #                         fc = 0, vc = 0 )
  # 
  # Obj <- DAARC::BBTEWMA$new()
  # Obj$setCtrl( universe = "dm",
  #              method = "base" )
  # Obj$spec$name <- "bbtewma_base_dm_eqw"
  # Obj$data <- list( X = Object$data$X[rownames(bm_eqw), ],
  #                   BBS = Object$data$BBS,
  #                   X_bm = bm_eqw,
  #                   capw = eqw )
  # Obj$computeSignalInsample()
  # # debugonce( Obj$computeSignal )
  # Obj$updateSignal()
  # Obj$save()
  
  
  bbtewma_base_eqw <- loadSensor( sensor_name = "bbtewma_base_dm_eqw" )
  bbtewma_base_eqw$update()
  bbtewma_base_eqw$save()
  
  
  #// No need to compute scl or scl2 because there is no difference to base 
  #// since base already assumes equal weighting.
  
  sig_tmp <- cbind( bbtewma_base_eqw$signal$insample[ ,"delta"],
                    bbtewma_base_eqw$getSignal() )
  test <- signalTesting.byTrading( X = bbtewma_base_eqw$data$X_bm,
                                   sig = (1 + sign(sig_tmp)) / 2,
                                   n_lag = 1, 
                                   tc = 0 )
  X_tmp <- na.omit(cbind(bbtewma_base_eqw$data$X_bm, test))
  plot( as.simTS(X_tmp) )
  
  stats_fields <- c("cumret", "sds", "sharpe", "maxDD")
  descStats( X_tmp )$stats[stats_fields, ]
  
  
  
  
  plot( sig_tmp, plot.type = "single" )
  abline( h = 0, col = "grey" )
  5
  
  plot( as.simTS(na.omit(cbind(capw = bbtewma_base$data$X_bm,
                               eqw = bbtewma_base_eqw$data$X_bm,
                               eqw2 = bm_eqw,
                               eqw3 = Obj$data$X_bm))) )
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # RANDOM WEIGHTED BENCHMARK
  # --------------------------------------------------------------------------
  
  Object <- loadSensor( sensor_name = "bbtewma_base_dm" )
  
  set.seed(1234)
  capw <- Object$data$capw
  S <- Simplex$new( d = ncol(capw) )
  samples <- S$runif( n_sim = 10 )
  barplot(t(samples))
  
  random_weights <- capw
  random_weights[1:nrow(capw), ] <- samples[rep(1, nrow(capw)), ]
  bm_rp1 <- simPortfolio( X = Object$data$X, 
                          wghts = random_weights, 
                          fc = 0, vc = 0 )
  
  Obj <- DAARC::BBTEWMA$new()
  Obj$setCtrl( universe = "dm", 
               method = "base" )
  Obj$spec$name <- "bbtewma_base_dm_rp1"
  Obj$data <- list( X = Object$data$X[rownames(bm_rp1), ],
                    BBS = Object$data$BBS,
                    X_bm = bm_rp1,
                    capw = random_weights )
  Obj$computeSignalInsample()
  Obj$updateSignal()
  Obj$save()
  
  Obj_scl <- Obj$copy()
  Obj_scl$spec$method <- "scl"
  Obj_scl$spec$name <- "bbtewma_scl_dm_rp1"
  Obj_scl$signal <- list()
  Obj_scl$computeSignalInsample()
  Obj_scl$updateSignal()
  Obj_scl$save()
  
  Obj_scl2 <- Obj$copy()
  Obj_scl2$spec$method <- "scl2"
  Obj_scl2$spec$name <- "bbtewma_scl2_dm_rp1"
  Obj_scl2$signal <- list()
  Obj_scl2$computeSignalInsample()
  Obj_scl2$updateSignal()
  Obj_scl2$save()
  
  
  
  
  
  
  bbtewma_base_rp1 <- loadSensor( sensor_name = "bbtewma_base_dm_rp1" )
  bbtewma_scl_rp1 <- loadSensor( sensor_name = "bbtewma_scl_dm_rp1" )
  bbtewma_scl2_rp1 <- loadSensor( sensor_name = "bbtewma_scl2_dm_rp1" )
  
  
  sig_tmp <- cbind( bbtewma_base_rp1$signal$insample[ ,"delta"],
                    bbtewma_scl_rp1$signal$insample[ ,"delta"],
                    bbtewma_scl2_rp1$signal$insample[ ,"delta"] )
  
  sig_tmp <- cbind( bbtewma_base_rp1$getSignal(),
                    bbtewma_scl_rp1$getSignal(),
                    bbtewma_scl2_rp1$getSignal() )
  
  
  plot( sig_tmp )
  plot(sig_tmp, plot.type = "single" )
  
  
  test <- signalTesting.byTrading( X = Obj$data$X_bm,
                                   sig = (1 + sign(sig_tmp)) / 2,
                                   n_lag = 1, 
                                   tc = 0.004 )
  X_tmp <- na.omit(cbind(Obj$data$X_bm, test))
  plot( as.simTS(X_tmp) )
  
  stats_fields <- c("cumret", "sds", "sharpe", "maxDD")
  descStats( X_tmp )$stats[stats_fields, ]
  
  
  
  
