  
  
  ############################################################################
  ### GOOD AND BAD TURBULENCE - OUT-OF-SAMPLE ANALYSIS - LAGGED VOLA
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     27.06.2022
  # First version:    27.06.2022
  # --------------------------------------------------------------------------
  
  
  # --------------------------------------------------------------------------
  # REQUIRE
  # --------------------------------------------------------------------------
  
  require(RP)
  require(BSS)
  require(garcholz)
  require(simolz)
  require(DAARC)
  require(visolz)
  
  wd <- "R:/Asset_Management/Research_Projects/Multi_Asset/DAA/DAARC/Sensor_Testing/Univariate/Good_And_Bad_Turbulence/"
  source( paste0(wd, "Source/source_functions.R") )
  
  
  stats_fields <- c("cumret", "means", "meansAr", "sds", "sharpe", "maxDD")
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Load VIX
  # --------------------------------------------------------------------------
  
  # VRP <- loadSensor( sensor_name = "vrp" )
  # VRP$update()
  # vix <- VRP$signal$base[ ,"iv"]
  
  X_sptr <-  rodbcGetOLZDBReturns( assetName = "SPTR",
                                   refCcy = "USD", 
                                   frqncy = "daily", 
                                   compounding = "discrete" )
  X_sptr <- X_sptr[isWeekday(time(X_sptr)), ]
  
  vix_ret <- rodbcGetOLZDBReturns( assetName = "VIX",
                                   refCcy = "USD", 
                                   frqncy = "daily", 
                                   compounding = "continuous" )
  
  vix_ret <- vix_ret[isWeekday(time(vix_ret)), ]
  vix_base <- 21.8999999999999  # this is the level of the vix as of 01.03.1990
  vix <- cumulated( vix_ret, compounding = "continuous" ) * vix_base
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Universe: MSCI World
  # Compute signal on lagged vix or lagged bm-returns
  # BBTEWMA base
  # --------------------------------------------------------------------------
  
  # On vix returns
  vix_lagged <- vix_ret
  for ( i in 1:30 ) {
    vix_lagged <- cbind( vix_lagged, lag(vix_ret, i) )
  }
  colnames(vix_lagged) <- paste0("vix_lag", 0:(ncol(vix_lagged)-1))
  vix_lagged <- na.omit(vix_lagged)
  
  bbt <- loadSensor( sensor_name = "bbtewma_base_dm" )
  bbt$data$capw <- NULL
  bbt$spec$name <- "bbtewma_base_laggedvix"
  bbt$signal <- list()
  bbt$data$X <- vix_lagged
  
  bbt$computeSignalInsample()
  bbt$updateSignal()
  bbt$save()
  
  
  # On vix levels  
  vix_lagged <- vix
  for ( i in 1:30 ) {
    vix_lagged <- cbind( vix_lagged, lag(vix, i) )
  }
  colnames(vix_lagged) <- paste0("vix_lag", 0:(ncol(vix_lagged)-1))
  vix_lagged <- na.omit(vix_lagged)
  
  bbt_level <- loadSensor( sensor_name = "bbtewma_base_dm" )
  bbt_level$data$capw <- NULL
  bbt_level$spec$name <- "bbtewma_base_laggedvixlevel"
  bbt_level$signal <- list()
  bbt_level$data$X <- vix_lagged
  
  bbt_level$computeSignalInsample()
  bbt$updateSignal()
  bbt$save()
  
  
  
  # On lagged bm-returns
  bbt_xlag <- loadSensor( sensor_name = "bbtewma_base_dm" )
  X_bm <- bbt_xlag$data$X_bm
  x_lagged <- X_bm
  for ( i in 1:30 ) {
    x_lagged <- cbind( x_lagged, lag(X_bm, i) )
  }
  colnames(x_lagged) <- paste0("x_lag", 0:(ncol(x_lagged)-1))
  x_lagged <- na.omit(x_lagged)
  
  bbt_xlag$data$capw <- NULL
  bbt_xlag$spec$name <- "bbtewma_base_xbmlagged"
  bbt_xlag$signal <- list()
  bbt_xlag$data$X <- x_lagged
  
  bbt_xlag$computeSignalInsample()
  bbt_xlag$updateSignal()
  bbt_xlag$save()
  
  
  
  
  
  # bbt <- loadSensor( sensor_name = "bbtewma_base_laggedvix" )
  # bbt_level <- loadSensor( sensor_name = "bbtewma_base_laggedvixlevel" )
  # bbt_xlag <- loadSensor( sensor_name = "bbtewma_base_xbmlagged" )
  # bbt$updateData()
  # bbt$computeSignalInsample()
  # bbt$updateSignal()
  # bbt$save()
  
  plot( bbt$signal$insample )
  plot( bbt_level$signal$insample )
  plot( bbt_xlag$signal$insample )
  plot( bbt$signal$base )
  tail( bbt$signal$base )
  
  
  
  # sig <- na.omit( cbind( bbt$signal$getSignal(),
  #                        bbt_level$getSignal(),
  #                        bbt_xlag$getSignal() ) )
  sig <- na.omit( cbind( bbt$signal$insample[ ,"delta"],
                         bbt_level$signal$insample[ ,"delta"],
                         bbt_xlag$signal$insample[ ,"delta"] ) )
  sig <- cbind( (sig[ ,1] + sig[ ,2]) / 2, sig )
  sig_binary <- (1 + sign(ema(sig, 1))) / 2
  X_bm <- bbt$data$X_bm
  n_lag <- 1
  test <- signalTesting.byTrading( X = X_bm,
                                   sig = sig_binary,
                                   n_lag = n_lag, 
                                   tc = 0.004 )
  test_nc <- signalTesting.byTrading( X =X_bm,
                                      sig = sig_binary,
                                      n_lag = n_lag, 
                                      tc = 0 )
  X_tmp <- na.omit( cbind( X_bm, test ) )
  X_tmp_nc <- na.omit( cbind( X_bm, test_nc ) )
  # plot( cumulated(X_tmp, "discrete"), plot.type = "single" )
  #
  # X_tmp <- X_tmp[rownames(X_tmp) > "2004-01-01", ]
  # X_tmp_nc <- X_tmp_nc[rownames(X_tmp_nc) > "2004-01-01", ]
  #
  plot( as.simTS(X_tmp), logscale = TRUE )
  plot( as.simTS(head(tail(X_tmp, 700), 200)), logscale = TRUE )
  t( descStats( X_tmp )$stats[stats_fields, ] )
  
  plot( as.simTS(X_tmp_nc), logscale = TRUE )
  t( descStats( X_tmp_nc )$stats[stats_fields, ] )
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Universe: MSCI World
  # Compute signal on lagged vix or lagged bm-returns
  # BBTGARCH base
  # --------------------------------------------------------------------------

  # On lagged bm-returns
  bbt_xlag <- loadSensor( sensor_name = "bbtgarch_base_dm" )
  X_bm <- bbt_xlag$data$X_bm
  
  x_lagged <- X_bm
  for ( i in 1:20 ) {
    x_lagged <- cbind( x_lagged, lag(X_bm, i) )
  }
  colnames(x_lagged) <- paste0("x_lag", 0:(ncol(x_lagged)-1))
  x_lagged <- na.omit(x_lagged)
  
  bbt_xlag$data$capw <- NULL
  bbt_xlag$spec$name <- "bbtgarch_base_xbmlagged"
  bbt_xlag$signal <- list()
  bbt_xlag$data$X <- x_lagged
  
  bbt_xlag$computeSignalInsample()
  # bbt_xlag$updateSignal()
  # bbt_xlag$save()
  
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Universe: MSCI World Country Indics (loop over indices)
  # Compute signal on lagged realized volas
  # BBTEWMA base
  # --------------------------------------------------------------------------
  
  
  bbt_tmp <- loadSensor( sensor_name = "bbtewma_base_dm" )
  bbt_tmp$data$capw <- NULL
  
  # Country indices
  X <- getMSCIData( universe = "dm", 
                    frqncy = "d", 
                    ccy = "USD" )
  X <- X[isWeekday(time(X)), ]
  
  # Realized volas
  vola <- applyRoll( Data = X,
                     Width = 21,
                     Gap = 0, 
                     By = 1, 
                     FUN = function(X) { apply(X, 2, sd) } )
  vola_ret <- returns( vola, "discrete" )[-1, ]
  
  lSignal <- list()
  lSignal_level <- list()
  lSignal_XBMlagged <- list()
  
  for ( j in 1:ncol(X) ) {
    
    # Benchmark (used to compute BBQ)
    X_bm <- X[ ,j]
    
    # # Vola returns
    # vola_lagged <- vola_ret[ ,j]
    # for ( i in 1:30 ) {
    #   vola_lagged <- cbind( vola_lagged, lag(vola_ret[ ,j], i) )
    # }
    # colnames(vola_lagged) <- paste0("vola_lag", 0:(ncol(vola_lagged)-1))
    # vola_lagged <- na.omit(vola_lagged)
    # 
    # bbt <- bbt_tmp$copy()
    # bbt$spec$name <- paste0("bbtewma_base_laggedvola_", colnames(X)[j])
    # bbt$data$X <- vola_lagged
    # bbt$data$X_bm <- X_bm
    # bbt$signal <- list()
    # 
    # bbt$computeSignalInsample()
    # bbt$updateSignal()
    # bbt$save()
    # lSignal[[j]] <- bbt$signal$insample
    # 
    # # Vola levels
    # vola_lagged <- vola
    # for ( i in 1:30 ) {
    #   vola_lagged <- cbind( vola_lagged, lag(vola[ ,j], i) )
    # }
    # colnames(vola_lagged) <- paste0("vola_lag", 0:(ncol(vola_lagged)-1))
    # vola_lagged <- na.omit(vola_lagged)
    # 
    # bbt <- bbt_tmp$copy()
    # bbt$spec$name <- paste0("bbtewma_base_laggedvolalevel_", colnames(X)[j])
    # bbt$data$X <- vola_lagged
    # bbt$data$X_bm <- X_bm
    # bbt$signal <- list()
    # 
    # bbt$computeSignalInsample()
    # bbt$updateSignal()
    # bbt$save()
    # lSignal_level[[j]] <- bbt$signal$insample
    
    # Lagged returns
    X_lagged <- X_bm
    for ( i in 1:30 ) {
      X_lagged <- cbind( X_lagged, lag(X_bm, i) )
    }
    colnames(X_lagged) <- paste0("X_lag", 0:(ncol(X_lagged)-1))
    X_lagged <- na.omit(X_lagged)
    
    bbt <- bbt_tmp$copy()
    bbt$spec$name <- paste0("bbtewma_base_laggedXbm_", colnames(X)[j])
    bbt$data$X <- X_lagged
    bbt$data$X_bm <- X_bm
    bbt$signal <- list()
    
    bbt$computeSignalInsample()
    bbt$updateSignal()
    bbt$save()
    lSignal_XBMlagged[[j]] <- bbt$signal$insample
    
  }
  names(lSignal) <- colnames(X)
  names(lSignal_level) <- colnames(X)
  
  
  lX_tmp <- list()
  for ( j in 1:length(lSignal) ) {
    X_bm <- X[ ,j]
    sig <- lSignal[[j]][ ,"delta"]
    sig_level <- lSignal_level[[j]][ ,"delta"]
    sig_xlagged <- lSignal_XBMlagged[[j]][ ,"delta"]
    sig <- na.omit( cbind( vola = sig, 
                           vola_level = sig_level, 
                           x_lagged = sig_xlagged) )
    sig <- cbind( sig, avg = (sig[ ,1] + sig[ ,2]) / 2 )
    colnames(sig) <- c("vola", "vola_level", "x_lagged", "vola_avg")
    sig_binary <- (1 + sign(sig)) / 2
    test <- signalTesting.byTrading( X = X_bm,
                                     sig = sig_binary,
                                     n_lag = 2, 
                                     tc = 0.004 )
    lX_tmp[[j]] <- na.omit( cbind( X_bm, test ) )
  }
  names(lX_tmp) <- names(lSignal)
  
  
  stats_fields <- c("cumret", "means", "meansAr", "sds", "sharpe", "maxDD")
  lStats <- lapply( lX_tmp, function(X) { descStats(X)$stats[stats_fields, ] } )
  stats <- do.call( cbind, lStats )
  stats
  barplot( stats["means", ], col = 1:5 )
  barplot( stats["sds", ], col = 1:5 )
  barplot( stats["maxDD", ], col = 1:5 )
  
  
  plot( as.simTS(lX_tmp[["GR"]]) )
  plot( lSignal[["GR"]] )
  
  
  
  
  
  
  
  sig <- bbt$signal$insample[ ,"delta"]
  sig_binary <- ( 1 + sign( ema( sig, 0.5 ) ) ) / 2
  
  X_bm <- bbt$data$X_bm
  test <- signalTesting.byTrading( X = X_bm,
                                   sig = sig_binary,
                                   n_lag = 2, 
                                   tc = 0.004 )
  test_nc <- signalTesting.byTrading( X =X_bm,
                                      sig = sig_binary,
                                      n_lag = 2, 
                                      tc = 0 )
  X_tmp <- na.omit( cbind( X_bm, test ) )
  X_tmp_nc <- na.omit( cbind( X_bm, test_nc ) )
  # plot( cumulated(X_tmp, "discrete"), plot.type = "single" )
  
  
  
  # #
  # X_tmp <- X_tmp[rownames(X_tmp) > "2004-01-01", ]
  # X_tmp_nc <- X_tmp_nc[rownames(X_tmp_nc) > "2004-01-01", ]
  # #
  plot( as.simTS(X_tmp), logscale = TRUE )
  plot( as.simTS(head(tail(X_tmp, 700), 200)), logscale = TRUE )
  t( descStats( X_tmp )$stats[stats_fields, ] )
  
  plot( as.simTS(X_tmp_nc), logscale = TRUE )
  t( descStats( X_tmp_nc )$stats[stats_fields, ] )
  
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Universe: S&P 500 - long history
  # Compute signal on lagged realized volas
  # BBTEWMA base
  # --------------------------------------------------------------------------
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Universe: Fama French 10 Industry Portfolios - long history
  # Compute signal on lagged realized volas
  # BBTEWMA base
  # --------------------------------------------------------------------------
  
  # Return series (FF industry portfolios)
  X <- get.file( wd = paste0(wd, "Data/"), filename = "FF_Industry_10", sep = ";" )
  X <- X[isWeekday(time(X)), ]

  # Realized volas
  vola <- applyRoll( Data = X,
                     Width = 21,
                     Gap = 0, 
                     By = 1, 
                     FUN = function(X) { apply(X, 2, sd) } )
  vola_ret <- returns( vola, "discrete" )[-1, ]
  
  # Load sensor (template)
  bbt_tmp <- loadSensor( sensor_name = "bbtewma_base_dm" )
  bbt_tmp$data$capw <- NULL
  
  lSignal <- list()
  lSignal_level <- list()
  lSignal_XBMlagged <- list()
  
  for ( j in 1:ncol(X) ) {
    
    # Benchmark (used to compute BBQ)
    X_bm <- X[ ,j]
    
    # Vola returns
    vola_lagged <- vola_ret[ ,j]
    for ( i in 1:30 ) {
      vola_lagged <- cbind( vola_lagged, lag(vola_ret[ ,j], i) )
    }
    colnames(vola_lagged) <- paste0("vola_lag", 0:(ncol(vola_lagged)-1))
    vola_lagged <- na.omit(vola_lagged)
    
    bbt <- bbt_tmp$copy()
    bbt$spec$name <- paste0("bbtewma_base_laggedvola_", colnames(X)[j])
    bbt$data$X <- vola_lagged
    bbt$data$X_bm <- X_bm
    bbt$signal <- list()
    
    bbt$computeSignalInsample()
    # debugonce( bbt$computeSignal )
    bbt$updateSignal()
    bbt$save()
    lSignal[[j]] <- bbt$signal$insample
    
    # Vola levels
    vola_lagged <- vola
    for ( i in 1:30 ) {
      vola_lagged <- cbind( vola_lagged, lag(vola[ ,j], i) )
    }
    colnames(vola_lagged) <- paste0("vola_lag", 0:(ncol(vola_lagged)-1))
    vola_lagged <- na.omit(vola_lagged)
    
    bbt <- bbt_tmp$copy()
    bbt$spec$name <- paste0("bbtewma_base_laggedvolalevel_", colnames(X)[j])
    bbt$data$X <- vola_lagged
    bbt$data$X_bm <- X_bm
    bbt$signal <- list()
    #
    bbt$data$BBS$spec$language <- "C"
    #
    bbt$computeSignalInsample()
    bbt$updateSignal()
    bbt$save()
    lSignal_level[[j]] <- bbt$signal$insample
    
    # Lagged returns
    X_lagged <- X_bm
    for ( i in 1:30 ) {
      X_lagged <- cbind( X_lagged, lag(X_bm, i) )
    }
    colnames(X_lagged) <- paste0("X_lag", 0:(ncol(X_lagged)-1))
    X_lagged <- na.omit(X_lagged)
    
    bbt <- bbt_tmp$copy()
    bbt$spec$name <- paste0("bbtewma_base_xbmlagged_", colnames(X)[j])
    bbt$data$X <- X_lagged
    bbt$data$X_bm <- X_bm
    bbt$signal <- list()
    #
    bbt$data$BBS$spec$language <- "C"
    #
    bbt$computeSignalInsample()
    bbt$updateSignal()
    bbt$save()
    lSignal_XBMlagged[[j]] <- bbt$signal$insample
    
    
  }
  names(lSignal) <- colnames(X)
  names(lSignal_level) <- colnames(X)
  
  
  lX_tmp <- list()
  for ( j in 1:length(lSignal) ) {
    X_bm <- X[ ,j]
    sig <- lSignal[[j]][ ,"delta"]
    sig_level <- lSignal_level[[j]][ ,"delta"]
    sig_xlagged <- lSignal_XBMlagged[[j]][ ,"delta"]
    sig <- na.omit( cbind( vola = sig, 
                           vola_level = sig_level, 
                           x_lagged = sig_xlagged) )
    sig <- cbind( sig, avg = (sig[ ,1] + sig[ ,2]) / 2 )
    colnames(sig) <- c("vola", "vola_level", "x_lagged", "vola_avg")
    sig_binary <- (1 + sign(sig)) / 2
    test <- signalTesting.byTrading( X = X_bm,
                                     sig = sig_binary,
                                     n_lag = 2, 
                                     tc = 0.004 )
    lX_tmp[[j]] <- na.omit( cbind( X_bm, test ) )
  }
  names(lX_tmp) <- names(lSignal)
  
  
  stats_fields <- c("cumret", "means", "meansAr", "sds", "sharpe", "maxDD")
  lStats <- lapply( lX_tmp, function(X) { descStats(X)$stats[stats_fields, ] } )
  stats <- do.call( cbind, lStats )
  stats
  barplot( stats["means", ], col = 1:5 )
  barplot( stats["sds", ], col = 1:5 )
  barplot( stats["maxDD", ], col = 1:5 )
  
  
  plot( as.simTS(lX_tmp[[1]]) )
  
  
  
  
  
  # bbt_tmp <- loadSensor( sensor_name = "bbtewma_base_laggedvola_NoDur" )
  # plot( bbt_tmp$signal$insample )
  # plot( bbt_tmp$signal$base )
  
  
  
  