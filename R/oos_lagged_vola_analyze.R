  
  
  ############################################################################
  ### GOOD AND BAD TURBULENCE - OUT-OF-SAMPLE ANALYSIS - LAGGED VOLA - ANALYZE
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
  require(simolz)
  require(DAARC)
  require(visolz)
  
  wd <- "R:/Asset_Management/Research_Projects/Multi_Asset/DAA/DAARC/Sensor_Testing/Univariate/Good_And_Bad_Turbulence/"
  source( paste0(wd, "Source/source_functions.R") )
  
  
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Compute signal on lagged vix - BBTEWMA base
  # --------------------------------------------------------------------------
  
  bbt <- loadSensor( sensor_name = "bbtewma_base_laggedvix" )
  bbt_level <- loadSensor( sensor_name = "bbtewma_base_laggedvixlevel" )
  
  
  sig_is <- bbt$signal$insample[ ,"delta"]
  sig_oos <- bbt$getSignal()
  sig_is_level <- bbt_level$signal$insample[ ,"delta"]
  sig_oos_level <- bbt_level$getSignal()
  sig_oos_avg <- (bbt$getSignal() + bbt_level$getSignal()) / 2
  sig <- na.omit( cbind( is = sig_is, 
                         oos = sig_oos,
                         is_level = sig_is_level,
                         oos_level = sig_oos_level,
                         oos_avg = sig_oos_avg ) )
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
  
  stats_fields <- c("cumret", "means", "meansAr", "sds", "sharpe", "maxDD")
  
  #
  X_tmp <- X_tmp[rownames(X_tmp) > "2004-01-01", ]
  X_tmp_nc <- X_tmp_nc[rownames(X_tmp_nc) > "2004-01-01", ]
  #
  plot( as.simTS(X_tmp), logscale = TRUE )
  plot( as.simTS(head(tail(X_tmp, 700), 200)), logscale = TRUE )
  t( descStats( X_tmp )$stats[stats_fields, ] )
  
  plot( as.simTS(X_tmp_nc), logscale = TRUE )
  t( descStats( X_tmp_nc )$stats[stats_fields, ] )
  
  
  
  HC <- hc_lineChart( X = cumulated(X_tmp, "discrete"), 
                      type = "stock",
                      decimals_y = 2,
                      spliner = FALSE,
                      title = "Cumulative Relative Performance",
                      compare = 'percent' )
  HC
  
  
  
  
  
  
  
  
  
  
  