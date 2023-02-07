  
  
  ############################################################################
  ### VIX VS. MAHALANOBIS DISTANCE
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     21.11.2022
  # First version:    21.11.2022
  # --------------------------------------------------------------------------
  
  
  
  require(DAARC)
  require(visolz)
  
  
  
  # --------------------------------------------------------------------------
  # Load sensors
  # --------------------------------------------------------------------------
  
  
  # VIX
  vrp_sensor <- loadSensor( sensor_name = "vrp" )
  
  # Mahalanobis distance based sensor
  bbt_base_sensor <- loadSensor( sensor_name = "bbt_base_dm" )
  bbt_scl_sensor <- loadSensor( sensor_name = "bbt_scl_dm" )
  
  # Mahalanobis distance based sensor with ewma smooting on the input data
  bbtewma_base_sensor <- loadSensor( sensor_name = "bbtewma_base_dm" )
  bbtewma_scl_sensor <- loadSensor( sensor_name = "bbtewma_scl_dm" )
  bbtewma_base_pca_sensor <- loadSensor( sensor_name = "bbtewma_base_dm_pca" )
  bbtewma_scl_pca_sensor <- loadSensor( sensor_name = "bbtewma_scl_dm_pca" )
  
  # Mahalanobis distance based sensor with garch smooting and prediction on the input data
  bbtgarch_base_sensor <- loadSensor( sensor_name = "bbtgarch_base_dm" )
  bbtgarch_scl_sensor <- loadSensor( sensor_name = "bbtgarch_scl_dm" )
  bbtgarch_base_pca_sensor <- loadSensor( sensor_name = "bbtgarch_base_dm_pca" )
  bbtgarch_scl_pca_sensor <- loadSensor( sensor_name = "bbtgarch_scl_dm_pca" )
  
  # Mahalanobis distance based sensor with stochastiv volatility smooting and prediction on the input data
  bbtsv_base_pca_sensor <- loadSensor( sensor_name = "bbtsv_base_dm_pca" )
  bbtsv_scl_pca_sensor <- loadSensor( sensor_name = "bbtsv_scl_dm_pca" )
  
  
  
    
  # --------------------------------------------------------------------------
  # Update sensors
  # --------------------------------------------------------------------------
  
  # vrp_sensor$update()
  # bbt_base_sensor$update()
  # bbt_scl_sensor$update()
  # bbtewma_base_sensor$update()
  # bbtewma_scl_sensor$update()
  # bbtewma_base_pca_sensor$update()
  # bbtewma_scl_pca_sensor$update()
  # bbtgarch_base_sensor$update()
  # bbtgarch_scl_sensor$update()
  # bbtgarch_base_pca_sensor$update()
  # bbtgarch_scl_pca_sensor$update()
  # bbtsv_base_pca_sensor$update()
  # bbtsv_scl_pca_sensor$update()
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Extract signals from sensors
  # --------------------------------------------------------------------------
  
  vix <- vrp_sensor$signal$base[ ,"iv"]
  vix_ewma <- ema( vix, 0.1 )
  sig_vix <- na.omit( cbind( vix, vix_ewma ) )
  colnames(sig_vix) <- c("vix", "vix_ewma")
  
  # Raw
  md_base <- bbt_base_sensor$signal$base[ ,"md_all"]
  md_scl <- bbt_scl_sensor$signal$scl[ ,"md_all"]
  sig_md <- cbind( md_base, md_scl )
  colnames(sig_md) <- paste0("raw_", c("base", "scl"))
  sig_md <- na.omit( sig_md )
  
  # Ewma
  md_base_ewma <- bbtewma_base_sensor$signal$base[ ,"md_all"]
  md_scl_ewma <- bbtewma_scl_sensor$signal$scl[ ,"md_all_scl"]
  md_scl_ewma_nostd <- bbtewma_scl_sensor$signal$scl[ ,"md_all"]
  md_base_pca_ewma <- bbtewma_base_pca_sensor$signal$base[ ,"md_all"]
  md_scl_pca_ewma <- bbtewma_scl_pca_sensor$signal$scl[ ,"md_all"]
  sig_md_ewma <- cbind( md_base_ewma, md_scl_ewma, md_scl_ewma_nostd, md_base_pca_ewma, md_scl_pca_ewma )
  colnames(sig_md_ewma) <- paste0("ewma_", c("base", "scl", "scl_nostd", "base_pca", "scl_pca"))
  sig_md_ewma <- na.omit( sig_md_ewma )
  
  # Garch
  md_base_garch <- bbtgarch_base_sensor$signal$base[ ,"md_all"]
  md_scl_garch <- bbtgarch_scl_sensor$signal$scl[ ,"md_all"]
  md_scl_garch_v2 <- (bbtgarch_scl_sensor$signal$scl[ ,"md_bear_scl"] + bbtgarch_scl_sensor$signal$scl[ ,"md_bull_scl"]) / 2
  md_scl_garch_v2_nostd <- (bbtgarch_scl_sensor$signal$scl[ ,"md_bear"] + bbtgarch_scl_sensor$signal$scl[ ,"md_bull"]) / 2
  md_base_pca_garch <- bbtgarch_base_pca_sensor$signal$base[ ,"md_all"]
  md_base_pca_pred_garch <- bbtgarch_base_pca_sensor$signal$base[ ,"md_all_pred5"]
  md_scl_pca_garch <- bbtgarch_scl_pca_sensor$signal$scl[ ,"md_all"]
  md_scl_pca_pred_garch <- bbtgarch_scl_pca_sensor$signal$scl[ ,"md_all_pred5"]
  sig_md_garch <- cbind( md_base_garch, md_scl_garch, md_scl_garch_v2, md_scl_garch_v2_nostd,
                         md_base_pca_garch, md_base_pca_pred_garch , md_scl_pca_garch, md_scl_pca_pred_garch )
  colnames(sig_md_garch) <- paste0("garch_", c("base", "scl", "scl_v2", "scl_v2_nostd", "base_pca",
                                               "base_pca_pred", "scl_pca", "scl_pca_pred"))
  sig_md_garch <- na.omit( sig_md_garch )

  # # SV
  # md_base_pca_sv <- bbtsv_base_pca_sensor$signal$base[ ,"md_all"]
  # md_scl_pca_sv <- bbtsv_scl_pca_sensor$signal$scl[ ,"md_all"]
  # sig_md_sv <- cbind( md_base_pca_sv, md_scl_pca_sv )
  # colnames(sig_md_sv) <- c("base_pca", "scl_pca")
  # sig_md_sv <- na.omit( sig_md_sv )
  
  # SV
  md_base_pca_sv <- (bbtsv_base_pca_sensor$signal$base[ ,"md_bull_scl"]) # + bbtsv_base_pca_sensor$signal$base[ ,"md_bear_scl"]) / 2
  md_scl_pca_sv <- (bbtsv_scl_pca_sensor$signal$scl[ ,"md_bull_scl"] + bbtsv_scl_pca_sensor$signal$scl[ ,"md_bear_scl"]) / 2
  md_scl_pca_sv <- md_scl_pca_sv / 100
  md_scl_pca_sv_nostd <- (bbtsv_scl_pca_sensor$signal$scl[ ,"md_bull"] + bbtsv_scl_pca_sensor$signal$scl[ ,"md_bear"]) / 2
  # md_scl_pca_sv_pred <- (bbtsv_scl_pca_sensor$signal$scl[ ,"md_bull_scl_pred5"] + bbtsv_scl_pca_sensor$signal$scl[ ,"md_bear_scl_pred5"]) / 2
  md_scl_pca_sv_pred <- (bbtsv_scl_pca_sensor$signal$scl[ ,"md_bull_scl_pred1"] + bbtsv_scl_pca_sensor$signal$scl[ ,"md_bear_scl_pred1"]) / 2
  sig_md_sv <- cbind( md_scl_pca_sv, md_scl_pca_sv_nostd, md_scl_pca_sv_pred )
  colnames(sig_md_sv) <- paste0("sv_", c("scl_pca", "scl_pca_nostd", "scl_pca_pred"))
  sig_md_sv <- na.omit( sig_md_sv )
    
 
  # p_base <- pchisq( q = bbt_base_sensor$signal$base[ ,"md_all"], df = ncol(bbt_base_sensor$data$X) )
  # s_base <- (p_base - 1) * (-1)
  
  
  
  
  
  
  plot( sig_vix )
  plot( sig_md )
  plot( sig_md_ewma )
  plot( sig_md_garch )
  plot( sig_md_sv )
  
  plot( scale( sig_md, FALSE, TRUE ), plot.type = "single" )
  plot( scale( sig_md_ewma, FALSE, TRUE ), plot.type = "single" )
  plot( scale( sig_md_garch, FALSE, TRUE ), plot.type = "single" )
  plot( scale( sig_md_sv, FALSE, TRUE ), plot.type = "single" )
  
  
  
  
  

  
  
  # --------------------------------------------------------------------------
  # Rolling quantiles
  # --------------------------------------------------------------------------
  
  
  quantile_th <- 0.95 
  
  
  # Vix
  sig_vix_q <- applyRoll( Data = sig_vix,
                          Width = 0, # expanding window
                          Gap = 252,
                          By = 1, 
                          FUN = function(X) { apply( X, 2, quantile, quantile_th ) } )
  sig_vix_binary <- do.call( cbind, lapply( 1:ncol(sig_vix), function(i) { sig_vix[rownames(sig_vix_q), i] < sig_vix_q[ ,i] } ) ) * 1
  
  
  # Raw
  sig_md_q <- applyRoll( Data = sig_md,
                         Width = 0, # expanding window
                         Gap = 252,
                         By = 1, 
                         FUN = function(X) { apply( X, 2, quantile, quantile_th ) } )
  sig_md_binary <- do.call( cbind, lapply( 1:ncol(sig_md), function(i) { sig_md[rownames(sig_md_q), i] < sig_md_q[ ,i] } ) ) * 1
  
  
  # Ewma
  sig_md_ewma_q <- applyRoll( Data = sig_md_ewma,
                              Width = 0, # expanding window
                              Gap = 252,
                              By = 1, 
                              FUN = function(X) { apply( X, 2, quantile, quantile_th ) } )
  sig_md_ewma_binary <- do.call( cbind, lapply( 1:ncol(sig_md_ewma), function(i) { sig_md_ewma[rownames(sig_md_ewma_q), i] < sig_md_ewma_q[ ,i] } ) ) * 1
  
  
  # Garch
  sig_md_garch_q <- applyRoll( Data = sig_md_garch,
                               Width = 0, # expanding window
                               Gap = 252,
                               By = 1, 
                               FUN = function(X) { apply( X, 2, quantile, quantile_th ) } )
  sig_md_garch_binary <- do.call( cbind, lapply( 1:ncol(sig_md_garch), function(i) { sig_md_garch[rownames(sig_md_garch_q), i] < sig_md_garch_q[ ,i] } ) ) * 1
  
  
  # SV
  sig_md_sv_q <- applyRoll( Data = sig_md_sv,
                            Width = 0, # expanding window
                            Gap = 252,
                            By = 1, 
                            FUN = function(X) { apply( X, 2, quantile, quantile_th ) } )
  sig_md_sv_binary <- do.call( cbind, lapply( 1:ncol(sig_md_sv), function(i) { sig_md_sv[rownames(sig_md_sv_q), i] < sig_md_sv_q[ ,i] } ) ) * 1
  
  
  
  
  
  plot( sig_md_q )
  plot( sig_md_binary )
  
  
  
  
  # --------------------------------------------------------------------------
  # Trading backtest
  # --------------------------------------------------------------------------
  
  
  
  sig_binary <- sig_vix_binary
  sig_binary <- sig_md_binary
  sig_binary <- sig_md_ewma_binary
  sig_binary <- sig_md_garch_binary
  sig_binary <- sig_md_sv_binary
  
  sig_binary <- cbind( sig_vix_binary, sig_md_binary, sig_md_ewma_binary, sig_md_garch_binary, sig_md_sv_binary)
  sig_binary <- cbind( sig_binary, round( apply( sig_binary, 1, mean) ) )
  
  
  
  y <- bbtewma_base_sensor$data$X_bm
  tc <- 0.001
  fc <- 0
  penalty <- 0.5
  test_1 <- signalTesting.byTrading( X = y,
                                     sig = sig_binary,
                                     n_lag = 1,
                                     tc = tc,
                                     fc = fc,
                                     penalty = penalty )
  colnames(test_1) <- paste0(colnames(sig_binary), "_lag1")
  test_1_nc <- signalTesting.byTrading( X = y,
                                        sig = sig_binary,
                                        n_lag = 1,
                                        tc = 0,
                                        fc = fc,
                                        penalty = penalty )
  colnames(test_1_nc) <- paste0(colnames(sig_binary), "_lag1_nc")
  test_2 <- signalTesting.byTrading( X = y,
                                     sig = sig_binary,
                                     n_lag = 2,
                                     tc = tc,
                                     fc = fc,
                                     penalty = penalty )
  colnames(test_2) <- paste0(colnames(sig_binary), "_lag2")
  X_tmp <- na.omit( cbind(y, test_1 ) ) #, test_1_nc ) ) #, test_2) )

  plotSimTS( as.simTS(X_tmp) )
  stats_fields <- c("cumret", "means", "sds", "sharpe", "maxDD")
  lStats <- descStats(X_tmp)
  t( lStats$stats[stats_fields, ] )
  statsolz:::plot.stats( lStats, sortby = "sharpe" )
  statsolz:::plot.stats( lStats, sortby = NULL )
  
  statsolz:::plot.stats( descStats(window(X_tmp, "2010-01-01", end(X_tmp))), sortby = NULL )
  
   
  
  # Without costs
  X_tmp <- na.omit( cbind(y, test_1_nc ) )
  
  plotSimTS( as.simTS(X_tmp) )
  stats_fields <- c("cumret", "means", "sds", "sharpe", "maxDD")
  lStats <- descStats(X_tmp)
  t( lStats$stats[stats_fields, ] )
  statsolz:::plot.stats( lStats, sortby = "sharpe" )
  statsolz:::plot.stats( lStats, sortby = NULL )
  
  
  
  
  dd <- drawDownStats(X_tmp)
  barplot( do.call( cbind, lapply( dd, FUN = function(x) { x[ ,"Depth"] } ) ), beside = TRUE )
  
  apply( na.omit(sig_binary), 2, mean )
  
  
  
  # Highcharter
  colors <- c("black", "red", "green", "blue", "orange", "purple")
  HC <- hc_lineChart( X = cumulated(X_tmp, "discrete"),
                      
                      type = "stock",
                      compare = "percent",
                      color = colors )
  HC

  
  
  
  # --------------------------------------------------------------------------
  # Regime analysis of backtests
  # --------------------------------------------------------------------------
  
  BBObj <- BBSRC$new()
  BBObj$setCtrl()
  BBObj$data <- list( X = X_tmp,
                      X_level = cumulated(bbt_base_sensor$data$X_bm, "discrete") )
  BBObj$runRobust()
  BBObj$phaseStats()
  
  
  
  
  # --------------------------------------------------------------------------
  # Performance after signal
  # --------------------------------------------------------------------------
  
  y <- bbt_base_sensor$data$X_bm
  n_lags <- 0:5
  
  
  # VIX
  # debugonce( pas )
  PAS_vix <- pas( X = y,
                  sig = sig_vix,
                  n_lags = n_lags )
  lreg_vix <- lapply( 1:ncol(sig_vix), FUN = function(j) {
                  lapply( PAS_vix, function(x) { regression( Y_train = x$Y_train, 
                                                             X_train = x$X_train[ ,j], 
                                                             type = "ols" ) } ) } )
  names(lreg_vix) <- colnames(sig_vix)
  i <- 2
  lapply( lreg_vix[[i]], function(x) { summary( x$reg ) } )
  lapply( lreg_vix[[i]], function(x) { x$coeffmat[2, ] } )
  do.call( rbind, lapply( lreg_vix[[i]], function(x) { x$coeffmat[2, ] } ) )
  
  
  # RAW
  PAS_raw <- pas( X = y,
                   sig = sig_md,
                   n_lags = n_lags )
  lreg_raw <- lapply( 1:ncol(sig_md), FUN = function(j) {
    lapply( PAS_raw, function(x) { regression( Y_train = x$Y_train, 
                                                X_train = x$X_train[ ,j], 
                                                type = "ols" ) } ) } )
  names(lreg_raw) <- colnames(sig_md)
  i <- 2
  lapply( lreg_raw[[i]], function(x) { summary( x$reg ) } )
  lapply( lreg_raw[[i]], function(x) { x$coeffmat[2, ] } )
  do.call( rbind, lapply( lreg_raw[[i]], function(x) { x$coeffmat[2, ] } ) )
  
  
  # EWMA
  sig_md_ewma <- bbtewma_scl_sensor$signal$scl[ ,c("md_all", "md_all_scl", "md_bear", "md_bear_scl")]
  tmp <- (bbtewma_scl_sensor$signal$scl[ ,c("md_bear")] + bbtewma_scl_sensor$signal$scl[ ,c("md_bull")]) / 2
  tmp_std <- (bbtewma_scl_sensor$signal$scl[ ,c("md_bear_scl")] + bbtewma_scl_sensor$signal$scl[ ,c("md_bull_scl")]) / 2
  sig_md_ewma <- cbind(tmp, tmp_std)
  PAS_ewma <- pas( X = y,
                  sig = sig_md_ewma,
                  n_lags = n_lags )
  lreg_ewma <- lapply( 1:ncol(sig_md_ewma), FUN = function(j) {
                  lapply( PAS_ewma, function(x) { regression( Y_train = x$Y_train, 
                                                              X_train = x$X_train[ ,j], 
                                                              type = "ols" ) } ) } )
  names(lreg_ewma) <- colnames(sig_md_ewma)
  i <- 4
  lapply( lreg_ewma[[i]], function(x) { summary( x$reg ) } )
  lapply( lreg_ewma[[i]], function(x) { x$coeffmat[2, ] } )
  do.call( rbind, lapply( lreg_ewma[[1]], function(x) { x$coeffmat[2, ] } ) )
  do.call( rbind, lapply( lreg_ewma[[2]], function(x) { x$coeffmat[2, ] } ) )
  do.call( rbind, lapply( lreg_ewma[[3]], function(x) { x$coeffmat[2, ] } ) )
  do.call( rbind, lapply( lreg_ewma[[4]], function(x) { x$coeffmat[2, ] } ) )
  
  
  
  
  # GARCH
  PAS_garch <- pas( X = y,
                    sig = sig_md_garch,
                    n_lags = n_lags )
  lreg_garch <- lapply( 1:ncol(sig_md_garch), FUN = function(j) {
    lapply( PAS_garch, function(x) { regression( Y_train = x$Y_train, 
                                              X_train = x$X_train[ ,j], 
                                              type = "ols" ) } ) } )
  names(lreg_garch) <- colnames(sig_md_garch)
  i <- 2
  lapply( lreg_garch[[i]], function(x) { summary( x$reg ) } )
  lapply( lreg_garch[[i]], function(x) { x$coeffmat[2, ] } )
  do.call( rbind, lapply( lreg_garch[[2]], function(x) { x$coeffmat[2, ] } ) )
  do.call( rbind, lapply( lreg_garch[[3]], function(x) { x$coeffmat[2, ] } ) )
  do.call( rbind, lapply( lreg_garch[[4]], function(x) { x$coeffmat[2, ] } ) )
  
  
  # SV
  PAS_sv <- pas( X = y,
                   sig = sig_md_sv,
                   n_lags = n_lags )
  lreg_sv <- lapply( 1:ncol(sig_md_sv), FUN = function(j) {
    lapply( PAS_sv, function(x) { regression( Y_train = x$Y_train, 
                                                X_train = x$X_train[ ,j], 
                                                type = "ols" ) } ) } )
  names(lreg_sv) <- colnames(sig_md_sv)
  i <- 1
  lapply( lreg_sv[[i]], function(x) { summary( x$reg ) } )
  lapply( lreg_sv[[i]], function(x) { x$coeffmat[2, ] } )
  do.call( rbind, lapply( lreg_sv[[1]], function(x) { x$coeffmat[2, ] } ) )
  do.call( rbind, lapply( lreg_sv[[2]], function(x) { x$coeffmat[2, ] } ) )
  do.call( rbind, lapply( lreg_sv[[3]], function(x) { x$coeffmat[2, ] } ) )
  
  
  
  
  
  # Combined
  sig_comb <- na.omit( cbind( vix = sig_vix[ ,"vix"],
                              md = sig_md[ ,"raw_base"],
                              md_sv = sig_md_sv[ ,"sv_scl_pca"] ) )
  PAS <- pas( X = y,
              sig = sig_comb,
              n_lags = n_lags )
  lreg <- lapply( 1:ncol(sig_comb), FUN = function(j) {
    lapply( PAS, function(x) { regression( Y_train = x$Y_train, 
                                           X_train = x$X_train[ ,j], 
                                           type = "ols" ) } ) } )
  names(lreg) <- colnames(sig_comb)
  do.call( rbind, lapply( lreg[[1]], function(x) { x$coeffmat[2, ] } ) )
  do.call( rbind, lapply( lreg[[2]], function(x) { x$coeffmat[2, ] } ) )
  do.call( rbind, lapply( lreg[[3]], function(x) { x$coeffmat[2, ] } ) )
  
  

  # Plot
  j <- 3
  i <- 1; plot( x = PAS[[i]]$X_train[ ,j], y = PAS[[i]]$Y_train  )
  i <- 2; plot( PAS[[i]]$X_train[ ,j], PAS[[i]]$Y_train  )
  i <- 3; plot( PAS[[i]]$X_train[ ,j], PAS[[i]]$Y_train  )
  i <- 4; plot( PAS[[i]]$X_train[ ,j], PAS[[i]]$Y_train  )
  i <- 5; plot( PAS[[i]]$X_train[ ,j], PAS[[i]]$Y_train  )
  unlist( lapply( PAS, function(x) { nrow(x$DF) } ) )
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Conditional densities
  # --------------------------------------------------------------------------
  
  y <- bbt_base_sensor$data$X_bm
  n_lag <- 1
  # weights_fun <- "kernel"
  weights_fun <- "cmeans"
  correct_bias <- TRUE
  sclfct = 1
  scl_by_entropy = FALSE
  apply_roll_y <- list( fun = NULL,
                        width = NULL,
                        by = NULL )
  # apply_roll_y <- list( fun = function(x) { exp( sum( log(1 + x) ) ) - 1 },
  #                       width = 5,
  #                       by = 1 )

  
  # VIX
  # debugonce( dirichletSampling )
  # debugonce( weightsFun.cmeans )
  DS_vix <- dirichletSampling( Y_train = y,
                               X_train = sig_vix_binary[ ,2],
                               weights_fun = weights_fun,
                               # alpha = 0.9,
                               # centers = 4,
                               n_lag = n_lag,
                               apply_roll_y = apply_roll_y,
                               correct_bias = correct_bias,
                               sclfct = sclfct,
                               scl_by_entropy = scl_by_entropy )
  ldens_vix <- lapply( DS_vix, density )
  colors <- rev( fBasics::divPalette( n = length(DS_vix), "RdYlGn" ) )
  slolz:::plot.ldensity( ldens_vix, fillin = FALSE, col = colors )
  barplot( unlist( lapply( DS_vix, mean ) ), col = colors )
  
  
  # GARCH
  DS_garch <- dirichletSampling( Y_train = y,
                                 X_train = sig_md_garch_binary[ ,5],
                                 weights_fun = weights_fun,
                                 alpha = 0.9,
                                 n_lag = n_lag,
                                 correct_bias = correct_bias,
                                 sclfct = sclfct,
                                 scl_by_entropy = scl_by_entropy )
  ldens_garch <- lapply( DS_garch, density )
  colors <- rev( fBasics::divPalette( n = length(DS_garch), "RdYlGn" ) )
  slolz:::plot.ldensity( ldens_garch, fillin = FALSE, col = colors )
  barplot( unlist( lapply( DS_garch, mean ) ), col = colors )
  
  
  # SV
  DS_sv <- dirichletSampling( Y_train = y,
                              X_train = sig_md_sv_binary[ ,2],
                              weights_fun = weights_fun,
                              # alpha = 0.9,
                              n_lag = n_lag,
                              correct_bias = correct_bias,
                              sclfct = sclfct,
                              scl_by_entropy = scl_by_entropy )
  ldens_sv <- lapply( DS_sv, density )
  colors <- rev( fBasics::divPalette( n = length(DS_sv), "RdYlGn" ) )
  slolz:::plot.ldensity( ldens_sv, fillin = FALSE, col = colors )
  barplot( unlist( lapply( DS_sv, mean ) ), col = colors )
  
  
  
  
  
  
  