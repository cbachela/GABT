  
  
  ############################################################################
  ### GOOD AND BAD TURBULENCE - REGIME CONDITIONAL DENSITIES
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     20.12.2020
  # First version:    20.12.2020
  # --------------------------------------------------------------------------

  
  
  # Content:
  
  # CRISP TURBULENCE REGIMES
  # CRISP VIX REGIMES
  # COMPARE VIX AND MD
  
  
  
  
  
  require(DAARC)
  
  
  # universes <- c("dm", "em", "europe")
  universe <- "dm"
  
  
  # --------------------------------------------------------------------------
  # LOAD SENSORS
  # --------------------------------------------------------------------------
  
  MD <- loadSensor( sensor_name = "turbulence_base_dm" )
  MD_scl2 <- loadSensor( sensor_name = "turbulence_scl2_dm" )
  BBT <- loadSensor( sensor_name = "bbturbulence_base_dm" )
  BBT_scl2 <- loadSensor( sensor_name = "bbturbulence_scl2_dm" )
  VRP <- loadSensor( sensor_name = "vrp" )
  
  # MD$update()
  # MD_scl2$update()
  # BBT$update()
  # BBT_scl2$update()
  # VRP$update()
  
  
  
  
  
  # --------------------------------------------------------------------------
  # CRISP TURBULENCE REGIMES
  # --------------------------------------------------------------------------
  
  # Calm, Turbulent, Good Turbulent, Bad Turbulent
  
  # md <- BBT$signal$insample[ ,c("all_scl", "bear_scl", "bull_scl", "delta")]
  # md <- BBT$signal$base[ ,c("md_all_scl", "md_bear_scl", "md_bull_scl", "md_delta")]
  md <- BBT_scl2$signal$scl2[ ,c("md_all_scl", "md_bear_scl", "md_bull_scl", "md_delta")]
  colnames(md) <- c("all", "bear", "bull", "delta")
  TD <- trainingData( Y_train = BBT$data$X_bm,
                      X_train = md )
  
  md_regimes <- md[ ,rep(1, 4)] * 0
  colnames(md_regimes) <- c("quiet", "turbulent", "bad_turbulent", "good_turbulent")
  md_regimes[md[ ,"all"] < quantile(md[ ,"all"], 0.9), 1] <- 1
  md_regimes[md[ ,"all"] > quantile(md[ ,"all"], 0.9), 2] <- 1
  md_regimes[md[ ,"delta"] < quantile(md[ ,"delta"], 0.1), 3] <- 1
  md_regimes[md[ ,"delta"] > quantile(md[ ,"delta"], 0.9), 4] <- 1
 
  possible_regimes <- rbind( c(1, 0, 0, 0), 
                             c(1, 0, 1, 0), 
                             c(1, 0, 0, 1),
                             c(0, 1, 0, 0),
                             c(0, 1, 1, 0),
                             c(0, 1, 0, 1) )
  rownames(possible_regimes) <-  c("quiet", "quiet-bad", "quiet-good",
                                   "turbulent", "turbulent-bad", "turbulent-good")
  
  # Create binary timeseries for each regime 
  regimes <- t( apply( md_regimes, 1,
                       function(s) { apply(possible_regimes, 1, 
                                           function(r) { as.numeric(all(s == r)) } ) } ) )
  regimes <- timeSeries(regimes, time(md_regimes) )
  # colnames(regimes) <- paste0("regime_", 1:ncol(regimes))
  colnames(regimes) <- rownames(possible_regimes)
  
  # Create categorical timeseries containing regimes 1 to 6
  market_regimes <- apply( regimes, 1, function(x) { which(x == 1) } )
  
  plot( regimes )
  plot( market_regimes )
  
  
  n_sim <- 10^3
  Y_train <- lag(TD$Y_train, 0)
  Y_train[is.na(Y_train)] <- 0
  y <- as.numeric(Y_train)
  lmu_regimes <- list()
  lcdens_regimes <- list()
  lcdens_regimes_fev <- list()
  lcdens_regimes_kernel <- list()
  
  for ( j in 1:ncol(regimes) ) {
    
    alpha1 <- regimes[ ,j] / sum(regimes[ ,j])
    P <- rdirichlet( n = n_sim,
                     alpha = alpha1 )
    lcdens_regimes[[j]] <- P %*% y
    
    alpha2 <- regimes[ ,j] / sum(regimes[ ,j]) * nrow(regimes)
    P <- rdirichlet( n = n_sim, 
                     alpha = alpha2 )
    lmu_regimes[[j]] <- P %*% y
    
    P_fev <- RP::fevBias( P, q = 100 )
    lcdens_regimes_fev[[j]] <- P_fev %*% y
    
    lcdens_regimes_kernel[[j]] <- density( y[ regimes[ ,j] == 1 ] )
  }
  
  # ldens <- lapply( lcdens_regimes, density )
  # ldens <- lapply( lcdens_regimes_fev, density )
  ldens <- lcdens_regimes_kernel
  names(ldens) <- colnames(regimes)
  ldens[["unconditional"]] <- density(TD$Y_train)

  # colors <- 1:ncol(regimes)
  colors <- c("black", "green", "red", "darkred", "darkgreen")
  plot.ldensity( ldens[ c("unconditional", "quiet", "turbulent", "turbulent-bad", "turbulent-good")], 
                 colors = colors, 
                 main = "Regime-Conditional Densities" )
  lines(density(TD$Y_train))
  abline(v = 0, col = "darkgrey")
  
  
  # colors <- c("black", "green", "red", "darkred", "darkgreen", "magenta", "cyan")
  # plot.ldensity( ldens[ c("unconditional", "quiet", "turbulent", "turbulent-bad", "turbulent-good",
  #                         "quiet_bad", "quiet-good")], 
  #                colors = colors, 
  #                main = "Regime-Conditional Densities" )
  # lines(density(TD$Y_train))
  # abline(v = 0, col = "darkgrey")
  
  
  
  ldens <- lapply( lmu_regimes, density )
  names(ldens) <- colnames(regimes)
  colors <- 1:ncol(regimes)
  plot.ldensity( ldens, colors = colors, 
                 main = "Posterior Conditional Means" )
  abline(v = 0, col = "darkgrey")
  
  
  
  # DS <- dirichletSampling( Y_train = TD$Y_train,
  #                          X_train = TD$X_train,
  #                          weights_mat = apply(regimes, 2, function(x) { x / sum(x) } ),
  #                          correct_bias = FALSE,
  #                          sclfct = NULL,
  #                          scl_by_entropy = TRUE )
  # ldens <- lapply( DS, density )
  # plot.ldensity( ldens )
  
  
  plot( md )
  
  
  
  
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # CRISP VIX REGIMES
  # --------------------------------------------------------------------------
  
  vix <- VRP$signal$base[ ,"iv"]
  X_bm <- BBT$data$X_bm
  
  vix_regimes <- vix[ ,rep(1, 2)] * 0
  colnames(vix_regimes) <- c("calm", "volatile")
  vix_regimes[vix < quantile(vix, 0.9), 1] <- 1
  vix_regimes[vix > quantile(vix, 0.9), 2] <- 1
  vix_regimes
  
  
  dates <- intersect( rownames(X_bm), rownames(vix_regimes) )
  dates <- intersect( dates, rownames(md) )
  DS_VIX <- dirichletSampling( Y_train = X_bm[dates, ],
                               X_train = X_bm[dates, ] * 0,
                               weights_mat = apply(vix_regimes[dates, ], 2, function(p) { p / sum(p) } ),
                               sclfct = 1 )
  ldens_vix <- lapply( DS_VIX, density )
  names( ldens_vix ) <- colnames(vix_regimes)
  ldens_vix[["unconditional"]] <- density(X_bm[dates, ])
  
  plot.ldensity( ldens_vix[ c("unconditional", "calm", "volatile")],
                 colors = c("black", "darkgreen", "darkred"),
                 main = "", 
                 fillin = TRUE )
  lines(density(X_bm[dates, ]))
  abline(v = 0, col = "darkgrey")
  
  
  # For comparison:
  ldens_vix2 <- list()
  ldens_vix2[[1]] <- density( X_bm[dates[vix_regimes[dates, 1] == 1], ] )
  ldens_vix2[[2]] <- density( X_bm[dates[vix_regimes[dates, 2] == 1], ] )
  names( ldens_vix2 ) <- colnames(vix_regimes)
  ldens_vix2[["unconditional"]] <- density(X_bm[dates, ])
  
  plot.ldensity( ldens_vix2[ c("unconditional", "calm", "volatile")],
                 colors = c("black", "darkgreen", "darkred"),
                 main = "", 
                 fillin = TRUE )
  lines(density(X_bm[dates, ]))
  abline(v = 0, col = "darkgrey")
  
  
  
  
  # --------------------------------------------------------------------------
  # COMPARE VIX AND MD
  # --------------------------------------------------------------------------
  
  ldens <- lcdens_regimes_kernel
  names(ldens) <- colnames(regimes)
  ldens[["unconditional"]] <- density(TD$Y_train)
  ldens_all <- c(ldens, ldens_vix)
  
  Names <- c("unconditional", "quiet", "calm", "turbulent", "volatile")
  colors <- c("black", "green", "darkgreen", "red", "tomato")
  
  plot.ldensity( ldens_all[ Names ],
                 colors = colors )
  
  mu <- setNames( unlist( lapply( ldens_all[Names], FUN = function(x) { sum(x$x * x$y ) } ) ), Names )
  barplot( mu, col = colors )   
  
  
  
  
  