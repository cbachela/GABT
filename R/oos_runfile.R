  
  
  ############################################################################
  ### GOOD AND BAD TURBULENCE - OUT-OF-SAMPLE ANALYSIS
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     08.01.2021
  # First version:    10.10.2020
  # --------------------------------------------------------------------------
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # REQUIRE
  # --------------------------------------------------------------------------
  
  require(RP)
  require(DAARC)
  
  wd <- "R:/Asset_Management/Research_Projects/Multi_Asset/DAA/DAARC/Sensor_Testing/Univariate/Good_And_Bad_Turbulence/"
  source( paste0(wd, "Source/source_functions.R") )
  
  
  
  # --------------------------------------------------------------------------
  # PARAMETERS
  # --------------------------------------------------------------------------
  
 
  
  
  
  # --------------------------------------------------------------------------
  # 1) TURBULENCE MEASURES
  # --------------------------------------------------------------------------
  
  # Absolute turbulence
  Turb <- loadSensor( sensor_name = "turbulence_base_dm", b_new = TRUE )
  # TurbScl <- loadSensor( sensor_name = "turbulence_scl_dm", b_new = TRUE )
  TurbScl2 <- loadSensor( sensor_name = "turbulence_scl2_dm", b_new = TRUE )
  TurbMTEWMA <- loadSensor( sensor_name = "turbulence_mtewma_dm", b_new = TRUE )
  # TurbMTEWMAScl2 <- loadSensor( sensor_name = "turbulence_mtewma_scl2_dm", b_new = TRUE )
  

  
  # TurbMTEWMAScl2 <- TurbMTEWMA$copy()
  # TurbMTEWMAScl2$data <- TurbScl2$data
  # TurbMTEWMAScl2$signal <- list()
  # TurbMTEWMAScl2$computeSignalInsample()
  # TurbMTEWMAScl2$updateSignal()
  # TurbMTEWMAScl2$spec$name <- "turbulence_mtewma_scl2_dm"
  # TurbMTEWMAScl2$save()
  
  Turb$update()
  TurbScl2$update()
  TurbMTewma$update()
  
  
  
  
  
  sig_md <- cbind( base = Turb$signal$insample[ ,"md"],
                   scl2 = TurbScl2$signal$insample[ ,"md"],
                   mtewma = TurbMTEWMA$signal$insample[ ,"md"],
                   mtewmascl2 = TurbMTEWMAScl2$signal$insample[ ,"md"] )
  sig_md <- scale( sig_md, FALSE, TRUE )
  plot( sig_md, plot.type = "single" )
  
  sig_md <- cbind( base = Turb$signal$base[ ,"md"],
                   scl2 = TurbScl2$signal$scl2[ ,"md"],
                   mtewma = TurbMTEWMA$signal$mtewma[ ,"md"],
                   mtewma = TurbMTEWMAScl2$signal$insample[ ,"md"] )
  sig_md <- scale( sig_md, FALSE, TRUE )
  plot( sig_md, plot.type = "single" )
  
  
  
  
  
  
  
  # Relative turbulence
  
  BBT <- loadSensor(sensor_name = "bbturbulence_base_dm" )
  BBTScl2 <- loadSensor(sensor_name = "bbturbulence_scl2_dm" )
  BBTewma <- loadSensor(sensor_name = "bbtewma_base_dm")
  BBTewmaScl2 <- loadSensor(sensor_name = "bbtewma_scl2_dm")
  ST <- loadSensor(sensor_name = "bbturbulence_1_base_dm" )
  STScl2 <- loadSensor(sensor_name = "bbturbulence_1_scl2_dm" )
  STewma <- loadSensor(sensor_name = "bbtewma_1_base_dm" )
  STewmaScl2 <- loadSensor(sensor_name = "bbtewma_1_scl2_dm" )
  
  
  
  
  BBTewmaScl2 <- DAARC::BBTEWMA$new()
  BBTewmaScl2$setCtrl( universe = "dm", method = "scl2" )
  # BBTewmaScl2$updateData()
  # lData <- BBTewmaScl2$data
  BBTewmaScl2$data <- lData
  BBTewmaScl2$data$wmat <- BBTScl2$data$wmat
  
  
  
  BBT$computeSignalInsample()
  BBTScl2$computeSignalInsample()
  BBTewma$computeSignalInsample()
  BBTewmaScl2$computeSignalInsample()
  ST$computeSignalInsample()
  STScl2$computeSignalInsample()
  STewma$computeSignalInsample()
  STewmaScl2$computeSignalInsample()
  
  
  
  plot( BBTewmaScl2$signal$insample )

  

  sig_mdrel <- cbind( base = BBT$signal$insample[ ,"delta"],
                      scl2 = BBTScl2$signal$insample[ ,"delta"],
                      ewma = BBTewma$signal$insample[ ,"delta"],
                      ewmascl2 = BBTewmaScl2$signal$insample[ ,"delta"],
                      st = ST$signal$insample[ ,"delta"],
                      stscl2 = STScl2$signal$insample[ ,"delta"],
                      stewma = STewma$signal$insample[ ,"delta"],
                      stewmascl2 = STewmaScl2$signal$insample[ ,"delta"] )
  sig_mdrel <- scale( na.omit(sig_mdrel), FALSE, TRUE )
  plot( sig_mdrel, plot.type = "single" )
  
  
  
  
  ### BBT using sign of returns
  ST <- BBTEWMA$new()
  ST$setCtrl( universe = "dm", method = "scl2")
  ST$spec$k_peak <- ST$spec$l_peak <- ST$spec$k_trough <- ST$spec$l_trough <- 1
  ST$spec$minphase <- 1
  ST$spec$mincycle <- 1
  ST$spec$name <- "bbtewma_1_scl2_dm"
  # ST$updateData()
  # lData <- ST$data
  ST$data <- lData
  ST$computeSignalInsample()
  ST$update()
  ST$save()
  
  
  
  
  
  
  
  ### Checks
  X <- TurbScl2$data$X
  wmat <- TurbScl2$data$wmat
  MT <- mtca( X )
  S <- getter(MT, "sources")
  tmat <- getter(MT, "torsion")
  Sw <- timeSeries( (scale(X, TRUE, FALSE) * wmat) %*% tmat, time(X) )
  # Sw <- timeSeries( (scale(X, TRUE, FALSE) * (wmat * 0 + 1)) %*% tmat, time(X) )
  
  
  md1 <- mahalanobisTS(X = X, wmat = NULL, scl = FALSE)
  md2 <- mahalanobisTS(X = X, wmat = wmat, scl = FALSE)
  md3 <- mahalanobisTS(X = S, wmat = NULL, scl = FALSE)
  md4 <- mahalanobisTS(X = Sw, 
                       center = apply(S, 2, mean),
                       covmat = cov(S), 
                       wmat = NULL, 
                       scl = FALSE)
  md5 <- timeSeries( apply( scale(Sw^2, FALSE, apply(S, 2, var)), 1, sum ),
                     time(S) )
  
  plot( md1 - md3 )
  plot( md2 - md4 )
  plot( md4 - md5 )
  
  
  
  
  
  
  
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # 2) DESCRIPTIVE STATISTICS
  # --------------------------------------------------------------------------
  
  # --------------------------------------------------------------------------
  # Posterior means
  # --------------------------------------------------------------------------
  
  TD <- trainingData( Y_train = Turb$data$X_bm,
                      X_train = sig_md )
  
  TD <- trainingData( Y_train = BBT$signal$insample[ ,"states"],
                      X_train = sig_mdrel )
  
 
  lDS <- list()
  for ( j in 1:ncol(TD$X_train) ) {
    
    lDS[[j]] <- dirichletSampling( Y_train = TD$Y_train,
                                   X_train = TD$X_train[ ,j],
                                   # X_eval = X_eval,
                                   weights_fun = "l1",
                                   correct_bias = FALSE,
                                   sclfct = NULL )
    # lDSW[[j]] <- dirichletSampling( Y_train = TDW$Y_train,
    #                                 X_train = TDW$X_train[ ,j],
    #                                 # X_eval = X_eval,
    #                                 weights_fun = "kernel",
    #                                 # centers = 5,
    #                                 correct_bias = FALSE,
    #                                 sclfct = 1 )
    # lDSM[[j]] <- dirichletSampling( Y_train = TDM$Y_train,
    #                                 X_train = TDM$X_train[ ,j],
    #                                 weights_fun = "kernel",
    #                                 correct_bias = FALSE,
    #                                 sclfct = NULL )
  }
  
  names(lDS) <- colnames(TD$X_train)
  
  
  # from <- quantile( unlist(lDS), 0 )
  # to <- quantile( unlist(lDS), 1 )
  from <- abs(min( unlist(lDS) )) * 0.7 * sign(min( unlist(lDS)))
  to <- max( unlist(lDS) ) * 1.3
  n <- 999
  lldens <- lapply( lDS, FUN = function(x) { lapply(x, density, from = from, to = to, n = n) } )
  for ( i in 1:length(lldens) ) {
    plot.ldensity( lldens[[i]], main = paste0("Posterior means conditional on X_train[ ,", i, "]" ) )
  }
  

  lMu <- lapply( lDS, FUN = function(x) { unlist(lapply(x, FUN = mean)) } )
  Mu <- do.call( cbind, lMu )
  
  colors <- rev(fBasics::divPalette(n = nrow(Mu), "RdYlGn"))
  colors <- rev(fBasics::divPalette(n = nrow(Mu), "Spectral"))
  barplot( Mu, beside = TRUE, col = colors )
  names(lMu)
  
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # BACKTESTS
  # --------------------------------------------------------------------------

  tc <- 0.004
  n_lag <- 2
  
  # Binary signals
  
  TD <- trainingData( Y_train = BBT$data$X_bm,
                      X_train = (1 + sign(ema(sig_mdrel, 1))) / 2 )
  
  test <- signalTesting.byTrading( X = TD$Y_train,
                                   sig = TD$X_train,
                                   n_lag = n_lag,
                                   tc = tc )
  colnames(test) <- colnames(TD$X_train)
  test_nc <- signalTesting.byTrading( X = TD$Y_train,
                                      sig = TD$X_train,
                                      n_lag = n_lag,
                                      tc = 0 )
  colnames(test_nc) <- colnames(TD$X_train)
  
  X_tmp <- na.omit( cbind( bm = TD$Y_train, test ) )
  X_tmp_nc <- na.omit( cbind( bm = TD$Y_train, test_nc ) )
  
  descStats( X_tmp )
  descStats( X_tmp_nc )
  
  colors <- c(1, fBasics::rainbowPalette(n = ncol(X_tmp)-1))
  plot( as.simTS(X_tmp), col = colors )
  plot( as.simTS(X_tmp_nc), col = colors )
  
  
  plot( as.simTS(X_tmp[ ,c("bm", "ewma", "stewma", "ewmascl2", "stewmascl2")]) )
  plot( as.simTS(X_tmp_nc[ ,c("bm", "ewma", "stewma", "ewmascl2", "stewmascl2")]) )
  
  plot( as.simTS(X_tmp[ ,c("bm", "base", "st", "scl2", "stscl2")]) )
  plot( as.simTS(X_tmp_nc[ ,c("bm", "base", "st", "scl2", "stscl2")]) )
  
  
      
  
  # --------------------------------------------------------------------------
  # Kernel conditional CDF
  # --------------------------------------------------------------------------

  require(np)
  require(volesti)
  require(RP)
  
  states <- (1 + BBT$signal$insample[ ,"states"]) / 2
  z <- as.numeric(states)
  dens <- density(z)
  plot(dens)  
  CDF <- cumsum( dens$y / sum(dens$y) )
  
  plot( x = dens$x, y = CDF )  

  
  # Varsi
  q_vec <- seq(from = min(z), to = max(z), length.out = 500)
  # FUN <- function(z0) { tail( as.numeric(varsi( mu = z, b = z0 )), 1) }
  FUN <- function(z0) { frustum_of_simplex( a = z, z0 = z0 ) }
  p_vec <- unlist( lapply( q_vec, FUN ) )  
  
  plot( x = q_vec, y = p_vec )
  abline( v = mean(z) )
  
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # ABSOLUTE AND RELATIVE TURBULENCE
  # --------------------------------------------------------------------------
  
  require(RP)
  require(DAARC)
  
  wd <- "R:/Asset_Management/Research_Projects/Multi_Asset/DAA/DAARC/Sensor_Testing/Univariate/Good_And_Bad_Turbulence/"
  source( paste0(wd, "Source/source_functions.R") )
  

  # Absolute turbulence 
  MD <- loadSensor( sensor_name = "turbulence_scl2_dm" )
  MD$update()
  
  # Relative turbulence
  BBTEwma <- loadSensor( sensor_name = "bbtewma_base_dm_a0.1" )
  BBTEwma$update()
  
  
  sig_abs <- MD$signal$scl2[ ,c("md_ema0.04")]
  sig_rel <- ema( BBTEwma$signal$base[ ,c("signal")], 0.1 )
  
  signals <- cbind(abs = sig_abs, rel = sig_rel)
  plot(signals)
  
  
  
  quant <- quantile(sig_abs, 0.7)
  bsig_abs <- sig_abs * 0 + 1
  bsig_abs[ sig_abs >= quant, ] <- 0
  bsig_rel <- sig_rel * 0 + 1
  bsig_rel[ sig_rel < 0, ] <- 0
  
  
  bsignals <- na.omit( cbind(bsig_abs, bsig_rel) )
  colnames(bsignals) <- c("abs", "rel")
  plot(bsignals)
  
  
  
  test <- signalTesting.byTrading( X = MD$data$X_bm,
                                   sig = bsignals,
                                   n_lag = 1,
                                   tc = 0.004 )
  X_tmp <- na.omit( cbind( MD$data$X_bm, test ) )
  plot( as.simTS(X_tmp) )
  
  

  
  states <- bsignals[ ,1] * NA
  states[1, ] <- 1
  for ( i in 2:nrow(bsignals) ) {
  
    if ( bsignals[i, 1] == bsignals[i, 2] ) {
      
      states[i, ] <- bsignals[i, 1]
      
    } else {
      
      if ( states[i-1, ] == 1 ) {
        
        # Decision rule to sell
        if ( bsignals[i, "abs"] == 1 && bsignals[i, "rel"] == 0 ) {
          if ( bsignals[i-1, "rel"] == 1 ) {
            states[i, ] <- 0
          } else {
            states[i, ] <- 1
          }
        }
        if ( bsignals[i, "abs"] == 0 && bsignals[i, "rel"] == 1 ) {
          if ( bsignals[i-1, "abs"] == 1 ) {
            states[i, ] <- 0
          } else {
            states[i, ] <- 1
          }
        }
        
      } else {
        
        # Decision rule to buy
        if ( bsignals[i, "abs"] == 1 && bsignals[i, "rel"] == 0 ) {
          if ( bsignals[i-1, "rel"] == 1 ) {
            states[i, ] <- 1
          } else {
            states[i, ] <- 0
          }
        }
        if ( bsignals[i, "abs"] == 0 && bsignals[i, "rel"] == 1 ) {
          if ( bsignals[i-1, "abs"] == 1 ) {
            states[i, ] <- 1
          } else {
            states[i, ] <- 0
          }
        }
      }
    }
  }
  
  
  plot(states)
  
  
  par( mfrow = c(2, 1) ) 
  plot( signals[ ,1] )
  abline( h = quant )
  plot( signals[ ,2] )
  abline( h = 0 )
  
  
  
  states <- bsignals[ ,"abs"]
  idx <- apply( bsignals, 1, function(x) { x[1] == 0 && x[2] == 1 } )
  states[idx, ] <- 1
  
      
      
  
  test <- signalTesting.byTrading( X = MD$data$X_bm,
                                   sig = cbind(bsignals, states),
                                   n_lag = 1,
                                   tc = 0.00 )
  X_tmp <- na.omit( cbind( MD$data$X_bm, test ) )
  plot( as.simTS(X_tmp) )
  
    
  
  plot( cbind(bsignals, states) )
  
  
  
  
  
  possible_regimes <- rbind(c(1, 1),
                            c(1, 0),
                            c(0, 1),
                            c(0, 0))
  rownames(possible_regimes) <- 1:nrow(possible_regimes)
  colnames(possible_regimes) <- c("Absolute", "Relative")
  # Create binary timeseries for each regime
  FUN <- function ( s ) { apply(possible_regimes, 1,
                                FUN = function ( r ) { as.numeric(all(s == r)) }) }
  regimes <- timeSeries(t(apply(bsignals, 1, FUN = FUN)), time(bsignals))
  colnames(regimes) <- paste0("Regime ", 1:ncol(regimes))
  # Create categorical timeseries containing regimes 1 to 4
  regimes_flag <- apply( regimes, 1, function(x) { which(x == 1) } )
  
 
  # Regime lengths
  regimes_lengths <- apply( regimes, 2, sum )
  regimes_lengths
  
  
  
  # Geometric Return by Regime
  P <- apply(regimes, 2, function(x) { x / sum(x) } )
  Y <- MD$data$X_bm[rownames(regimes), ]
  mu_geo <- exp( t(P) %*% log(1 + Y) ) - 1
  mu_geo_pm <- exp( t(P) %*% log(1 + Y) * 21) - 1
  mu_geo_pa <- exp( t(P) %*% log(1 + Y) * 252) - 1
  # Volatility by Regime
  sds <- sqrt( t(P) %*% Y^2 - (t(P) %*% Y)^2 )
  sds_pm <- sds * sqrt(21)
  sds_pa <- sds * sqrt(252)
  # Statistics by Regime
  cond_stats <- cbind(mu_geo_pm, sds_pm)
  colnames(cond_stats) <- c(paste0("Mean Geo ", colnames(mu_geo_pm)),
                            paste0("SD ", colnames(sds_pm)))
  
  
  plot( regimes )
  plot( regimes_flag )  
  cond_stats    
  
  dev.off()
  plot( log(cumulated(Y, "discrete")) )
  abline( v = time(regimes_flag)[which(regimes_flag == 1)], col = "darkgreen" )  
  abline( v = time(regimes_flag)[which(regimes_flag == 2)], col = "yellow" )  
  abline( v = time(regimes_flag)[which(regimes_flag == 3)], col = "orange" )  
  abline( v = time(regimes_flag)[which(regimes_flag == 4)], col = "darkred" )  
  lines( log(cumulated(Y, "discrete")) )
  
  
  
  
  
  
