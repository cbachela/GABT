  

  ############################################################################
  ### BBTurbulence - ILLUSTRATIONS
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     10.10.2020
  # First version:    10.10.2020
  # --------------------------------------------------------------------------
  
  
  # Content:
  
  # 2-D ILLUSTRATION OF SYMMETRIC AND ASYMMETRIC MAHALANOBIS DISTANCE
  # 3-D ILLUSTRATION OF SYMMETRIC AND ASYMMETRIC MAHALANOBIS DISTANCE
  # COMPARE DIFFERENCE IN MD TO ANALYTIC SOLUTION
  # PLOT CONTINUOUS REGIMES COLORCODED (IN-SAMPLE)
  # DISCRIMINATE ON UP- OR DOWNWARD TRENDING TURBULENCE LEVELS
  # MD BY SUMMATION
  
  
  
  # --------------------------------------------------------------------------
  # REQUIRE
  # --------------------------------------------------------------------------
  
  require(RP)
  require(DAARC)
  
  # wd <- "H:/Papers Cyril/PhD_Papers/Good_And_Bad_Turbulence/R/"
  wd <- "R:/Asset_Management/Research_Projects/Multi_Asset/DAA/DAARC/Sensor_Testing/Univariate/Good_And_Bad_Turbulence/"
  source( paste0(wd, "Source/source_functions.R") )
  
  
  BBSObj <- loadSensor( sensor_name = "bbs", b_new = TRUE )
  BBTObj <- loadSensor( sensor_name = "bbturbulence_scl2_dm", b_new = TRUE )  
  X <- BBTObj$data$X[isWeekday(time(BBTObj$data$X)), ]
  X_bm <- BBTObj$data$X_bm[isWeekday(time(BBTObj$data$X_bm)), ]
  states_is <- BBSObj$signal$insample
  
  
  
  # minphase <- 5
  # mincycle <- 21
  theta <- 0.15
  minphase <- 5 * 4
  mincycle <- 21 * 4
  BBS <- bbs(X = cumulated(X_bm, "discrete"), 
             mincycle = mincycle, 
             minphase = minphase, 
             k.peak = minphase, 
             l.peak = minphase, 
             k.trough = minphase, 
             l.trough = minphase, 
             logarithmic = FALSE, 
             theta = theta, 
             e = 0)
  
  plot( cumulated(X_bm, "discrete") )
  abline( v = time(BBS)[which(BBS == -1)], col = "red" )
  abline( v = time(BBS)[which(BBS == 1)], col = "green" )
  lines( cumulated(X_bm, "discrete") )
  
  
  # --------------------------------------------------------------------------
  # 2-D ILLUSTRATION OF SYMMETRIC AND ASYMMETRIC MAHALANOBIS DISTANCE
  # --------------------------------------------------------------------------
  
  # Names <- c("IL", "GR")
  # Names <- c("CA", "CH")
  # Names <- c("IT", "FR")
  # Names <- c("DE", "FR")
  Names <- c("US", "GB")
  
  X_train <- X[ ,Names]
  
  # Mahalanobis distance
  mu <- apply(X_train , 2, mean)
  sigma <- cov(X_train)
  md_all <- mahalanobisTS(X = X)
  b <- quantile( md_all, 0.8 ) * 5
  E_all <- ellipsoid(q = mu, 
                      Q = solve(sigma), 
                      b = b)
  RPE <- rp( object = E_all, 
             spec = rpCtrl(n_sim = 10^4, 
                           algo = "sphere",
                           b_pushy = FALSE) )
  S_all <- getSamples(RPE)
  
  # Bear 
  p_bear <- (1 - BBS) / 2
  TD <- trainingData(Y_train = p_bear, 
                     X_train = X_train)
  mu_bear <- t(TD$X_train) %*% (TD$Y_train / sum(TD$Y_train))
  scnd_mom = t(TD$X_train) %*% (TD$X_train * ((TD$Y_train/sum(TD$Y_train)) %*% 
                                                matrix(1, 1, ncol(TD$X_train))))
  scnd_mom = (scnd_mom + t(scnd_mom))/2
  sigma_bear = scnd_mom - mu_bear %*% t(mu_bear)
  md_bear <- mahalanobisTS(X = X[ ,Names], 
                           center = as.numeric(mu_bear), 
                           covmat = sigma_bear,
                           scl = FALSE)
  E_bear <- ellipsoid(q = mu_bear, 
                      # Q = solve(sigma_bear), 
                      # Q = solve( sigma_bear * sd(md_bear) ),
                      Q = solve(sigma_bear) / sd(md_bear),
                      b = b)
  RPE <- rp( object = E_bear, 
             spec = rpCtrl(n_sim = 10^4, 
                           algo = "sphere",
                           b_pushy = FALSE) )
  S_bear <- getSamples(RPE)
  
  # Bull
  p_bull <- (1 + BBS) / 2
  TD <- trainingData(Y_train = p_bull, 
                     X_train = X_train)
  mu_bull <- t(TD$X_train) %*% (TD$Y_train / sum(TD$Y_train))
  scnd_mom = t(TD$X_train) %*% (TD$X_train * ((TD$Y_train/sum(TD$Y_train)) %*% 
                                                matrix(1, 1, ncol(TD$X_train))))
  scnd_mom = (scnd_mom + t(scnd_mom))/2
  sigma_bull = scnd_mom - mu_bull %*% t(mu_bull)
  md_bull <- mahalanobisTS(X = X[ ,Names], 
                           center = as.numeric(mu_bull), 
                           covmat = sigma_bull,
                           scl = FALSE)
  E_bull <- ellipsoid(q = mu_bull, 
                      # Q = solve(sigma_bull), 
                      # Q = solve( sigma_bull * sd(md_bull) ),
                      Q = solve(sigma_bull) / sd(md_bull), # same same,
                      b = b)
  # The inverse of a matrix that has been multiplied by a non-zero scalar c
  # is equal to the inverse of the scalar multiplied by the inverse of the matrix.
  # (cA)^-1 = c^-1 A^-1
  
  # debugonce(.rpE)
  # debugonce(runifE)
  RPE <- rp( object = E_bull, 
             spec = rpCtrl(n_sim = 10^4, 
                           algo = "sphere",
                           b_pushy = FALSE) )
  S_bull <- getSamples(RPE)
  
  # Delta
  md_delta <- md_bear / sd(md_bear) - md_bull / sd(md_bull)
  plot(md_delta)
  
  
  plot(as.CovMatrix(sigma_bear))
  plot(as.CovMatrix(sigma_bull))
  
  
  
  
  # Plot ellipses
  pch <- rep(19, nrow(md_delta))
  pch[ which(md_delta <= 0) ] <- 10
  colors <- rep("darkgreen", nrow(md_delta))
  colors[ which(md_delta <= 0) ] <- "darkred"
  plot( x = X_train[ ,1], y = X_train[ ,2], pch = pch, col = colors,
        xlim = range(c(X_train[ ,1], S_bear[ ,1], S_bull[ ,1])),
        ylim = range(c(X_train[ ,2], S_bear[ ,2], S_bull[ ,2])),
        xlab = Names[1], ylab = Names[2])
  points( x = as.numeric(X_train[ ,1]), y = as.numeric(X_train[ ,2]), 
          pch = pch, col = colors )
  points( x = S_bear[ ,1], y = S_bear[ ,2], col = "red", pch = 20, cex = 0.5 ) 
  points( x = S_bull[ ,1], y = S_bull[ ,2], col = "green", pch = 20, cex = 0.5 )
  # points( x = S_all[ ,1], y = S_all[ ,2], col = "blue", pch = 19, cex = 0.5 ) 
  points( x = mu_bear[1], y = mu_bear[2], pch = 19, cex = 1.2, col = "red" )
  points( x = mu_bull[1], y = mu_bull[2], pch = 19, cex = 1.2, col = "green" )
  abline(h = 0, col = "grey")
  abline(v = 0, col = "grey")
  
    
  
  
  # plot( x = S_bear[ ,1], y = S_bear[ ,2] )
  # plot( x = S_bull[ ,1], y = S_bull[ ,2] )
  
  
  # Scale and take difference
  MD <- cbind(bear = md_bear, 
              bull = md_bull)
  md_delta <- MD[, "bear"] - MD[, "bull"]
  MD_scl <- scale(MD, FALSE, TRUE)
  md_scl_delta <- MD_scl[, "bear"] - MD_scl[, "bull"]
  
  plot( md_scl_delta )
  plot( cumsum(md_scl_delta) )
  
  plot( ema( md_scl_delta, 0.1) )
  plot( MD_scl, plot.type = "single" )
  
  
  
  md_delta_tmp <- apply(X[ ,Names],
                        1, 
                        FUN = mahalanobisDelta,
                        center_1 = mu_bear,
                        covmat_1 = sigma_bear,
                        center_2 = mu_bull,
                        covmat_2 = sigma_bull,
                        # wghts = rep(1/length(Names), length(Names)),
                        wghts = NULL,
                        b_scl = FALSE )
  
  plot( cbind(md_delta_tmp, MD_scl[ ,"bear"] - MD_scl[ ,"bull"]),
        plot.type = "multiple" )
  
  plot( cbind(md_delta_tmp, MD[ ,"bear"] - MD[ ,"bull"]),
        plot.type = "single" )
   
  plot( md_delta_tmp - (MD[ ,"bear"] - MD[ ,"bull"]) )  # same same
  
  
  

  
  # --------------------------------------------------------------------------
  # 3-D ILLUSTRATION OF SYMMETRIC AND ASYMMETRIC MAHALANOBIS DISTANCE
  # --------------------------------------------------------------------------
  
  
  Names <- c("IL", "GR", "IT")
  # require(MAE)
  # MAE <- mae( b_update = FALSE )
  # X <- MAE$X_sim
  # Names <- c("FI_Corp", "FI_Govt_MT", "FI_Govt_LT")
  
  X_train <- X[ ,Names]
  
  # Mahalanobis distance
  md_all <- mahalanobisTS(X = X)
  b <- quantile( md_all, 0.8 ) * 5
  
  # Bear 
  p_bear <- (1 - BBS) / 2
  TD <- trainingData(Y_train = p_bear, 
                     X_train = X_train)
  mu_bear <- t(TD$X_train) %*% (TD$Y_train / sum(TD$Y_train))
  scnd_mom = t(TD$X_train) %*% (TD$X_train * ((TD$Y_train/sum(TD$Y_train)) %*% 
                                                matrix(1, 1, ncol(TD$X_train))))
  scnd_mom = (scnd_mom + t(scnd_mom))/2
  sigma_bear = scnd_mom - mu_bear %*% t(mu_bear)
  md_bear <- mahalanobisTS(X = X[ ,Names], 
                           center = as.numeric(mu_bear), 
                           covmat = sigma_bear,
                           scl = FALSE)
  E_bear <- ellipsoid(q = mu_bear, 
                      # Q = solve(sigma_bear), 
                      # b = b * 0.7)
                      Q = solve( sigma_bear * sd(md_bear) ),
                      b = b)
  RPE <- rp( object = E_bear, 
             spec = rpCtrl(n_sim = 10^5, 
                           algo = "sphere",
                           b_pushy = FALSE) )
  S_bear <- getSamples(RPE)
  
  # Bull
  p_bull <- (1 + BBS) / 2
  TD <- trainingData(Y_train = p_bull, 
                     X_train = X_train)
  mu_bull <- t(TD$X_train) %*% (TD$Y_train / sum(TD$Y_train))
  scnd_mom = t(TD$X_train) %*% (TD$X_train * ((TD$Y_train/sum(TD$Y_train)) %*% 
                                                matrix(1, 1, ncol(TD$X_train))))
  scnd_mom = (scnd_mom + t(scnd_mom))/2
  sigma_bull = scnd_mom - mu_bull %*% t(mu_bull)
  md_bull <- mahalanobisTS(X = X[ ,Names], 
                           center = as.numeric(mu_bull), 
                           covmat = sigma_bull,
                           scl = FALSE)
  E_bull <- ellipsoid(q = mu_bull, 
                      # Q = solve(sigma_bull), 
                      Q = solve( sigma_bull * sd(md_bull) ),
                      b = b)
  # debugonce(.rpE)
  # debugonce(runifE)
  RPE <- rp( object = E_bull, 
             spec = rpCtrl(n_sim = 10^5, 
                           algo = "sphere",
                           b_pushy = FALSE) )
  S_bull <- getSamples(RPE)
  
  
  plot(as.CovMatrix(sigma_bear))
  plot(as.CovMatrix(sigma_bull))
  
  
  # Plot ellipses
  require(rgl)
  plot3d( x = X_train[ ,1], y = X_train[ ,2], z = X_train[ ,3],
          pch = 4, col = 1,
          # xlim = range(c(X[ ,,Names[1]], S_bear[ ,1], S_bull[ ,1])),
          # ylim = range(c(X[ ,,Names[2]], S_bear[ ,2], S_bull[ ,2])),
          # zlim = range(c(X[ ,,Names[3]], S_bear[ ,3], S_bull[ ,3])),
          xlab = Names[1], ylab = Names[2], zlab = Names[3])
  points3d( x = S_bear[ ,1], y = S_bear[ ,2], z = S_bear[ ,3], col = "red", pch = 19 ) 
  points3d( x = S_bull[ ,1], y = S_bull[ ,2], z = S_bull[ ,3], col = "green", pch = 19 )
  points( x = mu_bear[1], y = mu_bear[2], z = mu_bear[3], size = 2, pch = 19, col = "red" )
  points( x = mu_bull[1], y = mu_bull[2], z = mu_bull[3], size = 2, pch = 19, col = "green" )
  
  

  
  # --------------------------------------------------------------------------
  # COMPARE DIFFERENCE IN MD TO ANALYTIC SOLUTION
  # --------------------------------------------------------------------------
  
  
  center_1 = mu_bear
  covmat_1 = sigma_bear
  center_2 = mu_bull
  covmat_2 = sigma_bull
  wghts = rep(1, length(Names))
  W <- diag(wghts)
  covmat_1_inv <- solve(covmat_1)
  covmat_2_inv <- solve(covmat_2)
  Sigma <- t(W) %*% ( covmat_1_inv - covmat_2_inv ) %*% W
  a_1 <- t(W) %*% covmat_1_inv %*% W %*% center_1
  a_2 <- t(W) %*% covmat_2_inv %*% W %*% center_2
  a <- 2 * (a_2 - a_1)
  g_1 <- t(center_1) %*% t(W) %*% covmat_1_inv %*% W %*% center_1
  g_2 <- t(center_2) %*% t(W) %*% covmat_2_inv %*% W %*% center_2
  constant <- g_1 - g_2

  
  tmp1 <- t(x - center_1) %*% covmat_1_inv %*% (x - center_1)
  tmp2 <- t(x - center_2) %*% covmat_2_inv %*% (x - center_2)
  tmp1 - tmp2
  mahalanobis(x = x, center = mu_bear, cov = sigma_bear)
  mahalanobis(x = x, center = mu_bull, cov = sigma_bull)
  
  
  
  
  today <- "2020-07-02"
  MD[today, ]
  md_delta[today, ]
  x <- as.numeric(X[today, Names])
  ( t(x) %*% Sigma %*% x + t(a) %*% x + constant )
  
 
  
  
  
  # Sample from new ellipse
  
  # x0 <- rep(0, length(Names))
  x0 <- as.numeric(X[today, Names])
  ( t(x0) %*% Sigma %*% x0 + t(a) %*% x0 + constant )
  
  
  rhs <- ( t(x0) %*% Sigma %*% x0 + t(a) %*% x0 ) * 1.1
  rhs <- as.numeric(rhs)
  
  E_delta <- Ellipsoid$new(centre = rep(0, length(Names)),
                           shape = Sigma,
                           # rhs = 0.2,
                           # linear = as.numeric(a),
                           # constant = as.numeric(constant),
                           rhs = rhs,
                           linear = rep(0, length(Names)),
                           constant = 0,
                           interior_point = x0)
  
  # E2$rhs <- as.numeric( sigma + constant + t(qvec) %*% w_mv )

  # debugonce(E_delta$sample)
  # debugonce(har)
  samples <- E_delta$sample( n_sim = 10^2, algo = "hitnrun" )[[2]]
  head(samples)
  plot( x = samples[ ,1], y = samples[ ,2] )
    
  

  
  # --------------------------------------------------------------------------
  # PLOT CONTINUOUS REGIMES COLORCODED (IN-SAMPLE)
  # --------------------------------------------------------------------------
  
  # Absolute turbulence
  signal <- ema( md_all, 0.1 )
  signal <- signal / max(signal)
  colors <- fBasics::divPalette(n = nrow(signal), "RdYlGn")
  color <- colors
  color[order(signal)] <- rev(colors)
  plot( log(cumulated(X_bm, "discrete")) )  
  abline( v = time(signal), col = color )
  abline( h = 0, col = "grey" )
  lines( log(cumulated(X_bm, "discrete")), col = "blue" ) 
  lines( signal )
  
  
  # Relative turbulence
  BBT <- loadSensor( sensor_name = "bbturbulence_scl2_dm", b_new = TRUE ) 
  # BBT$computeSignalInsample()
  
  signal <- ema( BBT$signal$scl2[ ,"md_delta"], 0.1 )
  # signal <- ema( BBT$signal$insample[ ,"delta"], 0.1 )
  sig_pos <- unique( signal[ signal > 0, ] )
  sig_neg <- unique( signal[ signal <= 0, ] )
  colors_pos <- fBasics::divPalette(n = length(sig_pos) * 2, "RdYlGn")
  colors_neg <- fBasics::divPalette(n = length(sig_neg) * 2, "RdYlGn")
  colors <- c( head(colors_neg, length(colors_neg) / 2),
               tail(colors_pos, length(colors_pos) / 2) )
  color <- colors
  color[order(signal)] <- colors
  plot( log(cumulated(X_bm, "discrete")) )  
  abline( v = time(signal), col = color )
  abline( h = 0, col = "grey" )
  lines( log(cumulated(X_bm, "discrete")), col = "blue" ) 
  # lines( signal )
  
  
  
  
  
  
  plot( x = BBT$signal$scl2[ ,"md_delta"],
        y = BBT$signal$scl2[ ,"md_all_scl"] )
  plot( x = abs(BBT$signal$scl2[ ,"md_delta"]),
        y = BBT$signal$scl2[ ,"md_all_scl"] )
  
      
  # BBT$computeSignalInsample()
  plot( x = BBT$signal$insample[ ,"delta"],
        y = BBT$signal$insample[ ,"all"] )
  
  plot( x = BBT$signal$insample[ ,"delta"],
        y = BBT$signal$insample[ ,"bear"] )
  plot( x = BBT$signal$insample[ ,"delta"],
        y = BBT$signal$insample[ ,"bull"] )
  plot( x = BBT$signal$insample[ ,"delta"],
        y = BBT$signal$insample[ ,"all"] )
  
  plot( x = ema(BBT$signal$insample[ ,"bull_scl"], 1),
        y = ema(BBT$signal$insample[ ,"bear_scl"], 1) )
  
  cor( ema(BBT$signal$insample[ ,c("bear_scl", "bull_scl", "all_scl")], 1) )
             
  
  
  cor( x = abs(BBT$signal$insample[ ,"delta"]),
       y = BBT$signal$insample[ ,"bear"] )
  cor( x = abs(BBT$signal$insample[ ,"delta"]),
       y = BBT$signal$insample[ ,"bull"] )
  cor( x = abs(BBT$signal$insample[ ,"delta"]),
       y = BBT$signal$insample[ ,"all"] )
 
  cor( x = BBT$signal$insample[ ,"delta"],
       y = BBT$signal$insample[ ,"bear"] )
  cor( x = BBT$signal$insample[ ,"delta"],
       y = BBT$signal$insample[ ,"bull"] )
  cor( x = BBT$signal$insample[ ,"delta"],
       y = BBT$signal$insample[ ,"all"] )
  
  
  
  
  
  

  
  # --------------------------------------------------------------------------
  # DISCRIMINATE ON UP- OR DOWNWARD TRENDING TURBULENCE LEVELS
  # --------------------------------------------------------------------------
  
  BBT <- loadSensor( sensor_name = "bbturbulence_scl2_dm", b_new = TRUE ) 
  # BBT$computeSignalInsample()
  
  md <- ema( BBT$signal$scl2[ ,"md_all"], 0.05 )
  
  minphase <- 21 * 1
  mincycle <- 21 * 3
  k_peak <- k_trough <- l_peak <- l_trough <- minphase
  theta <- 0.2
  language <- "C"
  logarithmic <- FALSE
  BBS <- bbs( X = md,
              minphase = minphase,
              mincycle = mincycle,
              k.peak = k_peak,
              k.trough = k_trough,
              l.peak = l_peak,
              l.trough = l_trough,
              theta = theta,
              language = language,
              logarithmic = logarithmic)
  
  plot( md )
  abline( v = time(BBS)[ which(BBS == 1)], col = 2 )
  abline( v = time(BBS)[ which(BBS == -1)], col = 3 )
  lines( md, col = 4 )
  
  X_bm <- BBT$data$X_bm
  mean(X_bm[which(BBS == 1), ])
  mean(X_bm[which(BBS == -1), ])
  

  # # LT algo
  # setpar_filtering_alg(tr_bull = 200,
  #                      tr_bear = 0)
  # bull_lt <- run_filtering_alg(index =  md )
  # LT <- timeSeries( matrix(sign(as.numeric(bull_lt)-0.5), ncol = 1), time(md) ) 
  # 
  # plot( md )
  # abline( v = time(LT)[ which(LT == 1)], col = 2 )
  # abline( v = time(LT)[ which(LT == -1)], col = 3 )
  # lines( md, col = 4 )
  # 
  # mean(X_bm[which(LT == 1), ])
  # mean(X_bm[which(LT == -1), ])
  
  
 
  
  BBT2 <- BBT$copy()
  BBT2$setCtrl( method = "base" )
  BBT2$spec$minphase = minphase
  BBT2$spec$mincycle = mincycle
  BBT2$spec$k_peak = k_peak
  BBT2$spec$k_trough = k_trough
  BBT2$spec$l_peak = l_peak
  BBT2$spec$l_trough = l_trough
  BBT2$spec$theta = theta
  BBT2$spec$language = language
  BBT2$spec$logarithmic = logarithmic
  BBT2$data$X_bm <- returns(md, "discrete")
  BBT2$computeSignalInsample()
  
  plot( BBT2$signal$insample )
  
  
  X_bm <- BBT$data$X_bm
  sig <- cbind( (1 - BBT2$signal$insample[ ,"states"]) / 2,
                (1 - sign(BBT2$signal$insample[ ,"delta"])) / 2 )
  sig <- 1 / (1 + exp(BBT2$signal$insample[ ,"delta"]/2) )

  
  MD <- loadSensor( "turbulence" )
  sig_tmp <- MD$signal$base[ ,"sig_md_ema0.157"]
  sig <- sig_tmp * 0 + 1
  
  a <- sig_tmp
  b <- ema( BBT$signal$scl2[ ,"md_delta"], 0.05 )
  plot( cbind(a, b, 1 / (1 + exp(-b * 100))) )
  

  a <- sig_tmp * 0 + 1
  a[ which(sig_tmp < 0.0001), ] <- 0
  sig <- a
  plot(sig)
  
  b <- 1 / (1 + exp(-ema( BBT$signal$scl2[ ,"md_delta"], 0.01 ) * 1000))
  sig <- b
  
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
  plot( as.simTS(X_tmp), col = colors )
  descStats( X_tmp )
  
  plot( as.simTS(X_tmp_nc), col = colors )
  descStats( X_tmp_nc )
  
  tmp <- as.simTS(cbind(X_tmp_nc, X_tmp_nc[ ,1] * 0.5))
  plot( tmp, col = 1:3 )
  descStats( tmp )
  
  
  tmp <- as.simTS(cbind(X_tmp, X_tmp[ ,1] * 0.5))
  plot( tmp, col = 1:3 )
  descStats( tmp )
  
  

  
  ###############################################
  
  BBT <- loadSensor( "bbturbulence_base_dm_lag0to21" )
  # BBT <- loadSensor( "bbturbulence_base_dm" )
  BBT$computeSignalInsample()
  signals <- BBT$signal$base[ ,c("md_bear_scl", "md_bull_scl", "md_all_scl")]
  signals <- BBT$signal$insample[ ,c("bear_scl", "bull_scl", "all_scl")]
  colnames(signals) <- c("bear", "bull", "all")
    
  md_scl <- apply( na.omit(signals[ ,c("bear", "bull")]), 1, 
                   function(x) { x / sum(x) } )
  md_scl <- timeSeries( t(md_scl), time(na.omit(signals)) )
  
  plot(md_scl, plot.type = "single" )
  abline(h = 0.5, col = "grey")
  plot(cbind(md_scl, apply(md_scl, 1, sum)))
  
  
  
  
  colors <- rep(1, nrow(md_scl))
  colors[ which(md_scl[ ,1] > 0.5) ] <- "darkgreen"
  colors[ which(md_scl[ ,1] <= 0.5) ] <- "darkred"
  # dates_good <- rownames(md_scl)[ which(md_scl[ ,1] > 0.5) ]
  # dates_bad <- rownames(md_scl)[ which(md_scl[ ,1] <= 0.5) ]
  md <- signals[rownames(md_scl), "all"] / max(signals[rownames(md_scl), "all"])
  plot( md, col = "white", ylim = c(0, 2.2) )
  points( md, col = colors, pch = 19)  
  lines( 1 + log(cumulated(BBT$data$X_bm[rownames(md_scl), ], "discrete")), col = 4 )
  
  
  
  a <- md_scl[ ,1] * 0
  a[ which(md_scl[ ,1] > 0.5), ] <- 1
  b <- (1 + sign(BBT$signal$insample[rownames(a), "delta"])) / 2
  
  plot( a - b ) # same same
  
  
  
  
  
  # Influence of variance and covariance on MD
  
  BBT <- loadSensor( "bbturbulence_base_dm" )
  X <- BBT$data$X
  Y <- scale( X, TRUE, FALSE )
  # Y <- rmvnorm( n = nrow(X), mean = rep(0, ncol(X)), sigma = cov(X) )
  # Y <- timeSeries(Y, time(X))
  covmat <- cov(Y)
  covmat_inv <- solve(covmat)
  mu <- apply(Y, 2, mean)
  MD <- mahalanobisTS(Y, scl = FALSE)
  
  y_max <- Y[rownames(MD)[which(MD == max(MD))], ]
  y_plus <- Y[which(apply(Y, 1, sum) == max(apply(Y, 1, sum))), ]
  y_minus <- Y[which(apply(Y, 1, sum) == min(apply(Y, 1, sum))), ]
  md_max <- mahalanobis(x = y_max, center = mu, cov = covmat)
  md_plus <- mahalanobis(x = y_plus, center = mu, cov = covmat)
  md_minus <- mahalanobis(x = y_minus, center = mu, cov = covmat)
  
  plot(MD)
  abline(h = md_plus, col = 3)
  abline(h = md_minus, col = 2)
  abline(h = md_max, col = 4)
  
  
  var_max <- y_max^2 * diag(covmat_inv)
  var_plus <- y_plus^2 * diag(covmat_inv)
  var_minus <- y_minus^2 * diag(covmat_inv)
  
  cov_max <- cov_plus <- cov_minus <- NULL
  for (i in 1:nrow(covmat_inv)) {
    for(j in 1:ncol(covmat_inv)) {
      if ( j > i ) {
        cov_max <- c(cov_max, y_max[i] * y_max[j] * covmat_inv[i, j])
        cov_plus <- c(cov_plus, y_plus[i] * y_plus[j] * covmat_inv[i, j])
        cov_minus <- c(cov_minus, y_minus[i] * y_minus[j] * covmat_inv[i, j])
      }
    }
  }
  
  md_max; sum(var_max) + sum(cov_max) * 2
  md_plus; sum(var_plus) + sum(cov_plus) * 2
  md_minus; sum(var_minus) + sum(cov_minus) * 2
  
  
  
  var_max
  cov_max
  mean(var_max); mean(var_minus); mean(var_plus)
  mean(cov_max); mean(cov_minus); mean(cov_plus)
  
  sum(var_max); sum(var_minus); sum(var_plus)
  sum(cov_max); sum(cov_minus); sum(cov_plus)
  
  
  
  
    
  ###############################################
  
  # MD in PC space. GARCH to predict MD
  
  BBT <- loadSensor( "bbturbulence_base_dm" )
  X <- BBT$data$X
  # X <- window(X, start(X_tmp), "2020-03-26") 
  # BBS <- pca( X )
  BBS <- mtca( X )   # do mtca and not pca because of non-convergence issues on pc's
  S <- getter( BBS, "sources" )
  evalues <- apply(S, 2, var)
  fit <- garch(S, garchCtrl(steps = 10) )
  cvol <- getCondVar(fit)
  mforecast <- MultiForecast(fit)
  cvol_fcast <- getCondVar(mforecast)
 
  res_synt <- S * 0 + rnorm(length(S))
  # res_synt <- scale(res_synt, TRUE, TRUE)
  # res_synt <- scale(res_synt, FALSE, 1 / evalues^0.5)
  S_synt <- res_synt * cvol^0.5
  # S_synt <- scale( scale(S_synt, FALSE, TRUE), FALSE, 1 / evalues^0.5 )
  tmp <- cbind( apply(S_synt, 2, var), evalues)
  barplot( t(tmp), beside = TRUE )
  
  
  require(stochvol)
  sv <- NULL
  for ( j in 1:ncol(S) ) {
    draws <- svsample( S[ ,j], draws = 1000, burnin = 100 )
    sv <- cbind(sv, timeSeries( apply(exp(draws$latent), 2, mean), time(S) ))
  }
  plot( sv[ ,1:10] )
  
  
  
  
  md_s <- apply( scale( S^2, FALSE, evalues ), 1, sum )
  md_cvol <- apply( scale( cvol, FALSE, evalues ), 1, sum )
  md <- mahalanobisTS(X, scl = FALSE)
  md_S_synt <- apply( scale(S_synt^2, FALSE, evalues ), 1, sum )
  md_sv <- apply( scale( sv, FALSE, evalues ), 1, sum )
  
  plot(md - md_s)
  plot(md_cvol - md_sv)
  
  
  MD <- cbind( md = md, 
               md_s = md_s,
               md_cvol = md_cvol,
               md_S_synt = md_S_synt,
               md_sv = md_sv )
  plot( MD )
  plot( tail(MD, 200) )
  apply( MD, 2, mean )
  
  # Map to signals
  P <- apply( MD, 2, function(x) { pchisq(x, df = ncol(X)) } )
  sig <- (P - 1) * (-1)
  sig <- timeSeries( sig, time(MD) )
  colnames(sig) <- paste0("sig_", colnames(MD))
  
  plot( sig, plot.type = "single" )
  lines( sig[ ,5], col = 5 )
  
  
  # Forecasts
  md_fcast <- apply( scale( cvol_fcast, FALSE, evalues), 1, sum )
  md_fcast <- timeSeries( md_fcast, as.Date(rownames(md_fcast))+1 )
  plot( tail(rbind(md_cvol, md_fcast), 200) )
  lines( md_fcast, col = 2 )
  
  
  
  # Signal testing
  X_bm <- BBT$data$X_bm
  test <- signalTesting.byTrading( X = X_bm,
                                   sig = sig,
                                   n_lag = 1,
                                   tc = 0.004 )
  colnames(test) <- colnames(sig)
  test_nc <- signalTesting.byTrading( X = X_bm,
                                      sig = sig,
                                      n_lag = 1,  # ~~~
                                      tc = 0 )
  colnames(test_nc) <- colnames(sig)
  X_tmp <- na.omit(cbind(X_bm, test))
  X_tmp_nc <- na.omit(cbind(X_bm, test_nc))
  
  colors <- c(1, fBasics::rainbowPalette(ncol(test)))
  stats_fields <- c("cumret", "means", "sds", "sharpe", "maxDD")
  
  plot( as.simTS(X_tmp), col = colors )
  t( descStats( X_tmp )$stats[stats_fields, ] )
  drawDownStats(X_tmp)
  
  plot( as.simTS(X_tmp_nc), col = colors )
  t( descStats( X_tmp_nc )$stats[stats_fields, ] )
  drawDownStats(X_tmp_nc)
  
  
  
  
  

  
  
  # --------------------------------------------------------------------------
  # MD BY SUMMATION
  # --------------------------------------------------------------------------
  
  
  require(DAARC)
  
  X <- Data[ ,1:5]
  mu <- apply(X, 2, mean)
  covmat <- cov(X)
  E <- eigen(covmat)
  tmat <- E$vectors
  S <- X %*% tmat
  cov(S)
  tmat %*% covmat %*% t(tmat)
  t(tmat) %*% covmat %*% tmat
  apply(tmat^2, 2, sum)
  tmat %*% cov(S) %*% t(tmat) - covmat
  
  covmat_inv <- solve(covmat)
  tmat %*% solve(cov(S)) %*% t(tmat) - covmat_inv
  
  md <- mahalanobisTS(X, scl = FALSE)
  md <- mahalanobis(X, center = apply(X, 2, mean), cov = covmat)
  
  Z <- scale(X, TRUE, FALSE)
  FUN <- function(z)
  {
    ans <- NULL
    for ( i in 1:length(z) ) {
      if ( i < length(z) ) {
        for ( j in (i+1):length(z) ) {
          ans <- c(ans, z[i] * z[j] * covmat_inv[i, j])
        }
      }
    }
    return( ans )
  }
  B <- t(apply(Z, 1, FUN ))
  b <- apply(B, 1, sum)
  a <- apply(Z, 1, function(z) { sum(z^2 * diag(covmat_inv)) } ) 
  MD <- cbind(a, b, a + 2 * b, md)
  plot( MD, plot.type = "single")
  tail(MD)    
  
  x <- as.numeric(X[1, ])
  mahalanobis(x, center = mu, cov = covmat)
  z <- as.numeric(Z[1, ])  
  sum(z^2 * diag(covmat_inv))
  sum(FUN(z))
  
  mean(b)  
  range(b)  
  
  
  
  # --------------------------------------------------------------------------
  # MD by summation
  # --------------------------------------------------------------------------
  
  X <- Data[ ,1:20]
  mu <- apply(X, 2, mean)
  covmat <- cov(X)
  covmat_inv <- solve(covmat)
  
  # Generate multivariate normal data
  Y_synt <- rmvnorm(n = 10^4,
                    mean = mu,
                    sigma = covmat,
                    method = "eigen")
  Y_synt <- asTimeSeries(Y_synt)
  
  md_synt <- mahalanobisTS(Y_synt, scl = FALSE)
  
  Z <- scale(Y_synt, TRUE, FALSE)
  B <- t(apply(Z, 1, FUN ))
  b <- apply(B, 1, sum)
  a <- apply(Z, 1, function(z) { sum(z^2 * diag(covmat_inv)) } ) 
  
  MD <- cbind(a, b, a + 2 * b, md_synt)
  plot( MD, plot.type = "single")
  tail(MD)    
  
  
  mean(2 * b)
  mean(a)
  mean(a) + 2 * mean(b)
  
  mean(md_synt)
  
  
  
  
  
  
  
  
  
  ###############################################
  
  
  # install.packages("biotools")
  # require(biotools)
  # 
  # D2 <- D2.disc( data = X,
  #                grouping = BBS )
  
  
  
  
  
  
  

  
  
  
  
  
  
  
  
  
  
  
  
  
  