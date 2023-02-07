
  
  ############################################################################
  ### BBTurbulence - DISTRIBUTION OF DISTANCES
  ############################################################################

  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     30.10.2020
  # First version:    30.10.2020
  # --------------------------------------------------------------------------
  
  
  # --------------------------------------------------------------------------
  # REQUIRE
  # --------------------------------------------------------------------------
  
  require(RP)
  require(DAARC)
  
  # wd <- "H:/Papers Cyril/PhD_Papers/Good_And_Bad_Turbulence/R/"
  wd <- "R:/Asset_Management/Research_Projects/Multi_Asset/DAA/DAARC/Sensor_Testing/Univariate/Good_And_Bad_Turbulence/"
  source( paste0(wd, "Source/source_functions.R") )
  




  # --------------------------------------------------------------------------
  # DISTRIBUTION OF DISTANCES
  # --------------------------------------------------------------------------
  
  # Compare the distribution of empirical distances with theoretical ones.
  
  Turb <- loadSensor( sensor_name = "turbulence" )
  Turb_scl2 <- loadSensor( sensor_name = "turbulence_scl2" )
  
  X <- Turb$data$X
  mu <- apply(X, 2, mean)
  covmat <- cov(X)
  md_is <- mahalanobisTS( X, center = mu, covmat = covmat, scl = FALSE )
  md_chisq <- timeSeries( rchisq(n = nrow(md_is), df = ncol(X)), time(X) )  
  
 
  # Synthetic
  # Calculate MD on synthetic data set to proxy theoretical distribution
  Y_synt <- mvtnorm::rmvnorm(n = nrow(X),
                             mean = mu,
                             sigma = covmat,
                             method = "eigen")
  Y_synt <- timeSeries(Y_synt, time(X))
  # Y_synt <- scale(Y_synt, FALSE, FALSE)
  
  cov(X) - cov(Y_synt)
  mu - apply(Y_synt, 2, mean)
  
  
  # md_synt_norm <- mahalanobisTS(Y_synt, center = mu, covmat = covmat, scl = FALSE)
  md_synt_norm <- mahalanobisTS(Y_synt, scl = FALSE )

  
  
  signals <- na.omit( cbind( md_oos = Turb$signal$base[ ,"md"],
                             # md_oos_scl2 = MD_scl2$signal$scl2[ ,"md"],
                             md_is = md_is,
                             md_synt_norm = md_synt_norm,
                             md_chisq = md_chisq ) )
  ldens <- apply( signals, 2, density )
  plot.ldensity( ldens, fillin = FALSE )

  plot( ldens[[1]] )  

      
  # fitdistr( x = signals[ ,2], 
  #           densfun = "chi-squared", 
  #           start = list(df = 1) )
  
  
  
  # --------------------------------------------------------------------------
  # COMPARE EMPRICAL VS. CHI-SQUARE
  # --------------------------------------------------------------------------
  
  signals <- na.omit( cbind( empirical = md_is,
                             md_chisq = md_chisq,
                             theoretical = md_synt_norm ) )
  ldens <- apply( signals, 2, density )
  
  
  colors <- c("darkgreen", "darkred", "blue")
  par( mfrow = c(1, 2) )
  plot( signals, plot.type = "single", col = colors,
        main = "Squared Mahalanobis Distances", ylab = "d_t" )
  legend("topright", colnames(signals), lwd = 2, 
         col = colors, text.col = colors, bty = "n", cex = 1)
  plot.ldensity( ldens, fillin = FALSE, col = colors,
                 main = "Density of Squared Mahalanobis Distances",
                 legend = NULL,
                 xlab = "d_t")
  legend("topright", colnames(signals), lwd = 2, 
         col = colors, text.col = colors, bty = "n", cex = 1)
 
  # dev.off()
  
  
  # --------------------------------------------------------------------------
  # DISTRIBUTION OF WEIGHTED DISTANCES
  # --------------------------------------------------------------------------
  
  # Generate random portfolio weights
  W <- rdirichlet(n = 10^3, alpha = rep(1/ncol(x), ncol(X)) )
  
  mu <- apply(X, 2, mean)
  covmat <- cov(X)
  today <- "2020-03-03"
  x <- X[today, ]
  
  FUN <- function(w) 
  { 
    weightedMahalanobis( wghts = w, 
                         x = x, 
                         center = mu, 
                         covmat = covmat, 
                         scl = FALSE )
  } 
  wdist <- apply( W, 1, FUN )
  mdist <- weightedMahalanobis( wghts = rep(1/ncol(X), ncol(x)), 
                                x = x, 
                                center = mu, 
                                cov = covmat, 
                                scl = FALSE )
  
  plot( density(wdist) )
  abline( v = mdist )  
  
  range(wdist)
  mdist  
  
  
  
  
  # --------------------------------------------------------------------------
  # DISTRIBUTION OF DISTANCES UNDER A GARCH MODEL
  # --------------------------------------------------------------------------

  md <- mahalanobisTS( X, scl = FALSE )
  covmat <- cov(X)
  PCA <- pca( X )
  S <- getter( PCA, "sources" )[ ,-21]
  evalues <- apply(S, 2, var)
  S_scl <- scale( S^2, FALSE, evalues )
  md_S <- apply( S_scl, 1, sum )
  fit <- garch( S )
  cvol <- getCondVar(fit)
  # lRes <- lapply( fit@fit, FUN = function(x) { x@fit$residuals } )
  # res <- timeSeries( do.call( cbind, lRes ), time(X) )
  # res_std <- res / cvol
  res_std <- S / cvol
  md_res <- mahalanobisTS( X = res_std, scl = FALSE )  

  # Generate synthetic data
  set.seed(1234)
  # lRes_synt <- lapply( apply(S, 2, sd), FUN = function(x) { rnorm(nrow(S), 0, sd = x) } )
  # res_synt <- timeSeries( do.call( cbind, lRes_synt ), time(S) )
  res_synt <- S * 0 + rnorm(length(S))  # leads to same result
  S_synt <- res_synt * cvol
  md_S_synt <- apply( scale(S_synt^2, FALSE, apply(S_synt, 2, var)), 1, sum )
  # S_synt2 <- scale( scale( S_synt, FALSE, TRUE), FALSE, 1 / evalues^0.5 )
  # md_S_synt2 <- apply( scale(S_synt2^2, FALSE, evalues), 1, sum ) # same same

  
  MD <- cbind( md = md, 
               md_S = md_S, 
               md_res = md_res, 
               md_chisq = md_chisq,
               md_S_synt = md_S_synt )
  ldens <- apply( MD, 2, density )
  plot.ldensity( ldens, col = 1:ncol(MD), fillin = FALSE )
  plot( MD )
  
  plot( cbind(md, md_S, md_S_synt), plot.type = "single" )
  
  
  # Using the rugarch package
  require(rugarch)
  spec <- ugarchspec( variance.model = list(model = "sGARCH", 
                                            garchOrder = c(4, 4), 
                                            submodel = NULL, 
                                            external.regressors = NULL, 
                                            variance.targeting = FALSE), 
                      mean.model = list(armaOrder = c(0, 0), 
                                        include.mean = FALSE, 
                                        archm = FALSE, 
                                        archpow = 1, 
                                        arfima = FALSE, 
                                        external.regressors = NULL,
                                        archex = FALSE),
                      distribution.model = "norm" )
  cvol2 <- cvol * NA
  for ( j in 1:ncol(S) ) {
    fit <- ugarchfit( spec = spec,
                      data = S[ ,j] )
    cvol2[ ,j] <- fit@fit$sigma^2
  }

  
  res2_std <- S / cvol2
  md_res2 <- mahalanobisTS( X = res2_std, scl = FALSE )
  
  set.seed(1234)
  res2_synt <- S * 0 + rnorm(length(S))
  S2_synt <- res2_synt * cvol2
  md_S2_synt <- apply( scale(S2_synt^2, FALSE, apply(S2_synt, 2, var)), 1, sum )
  
  
  
  MD <- cbind( md = md, 
               # md_S = md_S, 
               # md_res = md_res, 
               md_res2 = md_res2,
               md_chisq = md_chisq,
               md_S2_synt = md_S2_synt )
  ldens <- apply( MD, 2, density )
  plot.ldensity( ldens, col = 1:ncol(MD), fillin = FALSE )
  plot( MD )
  
  plot( md_res - md_res2 )
  plot( cbind(md_S_synt, md_S2_synt), plot.type = "single" )
  
  
  plot(S[ ,1:10])
  
  
  
  # Generate Chi-square signal with synthetic md series
  MD <- cbind( md = md, 
               md_res2 = md_res2, 
               md_chisq = md_chisq,
               md_S2_synt = md_S2_synt )
  P <- apply(MD, 2, function(x) {
    pchisq(x, df = ncol(X))
  })
  sig <- (P - 1) * (-1)
  
  plot( sig )     
  plot( MD )
  tail(sig, 200)
  
  
  # Signal testing
  X_bm <- Turb$data$X_bm
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
  
  
  
  
  
  
  
  
  
  
  