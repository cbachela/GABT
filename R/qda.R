  
  
  ############################################################################
  ### BBTurbulence - QDA
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     29.10.2020
  # First version:    10.10.2020
  # --------------------------------------------------------------------------
  
  
  # --------------------------------------------------------------------------
  # REQUIRE
  # --------------------------------------------------------------------------
  
  require(MASS)
  require(RP)
  require(DAARC)
  require(GPO)
  
  # wd <- "H:/Papers Cyril/PhD_Papers/Good_And_Bad_Turbulence/R/"
  wd <- "R:/Asset_Management/Research_Projects/Multi_Asset/DAA/DAARC/Sensor_Testing/Univariate/Good_And_Bad_Turbulence/"
  source( paste0(wd, "Source/source_functions.R") )
  
  
  BBSObj <- loadSensor( sensor_name = "bbs", b_new = TRUE )
  BBTObj <- loadSensor( sensor_name = "bbturbulence", b_new = TRUE )  
  X <- BBTObj$data$X[isWeekday(time(BBTObj$data$X)), ]
  X_bm <- BBTObj$data$X_bm[isWeekday(time(BBTObj$data$X_bm)), ]
  X_m <- aggMonthly(X, "discrete")
  X_bm_m <- aggMonthly(X_bm, "discrete")
  
  
  # minphase <- 3
  # mincycle <- 7
  minphase <- 3 * 4 * 5
  mincycle <- 7 * 4 * 5
  theta <- 0.15
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
  # QUADRATIC DISCRIMINANT ANALYSIS
  # --------------------------------------------------------------------------
  
  # debugonce( MASS:::qda.default )
  QDA <- qda( x = X,
              grouping = BBS,
              prior = c(0.5, 0.5) )
  QDA
  ls(QDA)
  QDA$ldet
  
  # debugonce( MASS:::predict.qda )
  pred <- predict( QDA, X )
  ls(pred)
  pred$class
  tailleft(pred$posterior)
  plot( as.timeSeries(pred$posterior) )
  plot( timeSeries(as.numeric(pred$class), rownames(pred$posterior)) )
  
  
  QDA2 <- qda( x = X,
              grouping = BBS,
              prior = c(0.2, 0.8) )
  pred2 <- predict(QDA2, X)
  
  
  tmp <- pred$posterior[ ,1] - pred2$posterior[ ,1]
  tmp <- as.timeSeries(tmp)
  plot(tmp)
  
  tmp <- as.numeric(pred$class) - as.numeric(pred2$class)
  tmp <- as.timeSeries(tmp)
  plot(tmp)
  
  
  
  
  
  QDA <- qda( x = X,
              CV = TRUE,
              grouping = BBS,
              prior = c(0.5, 0.5) )
  plot( as.timeSeries(QDA$posterior) )
  ls(QDA)
  
  
  
  dates <- intersect(rownames(X), rownames(BBS))
  X <- X[dates, ]
  BBS <- BBS[dates, ]
  FUN <- function(k) { apply( X[which(BBS == k), ], 2, mean) }
  lMu <- lapply( c(-1, 1), FUN = FUN )
  mu <- do.call( rbind, lMu )
  
  # Check
  apply( cbind( mu[1, ] - QDA$means[1, ],
                mu[2, ] - QDA$means[2, ] ), 2, range )
  
  
  FUN <- function(k) { cov(X[which(BBS == k), ]) }
  lCovmat <- lapply( c(-1, 1), FUN = FUN )
  # debugonce( mahalanobisDelta )
  md_delta <- apply( X_m, 
                     1, 
                     FUN = mahalanobisDelta,
                     center_1 = lMu[[1]] * 21, 
                     covmat_1 = lCovmat[[1]] * sqrt(21),
                     center_2 = lMu[[2]] * 21,
                     covmat_2 = lCovmat[[2]] * sqrt(21),
                     b_scl = FALSE,
                     b_det = TRUE,
                     prior_1 = 0.5 )
  md_delta2 <- apply( X_m, 
                     1, 
                     FUN = mahalanobisDelta,
                     center_1 = lMu[[1]] * 21, 
                     covmat_1 = lCovmat[[1]] * sqrt(21),
                     center_2 = lMu[[2]] * 21,
                     covmat_2 = lCovmat[[2]] * sqrt(21),
                     b_scl = FALSE,
                     b_det = FALSE,
                     prior_1 = 0.5 )
  md_delta_prior <- apply( X_m, 
                           1, 
                           FUN = mahalanobisDelta,
                           center_1 = lMu[[1]] * 21, 
                           covmat_1 = lCovmat[[1]] * sqrt(21),
                           center_2 = lMu[[2]] * 21,
                           covmat_2 = lCovmat[[2]] * sqrt(21),
                           b_scl = FALSE,
                           prior_1 = 0.3 )
  
  
  md_1 <- mahalanobisTS( X = X_m,
                         center = lMu[[1]],
                         covmat = lCovmat[[1]] )
  md_2 <- mahalanobisTS( X = X_m,
                         center = lMu[[2]],
                         covmat = lCovmat[[2]] )
  MD <- cbind(md_1, md_2)
  MD_scl <- scale(MD, FALSE, TRUE)
  md_delta_scl <- MD_scl[ ,1] - MD_scl[ ,2]
  md_delta_non_scl <- MD[ ,1] - MD[ ,2]
  
  MDDelta <- cbind(md_delta, 
                   md_delta2,
                   md_delta_prior,
                   md_delta_scl * 50,
                   md_delta_non_scl)
  plot( MDDelta, plot.type = "single" )
  
  
  t(QDA$scaling[ ,,1]) %*% lCovmat[[1]] %*% QDA$scaling[ ,,1]
  
  
  
  
  
  
  # Compare QDA and own implementation
  tmp <- cbind( qda = timeSeries(as.numeric(pred$class), 
                                 rownames(md_delta)),
                own = 1 / (1 + exp(-md_delta * 10^4)) + 0.9 )
  

  plot( tmp, ylim = c(0.8, 2.1), plot.type = "single")
  abline( v = time(BBS)[which(BBS == -1)], col = "red" )
  abline( v = time(BBS)[which(BBS == 1)], col = "green" )
  lines( tmp[ ,1], col = 4, lwd = 3 )
  lines( tmp[ ,2], col = 6, lwd = 2 )
  
  
  
  mu <- X * NA
  for ( i in 1:nrow(X) ) {
    
    x_eval <- as.numeric(X[i, ])
    md <- mahalanobisTS(X = X, center = x_eval)
    eps <- density(md)$bw
    p <- exp( -1 / (2 * eps) * md )
    p <- p / sum(p)
    mu[i, ] <- t(X) %*% p
      
  }
  
  md_delta_mu <- apply( mu, 1, FUN = mahalanobisDelta,
                        center_1 = lMu[[1]], 
                        covmat_1 = lCovmat[[1]],
                        center_2 = lMu[[2]],
                        covmat_2 = lCovmat[[2]],
                        b_scl = FALSE )
  
  plot(md_delta_mu)
  
  # Compare QDA and own implementation
  tmp <- cbind( qda = timeSeries(as.numeric(pred$class) - 1, 
                                 rownames(md_delta)),
                own = 1 / (1 + exp(-md_delta * 10^4)) - 0.1,
                own2 = 1 / (1 + exp(-md_delta_mu * 10^4)) - 0.2 )
  
  
  plot( tmp, ylim = c(-0.3, 1.3), plot.type = "single")
  abline( v = time(BBS)[which(BBS == -1)], col = "red" )
  abline( v = time(BBS)[which(BBS == 1)], col = "green" )
  lines( tmp[ ,1], col = 4, lwd = 3 )
  lines( tmp[ ,2], col = 5, lwd = 2 )
  lines( tmp[ ,3], col = 6, lwd = 2 )
  lines( log(cumulated(X_bm, "discrete")), lwd = 2, col = "orange" )
  
  
  
  # Signal testing
  
  sig <- cbind( qda = timeSeries(as.numeric(pred$class) - 1, 
                                 rownames(md_delta)),
                own = 1 / (1 + exp(-md_delta * 10^4)),
                own2 = 1 / (1 + exp(-md_delta_mu * 10^4)),
                bbs = (1 + BBS) / 2 )
  test <- signalTesting.byTrading( X = X_bm,
                                   sig = sig,
                                   n_lag = 1,
                                   tc = 0.004)
  X_tmp <- na.omit( cbind(X_bm, test) )
  
  plot( cumulated(X_tmp, "discrete"), plot.type = "single" )
  
  plot( cumulated(X_tmp[rownames(X_tmp) > "2009-01-01", ], "discrete"), plot.type = "single" )
  
  
  
  # QDA scales the difference in Mahalanobis distances 
  # by the log ratio of the determinants of thecovariances and multiplies
  # everyting by the log ratio of prior probabilities
  
  lDet <- lapply( lapply( lCovmat, FUN = function(x){ x * 10^0 } ), 
                  FUN = det )
  lDet
  log( lDet[[1]] / lDet[[1]] )
  lapply( lDet, FUN = log )
  
  
  
  Y1 <- scale(X, center = lMu[[1]], FALSE) / sqrt(sum(BBS == -1))
  qx <- qr(Y1)
  det_1 <- 2 * sum(log(abs(diag(qx$qr))))
  
  Y2 <- scale(X, center = lMu[[2]], FALSE) / sqrt(sum(BBS == 1))
  qx <- qr(Y2)
  det_2 <- 2 * sum(log(abs(diag(qx$qr))))
  
  log( det_1 / det_2 )
  
  QDA$ldet
  
  
  determinant( x = cov(X[BBS = -1, ]), logarithm = TRUE )
  log( det( cov(X[BBS = -1, ]) ) )
  
  
  
  
  
  # Idea: 
  # Use raw Mahalanobis distance to compute conditional distribution 
  # of class probability. Use that as prior
  
  md_all <- mahalanobisTS( X = X )
  FUN <- function( x_eval )
  {
    p <- weightsFun.kernel( data = list(X_train = md_all), 
                            x_eval = x_eval )
    mu <- t( (1 + BBS) / 2 ) %*% p
    return( mu )
  }
  prior <- apply( md_all, 1, FUN = FUN )
  
  
  md_delta_dynprior <- md_all * NA
  for ( i in 1:nrow(X) ) {
    tmp <- mahalanobisDelta(x = as.numeric(X[i, ]),
                            center_1 = lMu[[1]], 
                            covmat_1 = lCovmat[[1]],
                            center_2 = lMu[[2]],
                            covmat_2 = lCovmat[[2]],
                            b_scl = FALSE,
                            prior_1 = (1 - prior[i]) )
    md_delta_dynprior[i, ] <- tmp
  }
  
  MDDelta <- cbind(md_delta, 
                   md_delta_prior,
                   md_delta_scl * 50,
                   md_delta_dynprior)
  MDDelta[ MDDelta == -Inf ] <- NA
  plot( MDDelta, plot.type = "single" )
  abline(h = 0)
  
  
  ###
  qda_posterior <- md_all * NA
  qda_class <- md_all * NA
  for ( i in 1:nrow(X) ) {
    p <- as.numeric( prior[rownames(X)[i], ] )
    QDA <- qda( x = X,
                grouping = BBS,
                prior = c(1 - p, p) )
    pred <- predict( QDA, X[i, ] )
    qda_class[i, ] <- as.numeric( pred$class )
    qda_posterior[i, ] <- as.numeric( pred$posterior[2] )
  }
  
  plot( qda_posterior )
  plot( qda_class )
  
  
  QDA <- qda( x = X,
              grouping = BBS,
              prior = c(0.5, 0.5) )
  pred <- predict(QDA, X)
  qda_posterior_5050 <- timeSeries( pred$posterior[ ,2],
                                    rownames(pred$posterior) )
  qda_class_5050 <- timeSeries( as.numeric(pred$class),
                                rownames(pred$posterior) )
  QDACV <- qda( x = X,
                CV = TRUE,
                grouping = BBS,
                prior = c(0.5, 0.5) )
  qda_posterior_5050_cv <- timeSeries( QDACV$posterior[ ,2],
                                       rownames(QDACV$posterior) )
  qda_class_5050_cv <- timeSeries( as.numeric(QDACV$class),
                                   rownames(QDACV$posterior) )
  
  
  
  # sig <- cbind( c1 = qda_class - 1,
  #               c2 = qda_class_5050 - 1,
  #               c3 = qda_class_5050_cv - 1 )
  sig <- cbind( c1 = qda_posterior,
                c2 = qda_posterior_5050,
                c3 = qda_posterior_5050_cv )
  sig <- round( ema(sig, 0.02) )
  test <- signalTesting.byTrading( X = X_bm,
                                   sig = sig,
                                   n_lag = 1,
                                   tc = 0.004)
  X_tmp <- na.omit( cbind(X_bm, test) )
  plot( cumulated(X_tmp, "discrete"), plot.type = "single" )
  
  plot( sig )
  
  
  
  
  # --------------------------------------------------------------------------
  # COMPARE WITH NONPARAMETRIC NAIVE BAYES
  # --------------------------------------------------------------------------
  
  # SYNTHETIC DATA
  
  # Generate iidN
  
  n_sim <- 10^4
  d <- 2
  mu_1 <- 1
  sd_1 <- 1
  mu_2 <- -1
  sd_2 <- 1
  X1 <- matrix( rnorm(n_sim * d, mu_1, sd_1), ncol = d )
  X2 <- matrix( rnorm(n_sim * d, mu_2, sd_2), ncol = d )
  X <- rbind(X1, X2)
  states <- c( rep(1, n_sim), rep(-1, n_sim) )
  
  QDA <- qda( x = X,
              grouping = states,
              prior = c(0.5, 0.5) )
  x <- c(-2, -1)
  pred <- predict(QDA, x)
  ls(pred)
  pred
  
  pred <- predict(QDA, X)
  plot(pred$posterior[ ,1])
  
  
  # Naive Bayes
  x <- c(-2, -1)
  d1 <- dnorm( x, mu_1, sd_1 )
  d2 <- dnorm( x, mu_2, sd_2 )
  
  prod(d1) / ( prod(d1) + prod(d2))
  prod(d2) / ( prod(d1) + prod(d2))
  pred$posterior
  

  # Nonparametric naive Bayes
  ldens1 <- apply( X1, 2, density )
  ldens2 <- apply( X2, 2, density )
  
  FUN <- function( dens, a )
  {
    idx <- which(dens$x >= a)[1]
    prob <- (dens$y[idx] + dens$y[idx-1]) / 2
    return( prob )
  }
  npd1 <- c( FUN( ldens1[[1]], a = x[1] ),
             FUN( ldens1[[2]], a = x[2] ) )
  npd2 <- c( FUN( ldens2[[1]], a = x[1] ),
             FUN( ldens2[[2]], a = x[2] ) )
  prod(npd1) / ( prod(npd1) + prod(npd2))
  prod(npd2) / ( prod(npd1) + prod(npd2))
  pred$posterior
  
  
  
  
  # REAL DATA
  
  MD <- loadSensor( "turbulence" )
  BBT <- BBTurbulence$new()
  BBT$setCtrl( method = "base" )
  # BBT$updateData()
  # lData <- BBT$data
  lData <- MD$data
  BBT$data <- lData
  BBT$computeSignalInsample()
  states <- BBT$signal$insample[ ,"states"]
  X <- BBT$data$X
  mu_1 <- as.numeric( BBT$signal$insample_stats$mu_bull )
  sd_1 <- sqrt( diag( BBT$signal$insample_stats$sigma_bull ) )
  mu_2 <- as.numeric( BBT$signal$insample_stats$mu_bear )
  sd_2 <- sqrt( diag( BBT$signal$insample_stats$sigma_bear ) )
  
  today <- "2020-03-05"
  # today <- "2020-04-07"
  x <- X[today, ]
  
  barplot(x)
 
  
  
  prior <- c(0.5, 0.5)
  QDA <- qda( x = X,
              grouping = states,
              prior = prior )
  pred <- predict( QDA, x )
  pred
  BBT$signal$insample[today, ]
  
  
  require(mvtnorm)
  
  mvd1 <- pmvnorm( as.numeric(x), mean = mu_1, sigma = cov(X[which(states == 1), ]) )
  mvd2 <- pmvnorm( as.numeric(x), mean = mu_2, sigma = cov(X[which(states == -1), ]) )
  as.numeric( (mvd1 * prior[1]) / sum(mvd1 * prior[1] + mvd2 * prior[2]) )
  as.numeric( (mvd2 * prior[2]) / sum(mvd1 * prior[1] + mvd2 * prior[2]) )
  
  
  
  # Naive Bayes
  d1 <- d2 <- rep(NA, ncol(X))
  for ( i in 1:length(d1) ) {
    d1[i] <- dnorm( x[i], mu_1[i], sd_1[i] )
    d2[i] <- dnorm( x[i], mu_2[i], sd_2[i] )
  }

  prod(d1) / ( prod(d1) + prod(d2))
  prod(d2) / ( prod(d1) + prod(d2))
  pred$posterior
  
  # Nonparametric naive Bayes
  ldens1 <- apply( X[which(states == 1), ], 2, density )
  ldens2 <- apply( X[which(states == -1), ], 2, density )
  
  FUN <- function( dens, a )
  {
    idx <- which(dens$x >= a)[1]
    prob <- (dens$y[idx] + dens$y[idx-1]) / 2
    return( prob )
  }
  npd1 <- npd2 <- rep(NA, ncol(X))
  for ( i in 1:ncol(x) ) {
    npd1[i] <- FUN( ldens1[[i]], a = x[i] )
    npd2[i] <- FUN( ldens2[[i]], a = x[i] )
  }
  prod(npd1) / ( prod(npd1) + prod(npd2))
  prod(npd2) / ( prod(npd1) + prod(npd2))
  pred$posterior
  
  
  
  # ICA
  
  # Bear regime
  X_bear <- X[ which(states == -1), ]
  ICA_bear <- BSS::ica( X = X_bear )
  Y_bear <- getter( ICA_bear, "sources" )
  ldens_bear <- apply( Y_bear, 2, density )
  
  # Bull regime
  X_bull <- X[ which(states == 1), ]
  ICA_bull <- BSS::ica( X = X_bull )
  Y_bull <- getter( ICA_bull, "sources" )
  ldens_bull <- apply( Y_bull, 2, density )
  
  
  ldens1 <- apply( Y_bull, 2, density )
  ldens2 <- apply( Y_bear, 2, density )
  FUN <- function( dens, a )
  {
    idx <- which(dens$x >= a)[1]
    prob <- (dens$y[idx] + dens$y[idx-1]) / 2
    return( prob )
  }
  y_bear <- X[today, ] %*% getter(ICA_bear, "torsion")
  y_bull <- X[today, ] %*% getter(ICA_bull, "torsion")
  npd1 <- npd2 <- rep(NA, ncol(X))
  for ( i in 1:ncol(x) ) {
    npd1[i] <- FUN( ldens1[[i]], a = y_bull[i] )
    npd2[i] <- FUN( ldens2[[i]], a = y_bear[i] )
  }
  prod(npd1) / ( prod(npd1) + prod(npd2))
  prod(npd2) / ( prod(npd1) + prod(npd2))
  pred$posterior
  
  
  
  
  # Loop over all test dates (still in-sample)
  
  # --------------------------------------------------------------------------
  normalNaiveBayes <- function( X_test, mu_1, sd_1, mu_2, sd_2 )
  {
    
    d1 <- d2 <- X_test * NA
    post1 <- post2 <- X_test[ ,1] * NA
    for ( i in 1:nrow(X_test) ) {
      # for ( j in 1:ncol(X_test) ) {
      #   d1[i, j] <- dnorm( X_test[i, j], mu_1[j], sd_1[j] )
      #   d2[i, j] <- dnorm( X_test[i, j], mu_2[j], sd_2[j] )
      # }
      d1[i, ] <- dnorm( X_test[i, ], mu_1, sd_1 )
      d2[i, ] <- dnorm( X_test[i, ], mu_2, sd_2 )
      denom <- prod(d1[i, ]) + prod(d2[i, ])
      post1[i, ] <- prod(d1[i, ]) / denom
      post2[i, ] <- prod(d2[i, ]) / denom
    }
    post <- cbind( posterior1 = post1, posterior2 = post2 )
    if ( is.timeSeries(X_test) ) {
      post <- timeSeries( post, time(X_test) )
    }
    return( post )
  }
  
  # --------------------------------------------------------------------------
  npNaiveBayes <- function( X_train , X_test, states )
  {
    # Nonparametric naive Bayes
    ldens1 <- apply( X_train[which(states == 1), ], 2, density )
    ldens2 <- apply( X_train[which(states == -1), ], 2, density )
    
    densFUN <- function( a, dens )
    {
      idx <- which(dens$x >= a)[1]
      if ( !is.na(idx) ) {
        if ( idx == 1 ) {
          prob <- dens$y[1]
        } else if ( idx == length(dens$x) ) {
          prob <- dens$y[length(dens$y)]
        } else {
          prob <- (dens$y[idx] + dens$y[idx-1]) / 2
        }
      } else {
        prob <- dens$y[1]
      }
      return( prob )
    }
    npd1 <- npd2 <- X_test * NA
    for ( i in 1:nrow(X_test) ) {
      # tmp <- lapply( 1:ncol(X_test), 
      #                FUN = function(k) { densFUN( a = as.numeric(X_test[i, k]), 
      #                                             dens = ldens1[[k]] ) } )
      # npd1[i, ] <- unlist(tmp)
      # tmp <- lapply( 1:ncol(X_test), 
      #                FUN = function(k) { densFUN( a = as.numeric(X_test[i, k]), 
      #                                             dens = ldens2[[k]] ) } )
      # npd2[i, ] <- unlist(tmp)
      for ( j in 1:ncol(X_test) ) {
        # debugonce( densFUN )
        npd1[i, j] <- densFUN( a = as.numeric(X_test[i, j]), dens = ldens1[[j]] )
        npd2[i, j] <- densFUN( a = as.numeric(X_test[i, j]), dens = ldens2[[j]] )
      }
    }
    
    denom <- apply( npd1, 1, prod ) + apply( npd2, 1, prod )
    post1 <- apply( npd1, 1, prod ) / denom
    post2 <- apply( npd2, 1, prod ) / denom
    post <- cbind( posterior1 = post1, posterior2 = post2 )
    if ( is.timeSeries(X_test) ) {
      post <- timeSeries( post, time(X_test) )
    }
    return( post )
  }
  
  
  
  # --------------------------------------------------------------------------
  icnpNaiveBayes <- function( Y1_train, Y2_train, ICA1, ICA2, X_test )
  {
    # ICA Nonparametric naive Bayes
    
    ldens1 <- apply( Y1_train, 2, density )
    ldens2 <- apply( Y2_train, 2, density )
    Y1_test <- X_test %*% getter( ICA1, "torsion" )
    Y2_test <- X_test %*% getter( ICA1, "torsion" )
    
    densFUN <- function( a, dens )
    {
      idx <- which(dens$x >= a)[1]
      if ( !is.na(idx) ) {
        if ( idx == 1 ) {
          prob <- dens$y[1]
        } else if ( idx == length(dens$x) ) {
          prob <- dens$y[length(dens$y)]
        } else {
          prob <- (dens$y[idx] + dens$y[idx-1]) / 2
        }
      } else {
        prob <- dens$y[1]
      }
      return( prob )
    }
    npd1 <- npd2 <- X_test * NA
    for ( i in 1:nrow(X_test) ) {
      # tmp <- lapply( 1:ncol(X_test), 
      #                FUN = function(k) { densFUN( a = as.numeric(X_test[i, k]), 
      #                                             dens = ldens1[[k]] ) } )
      # npd1[i, ] <- unlist(tmp)
      # tmp <- lapply( 1:ncol(X_test), 
      #                FUN = function(k) { densFUN( a = as.numeric(X_test[i, k]), 
      #                                             dens = ldens2[[k]] ) } )
      # npd2[i, ] <- unlist(tmp)
      for ( j in 1:ncol(X_test) ) {
        # debugonce( densFUN )
        npd1[i, j] <- densFUN( a = as.numeric(Y1_test[i, j]), dens = ldens1[[j]] )
        npd2[i, j] <- densFUN( a = as.numeric(Y2_test[i, j]), dens = ldens2[[j]] )
      }
    }
    
    denom <- apply( npd1, 1, prod ) + apply( npd2, 1, prod )
    post1 <- apply( npd1, 1, prod ) / denom
    post2 <- apply( npd2, 1, prod ) / denom
    post <- cbind( posterior1 = post1, posterior2 = post2 )
    if ( is.timeSeries(X_test) ) {
      post <- timeSeries( post, time(X_test) )
    }
    return( post )
  }
  
  
  
  post_nnb <- normalNaiveBayes( X_test = X, 
                                mu_1 = mu_1, 
                                sd_1 = sd_1, 
                                mu_2 = mu_2, 
                                sd_2 = sd_2 )
  # debugonce( npNaiveBayes )
  post_npnb <- npNaiveBayes( X_train = X,
                             X_test = X,
                             states = states )
  # debugonce( icnpNaiveBayes )
  ICA1 <- BSS::ica( X[which(states == 1), ] )
  ICA2 <- BSS::ica( X[which(states == -1), ] )
  Y1_train <- getter(ICA1, "sources")
  Y2_train <- getter(ICA2, "sources")
  post_icnpnb <- icnpNaiveBayes( Y1_train = Y1_train,
                                 Y2_train = Y2_train,
                                 X_test = X,
                                 ICA1 = ICA1,
                                 ICA2 = ICA2 )
  # Alternatively:
  ICA1 <- ICA2 <- BSS::ica( X )
  Y1_train <- getter( ICA1, "sources" )[which(states == 1), ]
  Y2_train <- getter( ICA1, "sources" )[which(states == -1), ]
  post_icnpnb2 <- icnpNaiveBayes( Y1_train = Y1_train,
                                  Y2_train = Y2_train,
                                  X_test = X,
                                  ICA1 = ICA1,
                                  ICA2 = ICA2 )
  pred <- predict( QDA, X )
  post_qda <- timeSeries( pred$posterior, time(X) )
  
  
  sig_nnb <- apply(post_nnb, 1, function(x) { ifelse(x[1] > x[2], 1, 0)})
  sig_npnb <- apply(post_npnb, 1, function(x) { ifelse(x[1] > x[2], 1, 0)})
  sig_icnpnb <- apply(post_icnpnb, 1, function(x) { ifelse(x[1] > x[2], 1, 0)})
  sig_icnpnb2 <- apply(post_icnpnb2, 1, function(x) { ifelse(x[1] > x[2], 1, 0)})
  sig_qda <- apply(post_qda, 1, function(x) { ifelse(x[1] < x[2], 1, 0)})
  sig_mmd_base_is <- (1 + sign(BBT$signal$insample[ ,"delta"])) / 2
  sig_mmd_scl2_oos <- (1 + sign(loadSensor("bbturbulence")$signal$scl2[ ,"md_delta"])) / 2
  sig_states_is = (1 + BBT$signal$insample[ ,"states"]) / 2
  sig <- cbind( nnb = sig_nnb, 
                npnb = sig_npnb,
                icnpnb = sig_icnpnb, 
                icnpnb2 = sig_icnpnb2, 
                qda = sig_qda, 
                mmd_base_is = sig_mmd_base_is,
                mmd_scl2_oos = sig_mmd_scl2_oos,
                states_is = sig_states_is )
  colors <- c(1, fBasics::rainbowPalette(n = ncol(sig) ) )
  test <- signalTesting.byTrading( X = BBT$data$X_bm,
                                   sig = sig,
                                   n_lag = 1,
                                   tc = 0.00 )
  colnames(test) <- colnames(sig)
  X_tmp <- na.omit(cbind(BBT$data$X_bm, test))
  plot( as.simTS(X_tmp), col = colors )
  
  
  
  
  
  # Signal based nonparametric naive Bayes - Prior on Prior
  
  signal <- MD$signal$base[ ,"md_ema0.157"]
  # signal <- MD$signal$base[ ,"signal"]
  DS <- dirichletSampling( Y_train = (1 + states) / 2,
                           X_train = signal,
                           weights_fun = "l1",
                           # weights_fun = "kernel",
                           # scl_by_entropy = TRUE )
                           sclfct = 1 )
  lprior <- lapply( DS, density )
  plot.ldensity( lprior )
  
  today <- "2020-04-07"
  
  ldens1 <- apply( X[which(states == 1), ], 2, density )
  ldens2 <- apply( X[which(states == -1), ], 2, density )
  FUN <- function( dens, a )
  {
    idx <- which(dens$x >= a)[1]
    prob <- (dens$y[idx] + dens$y[idx-1]) / 2
    return( prob )
  }
  npd1 <- npd2 <- rep(NA, ncol(X))
  for ( j in 1:ncol(X_test) ) {
    npd1[j] <- FUN( ldens1[[j]], a = as.numeric(X_test[today, j]) )
    npd2[j] <- FUN( ldens2[[j]], a = as.numeric(X_test[today, j]) )
  }
  lPost <- list()
  lPrior <- list()
  for ( k in 1:length(DS) ) {
    post1 <- post2 <- NULL
    prior <- cbind( DS[[k]], 1 - DS[[k]] )
    lPrior[[k]] <- prior
    for ( i in 1:nrow(prior) ) {
      denom <-  prod(npd1) * prior[i, 1] + prod(npd2) * prior[i, 2]
      post1 <- c(post1, prod(npd1) * prior[i, 1] / denom)
      post2 <- c(post2, prod(npd2) * prior[i, 2] / denom)
    }
    lPost[[k]] <- cbind(post1, post2)
  }

  plot( density(lPost[[1]][ ,1]) )
  pred$posterior
  
  Post <- t( do.call( cbind, lapply( lPost, FUN = function(x) { apply(x, 2, mean) } ) ) )
  rownames(Post) <- names(DS)
  rownames(Post) <- round( unlist( lapply( lPrior, FUN = function(x) { mean(x[ ,1]) } ) ), 3 ) 
  Post
  
  barplot( t(Post), beside = TRUE, col = c(3, 2) )
  boxplot( do.call( cbind, lPost ), col = c(3, 2) )
  
  
  
  
  
  # --------------------------------------------------------------------------
  # BBT MINIMUM TORSION (MT) GARCH OR STOCHASTIC VOLATILITY (SV)
  # --------------------------------------------------------------------------
  
  # Use minimum torsion transformation because of non-convergence issues in PC's
  
  # Data
  BBT <- loadSensor( sensor_name = "bbturbulence_scl2_dm", b_new = TRUE )  
  BBT$computeSignalInsample()
  X <- BBT$data$X[isWeekday(time(BBT$data$X)), ]
  wmat <- BBT$data$wmat
  states <- BBT$signal$insample[ ,"states"]
  dates <- intersect(rownames(wmat), rownames(X))
  X <- X[dates, ]
  wmat <- wmat[dates, ]
  states <- states[dates, ]
  
  # Standard MD's
  X_bear <- X[which(states == -1), ]
  X_bull <- X[which(states == 1), ]
  covmat_bear <- cov(X_bear)
  covmat_bull <- cov(X_bull)
  md_all <- mahalanobisTS(X = X, scl = FALSE)
  md_bear <- mahalanobisTS(X = X, center = apply(X_bear, 2, mean), covmat = covmat_bear, scl = FALSE)
  md_bull <- mahalanobisTS(X = X, center = apply(X_bull, 2, mean), covmat = covmat_bull, scl = FALSE)
  
  # Minimum torsion series
  S_all <- getter( mtca(X), "sources" )
  tmat_bear <- torsion(covmat_bear, method = "mt")
  tmat_bull <- torsion(covmat_bull, method = "mt")
  S_bear <- timeSeries( scale(X, apply(X_bear, 2, mean), FALSE) %*% tmat_bear, time(X) )
  S_bull <- timeSeries( scale(X, apply(X_bull, 2, mean), FALSE) %*% tmat_bull, time(X) )
 
  
  # tmat_all <- torsion(cov(X), method = "mt")
  # S_all <- timeSeries( scale(X, TRUE, FALSE) %*% tmat_all, time(X) )
  # a <- apply( tmat_all, 2, function(x) { sqrt( x^2 / sum(x^2) ) } )
  # S_all_a <- timeSeries( scale(X, TRUE, FALSE) %*% a, time(X) )
  # cor(S_all)
  # cor(S_all_a)
  
  
  # Distances in MT space
  S_all_scl <- scale( S_all^2, FALSE, apply(S_all, 2, var) )
  S_bear_scl <- scale( S_bear^2, FALSE, apply(S_bear, 2, var) )
  S_bull_scl <- scale( S_bull^2, FALSE, apply(S_bull, 2, var) )
  md_all_s <- apply( S_all_scl, 1, sum )
  md_bear_s <- apply( S_bear_scl, 1, sum )
  md_bull_s <- apply( S_bull_scl, 1, sum )
  
  
  # GARCH
  fit_all <- garch(S_all)
  cvol_all <- getCondVar(fit_all)
  # md_cvol2 <- mahalanobisTS( cvol^0.5,
  #                            center = apply(S, 2, mean),
  #                            covmat = cov(S),
  #                            scl = FALSE )
  md_cvol <- apply( scale( cvol_all, FALSE, apply(S_all, 2, var) ), 1, sum )
  
  fit_bear <- garch( S_bear, garchCtrl(armapq = c(0, 0, 1, 1)) )
  cvol_bear <- getCondVar(fit_bear)
  idx <- apply( cvol_bear, 2, function(x) { all(x == 0) } )
  if ( length(idx) > 0 ) {
    cvol_bear[ ,idx] <- S_bear[ ,idx]^2
  }
  fit_bull <- garch(S_bull, garchCtrl(armapq = c(0, 0, 1, 1)))
  cvol_bull <- getCondVar(fit_bull)
  idx <- apply( cvol_bull, 2, function(x) { all(x == 0) } )
  if ( length(idx) > 0 ) {
    cvol_bull[ ,idx] <- S_bull[ ,idx]^2
  }
  fit_delta <- garch(S_bear - S_bull, garchCtrl(armapq = c(0, 0, 1, 1)))
  cvol_delta <- getCondVar(fit_delta)
  
  
  md_all_s_cvol <- apply( scale(cvol_all, FALSE, apply(S_all, 2, var)), 1, sum)
  md_bear_s_cvol <- apply( scale(cvol_bear, FALSE, apply(S_bear, 2, var)), 1, sum)
  md_bull_s_cvol <- apply( scale(cvol_bull, FALSE, apply(S_bull, 2, var)), 1, sum)
  md_delta_s_cvol <- apply( scale(cvol_delta, FALSE, apply(S_bear - S_bull, 2, var)), 1, sum)
  
  
  
  # EWMA
  S_all_sq_ewma <- ema( S_all^2, alpha = 0.1 )
  S_bear_sq_ewma <- ema( S_bear^2, alpha = 0.1 )
  S_bull_sq_ewma <- ema( S_bull^2, alpha = 0.1 )
  md_all_s_ewma <- apply( scale(S_all_sq_ewma, FALSE, apply(S_all, 2, var)), 1, sum)
  md_bear_s_ewma <- apply( scale(S_bear_sq_ewma, FALSE, apply(S_bear, 2, var)), 1, sum)
  md_bull_s_ewma <- apply( scale(S_bull_sq_ewma, FALSE, apply(S_bull, 2, var)), 1, sum)
  
  
  
  # Stochastic volatility
  require(stochvol)
  sv_all <- sv_bear <- sv_bull <- S_all * NA
  sv_delta <- sv_all
  for ( j in 1:ncol(S_all) ) {
    draws <- NULL
    draws <- svsample( S_all[ ,j], draws = 1000, burnin = 100 )
    sv_all[ ,j] <- apply( exp(draws$latent), 2, mean )
    draws <- NULL
    draws <- svsample( S_bear[ ,j], draws = 1000, burnin = 100 )
    sv_bear[ ,j] <- apply( exp(draws$latent), 2, mean )
    draws <- NULL
    draws <- svsample( S_bull[ ,j], draws = 1000, burnin = 100 )
    sv_bull[ ,j] <- apply( exp(draws$latent), 2, mean )
    draws <- NULL
    draws <- svsample( S_bear[ ,j] - S_bull[ ,j], draws = 1000, burnin = 100 )
    sv_delta[ ,j] <- apply( exp(draws$latent), 2, mean )
  }
  
  md_all_s_sv <- apply( scale(sv_all, FALSE, apply(S_all, 2, var)), 1, sum )
  md_bear_s_sv <- apply( scale(sv_bear, FALSE, apply(S_bear, 2, var)), 1, sum )
  md_bull_s_sv <- apply( scale(sv_bull, FALSE, apply(S_bull, 2, var)), 1, sum )
  md_delta_s_sv <- apply( scale(sv_delta, FALSE, apply( S_bear - S_bull, 2, var)), 1, sum )
  
  
  
  # Deltas
  # MD_raw <- cbind( bear = md_bear, bull = md_bull, all = md_all )
  # MD_raw_scl <- scale( MD_raw, FALSE, TRUE )
  # md_raw_delta <- MD_raw_scl[ ,"bear"] - MD_raw_scl[ ,"bull"]
  # md_raw_delta <- BBT$signal$scl2[ ,"signal"]
  md_raw_delta <- BBT$signal$insample[ ,"delta"]
  
  
  MD_s <- cbind( bear = md_bear_s, bull = md_bull_s, all = md_all_s )
  MD_s_scl <- scale( MD_s, FALSE, TRUE )
  md_s_delta <- MD_s_scl[ ,"bear"] - MD_s_scl[ ,"bull"]

  # GARCH
  MD_cvol <- cbind( bear = md_bear_s_cvol, bull = md_bull_s_cvol, all = md_all_s_cvol )
  MD_cvol_scl <- scale( MD_cvol, FALSE, TRUE )
  md_cvol_delta <- MD_cvol_scl[ ,"bear"] - MD_cvol_scl[ ,"bull"]

  # EWMA
  MD_ewma <- cbind( bear = md_bear_s_ewma, bull = md_bull_s_ewma, all = md_all_s_ewma )
  MD_ewma_scl <- scale( MD_ewma, FALSE, TRUE )
  md_ewma_delta <- MD_ewma_scl[ ,"bear"] - MD_ewma_scl[ ,"bull"]
  
  # Stochvol
  MD_sv <- cbind( bear = md_bear_s_sv, bull = md_bull_s_sv, all = md_all_s_sv )
  MD_sv_scl <- scale( MD_sv, FALSE, TRUE )
  md_sv_delta <- MD_sv_scl[ ,"bear"] - MD_sv_scl[ ,"bull"]

    
  # # Synthetic series
  # S_bear_synt <- cvol_bear^0.5 * (S_bear * 0 + rnorm(length(cvol_bear)))
  # S_bear_synt <- scale( S_bear_synt, TRUE, 1 / eigen(covmat_bear)$values )
  # S_bull_synt <- cvol_bull^0.5 * (S_bull * 0 + rnorm(length(cvol_bull)))
  # S_bull_synt <- scale( S_bull_synt, TRUE, 1 / eigen(covmat_bull)$values )
  # md_bear_mt_synt <- apply( scale(S_bear_synt^2, FALSE, apply(S_bear_synt, 2, var)), 1, sum)
  # md_bull_mt_synt <- apply( scale(S_bull_synt^2, FALSE, apply(S_bull_synt, 2, var)), 1, sum)
  # # md_bear_mt_synt <- apply( scale(S_bear_synt^2, FALSE, eigen(covmat_bear)$values), 1, sum)
  # # md_bull_mt_synt <- apply( scale(S_bull_synt^2, FALSE, eigen(covmat_bull)$values), 1, sum)
  # 
  # plot( cbind(md_bear_mt_synt, md_bull_mt_synt) )
  # plot( cbind(md_bear_mt_synt, md_bull_mt_synt), plot.type = "single" )
  # 
  # MD <- cbind(bear = md_bear_mt_synt,
  #             bull = md_bull_mt_synt)
  # MD_scl <- scale(MD, FALSE, TRUE)
  # md_delta <- MD_scl[ ,"bear"] - MD_scl[ ,"bull"]
  # 
  # plot(ema(md_delta, 0.05))
  # abline(h = 0, col = "grey")
  
  
  
  
  # Signal testsing
  # MD_delta <- cbind( raw = md_raw_delta,
  #                    cvol = md_cvol_delta,
  #                    ewma = md_ewma_delta,
  #                    sv = md_sv_delta )
  # colnames(MD_delta) <- c("raw", "cvol", "ewma", "sv")
  MD_delta <- cbind( raw = md_raw_delta,
                     sv = md_sv_delta )
  colnames(MD_delta) <- c("raw", "sv")
  tmp <- apply( MD_delta, 2, function(x) { (1 + sign(x)) / 2 } )
  # sig <- timeSeries( t(tmp), time(MD_delta) )
  sig <- tmp
  X_bm <- BBT$data$X_bm
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
  
  
  require(visolz)
  HC <- hc_lineChart( )
  
  
  
  # Using PCA
  
  # Minimum torsion series
  PC_all <- getter( pca(X), "sources" )
  tmat_bear <- torsion(covmat_bear, method = "pc")
  tmat_bull <- torsion(covmat_bull, method = "pc")
  PC_bear <- timeSeries( scale(X, apply(X_bear, 2, mean), FALSE) %*% tmat_bear, time(X) )
  PC_bull <- timeSeries( scale(X, apply(X_bull, 2, mean), FALSE) %*% tmat_bull, time(X) )
  
  # Stochastic volatility
  require(stochvol)
  svpc_all <- svpc_bear <- svpc_bull <- PC_all * NA
  for ( j in 1:ncol(PC_all) ) {
    draws <- NULL
    draws <- svsample( PC_all[ ,j], draws = 1000, burnin = 100 )
    svpc_all[ ,j] <- apply( exp(draws$latent), 2, mean )
    draws <- NULL
    draws <- svsample( PC_bear[ ,j], draws = 1000, burnin = 100 )
    svpc_bear[ ,j] <- apply( exp(draws$latent), 2, mean )
    draws <- NULL
    draws <- svsample( PC_bull[ ,j], draws = 1000, burnin = 100 )
    svpc_bull[ ,j] <- apply( exp(draws$latent), 2, mean )
  }
  
  md_all_pc_sv <- apply( scale(sv_all, FALSE, apply(PC_all, 2, var)), 1, sum )
  md_bear_pc_sv <- apply( scale(sv_bear, FALSE, apply(PC_bear, 2, var)), 1, sum )
  md_bull_pc_sv <- apply( scale(sv_bull, FALSE, apply(PC_bull, 2, var)), 1, sum )
  MD_svpc <- cbind( bear = md_bear_pc_sv, bull = md_bull_pc_sv, all = md_all_pc_sv )
  MD_svpc_scl <- scale( MD_svpc, FALSE, TRUE )
  md_svpc_delta <- MD_svpc_scl[ ,"bear"] - MD_svpc_scl[ ,"bull"]
  
  
  MD_delta <- cbind( raw = md_raw_delta,
                     sv = md_sv_delta,
                     svpc = md_svpc_delta )
  colnames(MD_delta) <- c("raw", "sv_mt", "sv_pc")
  tmp <- apply( MD_delta, 2, function(x) { (1 + sign(x)) / 2 } )
  # sig <- timeSeries( t(tmp), time(MD_delta) )
  sig <- tmp
  X_bm <- BBT$data$X_bm
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
  
  
  
  
  
  ##########################
  
  # Example from MASS package
  
  tr <- sample(1:50, 25)
  train <- rbind(iris3[tr,,1], iris3[tr,,2], iris3[tr,,3])
  test <- rbind(iris3[-tr,,1], iris3[-tr,,2], iris3[-tr,,3])
  cl <- factor(c(rep("s",25), rep("c",25), rep("v",25)))
  z <- qda(train, cl)
  predict(z,test)$class
  
  
  
  
  
  
  
