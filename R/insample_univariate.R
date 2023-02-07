  
  
  ############################################################################
  ### GOOD AND BAD TURBULENCE - UNIVARIATE
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     22.01.2021
  # First version:    10.10.2020
  # --------------------------------------------------------------------------

  
  # --------------------------------------------------------------------------
  # REQUIRE
  # --------------------------------------------------------------------------
  
  require(DAARC)
  wd <- "R:/Asset_Management/Research_Projects/Multi_Asset/DAA/DAARC/Sensor_Testing/Univariate/Good_And_Bad_Turbulence/"

  
  

  # --------------------------------------------------------------------------
  # UNIVARIATE
  # --------------------------------------------------------------------------
  
  # Return series
  X_bm <- GPO::getMSCIData( universe = "bm",
                            frqncy = "d",
                            ccy = "Local" )
  X_bm <- X_bm[isWeekday(time(X_bm)), ]
  
  # Bear Bull states
  BBS <- BBSRC$new()
  BBS$setCtrl( minphase = 21,
               mincycle = 21 * 3,
               k_peak = 10,
               k_trough = 10,
               l_peak = 10,
               l_trough = 10 )
  BBS$data <- list( X_level = cumulated(X_bm, "discrete") )
  BBS$run()
  BBS$plotStates()
  
  bbstates <- BBS$output$states
  
  # Conditional statistics
  dates_bear <- rownames(bbstates)[ which(bbstates == -1) ]
  dates_bull <- rownames(bbstates)[ which(bbstates == 1) ]
  mu_bear <- meanGeo( X_bm[dates_bear, ] )
  sigma_bear <- sd( X_bm[dates_bear, ] )
  mu_bull <- meanGeo( X_bm[dates_bull, ] )
  sigma_bull <- sd( X_bm[dates_bull, ] )
  
  mu_bear; mu_bull
  sigma_bear; sigma_bull
  
    
  d_bear <- dnorm( x = X_bm, mean = mu_bear, sd = sigma_bear )
  d_bull <- dnorm( x = X_bm, mean = mu_bull, sd = sigma_bull )
  
  plot( cbind(d_bear, d_bull) )  
  plot( log( cbind(d_bear, d_bull) ) )  
  
  
  plot( ema( d_bear - d_bull, 0.1 ) )
  abline( h = 0 )
  
  
  
  fun <- function(x, mu, sigma) { 1 / sqrt(2 * pi * sigma^2) * exp( -(x - mu)^2 / (2 * sigma^2) ) }
  fun(x = 0, mu = 0, sigma = 0.05)
  dnorm( x = 0, mean = 0, sd = 0.05)  

  
  
  fun2 <- function(x, sigma) { 1 / sqrt(2 * pi * sigma^2) * exp( -x / (2 * sigma^2) ) }
  X_bear <- scale( X_bm, mu_bear, FALSE )
  d_bear_ema <- fun2(x = ema( X_bear^2, 1 ), sigma = sigma_bear )
  X_bull <- scale( X_bm, mu_bull, FALSE )
  d_bull_ema <- fun2(x = ema( X_bull^2, 1 ), sigma = sigma_bull )
  
  
  
  logLikFun <- function(x, sigma) { -log(sigma) - 1/2 * log(2 * pi) - 1/2 * (x / sigma)^2 }
  
  
  discrFun <- function(x, mu_1, sigma_1, mu_2, sigma_2)
  {
    th <- log( sigma_1 / sigma_2 )
    delta <- (x - mu_2)^2 / (2 * sigma_2^2) - (x - mu_1)^2 / (2 * sigma_1^2)
    ans <- c(delta, th)
    return( ans )
  }
  
  
  ll_bear <- logLikFun(x = (X_bm - mu_bear), sigma = sigma_bear )
  ll_bull <- logLikFun(x = (X_bm - mu_bull), sigma = sigma_bull )
  plot( cbind(ll_bear, ll_bull) )
  plot( ll_bear - ll_bull )
  
  
  
  
  tmp <- apply( X_bm, 
                1, 
                FUN = discrFun, 
                mu_1 = mu_bear, 
                sigma_1 = sigma_bear,
                mu_2 = mu_bull, 
                sigma_2 = sigma_bull )
  delta <- timeSeries( t(tmp), time(X_bm) )
  sig <- (1 + sign( delta[ ,2] - delta[ ,1] )) / 2
  
  plot( delta, plot.type = "single" )
  
  
  
  
  # EMA on squared returns
  th <- log( sigma_bear / sigma_bull )
  X_bear <- X_bm - mu_bear
  X_bull <- X_bm - mu_bull
  delta <- ema(X_bull^2, 0.5) / (2 * sigma_bull^2) - ema(X_bear^2, 0.5) / (2 * sigma_bear^2)
  
  X_bear <- scale( X_bm, mu_bear, sigma_bear )
  X_bull <- scale( X_bm, mu_bull, sigma_bull )
  delta <- ema( X_bull^2, 0.5 ) - ema( X_bear^2, 0.5 )
  sig <- (1 + sign( delta < th )) / 2
  
  
  # GARCH
  cvar_bear <- getCondVar( garch(X_bear) )
  cvar_bull <- getCondVar( garch(X_bull) )
  # delta <- cvar_bull / (2 * sigma_bull^2) - cvar_bear / (2 * sigma_bear^2)
  
  sig <- (1 + sign( delta < th )) / 2
  
  plot(delta)
  abline(h = th, col = 2)
  
  plot( cbind( cvar_bull / (2 * sigma_bull^2), cvar_bear / (2 * sigma_bear^2), delta) )
  plot( cbind( (X_bm - mu_bull)^2 / (2 * sigma_bull^2),
               (X_bm - mu_bear)^2 / (2 * sigma_bear^2) ), plot.type = "single" )
  plot( cbind( (X_bm - mu_bull)^2,
               (X_bm - mu_bear)^2 ), plot.type = "single" )
  
  
  
  
  TD <- trainingData( Y_train = X_bm,
                      X_train = sig )
  n_lag <- 2
  penalty <- 0
  tc <- 0.004
  
  # debugonce(signalTesting.byTrading)
  test <- signalTesting.byTrading( X = TD$Y_train,
                                   sig = TD$X_train,
                                   n_lag = n_lag,
                                   tc = tc,
                                   penalty = penalty )
  colnames(test) <- colnames(TD$X_train)
  test_nc <- signalTesting.byTrading( X = TD$Y_train,
                                      sig = TD$X_train,
                                      n_lag = n_lag,
                                      tc = 0,
                                      penalty = penalty )
  colnames(test_nc) <- colnames(TD$X_train)
  
  X_tmp <- na.omit( cbind( bm = TD$Y_train, test ) )
  X_tmp_nc <- na.omit( cbind( bm = TD$Y_train, test_nc ) )
  
  lStats <- descStats( X_tmp )
  lStats_nc <- descStats( X_tmp_nc )
  
  stats_fields <- c("cumret", "sds", "maxDD")
  t(lStats_nc$stats[stats_fields, ])    
  t(lStats$stats[stats_fields, ]) 
  
  plot( as.simTS(X_tmp) ) 
  plot( as.simTS(X_tmp_nc) ) 
  
  
  
  
  
  # --------------------------------------------------------------------------
  # NONPARAMETRIC FILTERING
  # --------------------------------------------------------------------------
  
  # Return series
  X_bm <- GPO::getMSCIData( universe = "bm",
                            frqncy = "d",
                            ccy = "Local" )
  X_bm <- X_bm[isWeekday(time(X_bm)), ]
  
  # Fuzzy Bear Bull states
  BBS <- BBSRC$new()
  BBS$setCtrl( minphase = 21,
               mincycle = 21 * 3,
               k_peak = 10,
               k_trough = 10,
               l_peak = 10,
               l_trough = 10,
               language = "C" )
  BBS$data <- list( X_level = cumulated(X_bm, "discrete") )
  # BBS$runRobust( minphase_vec = floor(21 * seq(from = 0.1, to = 4, length.out = 10)) )
  BBS$runRobust()
  BBS$plotStates()
  
  states_fuzzy <- (1 + BBS$output$states_fuzzy) / 2
  plot(states_fuzzy)
  
  
  today <- "2020-10-30"
  X_eval <- states_fuzzy[today, ]
  
  # Likelihood (nonparametric)
  LH <- dirichletSampling( Y_train = X_bm,
                           X_train = states_fuzzy,
                           # X_eval = X_eval,
                           weights_fun = "l1",
                           scl_by_entropy = FALSE,
                           sclfct = 1 )
  ldens_lh <- lapply( LH, FUN = density )
  plot.ldensity( ldens_lh )
  
  mu <- setNames( unlist( lapply( LH, FUN = mean) ), names(LH) )
  barplot( mu ) 
  
  
  
  # Prior (nonparametric)
  Prior <- dirichletSampling( Y_train = states_fuzzy,
                              X_train = states_fuzzy,
                              # X_eval = X_eval,
                              weights_fun = "l1",
                              scl_by_entropy = FALSE,
                              sclfct = 1 )
  ldens_prior <- lapply( Prior, FUN = density )
  plot.ldensity( ldens_prior )
  

  
  
 
    
  # Marginal (nonparametric)
  Marginal <- dirichletSampling( Y_train = X_bm,
                                 X_train = X_bm,
                                 weights_mat = X_bm[ ,1] * 0 + 1 / nrow(X_bm),
                                 scl_by_entropy = FALSE,
                                 sclfct = NULL )
  ldens_marginal <- lapply( Marginal, FUN = densFUN )
  plot.ldensity( ldens_marginal )
  lines( density(X_bm) )
  
  
  
  
  
  FUN = function(i) { ldens_lh[[i]]$y * ldens_prior[[i]]$y }
  # FUN = function(i) { LH[[i]] * Prior[[i]] / Marginal[[1]] }
  Post <- lapply( 1:length(LH), FUN = FUN )
  names(Post) <- names(LH)
  
  ldens_post <- ldens_lh
  for ( i in 1:length(ldens_post) ) {
    ldens_post[[i]]$y <- ldens_lh[[i]]$y * ldens_prior[[i]]$y
  }
  
  plot.ldensity( ldens_post )
  
  
  
  mu <- setNames( unlist( lapply( Post, FUN = mean) ), names(Post) )
  barplot( mu ) 

  
  plot( x = LH[[1]], y = Prior[[1]] )
  points( x = LH[[1]][1], y = Prior[[1]][1], col = 2, pch = 19, cex = 2 )

  head( cbind( LH[[1]], Prior[[1]] ) )
  
  
  
  idx <- 6
  lh <- LH[[idx]]
  prior <- Prior[[idx]]
  dens_lh <- density(lh)
  dens_prior <- density(prior)

  d_lh <- lh * NA
  d_prior <- d_lh
  for ( i in seq(along = lh) ) {
    idx_leq <- which( dens_lh$x <= lh[i] )
    idx_geq <- which( dens_lh$x >= lh[i] )
    d_lh[i] <- (dens_lh$y[ idx_leq[length(idx_leq)]] + dens_lh$y[ idx_geq[1] ]) / 2
    idx_leq <- which( dens_prior$x <= prior[i] )
    idx_geq <- which( dens_prior$x >= prior[i] )
    d_prior[i] <- (dens_prior$y[ idx_leq[length(idx_leq)]] + dens_prior$y[ idx_geq[1] ]) / 2
  }

  dens_post <- dens_prior
  ordering <- order(prior)
  dens_post$y <- (d_lh * d_prior)[ordering]
  dens_post$x <- prior[ordering]
  
  plot( dens_post )
  lines( dens_prior, col = 2 )
  
  
  
  
  
  # --------------------------------------------------------------------------
  # TOY EXAMPLE: POSTERIOR MEAN FOR NORMAL-NORMAL MODEL WITH ONE OBSERVATION
  # --------------------------------------------------------------------------
  
  # X|\mu \sim N(\mu, \sigma^2); \sigma^2 known.
  # \mu \sim N(\mu_0, \sigma_0^2)
  
  n <- 10^4
  mu_0 <- -0.05
  sigma_0 <- 0.1
  sigma <- 2
  prior <- rnorm( n, mu_0, sigma_0 )
  sigma_1 <- 1 / (sigma^(-2) + sigma_0^(-2))
  # or:
  # sigma_1 <- (sigma^2 * sigma_0^2) / (sigma^2 + sigma_0^2)
  x <- 0.2
  mu_1 <- (mu_0 * sigma^2 + x * sigma_0^2) / (sigma^2 + sigma_0^2)
  # or:
  # mu_1 <- (mu_0 * sigma_0^(-2) + x * sigma^(-2)) / (sigma^(-2) + sigma_0^(-2))
  # mu_1 <- sigma_1^2 * (mu_0 * sigma_0^(-2) + x * sigma^(-2))
  post <- rnorm( n, mu_1, sigma_1 )
    
  range(prior)
  range(post)  

  
  ldens <- apply( cbind(prior, post), 2, density )
  plot.ldensity( ldens ) 
  abline( v = x )
  
  
  
    
  
  