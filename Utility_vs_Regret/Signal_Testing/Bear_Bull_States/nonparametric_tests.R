  
  
  ############################################################################
  ### SHINY - BEAR AND BULL STATES - NONPARAMETRIC TESTS
  ############################################################################
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version      24.09.2020
  # First version     24.09.2020
  # --------------------------------------------------------------------------
  
  
  require(BBSolz)
  require(DAARC)
  
  wd <- "R:/Asset_Management/R/Shiny/Signal_Testing/Bear_Bull_States/"
  
  
  # BBS
  BBS <- loadSensor(sensor_name = "bbs")
  states_is <- BBS$signal$insample
  
  # Return series
  ticker <- c("WFEWEBC LX", "NDDUWI")
  X <- rodbcGetOLZDBReturns( assetName = ticker,
                             frqncy = "daily",
                             refCcy = "CHF",
                             na.rm = "r" )
  X <- X[ ,ticker]
  colnames(X) <- c("OLZ", "BM")
  X_delta <- (X[ ,"OLZ"] - X[ ,"BM"]) / (1 + X[ ,"BM"])
  colnames(X_delta) <- "Outperformance"
  
  
  plot( as.simTS(X_delta) ) 
  
  
  
  
  # --------------------------------------------------------------------------
  # Wilcoxon rank test
  # H0: Returns are the same in both funds
  # Ha: Returns are higher in the BM
  # --------------------------------------------------------------------------
  
  Y_delta <- log(1 + X[ ,"OLZ"]) - log(1 + X[ ,"BM"])
  wtest <- wilcox.test( Y_delta, alternative = "less" )
  wtest
  
  
  
  

  # --------------------------------------------------------------------------
  # Permutation test for regimes
  # H0: Excess returns are the same in both regimes
  # Ha: Excess returns are higher in bear regimes
  # --------------------------------------------------------------------------
  
  
  n_sim <- 10^3
  dates <- intersect(rownames(BBS$signal$insample), rownames(X_delta))
  states <- as.numeric( BBS$signal$insample[dates, ] )
  y <- as.numeric(X_delta[dates, ])
  
  M1 <- numeric( n_sim )
  for ( i in 1:n_sim ) {
    
    states_tmp <- sample(states)
    S <- cbind( (1 + states_tmp) / 2,
                (1 - states_tmp) / 2 )
    P <- apply( S, 2, function(x) { x / sum(x) } )
    mu <- t(P) %*% y
    M1[i] <- mu[1] - mu[2]
    
  }
  

  
  
  
  X_bm <- BBS$data$X_bm
  M2 <- numeric( n_sim )
  ratio <- M2
  for ( i in 1:n_sim ) {
    
    # Permute days in the return series (i.e., generate synthetic series)
    Y <- as.timeSeries(sample(as.numeric(X_bm)), time(X_bm))
    BBS_tmp <- bbs( cumulated(Y, "discrete"),
                    minphase = BBS$spec$minphase,
                    mincycle = BBS$spec$mincycle,
                    e = BBS$spec$e,
                    language = "C" )
    states_tmp <- sample(as.numeric(BBS_tmp[dates, ]))
    S <- cbind( (1 + states_tmp) / 2,
                (1 - states_tmp) / 2 )
    P <- apply( S, 2, function(x) { x / sum(x) } )
    mu <- t(P) %*% y
    ratio[i] <- mean( (1 + states_tmp) / 2)
    M2[i] <- mu[1] - mu[2]
    
  }
  
  plot(ratio, ylim = c(0, 1))
  abline( h = mean( (1 + BBS$signal$insample) / 2 ), col = 2 )
  abline( h = 0.5 )
  
  M <- cbind(M1 = M1, M2 = M2)
  ldens <- lapply( colnames(M), function(x) { densFUN(M[ ,x]) } )
  names(ldens) <- colnames(M)
  plot.ldensity( ldens )
  
  apply(M, 2, mean)
  apply(M, 2, sd)

  # "greater" is the alternative that x has a larger mean than y.
  t.test( x = M[ ,1], 
          y = M[ ,2], 
          alternative = "less",
          paired = FALSE )  
  
  
  
  # Dirichlet sampling
  P <- rdirichlet( n = n_sim, 
                   alpha = (1 + BBS$signal$insample[dates, ]) / 2 )
  M_bull <- as.numeric( P %*% y )
  P <- rdirichlet( n = n_sim, 
                   (1 - BBS$signal$insample[dates, ]) / 2 )
  M_bear <- as.numeric( P %*% y )
  
  M <- cbind(M1 = M1,
             M2 = M2, 
             M_bull = M_bull, 
             M_bear = M_bear, 
             M_delta = M_bull - M_bear)
  ldens <- lapply( colnames(M), function(x) { densFUN(M[ ,x]) } )
  names(ldens) <- colnames(M)
  plot.ldensity( ldens )
  abline(v = 0)
  abline(v = mean(M_bull))
  abline(v = mean(M_bear))
  
  
  # Alternatively,
  alpha_mat <- cbind( (1 + BBS$signal$insample[dates, ]) / 2,
                      (1 - BBS$signal$insample[dates, ]) / 2,
                      rep(1, length(dates)) )
  M <- matrix(NA, nrow = n_sim, ncol = ncol(alpha_mat))
  for ( j in 1:ncol(alpha_mat) ) {
    P <- rdirichlet( n = n_sim, 
                     alpha = alpha_mat[ ,j] )
    m <- as.numeric( P %*% y )
    M[ ,j] <- m
  }  
  
  ldens <- lapply( 1:ncol(M), function(x) { densFUN(M[ ,x]) } )
  # names(ldens) <- colnames(M)
  plot.ldensity( ldens )
  
  
  
  
  # --------------------------------------------------------------------------
  # Alternatively,
  # using fuzzy states
  # --------------------------------------------------------------------------
  
  BBSF <- loadSensor( sensor_name = "bbsfuzzy" )
  dates <- intersect(rownames(X_delta), rownames(BBSF$signal$insample))
  y <- as.numeric(X_delta[dates, ])
  alpha_mat <- cbind( (1 + BBSF$signal$insample[dates, ]) / 2,
                      (1 - BBSF$signal$insample[dates, ]) / 2,
                      rep(1, length(dates)) )
  M <- matrix(NA, nrow = n_sim, ncol = ncol(alpha_mat))
  for ( j in 1:ncol(alpha_mat) ) {
    P <- rdirichlet( n = n_sim, 
                     alpha = alpha_mat[ ,j] )
    m <- as.numeric( P %*% y )
    M[ ,j] <- m 
  }  
  
  ldens <- lapply( 1:ncol(M), function(x) { densFUN(M[ ,x]) } )
  # names(ldens) <- colnames(M)
  plot.ldensity( ldens )
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Synthetic example, iidN data
  # --------------------------------------------------------------------------
  
  set.seed(999)
  n <- 10^3
  x1 <- rnorm(n, 0, 1)
  x2 <- rnorm(n/2, -0.2, 1)
  x <- c(x1, x2)
  states <- c(rep(1, length(x1)), rep(0, length(x2)))
    
  
  t.test( x = x[states == 1], 
          y = x[states == 0], 
          # alternative = "less",
          paired = FALSE )
  
  
  ldens <- lapply( list(x1, x2), FUN = densFUN )
  plot.ldensity( ldens )
  
  
  # Via Dirichlet sampling
  n_sim <- 999
  P1 <- rdirichlet( n = n_sim, alpha = states )
  P2 <- rdirichlet( n = n_sim, alpha = (states - 1) * (-1) )
  
  mu1_dist <- as.numeric(P1 %*% x)
  mu2_dist <- as.numeric(P2 %*% x)
  mu_dist <- cbind(mu1_dist, mu2_dist)
  ldens <- apply( mu_dist, 2, densFUN )
  plot.ldensity( ldens )
  
  
  # Fuzzy states
  meanBeta <- function(a, b) { a / (a + b) }
  meanBeta(a = 1, b = 3)
  
  states_fuzzy <- states * 0
  states_fuzzy[ states == 1 ] <- rbeta(n = sum(states == 1), 
                                       shape1 = 3, 
                                       shape2 = 1)
  states_fuzzy[ states == 0 ] <- rbeta(n = sum(states == 0), 
                                       shape1 = 1, 
                                       shape2 = 3)
  
  
  P1f <- rdirichlet( n = n_sim, alpha = states_fuzzy )
  P2f <- rdirichlet( n = n_sim, alpha = (states_fuzzy - 1) * (-1) )
  mu1f_dist <- as.numeric(P1f %*% x)
  mu2f_dist <- as.numeric(P2f %*% x)
  muf_dist <- cbind(mu1f_dist, mu2f_dist)
  ldensf <- apply( muf_dist, 2, densFUN )
  plot.ldensity( ldensf )
  
  ldens <- apply( cbind(mu_dist, muf_dist), 2, densFUN )
  plot.ldensity( ldens )
  
  
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Influence of BBS parameters
  # Conjecture: using horter phases increases the return difference but also the 
  # variance of the mean so that significance is not increased.
  # --------------------------------------------------------------------------
  
  
  X_level <- cumulated(X[ ,"BM"], "discrete")
  b_logarithmic <- FALSE
  e <- 0
  theta <- 0.15
  # minphase_vec <- c(5, 15, 25, 35, 45)
  # mincycle_vec <- (21 * 1:7)
  minphase_vec <- c(1, 5, 15, 25, 35, 45)
  mincycle_vec <- c(5, 12, 21 * 1:7)
  grid <- expand.grid(minphase_vec, mincycle_vec)
  grid <- grid[apply(grid, 1, function(x) {x[1] < x[2] } ), ]
  Names <- apply( grid, 1, 
                  function(x) { paste0("minphase=", 
                                       x[1], "_mincycle=", x[2]) } )
  states <- X_level[ ,rep(1, length(Names))] * NA
  colnames(states) <- Names
  lBBS <- list()
  
  for ( k in 1:nrow(grid) ) {
    
    minphase <- grid[k, 1]
    mincycle <- grid[k, 2]
   
    BBS <- bbs( X = X_level,
                mincycle = mincycle,
                minphase = minphase,
                k.peak = minphase,
                l.peak = minphase,
                k.trough = minphase,
                l.trough = minphase,
                logarithmic = b_logarithmic,
                theta = theta,
                e = e,
                language = "C" )
    states[ ,k] <- BBS
  
  }
  states_bear <- (1 - states) / 2
  states_bull <- (1 + states) / 2
  
  states_bear_scl <- apply( states_bear, 2, function(x) { x / sum(x)} )
  states_bull_scl <- apply( states_bull, 2, function(x) { x / sum(x)} )
  ens_bear <- apply( states_bear_scl, 2, entropy, exponential = TRUE, eps = 1e-05)
  ens_bull <- apply( states_bull_scl, 2, entropy, exponential = TRUE, eps = 1e-05)
  
  
  n_sim <- 10^3
  y <- as.numeric(X[ ,"BM"])
  M_bear <- M_bull <- matrix(NA, nrow = n_sim, ncol = ncol(states))
  
  for ( k in 1:ncol(states) ) {
    
    alpha_bear <- (1 + states_scl[ ,k]) / 2 * ens_bear[k]
    alpha_bull <- (1 - states_scl[ ,k]) / 2 * ens_bull[k]
    P <- rdirichlet( n = n_sim, 
                     alpha = alpha_bear )
    M_bear[ ,k] <- as.numeric( P %*% y )
    P <- rdirichlet( n = n_sim, 
                     alpha = alpha_bull )
    M_bull[ ,k] <- as.numeric( P %*% y )
    
  }
  
  ldens <- apply( M_bear - M_bull, 2, densFUN )
  plot.ldensity( ldens )
  
  
  # idx <- grepl("mincycle=147", Names)
  idx <- grepl("minphase=1_", Names)
  ldens <- apply( (M_bear - M_bull)[ ,idx], 2, densFUN )
  names(ldens) <- Names[idx]
  plot.ldensity( ldens )
  
  
  idx <- grepl("mincycle=147", Names)
  idx <- grepl("minphase=1_", Names)
  ldens <- apply( M_bull[ ,idx], 2, densFUN )
  names(ldens) <- Names[idx]
  plot.ldensity( ldens )
  
  
  
  
  
  
  