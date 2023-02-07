  
  ############################################################################
  ### GOOD AND BAD TURBULENCE - SOURCE FUNCTIONS
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     10.10.2020
  # First version:    10.10.2020
  # --------------------------------------------------------------------------
  

  
  # mahalanobisDelta
  # factorIndices
  # rmvgarch2
  
  
  
  
  # --------------------------------------------------------------------------
  mahalanobisDelta <- function( x, 
                                center_1, 
                                covmat_1, 
                                center_2, 
                                covmat_2,
                                wghts = NULL,
                                b_scl = TRUE,
                                b_det = TRUE,
                                prior_1 = 0.5 )
  {
    
    covmat_1_inv <- solve(covmat_1)
    covmat_2_inv <- solve(covmat_2)
    if ( is.null(wghts) ) {
      wghts <- rep(1, ncol(covmat_1))
    }
    W <- diag(wghts)
    Sigma <- t(W) %*% ( covmat_1_inv - covmat_2_inv ) %*% W
    a_1 <- t(W) %*% covmat_1_inv %*% W %*% center_1
    a_2 <- t(W) %*% covmat_2_inv %*% W %*% center_2
    a <- 2 * (a_2 - a_1)
    g_1 <- t(center_1) %*% t(W) %*% covmat_1_inv %*% W %*% center_1
    g_2 <- t(center_2) %*% t(W) %*% covmat_2_inv %*% W %*% center_2
    constant <- g_1 - g_2
    delta <- t(x) %*% Sigma %*% x + t(a) %*% x + constant
    
    # prior <- log( prior_1 / (1 - prior_1) )
    prior_ratio <- log( (1 - prior_1) / prior_1 )
    
    if ( isTRUE(b_det) ) {
      det_ratio <- log(det(covmat_1)) - log(det(covmat_2))
    } else {
      det_ratio <- 0
    }
    
    if ( isTRUE(b_scl)) {
      # scl <- 1 / sum(wghts^2)
      scl <- 1 / sum(wghts)
    } else {
      scl <- 1
    }
    ans <- scl * delta + det_ratio + prior_ratio
      
    return( ans )
  }
  
  
  # --------------------------------------------------------------------------
  quadDisc <- function( x, 
                        center_1, 
                        covmat_1, 
                        center_2, 
                        covmat_2,
                        wghts = NULL,
                        b_scl = TRUE,
                        prior_1 = 0.5 )
  {
    
   
    return( ans )
  }
  
  
  
  
  
  # --------------------------------------------------------------------------
  factorIndices <- function( universe = c("europe_ex_ch", 
                                          "usa", 
                                          "dm") )
  {
    universe <- match.arg(universe)
    
    Names <- c("capw",
               # "minvar", 
               "value", 
               "quality", 
               "momentum",
               "lowvola")
    lTicker <- list(europe_ex_ch = c("NDDUE15", 
                                     # "M00IER$O", 
                                     "M7EUEV", 
                                     "M7EUQU", 
                                     "MAEUMMT",
                                     "M00IER$O"),
                    usa = c("NDDUUS", 
                            # "M00IMVSO", 
                            "M1USEV", 
                            "M1USQU", 
                            "M1USMMT",
                            "M00IMVSO"),
                    dm = c("NDDUWI", 
                           # "M00IWO$O",
                           "M1WOEV", 
                           "M1WOQU", 
                           "M1WOMOM",
                           "M00IWO$O"))
    
    ticker <- lTicker[[universe]]
    X <- rodbcGetOLZDBReturns( assetName = ticker,
                               refCcy = "USD",
                               frqncy = "daily" )
    X <- X[ ,ticker]
    colnames(X) <- Names
    
    return( X )
  }
  
  
 
  
  
  # --------------------------------------------------------------------------
  rmvgarch2 <- function( l_fit_or_spec = list(), 
                         M1 = NULL, 
                         M2 = NULL,
                         n_sample = 1000,
                         asTS = TRUE,
                         start_date = NULL,
                         frqncy = c('day', 'week', 'month', 'quarter', 'year'),
                         seed = NULL )
  {
    if ( is.null(l_fit_or_spec) ) {
      stop("list 'l_fit_or_spec' is empty.")
    }
    n_assets <- length(l_fit_or_spec)
    
    if ( is.null(M1) ) {
      M1 <- rep(0, n_assets)
    }
    if ( is.null(M2) ) {
      M2 <- diag(n_assets)
    }
    if ( !isSymmetric(M2,
                      tol = sqrt(.Machine$double.eps),
                      check.attributes = FALSE)) {
      stop("The covariance matrix must be a symmetric")
    }
    
    M2_tmp <- M2
    dimnames(M2_tmp) <- NULL
    if ( !isTRUE( all.equal( M2_tmp, t(M2_tmp)) ) ) {
      warning("The covariance matrix is numerically not symmetric")
    }
    
    E <- eigen( M2, symmetric = TRUE )

    if ( !all(E$values >= -sqrt(.Machine$double.eps) * abs(E$values[1])) ) {
      warning("The covariance matrix is numerically not positive definite")
    }
    # tmat <- torsion(M2, "mt")
    
    # Generate garch simulations based on fit
    frqncy <- match.arg(frqncy)
    if ( is.null(start_date) ) {
      start_date <- Sys.Date()
    }
    charvec <- seq.Date( from = as.Date(start_date),
                         by = frqncy,
                         length.out = n_sample )
    
    # Set seed
    if ( !is.null(seed) ) {
      set.seed( seed )
    }
    seed_vec <- round( runif(n_assets, 1, 10^6) )
    
    # Simulate 
    simFUN <- function(i)
    {
      sim <- ugarchsim( l_fit_or_spec[[i]], 
                        n.sim = n_sample, 
                        n.start = 1, 
                        m.sim = 1, 
                        startMethod = "sample",
                        rseed = seed_vec[i] )
      return( sim )
    }
    lSim <- lapply( 1:length(l_fit_or_spec), 
                    FUN = simFUN )
    PC <- do.call( cbind, lapply(lSim, function(x) { x@simulation$seriesSim } ) )
    PC <- timeSeries( PC, as.timeDate( charvec ) ) 
    
    # Transform simulation of PC to be uncorrelated
    PC_uncor <- mtca(PC)@sources
    
    # Scale PC's to reach desired variances
    PC_scl <- t( apply(scale(PC_uncor, TRUE, TRUE), 
                       1, 
                       function(x, y) { as.numeric(x) * as.numeric(y) }, 
                       y = sqrt(E$values) ) )
    
    # Transform PC's to correlated series
    X <- PC_scl %*% t(E$vectors)
    
    # Shift data to reach desired mean vector
    X <- sweep( scale(X, TRUE, FALSE), 2, M1, "+" )
    X <- exp(X) - 1  #~~~~~~~~~~~~~~~ 
    
    # Cast data to timeSeries
    if ( isTRUE(asTS) ) {
      X <- timeSeries( X, charvec )
    }
    colnames(X) <- colnames(M2)
    
    attr(X, "seed") <- seed
    attr(X, "seed_vec") <- seed_vec
    # attr(X, "PC") <- PC
    
    return( X ) 
  }
  
  
  