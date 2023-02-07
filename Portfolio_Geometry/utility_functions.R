

############################################################################
### AAA - UTILITY FUNCTIONS
############################################################################


# --------------------------------------------------------------------------
# Cyril Bachelard
# This version:     15.02.2018
# First version:    15.02.2018
# --------------------------------------------------------------------------
  
  
  
  
  # --------------------------------------------------------------------------
  # UNIVARIATE UTILITY FUNCTIONS
  # --------------------------------------------------------------------------
  
  # meanvarianceUtility
  # minvarianceUtility
  # mvUtility
  # entropyUtility
  # ktUtility
  # kinkedLogWealthUtility
  # kinkedPowerUtility
  # acUtility
  # lambda2b
  # caraUtility
  # crraUtility
  


  
  # --------------------------------------------------------------------------
  meanvarianceUtility <- function(x, lambda = 1 ) { meanGeo(x) - lambda * var(x) }
  # --------------------------------------------------------------------------
  minvarianceUtility <- function(x) { x^2 }
  # --------------------------------------------------------------------------
  mvUtility <- function(x, lambda = 1) { x - lambda/2 * (x - mean(x))^2 }
  # --------------------------------------------------------------------------
  entropyUtility <- function(x, lambda = 1) 
  { 
    ans <- -1 / lambda * log( mean(exp(-lambda * x)) )
    return( ans ) 
  }
  # --------------------------------------------------------------------------
  ktUtility <- function(mu, 
                        RP = 0, 
                        # a_plus = 0.88,
                        # a_minus = 0.88,
                        a = 0.88,
                        b = 2.25) 
  {
    ans <- mu
    idx_pos <- which(mu >= RP)
    if ( length(idx_pos) > 0 ) {
      ans[idx_pos] <- (mu[idx_pos] - RP)^a
    }
    idx_neg <- which(mu < RP)
    if ( length(idx_neg) > 0 ) {
      ans[idx_neg] <- -b * (RP - mu[idx_neg])^a
    }
    
    return( ans ) 
  }
  # --------------------------------------------------------------------------
  kinkedLogWealthUtility <- function(x, RP = 0, beta = 1)
  { 
    y <- x
    idx_neg <- which(x < RP)
    if ( length(idx_neg) > 0 ) {
      y[idx_neg] <- beta * (x[idx_neg] - RP) + log(1 + x[idx_neg])
    } else {
      y[idx_pos] <- log(1 + x[idx_pos])
    }
    return( y )
  }
  # --------------------------------------------------------------------------
  kinkedPowerUtility <- function(x, RP = 0, a = 1)
  { 
    y <- x
    idx_neg <- which(x < RP)
    if ( length(idx_neg) > 0 ) {
      y[idx_neg] <- RP - a * (RP - x[idx_neg])
    }
    if ( a == 1 ) {
      ans <- log(1 + y)
    } else {
      ans <- ( (1 + y)^(1 - RP) - 1) / (1 - RP)  
    }
    return( ans )
  }
  # --------------------------------------------------------------------------
  acUtility <- function(x, b1 = 0.5, b2 = 0.5)
  {
    ans <- x
    a <- 1 / b1
    idx_pos <- which(x >= 0)
    if ( length(idx_pos) > 0 ) {
      # ans[idx_pos] <- a * x[idx_pos]^b
      # ans[idx_pos] <- a * ((x[idx_pos] + 1)^b - 1)
      ans[idx_pos] <- a * ((x[idx_pos] + 1)^b1 - 1)
    }
    idx_neg <- which(x < 0)
    if ( length(idx_neg) > 0 ) {
      # ans[idx_neg] <- -a * (RP - x[idx_neg])^(1/b)
      # ans[idx_neg] <- -a * b^2 * ((-x[idx_neg] + 1)^(1/b) - 1)
      ans[idx_neg] <- -a * b1 * b2 * ((-x[idx_neg] + 1)^(1/b2) - 1)
    }
    
    return( ans ) 
  }
  # --------------------------------------------------------------------------
  lambda2b <- function ( lambda ) {
    ll <- ifelse(lambda < 0.5, log2(lambda), 0.5 * log2((0.5+lambda)^(0.5+lambda)))
    b2 <- 1 - 0.5 * (ll == 0) - as.numeric(substituteNA(((ll - 2 + sqrt(ll^2 + 4)) / (2 * ll)), method = "zeros"))
    b1 <- 1 - (1 - b2)^2
    if (lambda < 0.212264 ) {
      b2 <- b2 + ((1 - 5 * (2 * abs(b1 - b2 - 2 * 0.212264/3)))^2) / 5 - 0.0344605
    }
    return(list(b1 = b1, b2 = b2))
  }
  
  # --------------------------------------------------------------------------
  caraUtility <- function(x) { log(1 + x) }
  # --------------------------------------------------------------------------
  crraUtility <- function(x, g = 1)
  {
    if ( g == 1 ) {
      U <- log( 1 + x)
    } else {
      U <- ((1 + x)^(1 - g) - 1) / (1 - g)
    }
    return( U )
  }
  
  
  
  
  
  
  
  
  
  
  