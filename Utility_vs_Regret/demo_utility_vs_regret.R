  
  
  ############################################################################
  ### UTILITY - REGRET
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     02.06.2020
  # First version:    02.06.2020
  # --------------------------------------------------------------------------

  
  require(RP)
  require(DAARC)
  wd <- "R:/Asset_Management/Research_Projects/Multi_Asset/DAA/DAARC/"
  source( paste0(wd, "Source/custom_functions.R") )
  
  
  # --------------------------------------------------------------------------
  # STATES ARE -1 OR 1
  # --------------------------------------------------------------------------
  
  expected_state <- seq(from = -1, to = 1, length.out = 100)
  u_kt <- ktUtility(mu = expected_state, RP = 0)
  u_mv1 <- mvUtility(x = expected_state, lambda = 1)
  u_mv2 <- mvUtility(x = expected_state, lambda = 2)
  
  
  plot( x = expected_state, y = u_kt )
  points( x = expected_state, y = u_mv1, col = 2 )
  points( x = expected_state, y = u_mv2, col = 3 )
  
  plot( density(expected_state), ylim = c(-2, 1) )
  
  
  
  
  
  
  
  
  # Utility of expected states * portfolio weights
  wghts <- seq(from = 0, to = 1, length.out = 100)
  w_mat <- matrix( wghts, ncol = 1 )
  U_kt <- apply( w_mat, 1, function(w) { ktUtility( mu = expected_state * w, RP = 0 ) } )
  U_mv1 <- apply( w_mat, 1, function(w) { mvUtility( x = expected_state * w, lambda = 1 ) } )
  U_mv2 <- apply( w_mat, 1, function(w) { mvUtility( x = expected_state * w, lambda = 2 ) } )
  
  EU_kt <- apply(U_kt, 2, mean)
  ER_kt <- rev(EU_kt) * (-1) + max(EU_kt)
  EU_mv1 <- apply(U_mv1, 2, mean)
  ER_mv1 <- rev(EU_mv1) * (-1) + max(EU_mv1)
  EU_mv2 <- apply(U_mv2, 2, mean)
  ER_mv2 <- rev(EU_mv2) * (-1) + max(EU_mv2)
  
  
  EU <- EU_mv1
  ER <- ER_mv1
  lambda <- 0.7
  
  plot( x = wghts, y = EU - ER, ylim = range(c(EU, ER, EU - ER)) )
  abline(h = 0)
  points( x = wghts, y = EU, col = 3 )
  points( x = wghts, y = ER, col = 2 )
  points( x = wghts, y = EU - ER * lambda, col = 4 )
  
  
  
  
 
  
  plot( x = wghts, y = EU - ER, ylim = range(c(EU, ER, EU - ER)),
        main = rebdate )
  abline(h = 0)
  points( x = wghts, y = EU, col = 3 )
  points( x = wghts, y = ER, col = 2 )
  points( x = wghts, y = EU - ER * lambda, col = 4 )
  
    
  
  
  # --------------------------------------------------------------------------
  # STATES ARE 0 OR 1  -  BETA DISTRIBUTION
  # --------------------------------------------------------------------------
  
  
  
  expected_state <- rbeta(n = 10^5, shape1 = 1.2, shape2 = 0.6)
  mean(expected_state)

  u_kt <- ktUtility(mu = expected_state, RP = 0)
  u_mv1 <- mvUtility(x = expected_state, lambda = 1)
  u_mv2 <- mvUtility(x = expected_state, lambda = 2)
  
  plot( density(expected_state) )
  
  
  
  # Utility of expected states * portfolio weights
  wghts <- seq(from = 0, to = 1, length.out = 100)
  w_mat <- matrix( wghts, ncol = 1 )
  U_kt <- apply( w_mat, 1, function(w) { ktUtility( mu = expected_state * w, RP = 0 ) } )
  U_mv1 <- apply( w_mat, 1, function(w) { mvUtility( x = expected_state * w, lambda = 1 ) } )
  U_mv2 <- apply( w_mat, 1, function(w) { mvUtility( x = expected_state * w, lambda = 2 ) } )
  
  EU_kt <- apply(U_kt, 2, mean)
  ER_kt <- rev(EU_kt) * (-1) + max(EU_kt)
  EU_mv1 <- apply(U_mv1, 2, mean)
  ER_mv1 <- rev(EU_mv1) * (-1) + max(EU_mv1)
  EU_mv2 <- apply(U_mv2, 2, mean)
  ER_mv2 <- rev(EU_mv2) * (-1) + max(EU_mv2)
  
  
  EU <- EU_mv2
  ER <- ER_mv2
  es_dens <- densFUN( x = expected_state )
  lambda <- entropy( es_dens$y, exponential = TRUE ) / length(es_dens$y)
  lambda
  
  plot( x = wghts, y = EU - ER, ylim = range(c(EU, ER, EU - ER)),
        main = paste0("lambda = ", round(lambda, 3)) )
  abline(h = 0)
  points( x = wghts, y = EU, col = 3 )
  points( x = wghts, y = ER, col = 2 )
  points( x = wghts, y = EU - ER * lambda, col = 4 )
  idx <- which( (EU - ER * lambda) == max((EU - ER * lambda)) )
  points( x = wghts[idx], y = (EU - ER * lambda)[idx], col = 4, pch = 4, cex = 2 )
  idx
  
  
  
  
  
  plot( es_dens, ylab = "density of expected states", xlim = c(0, 1), 
        main = "Density and Utility", xlab = "" )
  par( new = TRUE )
  plot( x = sort(expected_state), y = sort(u_kt), type = "l", 
        ylab = "", xlab = "", xlim = c(0, 1), 
        ylim = range(c(u_kt, u_mv1, u_mv2)), 
        xaxt = "n", yaxt = "n", pos = 4)
  points( x = sort(expected_state), y = sort(u_mv1), type = "l", col = 2 )
  points( x = sort(expected_state), y = sort(u_mv2), type = "l", col = 3 )
  axis(side = 4)
  mtext("Utility", side = 4, line = 3)
  
  
  
  
  
  
  
  