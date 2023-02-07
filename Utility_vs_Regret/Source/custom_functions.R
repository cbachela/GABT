  
  
  ############################################################################
  ### BacktestDAA - CUSTOM FUNCTIONS
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     29.06.2020
  # First version:    26.05.2020
  # --------------------------------------------------------------------------
  
  
  # Description:
  # 
  # This file contains custom functions for testing and work in progress that are
  # not in the DAARC package. When deemed useful, functions are copied to the package.
  
  
 
  # BacktestDAA.rpm
  # BacktestDAA.rpm_mvu
  # BacktestDAA.rpm_mvu_disc
  # BacktestDAA.rpm_mvu_qnn
  # BacktestDAA.rpm_ada
  # BacktestDAA.rpmAsyTurb
  # BacktestDAA.rpmAsyTurb_2
  # BacktestDAA.mom_lo
  # BacktestDAA.mom_ls
  # BacktestDAA.momvariancePortfolio
  # BacktestDAA.runDAA
  
  
  
  
  # --------------------------------------------------------------------------
  BacktestDAA.rpm <- function( rebdate )
  {
    
    # BBS is in-sample. 
    # To do: write function to extract last reliable state
    # for oos backtest.
    
    Y_train <- data$TD$Y_train[rownames(data$TD$Y_train) <= rebdate, ]
    X_train <- data$TD$X_train[rownames(data$TD$X_train) <= rebdate, ]
    X_eval <- data$TD$X_train[rebdate, ]
    # md <- mahalanobis(X_train, 
    #                   center = X_eval, 
    #                   cov = covmat_X_train)
    md <- abs(X_train - as.numeric(X_eval))^2 
    eps <- density(md)$bw
    p <- exp( -1 / (2 * eps) * md )
    p <- p / sum(p)
    
    n_sim <- ifelse( is.null(spec$n_sim), 10^3, spec$n_sim )
    P <- rdirichlet(n = n_sim, alpha = p )
    # P <- rdirichlet(n = n_sim, alpha = p * length(p) )
    
    # Expected states as RPM weighted average of state series
    y_train <- as.numeric(Y_train)
    es <- apply( P, 1, function(p) { t(p) %*% y_train })
    es_dens <- densFUN(es)
    lambda <- entropy( es_dens$y, exponential = TRUE) / length(es_dens$y)
    
    # Utility of expected states * portfolio weights
    wghts <- seq(from = 0, to = 1, length.out = 100)
    w_mat <- matrix( wghts, ncol = 1 )
    U <- apply( w_mat, 1, function(w) { ktUtility( mu = es * w, RP = 0 ) } )
    # U <- apply( w_mat, 1, function(w) { mvUtility( x = es * w, lambda = 2 ) } )
    EU <- apply(U, 2, mean)
    # Define regret as the reversed inverted empirical expected utility
    ER <- rev(EU) * (-1) + max(EU)
    
    
    # plot( es, ktUtility(mu = es) )
    # points( es, mvUtility(x = es, lambda = 1), col = 2 )
    # points( es, mvUtility(x = es, lambda = 2), col = 3 )
    
    
    
    
    
    
    
    
    # matplot( x = wghts, y = t(U), type = "l" )
    # lines( x = wghts, y = apply(U, 2, mean), lwd = 3, col = "orange" )
    
    # plot( x = wghts, y = EU - ER, ylim = range(c(EU, ER, EU - ER)),
    #       main = rebdate )
    # abline(h = 0)
    # points( x = wghts, y = EU, col = 3 )
    # points( x = wghts, y = ER, col = 2 )
    # points( x = wghts, y = EU - ER * lambda, col = 4 )
    
    
    
    w_star <- which( EU - ER * lambda == max(EU - ER * lambda) )[1] / 100
    w_star_nr <- which( EU == max(EU) )[1] / 100
    
    mat <- matrix( c(w_star, w_star_nr), nrow = 1,
                   dimnames = list(rebdate, c("w_star", "w_star_nr")) )
    lOut <- list( weights = as.timeSeries(mat),
                  lambda = timeSeries(matrix(lambda), rebdate),
                  exp_state = timeSeries(matrix(es, nrow = 1), rebdate),
                  exp_util = timeSeries(matrix(EU, nrow = 1), rebdate),
                  exp_regret = timeSeries(matrix(EU, nrow = 1), rebdate) )
    
    .self$appendOutput( value = lOut )
    
    return( TRUE )
  }
  
  
  
  # --------------------------------------------------------------------------
  BacktestDAA.rpm_mvu <- function( rebdate )
  {
    
    # Same as BacktestDAA.rpm but with mvUtility function with lambda = 1
    # instead of ktUtility.
    # Also, one has to make sure that Y_train is \in {-1, 1}.

    
    Y_train <- data$TD$Y_train[rownames(data$TD$Y_train) <= rebdate, ]
    X_train <- data$TD$X_train[rownames(data$TD$X_train) <= rebdate, ]
    # X_eval <- data$TD$X_train[rebdate, ]
    X_eval <- data$signal[rebdate, ]
    # md <- mahalanobis(X_train, 
    #                   center = X_eval, 
    #                   cov = covmat_X_train)
    md <- abs(X_train - as.numeric(X_eval))^2 
    eps <- density(md)$bw
    p <- exp( -1 / (2 * eps) * md )
    p <- p / sum(p)
    sclfct <- entropy(p, exponential = TRUE)
    
  
    n_sim <- ifelse( is.null(spec$n_sim), 10^3, spec$n_sim )
    P <- rdirichlet(n = n_sim, alpha = p )
    # P <- rdirichlet(n = n_sim, alpha = p * length(p) )
    # P <- rdirichlet(n = n_sim, alpha = p * sclfct )
    
    
    # Expected states as RPM weighted average of state series
    y_train <- as.numeric(Y_train)
    es <- apply( P, 1, function(p) { t(p) %*% y_train })
    # Correct bias
    es <- scale(es, TRUE, FALSE) + as.numeric( t(p) %*% y_train )
    es_dens <- densFUN(es)
    h <- entropy( es_dens$y, exponential = TRUE) / length(es_dens$y)
    
    # Utility of expected states * portfolio weights
    wghts <- seq(from = 0, to = 1, by = 0.01)
    w_mat <- matrix( wghts, ncol = 1 )
    U <- apply( w_mat, 1, function(w) { mvUtility( x = es * w, lambda = spec$mvu_lambda ) } )
    # Average utility, i.e., expected utility of expected state E(U(E_p(S)))
    EU <- apply(U, 2, mean)
    # Define regret as the reversed inverted empirical expected utility
    ER <- rev(EU) * (-1) + max(EU)
    
    # U2 <- apply( (1 - w_mat), 1, function(w) { mvUtility( x = es * w, lambda = spec$mvu_lambda ) } )
    # EU2 <- apply(U2, 2, mean)
    # 
    # plot( wghts, EU )
    # lines( wghts, EU2 * -1 + max(EU), col = 2)
    
    
    # plot( es, ktUtility(mu = es) )
    # points( es, mvUtility(x = es, lambda = 1), col = 2 )
    # points( es, mvUtility(x = es, lambda = 2), col = 3 )
    
    # matplot( x = wghts, y = t(U), type = "l" )
    # lines( x = wghts, y = apply(U, 2, mean), lwd = 3, col = "orange" )
    
    matplot( x = wghts, 
             y = cbind(EU, ER, EU - ER, EU - h * ER ), 
             type = "l", 
             lwd = 2,
             main = paste0(rebdate, " (entropy = ", round(h, 3), ")") )
    abline(h = 0, col = "grey")
    
    
    # Optimal weights
    obj_value <- EU - ER * lambda
    idx_star <- which( obj_value == max(obj_value) )[1]
    idx_star_nr <- which( EU == max(EU) )[1]
    w_star <- wghts[ idx_star ]
    w_star_nr <- wghts[ idx_star_nr ]
    w_star_ntr <- w_star
    w_star_nr_ntr <- w_star_nr
    
    # Compare to initial weights
    idx <- which(spec$rebdates == rebdate)
    if ( idx > 1 ) {
      
      yesterday <- spec$rebdates[ idx - 1 ]
      w_init <- as.numeric( output$weights[yesterday, "w_star_ntr"] )
      w_init_nr <- as.numeric( output$weights[yesterday, "w_star_nr_ntr"] )
      util_w_init <- EU[ which(wghts %in% w_init) ]
      util_w_init_nr <- EU[ which(wghts %in% w_init_nr) ]
      util_w_star <- EU[ idx_star ]
      util_w_star_nr <- EU[ idx_star_nr ]
      
      cbind( c(w_star, w_init), 
             c(w_star_nr, w_init_nr) )
      
      cbind( c(util_w_star, util_w_init), 
             c(util_w_star_nr, util_w_init_nr) )
      
      # Only rebalance if gain in utility is larger than turnover
      if ( (util_w_star - util_w_init) > (abs(w_star - w_init) / 3) ) {
        w_star_ntr <- w_star
      } else {
        w_star_ntr <- w_init
      }
      if ( (util_w_star_nr - util_w_init_nr) > (abs(w_star_nr - w_init_nr) / 3) ) {
        w_star_nr_ntr <- w_star_nr
      } else {
        w_star_nr_ntr <- w_init_nr
      }
      
    }
    
    # idx <- c( which(wghts == w_init), which(wghts == w_star) )
    # idx_delta <- seq(from = min(idx), to = max(idx), by = 1)
    # w_delta <- wghts[ idx_delta ]
    # plot( x = abs(w_delta - w_init), y = (EU - ER * lambda)[ idx_delta ], type = "o" )
    
    
    mat <- matrix( c(w_star, w_star_nr, 
                     w_star_ntr, w_star_nr_ntr), nrow = 1,
                   dimnames = list(rebdate, c("w_star", "w_star_nr",
                                              "w_star_ntr", "w_star_nr_ntr")) )
    lOut <- list( weights = as.timeSeries(mat),
                  H = timeSeries(matrix(H), rebdate),
                  exp_state = timeSeries(matrix(es, nrow = 1), rebdate),
                  exp_util = timeSeries(matrix(EU, nrow = 1), rebdate),
                  exp_regret = timeSeries(matrix(EU, nrow = 1), rebdate) )
    
    .self$appendOutput( value = lOut )
    
    return( TRUE )
  }
  
  # --------------------------------------------------------------------------
  BacktestDAA.rpm_mvu_disc <- function( rebdate )
  {
    
    # Same as BacktestDAA.rpm but with mvUtility function with lambda = 1
    # instead of ktUtility.
    # Also, one has to make sure that Y_train is \in {-1, 1}.
    # Further, weights are discretized.
    
    Y_train <- data$TD$Y_train[rownames(data$TD$Y_train) <= rebdate, ]
    X_train <- data$TD$X_train[rownames(data$TD$X_train) <= rebdate, ]
    # X_eval <- data$TD$X_train[rebdate, ]
    X_eval <- data$signal[rebdate, ]
    # md <- mahalanobis(X_train, 
    #                   center = X_eval, 
    #                   cov = covmat_X_train)
    md <- abs(X_train - as.numeric(X_eval))^2 
    eps <- density(md)$bw
    p <- exp( -1 / (2 * eps) * md )
    p <- p / sum(p)
    sclfct <- entropy(p, exponential = TRUE)
    
    
    n_sim <- ifelse( is.null(spec$n_sim), 10^3, spec$n_sim )
    P <- rdirichlet(n = n_sim, alpha = p )
    # P <- rdirichlet(n = n_sim, alpha = p * length(p) )
    # P <- rdirichlet(n = n_sim, alpha = p * sclfct )
    

    # Expected states as RPM weighted average of state series
    y_train <- as.numeric(Y_train)
    es <- apply( P, 1, function(p) { t(p) %*% y_train })
    # Correct bias
    es <- scale(es, TRUE, FALSE) + as.numeric( t(p) %*% y_train )
    es_dens <- densFUN(es)
    lambda <- entropy( es_dens$y, exponential = TRUE) / length(es_dens$y)
    
    # Utility of expected states * portfolio weights
    wghts <- c(0, 0.25, 0.5, 0.75, 1)
    # wghts <- seq(from = 0, to = 1, by = 0.01)
    w_mat <- matrix( wghts, ncol = 1 )
    U <- apply( w_mat, 1, function(w) { mvUtility( x = es * w, lambda = spec$mvu_lambda ) } )
    # Average utility, i.e., expected utility of expected state E(U(E_p(S)))
    EU <- apply(U, 2, mean)
    # Define regret as the reversed inverted empirical expected utility
    ER <- rev(EU) * (-1) + max(EU)
    
  
    # plot( es, ktUtility(mu = es) )
    # points( es, mvUtility(x = es, lambda = 1), col = 2 )
    # points( es, mvUtility(x = es, lambda = 2), col = 3 )
    
    # matplot( x = wghts, y = t(U), type = "l" )
    # lines( x = wghts, y = apply(U, 2, mean), lwd = 3, col = "orange" )
    
    matplot( x = wghts, 
             y = cbind(EU, ER, EU-ER, EU-ER*lambda ), 
             type = "o", lwd = 2,
             main = rebdate )
    
   
    # Optimal weights
    obj_value <- EU - ER * lambda
    idx_star <- which( obj_value == max(obj_value) )[1]
    idx_star_nr <- which( EU == max(EU) )[1]
    w_star <- wghts[ idx_star ]
    w_star_nr <- wghts[ idx_star_nr ]
    w_star_ntr <- w_star
    w_star_nr_ntr <- w_star_nr
    
    # Compare to initial weights
    idx <- which(spec$rebdates == rebdate)
    if ( idx > 1 ) {
      
      yesterday <- spec$rebdates[ idx - 1 ]
      w_init <- as.numeric( output$weights[yesterday, "w_star_ntr"] )
      w_init_nr <- as.numeric( output$weights[yesterday, "w_star_nr_ntr"] )
      util_w_init <- EU[ which(wghts %in% w_init) ]
      util_w_init_nr <- EU[ which(wghts %in% w_init_nr) ]
      util_w_star <- EU[ idx_star ]
      util_w_star_nr <- EU[ idx_star_nr ]
      
      cbind( c(w_star, w_init), 
             c(w_star_nr, w_init_nr) )
      
      cbind( c(util_w_star, util_w_init), 
             c(util_w_star_nr, util_w_init_nr) )
      
      # Only rebalance if gain in utility is larger than turnover
      if ( (util_w_star - util_w_init) > (abs(w_star - w_init) / 3) ) {
        w_star_ntr <- w_star
      } else {
        w_star_ntr <- w_init
      }
      if ( (util_w_star_nr - util_w_init_nr) > (abs(w_star_nr - w_init_nr) / 3) ) {
        w_star_nr_ntr <- w_star_nr
      } else {
        w_star_nr_ntr <- w_init_nr
      }
      
    }
    
    # idx <- c( which(wghts == w_init), which(wghts == w_star) )
    # idx_delta <- seq(from = min(idx), to = max(idx), by = 1)
    # w_delta <- wghts[ idx_delta ]
    # plot( x = abs(w_delta - w_init), y = (EU - ER * lambda)[ idx_delta ], type = "o" )
    
    
    mat <- matrix( c(w_star, w_star_nr, 
                     w_star_ntr, w_star_nr_ntr), nrow = 1,
                   dimnames = list(rebdate, c("w_star", "w_star_nr",
                                              "w_star_ntr", "w_star_nr_ntr")) )
    lOut <- list( weights = as.timeSeries(mat),
                  lambda = timeSeries(matrix(lambda), rebdate),
                  exp_state = timeSeries(matrix(es, nrow = 1), rebdate),
                  exp_util = timeSeries(matrix(EU, nrow = 1), rebdate),
                  exp_regret = timeSeries(matrix(EU, nrow = 1), rebdate) )
    
    .self$appendOutput( value = lOut )
    
    return( TRUE )
  }
  
  
  # --------------------------------------------------------------------------
  BacktestDAA.rpm_mvu_qnn <- function( rebdate )
  {
    
    # Same as BacktestDAA.rpm_mvu but with an alternative distance function
    # to weight historical data.
    
    
    Y_train <- data$TD$Y_train[rownames(data$TD$Y_train) <= rebdate, ]
    X_train <- data$TD$X_train[rownames(data$TD$X_train) <= rebdate, ]
    # X_eval <- data$TD$X_train[rebdate, ]
    X_eval <- data$signal[rebdate, ]
    # # md <- mahalanobis(X_train, 
    # #                   center = X_eval, 
    # #                   cov = covmat_X_train)
    # md <- abs(X_train - as.numeric(X_eval))^2 
    # eps <- density(md)$bw
    # p <- exp( -1 / (2 * eps) * md )
    # p <- p / sum(p)
    # sclfct <- entropy(p, exponential = TRUE)
    
    # DS <- dirichletSampling( Y_train = Y_train,
    #                          X_train = X_train,
    #                          weights_fun = "cmeans",
    #                          centers = 2,
    #                          correct_bias = FALSE,
    #                          sclfct = 1 )
    # ldens <- lapply( DS, FUN = densFUN )
    # plot.ldensity( ldens )
    
    w <- weightsFun( data = trainingData(Y_train = Y_train,
                                         X_train = X_train),
                     x_eval = X_eval,
                     method = "cmeans", 
                     centers = 3 )
    w <- as.timeSeries(w, time(X_train))
    # plot( w)
    p <- w
    p <- p / sum(p)
    # sclfct <- entropy(p, exponential = TRUE)
    
    
    # plot( X_train )
    # abline( h = X_eval )
    # lines( p * 10^4, col = 3 )
    
    
    n_sim <- ifelse( is.null(spec$n_sim), 10^3, spec$n_sim )
    P <- rdirichlet(n = n_sim, alpha = p )
    # P <- rdirichlet(n = n_sim, alpha = p * length(p) )
    # P <- rdirichlet(n = n_sim, alpha = p * sclfct )
    
    
    # Expected states as RPM weighted average of state series
    y_train <- as.numeric(Y_train)
    es <- apply( P, 1, function(p) { t(p) %*% y_train })
    # Correct bias
    es <- scale(es, TRUE, FALSE) + as.numeric( t(p) %*% y_train )
    es_dens <- densFUN(es)
    lambda <- entropy( es_dens$y, exponential = TRUE) / length(es_dens$y)
    
    # Utility of expected states * portfolio weights
    wghts <- seq(from = 0, to = 1, by = 0.01)
    w_mat <- matrix( wghts, ncol = 1 )
    U <- apply( w_mat, 1, function(w) { mvUtility( x = es * w, lambda = spec$mvu_lambda ) } )
    # Average utility, i.e., expected utility of expected state E(U(E_p(S)))
    EU <- apply(U, 2, mean)
    # Define regret as the reversed inverted empirical expected utility
    ER <- rev(EU) * (-1) + max(EU)
    
    
    # plot( es, ktUtility(mu = es) )
    # points( es, mvUtility(x = es, lambda = 1), col = 2 )
    # points( es, mvUtility(x = es, lambda = 2), col = 3 )
    
    # matplot( x = wghts, y = t(U), type = "l" )
    # lines( x = wghts, y = apply(U, 2, mean), lwd = 3, col = "orange" )
    
    matplot( x = wghts, 
             y = cbind(EU, ER, EU-ER, EU-ER*lambda ), 
             type = "o", lwd = 2,
             main = rebdate )
    
    
    # Optimal weights
    obj_value <- EU - ER * lambda
    idx_star <- which( obj_value == max(obj_value) )[1]
    idx_star_nr <- which( EU == max(EU) )[1]
    w_star <- wghts[ idx_star ]
    w_star_nr <- wghts[ idx_star_nr ]
    w_star_ntr <- w_star
    w_star_nr_ntr <- w_star_nr
    
    # Compare to initial weights
    idx <- which(spec$rebdates == rebdate)
    if ( idx > 1 ) {
      
      yesterday <- spec$rebdates[ idx - 1 ]
      w_init <- as.numeric( output$weights[yesterday, "w_star_ntr"] )
      w_init_nr <- as.numeric( output$weights[yesterday, "w_star_nr_ntr"] )
      util_w_init <- EU[ which(wghts %in% w_init) ]
      util_w_init_nr <- EU[ which(wghts %in% w_init_nr) ]
      util_w_star <- EU[ idx_star ]
      util_w_star_nr <- EU[ idx_star_nr ]
      
      cbind( c(w_star, w_init), 
             c(w_star_nr, w_init_nr) )
      
      cbind( c(util_w_star, util_w_init), 
             c(util_w_star_nr, util_w_init_nr) )
      
      # Only rebalance if gain in utility is larger than turnover
      if ( (util_w_star - util_w_init) > (abs(w_star - w_init) / 3) ) {
        w_star_ntr <- w_star
      } else {
        w_star_ntr <- w_init
      }
      if ( (util_w_star_nr - util_w_init_nr) > (abs(w_star_nr - w_init_nr) / 3) ) {
        w_star_nr_ntr <- w_star_nr
      } else {
        w_star_nr_ntr <- w_init_nr
      }
      
    }
    
    # idx <- c( which(wghts == w_init), which(wghts == w_star) )
    # idx_delta <- seq(from = min(idx), to = max(idx), by = 1)
    # w_delta <- wghts[ idx_delta ]
    # plot( x = abs(w_delta - w_init), y = (EU - ER * lambda)[ idx_delta ], type = "o" )
    
    
    mat <- matrix( c(w_star, w_star_nr, 
                     w_star_ntr, w_star_nr_ntr), nrow = 1,
                   dimnames = list(rebdate, c("w_star", "w_star_nr",
                                              "w_star_ntr", "w_star_nr_ntr")) )
    lOut <- list( weights = as.timeSeries(mat),
                  lambda = timeSeries(matrix(lambda), rebdate),
                  exp_state = timeSeries(matrix(es, nrow = 1), rebdate),
                  exp_util = timeSeries(matrix(EU, nrow = 1), rebdate),
                  exp_regret = timeSeries(matrix(EU, nrow = 1), rebdate) )
    
    .self$appendOutput( value = lOut )
    
    return( TRUE )
  }
  
  
  
  # --------------------------------------------------------------------------
  BacktestDAA.rpm_ada <- function( rebdate )
  {
    
    # Same as BacktestDAA.rpm_mvu but with adaptive risk and regret aversion
    # parameter depending on current allocation and a no-trade-region.
    
    Y_train <- data$TD$Y_train[rownames(data$TD$Y_train) <= rebdate, ]
    X_train <- data$TD$X_train[rownames(data$TD$X_train) <= rebdate, ]
    X_eval <- data$TD$X_train[rebdate, ]
    # md <- mahalanobis(X_train, 
    #                   center = X_eval, 
    #                   cov = covmat_X_train)
    md <- abs(X_train - as.numeric(X_eval))^2 
    eps <- density(md)$bw
    p <- exp( -1 / (2 * eps) * md )
    p <- p / sum(p)
    
    n_sim <- ifelse( is.null(spec$n_sim), 10^3, spec$n_sim )
    P <- rdirichlet(n = n_sim, alpha = p )
    # P <- rdirichlet(n = n_sim, alpha = p * length(p) )
    
    # Expected states as RPM weighted average of state series
    y_train <- as.numeric(Y_train)
    es <- apply( P, 1, function(p) { t(p) %*% y_train } )
    es_dens <- densFUN(es)
    lambda <- entropy( es_dens$y, exponential = TRUE) / length(es_dens$y)
    
    # Initial weights
    # w_init <- .self$initialWeights( rebdate = rebdate )
    w_init <- 0
    if ( !is.null(output$weights) ) {
      idx <- which(rownames(output$weights) < rebdate )
      if ( length(idx) > 0 ) {
        w_init <- as.numeric( output$weights[idx[length(idx)], ] )
      }
    }
    
    # Potential new weights 
    wghts <- seq(from = 0, to = 1, length.out = 100)
    idx_ntr <- which(abs(wghts - w_init) < 0.2)
    
    # Utility of expected states * portfolio weights
    w_mat <- matrix( wghts, ncol = 1 )
    riskaversion <- 1 * (1 + w_init)
    U <- apply( w_mat, 1, function(w) { mvUtility( x = es * w, lambda = riskaversion ) } )
    EU <- apply(U, 2, mean)
    # Define regret as the reversed inverted empirical expected utility
    ER <- rev(EU) * (-1) + max(EU)
    
    
    # plot( es, ktUtility(mu = es) )
    # points( es, mvUtility(x = es, lambda = 1), col = 2 )
    # points( es, mvUtility(x = es, lambda = 2), col = 3 )
    
    # matplot( x = wghts, y = t(U), type = "l" )
    # lines( x = wghts, y = apply(U, 2, mean), lwd = 3, col = "orange" )
    
    matplot( x = wghts, 
             y = cbind(EU, ER, EU-ER, EU-ER*lambda ), 
             type = "o", lwd = 2,
             main = rebdate )
    
    
    idx_star <- which( EU - ER * lambda == max(EU - ER * lambda) )[1]
    w_star <- wghts[ idx_star ]
    if ( abs(w_star - w_init) < 0.2 ) {
      w_star <- w_init
    }
    idx_star_nr <- which( EU == max(EU) )[1]
    w_star_nr <- wghts[ idx_star_nr ]
    if ( abs(w_star_nr - w_init) < 0.2 ) {
      w_star_nr <- w_init
    }
    
    
    w_star <- matrix( w_star, nrow = 1, dimnames = list(rebdate, "w_star") )
    w_star_nr <- matrix( w_star_nr, nrow = 1, dimnames = list(rebdate, "w_star_nr") )
    lOut <- list( weights = as.timeSeries(w_star),
                  weights2 = as.timeSeries(w_star_nr),
                  lambda = timeSeries(matrix(lambda), rebdate),
                  exp_state = timeSeries(matrix(es, nrow = 1), rebdate),
                  exp_util = timeSeries(matrix(EU, nrow = 1), rebdate),
                  exp_regret = timeSeries(matrix(EU, nrow = 1), rebdate) )
    
    .self$appendOutput( value = lOut )
    
    return( TRUE )
  }
  
  
  
  
  
  # --------------------------------------------------------------------------
  Backtest.rpmAsyTurb <- function( rebdate )
  {
    
    
    Y_train <- data$TD$Y_train[rownames(data$TD$Y_train) <= rebdate, ]
    X_train <- data$TD$X_train[rownames(data$TD$X_train) <= rebdate, ]
    X_eval <- data$TD$X_train[rebdate, ]
    
    # Compute mahalanobis distance in signal space seperately 
    # for bear and bull specific coefficients.
    
    p_bear <- data$p_bear[rownames(X_train), ] / sum(data$p_bear[rownames(X_train), ])
    mu_bear <- t(X_train) %*% p_bear
    scnd_mom = t(X_train) %*% (X_train * (p_bear %*% matrix(1, 1, ncol(X_train))) )
    scnd_mom = (scnd_mom + t(scnd_mom)) / 2
    sigma_bear = scnd_mom - mu_bear %*% t(mu_bear)
    md_bear <- mahalanobis.ts(X = X_train, # scale(X_train, mu_bear, TRUE),
                              center = mu_bear, # X_eval,
                              covmat = sigma_bear)
    eps_bear <- density(md_bear)$bw
    w_bear <- exp( -1 / (2 * eps_bear) * md_bear )
    w_bear <- w_bear / sum(w_bear)
    
    
    
    p_bull <- data$p_bull[rownames(X_train), ] / sum(data$p_bull[rownames(X_train), ])
    mu_bull <- t(X_train) %*% p_bull
    scnd_mom = t(X_train) %*% (X_train * (p_bull %*% matrix(1, 1, ncol(X_train))) )
    scnd_mom = (scnd_mom + t(scnd_mom)) / 2
    sigma_bull = scnd_mom - mu_bull %*% t(mu_bull)
    md_bull <- mahalanobis.ts(X = X_train, # scale(X_train, mu_bull, TRUE),
                              center = mu_bull, # X_eval,
                              covmat = sigma_bull)
    eps_bull <- density(md_bull)$bw
    w_bull <- exp( -1 / (2 * eps_bull) * md_bull )
    w_bull <- w_bull / sum(w_bull)
    
    
    
    
    md_delta <- md_bull - md_bear
    d <- abs(md_delta - as.numeric(md_delta[rebdate, ]))^2 
    eps <- density(d)$bw
    p <- exp( -1 / (2 * eps) * d )
    p <- p / sum(p)
    n_sim <- ifelse( is.null(spec$n_sim), 10^3, spec$n_sim )
    P <- rdirichlet(n = n_sim, alpha = p )
    y_train <- as.numeric(Y_train)
    es <- apply( P, 1, function(p) { t(p) %*% y_train })
    es_dens <- densFUN(es)
    lambda <- entropy( es_dens$y, exponential = TRUE) / length(es_dens$y)
    wghts <- seq(from = 0, to = 1, length.out = 100)
    w_mat <- matrix( wghts, ncol = 1 )
    U <- apply( w_mat, 1, function(w) { ktUtility( mu = es * w, RP = 0 ) } )
    EU <- apply(U, 2, mean)
    ER <- rev(EU) * (-1) + max(EU)
    
    
    
    # matplot( x = wghts, y = t(U_bull), type = "l" )
    # lines( x = wghts, y = apply(U, 2, mean), lwd = 3, col = "orange" )
    # 
    # plot( x = wghts, y = EU - ER, ylim = range(c(EU, ER, EU - ER)),
    #       main = rebdate )
    # abline(h = 0)
    # points( x = wghts, y = EU, col = 3 )
    # points( x = wghts, y = ER, col = 2 )
    # points( x = wghts, y = EU - ER * lambda, col = 4 )
    
    
    
    w_star <- which( EU - ER * lambda == max(EU - ER * lambda) )[1] / 100
    w_star_nr <- which( EU == max(EU) )[1] / 100
    
    mat <- matrix( c(w_star, w_star_nr), nrow = 1,
                   dimnames = list(rebdate, c("w_star", "w_star_nr")) )
    weights <- list( weights = as.timeSeries(mat),
                     lambda = timeSeries(matrix(lambda), rebdate),
                     exp_util = timeSeries(matrix(EU, nrow = 1), rebdate),
                     exp_regret = timeSeries(matrix(EU, nrow = 1), rebdate) )
    
    .self$appendOutput( value = weights )
    
    return( TRUE )
  }
  
  
  
  # --------------------------------------------------------------------------
  Backtest.rpmAsyTurb_2 <- function( rebdate )
  {
    
    
    Y_train <- data$TD$Y_train[rownames(data$TD$Y_train) <= rebdate, ]
    X_train <- data$TD$X_train[rownames(data$TD$X_train) <= rebdate, ]
    X_eval <- data$TD$X_train[rebdate, ]
    
    # Compute mahalanobis distance in signal space seperately 
    # for bear and bull specific coefficients.
    
    p_bear <- data$p_bear[rownames(X_train), ] / sum(data$p_bear[rownames(X_train), ])
    mu_bear <- t(X_train) %*% p_bear
    scnd_mom = t(X_train) %*% (X_train * (p_bear %*% matrix(1, 1, ncol(X_train))) )
    scnd_mom = (scnd_mom + t(scnd_mom)) / 2
    sigma_bear = scnd_mom - mu_bear %*% t(mu_bear)
    # md_bear <- mahalanobis.ts(X = X_train, # scale(X_train, mu_bear, TRUE),
    #                           center = mu_bear, # X_eval,
    #                           covmat = sigma_bear)
    # eps_bear <- density(md_bear)$bw
    # w_bear <- exp( -1 / (2 * eps_bear) * md_bear )
    # w_bear <- w_bear / sum(w_bear)
    #
    md_bear <- mahalanobis.ts(X = scale(X_train, mu_bear, TRUE),
                              center = X_eval,
                              covmat = sigma_bear)
    eps_bear <- density(md_bear)$bw
    w_bear <- exp( -1 / (2 * eps_bear) * md_bear )
    w_bear <- w_bear / sum(w_bear)
    
    
    
    p_bull <- data$p_bull[rownames(X_train), ] / sum(data$p_bull[rownames(X_train), ])
    mu_bull <- t(X_train) %*% p_bull
    scnd_mom = t(X_train) %*% (X_train * (p_bull %*% matrix(1, 1, ncol(X_train))) )
    scnd_mom = (scnd_mom + t(scnd_mom)) / 2
    sigma_bull = scnd_mom - mu_bull %*% t(mu_bull)
    # md_bull <- mahalanobis.ts(X = X_train, # scale(X_train, mu_bull, TRUE),
    #                           center = mu_bull, # X_eval,
    #                           covmat = sigma_bull)
    # eps_bull <- density(md_bull)$bw
    # w_bull <- exp( -1 / (2 * eps_bull) * md_bull )
    # w_bull <- w_bull / sum(w_bull)
    # #
    md_bull <- mahalanobis.ts(X = scale(X_train, mu_bull, TRUE),
                              center = X_eval,
                              covmat = sigma_bull)
    eps_bull <- density(md_bull)$bw
    w_bull <- exp( -1 / (2 * eps_bull) * md_bull )
    w_bull <- w_bull / sum(w_bull)
    
    
    md_all <- mahalanobis.ts(X = scale(X_train, TRUE, TRUE),
                             center = X_eval)
    eps_all <- density(md_all)$bw
    w_all <- exp( -1 / (2 * eps_all) * md_all )
    w_all <- w_all / sum(w_all)
    
    
    
    plot( cbind(md_all, md_bear, md_bull), plot.type = "multiple" )
    plot( cbind(md_all, md_bear, md_bull), plot.type = "single" )
    plot( cbind(w_all, w_bear, w_bull), plot.type = "single" )
    
    
    
    IS <- importanceSampling(Y_train = Y_train,
                             X_train = md_bear)
    plot.is( IS )
    abline(v = mean(Y_train))
    
    
    
    
    n_sim <- ifelse( is.null(spec$n_sim), 10^3, spec$n_sim )
    P_all <- rdirichlet(n = n_sim, alpha = w_all )
    P_bear <- rdirichlet(n = n_sim, alpha = w_bear )
    P_bull <- rdirichlet(n = n_sim, alpha = w_bull )
    
    # P <- rdirichlet(n = n_sim, alpha = p * length(p) )
    
    # Expected states as RPM weighted average of state series
    y_train <- as.numeric(Y_train)
    es_all <- apply( P_all, 1, function(p) { t(p) %*% y_train })
    es_bear <- apply( P_bear, 1, function(p) { t(p) %*% y_train })
    es_bull <- apply( P_bull, 1, function(p) { t(p) %*% y_train })
    
    apply( cbind(es_all, es_bear, es_bull), 2, mean )
    
    
    es_all_dens <- densFUN(es_all)
    es_bear_dens <- densFUN(es_bear)
    es_bull_dens <- densFUN(es_bull)
    
    lambda_all <- entropy( es_all_dens$y, exponential = TRUE) / length(es_all_dens$y)
    lambda_bear <- entropy( es_bear_dens$y, exponential = TRUE) / length(es_bear_dens$y)
    lambda_bull <- entropy( es_bull_dens$y, exponential = TRUE) / length(es_bull_dens$y)
    
    lambda_all
    lambda_bear
    lambda_bull
    
    
    plot.ldensity( list(all = es_all_dens,
                        bear = es_bear_dens,
                        bull = es_bull_dens) )
    
    
    
    # Utility of expected states * portfolio weights
    wghts <- seq(from = 0, to = 1, length.out = 100)
    w_mat <- matrix( wghts, ncol = 1 )
    
    U_all <- apply( w_mat, 1, function(w) { ktUtility( mu = es_all * w, RP = 0 ) } )
    U_bear <- apply( w_mat, 1, function(w) { ktUtility( mu = es_bear * w, RP = 0 ) } )
    U_bull <- apply( w_mat, 1, function(w) { ktUtility( mu = es_bull * w, RP = 0 ) } )
    
    
    # U <- apply( w_mat, 1, function(w) { mvUtility( x = es * w ) } )
    EU_all <- apply(U_all, 2, mean)
    EU_bear <- apply(U_bear, 2, mean)
    EU_bull <- apply(U_bull, 2, mean)
    
    # Define regret as the reversed inverted empirical expected utility
    ER_all <- rev(EU_all) * (-1) + max(EU_all)
    ER_bear <- rev(EU_bear) * (-1) + max(EU_bear)
    ER_bull <- rev(EU_bull) * (-1) + max(EU_bull)
    
    
    # matplot( x = wghts, y = cbind(EU_all, EU_bear, EU_bull) )
    
    
    md_delta <- md_bull - md_bear
    d <- abs(md_delta - as.numeric(md_delta[rebdate, ]))^2 
    eps <- density(d)$bw
    p <- exp( -1 / (2 * eps) * d )
    p <- p / sum(p)
    n_sim <- ifelse( is.null(spec$n_sim), 10^3, spec$n_sim )
    P <- rdirichlet(n = n_sim, alpha = p )
    y_train <- as.numeric(Y_train)
    es <- apply( P, 1, function(p) { t(p) %*% y_train })
    es_dens <- densFUN(es)
    lambda <- entropy( es_dens$y, exponential = TRUE) / length(es_dens$y)
    wghts <- seq(from = 0, to = 1, length.out = 100)
    w_mat <- matrix( wghts, ncol = 1 )
    U <- apply( w_mat, 1, function(w) { ktUtility( mu = es * w, RP = 0 ) } )
    EU <- apply(U, 2, mean)
    ER <- rev(EU) * (-1) + max(EU)
    
    
    
    # matplot( x = wghts, y = t(U_bull), type = "l" )
    # lines( x = wghts, y = apply(U, 2, mean), lwd = 3, col = "orange" )
    # 
    # plot( x = wghts, y = EU - ER, ylim = range(c(EU, ER, EU - ER)),
    #       main = rebdate )
    # abline(h = 0)
    # points( x = wghts, y = EU, col = 3 )
    # points( x = wghts, y = ER, col = 2 )
    # points( x = wghts, y = EU - ER * lambda, col = 4 )
    
    
    
    w_star <- which( EU - ER * lambda == max(EU - ER * lambda) )[1] / 100
    w_star_nr <- which( EU == max(EU) )[1] / 100
    
    mat <- matrix( c(w_star, w_star_nr), nrow = 1,
                   dimnames = list(rebdate, c("w_star", "w_star_nr")) )
    weights <- list( weights = as.timeSeries(mat),
                     lambda = timeSeries(matrix(lambda), rebdate),
                     exp_util = timeSeries(matrix(EU, nrow = 1), rebdate),
                     exp_regret = timeSeries(matrix(EU, nrow = 1), rebdate) )
    
    .self$appendOutput( value = weights )
    
    return( TRUE )
  }
  
  
  
  # --------------------------------------------------------------------------
  BacktestDAA.mom_lo <- function( rebdate = NULL )
  {
    
    X <- data$X[(rownames(data$X) <= rebdate), "BM"]
    mom_spec <- momCtrl(method = "cumretEwma",
                        compounding = "discrete",
                        scale = NULL,
                        width = spec$mom_ctrl_width,
                        lag = 0,
                        zscore_flag = FALSE,
                        ellipsis = list(tau = spec$mom_ctrl_tau, stdz = FALSE))
    mu <- momentum(Data = X,
                   spec = mom_spec)
    
    lOut <- list( mu = as.timeSeries(matrix(mu, dimnames = list(rebdate, "mu"))) )
    
    .self$appendOutput( value = lOut )
    
    return( TRUE )
  }
  
  
  
  # --------------------------------------------------------------------------
  BacktestDAA.mom_ls <- function( rebdate = NULL )
  {
    
    X <- data$X[(rownames(data$X) <= rebdate), ]
    mom_spec <- momCtrl(method = "cumretEwma",
                        compounding = "discrete",
                        scale = NULL,
                        width = spec$mom_ctrl_width,
                        lag = 0,
                        zscore_flag = FALSE,
                        ellipsis = list(tau = spec$mom_ctrl_tau, stdz = FALSE))
    mu <- momentum(Data = X,
                   spec = mom_spec)
    
    lOut <- list( mu = as.timeSeries(matrix(mu, dimnames = list(rebdate, "mu"))) )
    
    .self$appendOutput( value = lOut )
    
    return( TRUE )
    
  }
  
  
  # --------------------------------------------------------------------------
  BacktestDAA.momvariancePortfolio <- function( rebdate = NULL )
  {
    
    # Return series
    X <- data$X[(rownames(data$X) <= rebdate), ]
    
    # Exponentially weighted mean and variance
    tau <- spec$mom_ctrl_tau
    lambda <- exp(-log(2) / tau)
    i <- (0:(nrow(X) - 1))
    wt <- lambda^i
    p_vec <- rev( wt / sum(wt) )
    mu = t(X) %*% p_vec
    scnd_mom = t(X) %*% (X * (p_vec %*% matrix(1, 1, ncol(X))) )
    scnd_mom = (scnd_mom + t(scnd_mom)) / 2
    covmat = scnd_mom - mu %*% t(mu)
    
    # Get initial weights
    # w_init <- .self$initialWeights( rebdate = rebdate )
    # w_init <- .self$floatingWeights( rebdate = rebdate )
    w_init <- NULL
    
    # Run mean-variance optimization
    alpha <- ifelse( is.null(spec$alpha), 0.05, spec$alpha )
    riskaversion <- ifelse( is.null(spec$riskaversion), 1, spec$riskaversion)
    
    GPS <- GPO::gps( Data = X,
                     Constraints = spec$Constraints,
                     Covariance = covmat,
                     Solver = solverCtrl( portfolio = "meanvariance",
                                          recursion = TRUE,
                                          recur_method = "mvtob_tree",
                                          w_init = w_init,
                                          utility = list(riskaversion = riskaversion,
                                                         alpha = alpha),
                                          obj_lin = mu * (-1) ) )
    if ( !is.null(w_init) )  {
      addConstraint(GPS) <- turnoverConstraint( rhs = Inf,
                                                w_init = w_init,
                                                linearize = TRUE )
    }
    GPO <- gpo(GPS = GPS)
    wghts <- getWeights(GPO)
    
    lOut <- list( weights = as.timeSeries( matrix(wghts, nrow = 1, 
                                                  dimnames = list(rebdate, names(wghts)))) )
    
    .self$appendOutput( value = lOut )
    
    return( TRUE )
  }
  
  
  
  
  # --------------------------------------------------------------------------
  BacktestDAA.runDAA <- function( rebdates = NULL, 
                                  method = "rpm" )
  {
    if ( is.null(rebdates) ) {
      rebdates <- spec$rebdates
    }
    
    # Check if some weights already exist.
    # If so, run backtest only for missing dates.
    dates <- rownames(output$weights) 
    if ( is.null(dates) ) {
      missing_dates <- rebdates
    } else {
      idx <- which(rebdates > dates[length(dates)])
      if ( length(idx) > 0 ) {
        missing_dates <- rebdates[idx]
      } else {
        if ( isTRUE(spec$verbose) ) {
          cat("Backtest is already up-to-date.\n")
        }
        return( TRUE )
      }
    }
    for ( today in missing_dates ) {
      
      if ( isTRUE(spec$verbose) ) {
        cat("Rebalancing date: ", today, "\n")
      }
      txt <- paste0( ".self$", method, "(rebdate = '", today, "')" )
      eval( parse( text = txt) )
      
      if ( isTRUE(spec$verbose) ) {
        cat("\n")
      }
    }
    return( TRUE )
  }
  
  
