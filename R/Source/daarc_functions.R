  
  
  ############################################################################
  ### DAARC TEMPORARY FUNCTIONS
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     13.03.2021
  # First version:    13.03.2021
  # --------------------------------------------------------------------------
  
  

  # BBTSV.setCtrl
  # BBTSV.computeSignalInsample
  # BBTSV.computeSignal
  # BBTSV
  
  # BBTGARCH.setCtrl
  # BBTGARCH.computeSignalInsample
  # BBTGARCH
  
  # BBTEWMA.setCtrl
  # BBTEWMA.computeSignalInsample
  # BBTEWMA
  
  
  # --------------------------------------------------------------------------
  #' @export
  # --------------------------------------------------------------------------
  BBTSV.setCtrl <- function( universe, method, verbose = FALSE )
  {
    BBT <- BBTurbulence$new()
    BBT$setCtrl( universe = universe,
                 method = method,
                 verbose = verbose )
    spec <<- BBT$spec
    spec$name <<- paste0( "bbtsv_", spec$method, "_", spec$universe )
    return( TRUE )
  }
  
 
  # --------------------------------------------------------------------------
  #' @export
  # --------------------------------------------------------------------------
  BBTSV.computeSignalInsample <- function( )
  {
    
    if ( isTRUE(spec$verbose) ) {
      cat( "computing", spec$name, "signal insample...\n")
    }
    
    # Bears n' bulls
    BBSObj <- spec$BBS$copy()
    BBSObj$data$X_level <- cumulated(data$X_bm, "discrete")
    if ( isTRUE(spec$robust) ) {
      BBSObj$runRobust( minphase_vec = spec$minphase_vec,
                        mincycle_vec = spec$mincycle_vec )
    } else {
      BBSObj$run()
    }
    data$BBS <<- BBSObj
    bbstates <- na.omit( BBSObj$output$states )
    # # Omit last phase
    # bbstates <- bbstates[1:(nrow(bbstates)-spec$minphase), ]
    
    dates <- intersect(rownames(data$X), rownames(bbstates))
    X_train <- data$X[dates, ]
    bbstates <- bbstates[dates, ]
    
    # Compute regime specific stats
    
    # All
    covmat <- cov(X_train)
    
    # Bear
    p_bear <- (1 - bbstates) / 2
    stats_bear <- cov.wt(x = X_train, 
                         wt = as.numeric( p_bear / sum(p_bear) ),
                         cor = FALSE, 
                         center = TRUE, 
                         method = "ML")
    mu_bear <- stats_bear$center
    covmat_bear <- make.positive.definite( stats_bear$cov )
    
    # Bull
    p_bull <- (1 + bbstates) / 2
    stats_bull <- cov.wt(x = X_train, 
                         wt = as.numeric( p_bull / sum(p_bull) ),
                         cor = FALSE, 
                         center = TRUE, 
                         method = "ML")
    mu_bull <- stats_bull$center
    covmat_bull <- make.positive.definite( stats_bull$cov )
    
    # Minimum torsion series
    S_all <- getter( mtca(X_train), "sources" )
    S_bear <- timeSeries( scale(X_train, mu_bear, FALSE) %*% 
                            torsion(covmat_bear, method = "mt"), time(X_train) )
    S_bull <- timeSeries( scale(X_train, mu_bull, FALSE) %*% 
                            torsion(covmat_bull, method = "mt"), time(X_train) )
    
    # Weighted minimum torsion series
    if ( !is.null(data$wmat) ) {
      wmat <- data$wmat[rownames(X_train), ]
      S_all_wght <- timeSeries( (scale(X_train, TRUE, FALSE) * wmat) %*% 
                                  getter(MT, "torsion"), time(X_train) )
      S_bear_wght <- timeSeries( (scale(X_train, mu_bear, FALSE) * wmat) %*% 
                                   torsion(covmat_bear, method = "mt"), time(X_train) )
      S_bull_wght <- timeSeries( (scale(X_train, mu_bull, FALSE) * wmat) %*% 
                                   torsion(covmat_bull, method = "mt"), time(X_train) )
    } else {
      S_all_wght <- S_all
      S_bear_wght <- S_bear
      S_bull_wght <- S_bull
    }
    
    # Stochastic volatility MD in MT space
    sv_all <- sv_bear <- sv_bull <- S_all * NA
    n_sim <- 1000
    sv_array_bear <- array( NA, dim = c(nrow(S_all), ncol(S_all), n_sim) )
    sv_array_bull <- sv_array_bear
    for ( j in 1:ncol(S_all) ) {
      
      # draws <- NULL
      # draws <- svsample( S_all[ ,j], draws = 1000, burnin = 100 )
      # sv_all[ ,j] <- apply( exp(draws$latent), 2, mean )
      
      draws <- NULL
      draws <- svsample( S_bear[ ,j], draws = n_sim, burnin = 100, quiet = TRUE  )
      h_bear <- t(exp(draws$latent))
      sv_bear[ ,j] <- apply( h_bear, 1, mean )
      sv_array_bear[ ,j, ] <- scale( h_bear, FALSE, rep(var(S_bear[ ,j]), n_sim) )
      
      draws <- NULL
      draws <- svsample( S_bull[ ,j], draws = n_sim, burnin = 100, quiet = TRUE  )
      h_bull <- t(exp(draws$latent))
      sv_bull[ ,j] <- apply( h_bull, 1, mean )
      sv_array_bull[ ,j, ] <- scale( h_bull, FALSE, rep(var(S_bull[ ,j]), n_sim) ) 
    }
    
    md_bear_mat <- do.call( cbind, lapply(1:n_sim, FUN = function(k) { apply( sv_array_bear[ ,,k], 1, sum ) }) )
    md_bear_mat <- timeSeries( md_bear_mat, time(sv_bear) )
    md_bear_mat_scl <- scale( md_bear_mat, FALSE, TRUE )
    
    md_bull_mat <- do.call( cbind, lapply(1:n_sim, FUN = function(k) { apply( sv_array_bull[ ,,k], 1, sum ) }) )
    md_bull_mat <- timeSeries( md_bull_mat, time(sv_bull) )
    md_bull_mat_scl <- scale( md_bull_mat, FALSE, TRUE )
    
    # md_all <- apply( scale(sv_all, FALSE, apply(S_all, 2, var)), 1, sum)
    md_bear <- apply( scale(sv_bear, FALSE, apply(S_bear, 2, var)), 1, sum)
    md_bear <- timeSeries( md_bear, time(sv_bear) )
    md_bull <- apply( scale(sv_bull, FALSE, apply(S_bull, 2, var)), 1, sum)
    md_bull <- timeSeries( md_bull, time(sv_bull) )
    
    MD <- cbind( bear = md_bear, 
                 bull = md_bull )
    MD_scl <- scale( MD, FALSE, TRUE )
    md_delta <- MD_scl[ ,"bear"] - MD_scl[ ,"bull"]
    
    
    # k <- n_sim
    # ldens <- list( bear = density(md_bear_mat_scl[k, ]),
    #                bull = density(md_bull_mat_scl[k, ]) )
    # plot.ldensity( ldens )
    # abline( v = MD_scl[k, "bull"], col = "darkgreen" )
    # abline( v = MD_scl[k, "bear"], col = "darkred" )
    # 
    # 
    # plot( cbind(md_bear_mat[ ,k], md_bull_mat[ ,k]), plot.type = "single" )
    # plot( cbind(md_bear, md_bull), plot.type = "single" )
    
    
    
    # plot( density(md_delta_mat) )
    # plot( apply( md_delta_mat, 1, mean) )
    # lines( md_delta, col = 2 )
    
    # ldelta <- lapply( 1:10^4, 
    #                  FUN = function(k) { md_bear_mat_scl[ ,round(runif(1, 1, n_sim))] - md_bull_mat_scl[ ,round(runif(1, 1, n_sim))] } )
    # md_delta_mat <- do.call( cbind, ldelta )
    # prob <- apply( md_delta_mat, 1, function(x) { sum( x > 0) / ncol(md_delta_mat) } )
    md_delta_mat <- md_bear_mat_scl - md_bull_mat_scl
    prob <- apply( md_delta_mat, 1, function(x) { sum( x > 0) / n_sim } )
    
    # plot( prob )
    
    
    ans <- cbind( states = bbstates,
                  bear = md_bear,
                  bull = md_bull,
                  bear_scl = MD_scl[ ,"bear"],
                  bull_scl = MD_scl[ ,"bull"],
                  signal = md_delta,
                  prob = prob )
    colnames(ans) <-  c("states", "bear", "bull","bear_scl", "bull_scl", "delta", "prob")
    signal$insample <<- ans
    
    if ( isTRUE(spec$verbose) ) {
      cat( " done.\n")
    }
    
    return( TRUE )    
  }
  
  
  # --------------------------------------------------------------------------
  #' @export
  # --------------------------------------------------------------------------
  BBTSV.computeSignal <- function( dates = NULL )
  {
   
    stopifnot( is.character(dates) )
    
    # Initialize outputs
    md_all <- data$X_bm[dates, ] * NA
    md_bear <- md_bull <- md_delta <- md_all
    md_all_scl <- md_bear_scl <- md_bull_scl <- md_all
    prob <- md_all
    states <- md_all
    
    # Append signals on exit
    on.exit( { 
      lSignal <- list( na.omit( cbind( md_bear = md_bear,
                                       md_bull = md_bull,
                                       md_bear_scl = md_bear_scl,
                                       md_bull_scl = md_bull_scl,
                                       signal = md_delta,
                                       prob = prob,
                                       states = states ) ) )
      names(lSignal) <- spec$method
      .self$appendSignal( value = lSignal )
    } )
    
    # Loop over dates
    for ( today in dates ) {
      
      if ( isTRUE(spec$verbose) ) {
        cat( "Running BBTSV.computeSignal: date =", today, "\n" )
      }
      
      # Bears n' bulls
      BBSObj <- spec$BBS$copy()
      BBSObj$data$X_level <- cumulated(data$X_bm[(rownames(data$X_bm) <= today), ], "discrete")
      if ( isTRUE(spec$robust) ) {
        BBSObj$runRobust( minphase_vec = spec$minphase_vec,
                          mincycle_vec = spec$mincycle_vec )
      } else {
        BBSObj$run()
      }
      data$BBS <<- BBSObj
      bbstates <- na.omit( BBSObj$output$states )
      # # Omit last phase
      # bbstates <- bbstates[1:(nrow(bbstates)-spec$minphase), ]
      
      X_train <- data$X[rownames(data$X) <= today, ]
      bbstates <- bbstates[rownames(X_train), ]
      
      
      # Compute regime specific stats
      
      # All
      covmat <- cov(X_train)
      
      # Bear
      p_bear <- (1 - bbstates) / 2
      stats_bear <- cov.wt(x = X_train, 
                           wt = as.numeric( p_bear / sum(p_bear) ),
                           cor = FALSE, 
                           center = TRUE, 
                           method = "ML")
      mu_bear <- stats_bear$center
      covmat_bear <- make.positive.definite( stats_bear$cov )
      
      # Bull
      p_bull <- (1 + bbstates) / 2
      stats_bull <- cov.wt(x = X_train, 
                           wt = as.numeric( p_bull / sum(p_bull) ),
                           cor = FALSE, 
                           center = TRUE, 
                           method = "ML")
      mu_bull <- stats_bull$center
      covmat_bull <- make.positive.definite( stats_bull$cov )
      
      # Minimum torsion series
      S_all <- getter( mtca(X_train), "sources" )
      S_bear <- timeSeries( scale(X_train, mu_bear, FALSE) %*% 
                              torsion(covmat_bear, method = "mt"), time(X_train) )
      S_bull <- timeSeries( scale(X_train, mu_bull, FALSE) %*% 
                              torsion(covmat_bull, method = "mt"), time(X_train) )
      
      # Weighted minimum torsion series
      if ( !is.null(data$wmat) ) {
        wmat <- data$wmat[rownames(X_train), ]
        S_all_wght <- timeSeries( (scale(X_train, TRUE, FALSE) * wmat) %*% 
                                    getter(MT, "torsion"), time(X_train) )
        S_bear_wght <- timeSeries( (scale(X_train, mu_bear, FALSE) * wmat) %*% 
                                     torsion(covmat_bear, method = "mt"), time(X_train) )
        S_bull_wght <- timeSeries( (scale(X_train, mu_bull, FALSE) * wmat) %*% 
                                     torsion(covmat_bull, method = "mt"), time(X_train) )
      } else {
        S_all_wght <- S_all
        S_bear_wght <- S_bear
        S_bull_wght <- S_bull
      }
      
      # Stochastic volatility MD in MT space
      n_sim <- 1000
      sv_all <- sv_bear <- sv_bull <- S_all * NA
      sv_array_bear <- array( NA, dim = c(nrow(S_all), ncol(S_all), n_sim) )
      sv_array_bull <- sv_array_bear
      for ( j in 1:ncol(S_all) ) {
        
        # draws <- NULL
        # draws <- svsample( S_all[ ,j], draws = 1000, burnin = 100 )
        # sv_all[ ,j] <- apply( exp(draws$latent), 2, mean )
        
        draws <- NULL
        draws <- svsample( S_bear[ ,j], draws = n_sim, burnin = 100, quiet = TRUE  )
        h_bear <- t(exp(draws$latent))
        sv_bear[ ,j] <- apply( h_bear, 1, mean )
        sv_array_bear[ ,j, ] <- scale( h_bear, FALSE, rep(var(S_bear[ ,j]), n_sim) )
        
        draws <- NULL
        draws <- svsample( S_bull[ ,j], draws = n_sim, burnin = 100, quiet = TRUE  )
        h_bull <- t(exp(draws$latent))
        sv_bull[ ,j] <- apply( h_bull, 1, mean )
        sv_array_bull[ ,j, ] <- scale( h_bull, FALSE, rep(var(S_bull[ ,j]), n_sim) ) 
      }
      
      md_bear_mat <- do.call( cbind, lapply(1:n_sim, FUN = function(k) { apply( sv_array_bear[ ,,k], 1, sum ) }) )
      md_bear_mat <- timeSeries( md_bear_mat, time(sv_bear) )
      md_bear_mat_scl <- scale( md_bear_mat, FALSE, TRUE )
      
      md_bull_mat <- do.call( cbind, lapply(1:n_sim, FUN = function(k) { apply( sv_array_bull[ ,,k], 1, sum ) }) )
      md_bull_mat <- timeSeries( md_bull_mat, time(sv_bull) )
      md_bull_mat_scl <- scale( md_bull_mat, FALSE, TRUE )
      
      md_delta_mat <- md_bear_mat_scl - md_bull_mat_scl
      prob_tmp <- apply( md_delta_mat, 1, function(x) { sum( x > 0) / n_sim } )
      
      
      # md_all_tmp <- apply( scale(sv_all, FALSE, apply(S_all, 2, var)), 1, sum)
      md_bear_tmp <- apply( scale(sv_bear, FALSE, apply(S_bear, 2, var)), 1, sum)
      md_bear_tmp <- timeSeries( md_bear_tmp, time(sv_bear) )
      md_bull_tmp <- apply( scale(sv_bull, FALSE, apply(S_bull, 2, var)), 1, sum)
      md_bull_tmp <- timeSeries( md_bull_tmp, time(sv_bull) )
      
      MD <- cbind( bear = md_bear_tmp, 
                   bull = md_bull_tmp )
      MD_scl <- scale( MD, FALSE, TRUE )
      md_delta_tmp <- MD_scl[ ,"bear"] - MD_scl[ ,"bull"]
      
      # Output
      md_bear[today, ] <- md_bear_tmp[today, ]
      md_bull[today, ] <- md_bull_tmp[today, ]
      md_bear_scl[today, ] <- MD_scl[today, "bear"]
      md_bull_scl[today, ] <- MD_scl[today, "bull"]
      md_delta[today, ] <- md_delta_tmp[today, ]
      prob[today, ] <- prob_tmp[today, ]
      states[today, ] <- bbstates[today, ]
      
    }
    
    return( TRUE ) 
  }
    
    
  
  # --------------------------------------------------------------------------
  #' @title Reference Class BBTSV
  #' @description setRefClass call for object of class BBTSV.
  #' @import methods
  #' @include class_sensor.R
  #' @export BBTSV
  # --------------------------------------------------------------------------
  BBTSV <- setRefClass( Class = "BBTSV", 
                        contains = "BBTurbulence",
                        methods = list( setCtrl = BBTSV.setCtrl,
                                        computeSignal = BBTSV.computeSignal,
                                        computeSignalInsample = BBTSV.computeSignalInsample) )
  
  
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  #' @export
  # --------------------------------------------------------------------------
  BBTGARCH.setCtrl <- function( universe, method, verbose = FALSE )
  {
    BBT <- BBTurbulence$new()
    BBT$setCtrl( universe = universe,
                 method = method,
                 verbose = verbose )
    spec <<- BBT$spec
    spec$name <<- paste0( "bbtgarch_", spec$method, "_", spec$universe )
    return( TRUE )
  }
  
  
  # --------------------------------------------------------------------------
  #' @export
  # --------------------------------------------------------------------------
  BBTGARCH.computeSignalInsample <- function( )
  {
    
    if ( isTRUE(spec$verbose) ) {
      cat( "computing", spec$name, "signal insample...\n")
    }
    
    # Bears n' bulls
    BBSObj <- spec$BBS$copy()
    BBSObj$data$X_level <- cumulated(data$X_bm, "discrete")
    if ( isTRUE(spec$robust) ) {
      BBSObj$runRobust( minphase_vec = spec$minphase_vec,
                        mincycle_vec = spec$mincycle_vec )
    } else {
      BBSObj$run()
    }
    data$BBS <<- BBSObj
    bbstates <- na.omit( BBSObj$output$states )
    # # Omit last phase
    # bbstates <- bbstates[1:(nrow(bbstates)-spec$minphase), ]
    
    dates <- intersect(rownames(data$X), rownames(bbstates))
    X_train <- data$X[dates, ]
    bbstates <- bbstates[dates, ]
    
    # Compute regime specific stats
    
    # All
    covmat <- cov(X_train)
    
    # Bear
    p_bear <- (1 - bbstates) / 2
    stats_bear <- cov.wt(x = X_train, 
                         wt = as.numeric( p_bear / sum(p_bear) ),
                         cor = FALSE, 
                         center = TRUE, 
                         method = "ML")
    mu_bear <- stats_bear$center
    covmat_bear <- make.positive.definite( stats_bear$cov )
    
    # Bull
    p_bull <- (1 + bbstates) / 2
    stats_bull <- cov.wt(x = X_train, 
                         wt = as.numeric( p_bull / sum(p_bull) ),
                         cor = FALSE, 
                         center = TRUE, 
                         method = "ML")
    mu_bull <- stats_bull$center
    covmat_bull <- make.positive.definite( stats_bull$cov )
    
    # Minimum torsion series
    S_all <- getter( mtca(X_train), "sources" )
    S_bear <- timeSeries( scale(X_train, mu_bear, FALSE) %*% 
                            torsion(covmat_bear, method = "mt"), time(X_train) )
    S_bull <- timeSeries( scale(X_train, mu_bull, FALSE) %*% 
                            torsion(covmat_bull, method = "mt"), time(X_train) )
    
    # Weighted minimum torsion series
    if ( !is.null(data$wmat) ) {
      wmat <- data$wmat[rownames(X_train), ]
      S_all_wght <- timeSeries( (scale(X_train, TRUE, FALSE) * wmat) %*% 
                                  getter(MT, "torsion"), time(X_train) )
      S_bear_wght <- timeSeries( (scale(X_train, mu_bear, FALSE) * wmat) %*% 
                                   torsion(covmat_bear, method = "mt"), time(X_train) )
      S_bull_wght <- timeSeries( (scale(X_train, mu_bull, FALSE) * wmat) %*% 
                                   torsion(covmat_bull, method = "mt"), time(X_train) )
    } else {
      S_all_wght <- S_all
      S_bear_wght <- S_bear
      S_bull_wght <- S_bull
    }
    
    # GARCH MD in MT space
    cvar_bear <- getCondVar( garch( S_bear ) )
    cvar_bull <- getCondVar( garch( S_bull ) ) 
    cvar_all <- getCondVar( garch( S_all ) )
    
    md_bear <- apply( scale(cvar_bear, FALSE, apply(S_bear, 2, var)), 1, sum )
    md_bear <- timeSeries( md_bear, time(cvar_bear) )
    md_bull <- apply( scale(cvar_bull, FALSE, apply(S_bull, 2, var)), 1, sum )
    md_bull <- timeSeries( md_bull, time(cvar_bull) )
    md_all <- apply( scale(cvar_all, FALSE, apply(S_all, 2, var)), 1, sum )
    md_all <- timeSeries( md_all, time(cvar_all) )
    
    MD <- cbind( bear = md_bear, 
                 bull = md_bull,
                 all = md_all )
    MD_scl <- scale( MD, FALSE, TRUE )
    md_delta <- MD_scl[ ,"bear"] - MD_scl[ ,"bull"]
    
    ans <- cbind( states = bbstates,
                  bear = md_bear,
                  bull = md_bull,
                  all = md_all,
                  bear_scl = MD_scl[ ,"bear"],
                  bull_scl = MD_scl[ ,"bull"],
                  all_scl = MD_scl[ ,"all"],
                  signal = md_delta )
    colnames(ans) <-  c("states", "bear", "bull", "all", "bear_scl", "bull_scl", "all_scl", "delta")
    signal$insample <<- ans
    
    if ( isTRUE(spec$verbose) ) {
      cat( " done.\n")
    }
    
    return( TRUE )
  }
  
  
  # --------------------------------------------------------------------------
  #' @export
  # --------------------------------------------------------------------------
  BBTGARCH.computeSignal <- function( dates = NULL )
  {
    if ( spec$garch_ctrl$steps > 0 ) {
      .self$computeSignalPrediction( dates = dates )
    } else {
      .self$computeSignalOOS( dates = dates )
    }
    return( TRUE )
  }
  
  # --------------------------------------------------------------------------
  #' @export
  # --------------------------------------------------------------------------
  BBTGARCH.computeSignalOOS <- function( dates = NULL )
  {
    
    stopifnot( is.character(dates) )
    
    # Initialize outputs
    md_all <- data$X_bm[dates, ] * NA
    md_bear <- md_bull <- md_delta <- md_all
    md_all_scl <- md_bear_scl <- md_bull_scl <- md_all
    states <- md_all
    
    # Append signals on exit
    on.exit( { 
      lSignal <- list( na.omit( cbind( md_bear = md_bear,
                                       md_bull = md_bull,
                                       md_bear_scl = md_bear_scl,
                                       md_bull_scl = md_bull_scl,
                                       signal = md_delta,
                                       states = states ) ) )
      names(lSignal) <- spec$method
      .self$appendSignal( value = lSignal )
    } )
    
    # Loop over dates
    for ( today in dates ) {
      
      if ( isTRUE(spec$verbose) ) {
        cat( "Running BBTGARCH.computeSignal: date =", today, "\n" )
      }
      
      # Bears n' bulls
      BBSObj <- spec$BBS$copy()
      BBSObj$data$X_level <- cumulated(data$X_bm[(rownames(data$X_bm) <= today), ], "discrete")
      if ( isTRUE(spec$robust) ) {
        BBSObj$runRobust( minphase_vec = spec$minphase_vec,
                          mincycle_vec = spec$mincycle_vec )
      } else {
        BBSObj$run()
      }
      data$BBS <<- BBSObj
      bbstates <- na.omit( BBSObj$output$states )
      # # Omit last phase
      # bbstates <- bbstates[1:(nrow(bbstates)-spec$minphase), ]
      
      # Remove series with more than 10% of zero returns
      idx <- apply(data$X[rownames(data$X) <= today, ], 2, function(x) { sum(x == 0) / length(x) <  0.1 } )
      X_train <- data$X[rownames(data$X) <= today, idx]
      bbstates <- bbstates[rownames(X_train), ]
      
      
      # Compute regime specific stats
      
      # All
      covmat <- cov(X_train)
      
      # Bear
      p_bear <- (1 - bbstates) / 2
      stats_bear <- cov.wt(x = X_train, 
                           wt = as.numeric( p_bear / sum(p_bear) ),
                           cor = FALSE, 
                           center = TRUE, 
                           method = "ML")
      mu_bear <- stats_bear$center
      covmat_bear <- make.positive.definite( stats_bear$cov )
      
      # Bull
      p_bull <- (1 + bbstates) / 2
      stats_bull <- cov.wt(x = X_train, 
                           wt = as.numeric( p_bull / sum(p_bull) ),
                           cor = FALSE, 
                           center = TRUE, 
                           method = "ML")
      mu_bull <- stats_bull$center
      covmat_bull <- make.positive.definite( stats_bull$cov )
      
      # Minimum torsion series
      S_all <- getter( mtca(X_train), "sources" )
      S_bear <- timeSeries( scale(X_train, mu_bear, FALSE) %*% 
                              torsion(covmat_bear, method = "mt"), time(X_train) )
      S_bull <- timeSeries( scale(X_train, mu_bull, FALSE) %*% 
                              torsion(covmat_bull, method = "mt"), time(X_train) )
      
      # Weighted minimum torsion series
      if ( !is.null(data$wmat) ) {
        wmat <- data$wmat[rownames(X_train), ]
        S_all_wght <- timeSeries( (scale(X_train, TRUE, FALSE) * wmat) %*% 
                                    getter(MT, "torsion"), time(X_train) )
        S_bear_wght <- timeSeries( (scale(X_train, mu_bear, FALSE) * wmat) %*% 
                                     torsion(covmat_bear, method = "mt"), time(X_train) )
        S_bull_wght <- timeSeries( (scale(X_train, mu_bull, FALSE) * wmat) %*% 
                                     torsion(covmat_bull, method = "mt"), time(X_train) )
      } else {
        S_all_wght <- S_all
        S_bear_wght <- S_bear
        S_bull_wght <- S_bull
      }
      
      # GARCH MD in MT space
      cvar_bear <- getCondVar( garch( S_bear ) )
      cvar_bull <- getCondVar( garch( S_bull ) ) 
      cvar_all <- getCondVar( garch( S_all ) )
      
      md_bear_tmp <- apply( scale(cvar_bear, FALSE, apply(S_bear, 2, var)), 1, sum)
      md_bear_tmp <- timeSeries( md_bear_tmp, time(cvar_bear) )
      md_bull_tmp <- apply( scale(cvar_bull, FALSE, apply(S_bull, 2, var)), 1, sum)
      md_bull_tmp <- timeSeries( md_bull_tmp, time(cvar_bull) )
      md_all_tmp <- apply( scale(cvar_all, FALSE, apply(S_all, 2, var)), 1, sum)
      md_all_tmp <- timeSeries( md_all_tmp, time(cvar_all) )
      
      MD <- cbind( bear = md_bear_tmp, 
                   bull = md_bull_tmp,
                   all = md_all_tmp )
      MD_scl <- scale( MD, FALSE, TRUE )
      md_delta_tmp <- MD_scl[ ,"bear"] - MD_scl[ ,"bull"]
      
      md_bear[today, ] <- md_bear_tmp[today, ]
      md_bull[today, ] <- md_bull_tmp[today, ]
      md_all[today, ] <- md_all_tmp[today, ]
      md_bear_scl[today, ] <- MD_scl[today, "bear"]
      md_bull_scl[today, ] <- MD_scl[today, "bull"]
      md_all_scl[today, ] <- MD_scl[today, "all"]
      md_delta[today, ] <- md_delta_tmp[today, ]
      states[today, ] <- bbstates[today, ]
      
    }

    return( TRUE ) 
  }
  
  
  # --------------------------------------------------------------------------
  #' @export
  # --------------------------------------------------------------------------
  BBTGARCH.computeSignalPrediction <- function( dates = NULL )
  {
    
    stopifnot( is.character(dates) )
    
    # Initialize outputs
    md_all <- data$X_bm[dates, rep(1, spec$garch_ctrl$steps + 1)] * NA
    md_bear <- md_bull <- md_delta <- md_all
    md_all_scl <- md_bear_scl <- md_bull_scl <- md_all
    states <- md_all[ ,1]
    
    # Append signals on exit
    on.exit( { 
      lSignal <- list()
      for ( j in 1:ncol(md_delta) ) {
        if ( j == 1 ) {
          lSignal[[j]] <- na.omit( cbind( md_bear = md_bear[ ,j],
                                          md_bull = md_bull[ ,j],
                                          md_bear_scl = md_bear_scl[ ,j],
                                          md_bull_scl = md_bull_scl[ ,j],
                                          signal = md_delta[ ,j],
                                          states = states ) )
        } else {
          lSignal[[j]] <- na.omit( cbind( md_bear_scl = md_bear_scl[ ,j],
                                          md_bull_scl = md_bull_scl[ ,j],
                                          signal = md_delta[ ,j] ) ) 
        }
      }
      names(lSignal) <- c( spec$method, paste0("pred", (2:ncol(md_delta))-1) )
      .self$appendSignal( value = lSignal )
    } )
    
    # Loop over dates
    for ( today in dates ) {
      
      if ( isTRUE(spec$verbose) ) {
        cat( "Running BBTGARCH.computeSignalPrediction: date =", today, "\n" )
      }
      
      # Bears n' bulls
      BBSObj <- spec$BBS$copy()
      BBSObj$data$X_level <- cumulated(data$X_bm[(rownames(data$X_bm) <= today), ], "discrete")
      if ( isTRUE(spec$robust) ) {
        BBSObj$runRobust( minphase_vec = spec$minphase_vec,
                          mincycle_vec = spec$mincycle_vec )
      } else {
        BBSObj$run()
      }
      data$BBS <<- BBSObj
      bbstates <- na.omit( BBSObj$output$states )
      # # Omit last phase
      # bbstates <- bbstates[1:(nrow(bbstates)-spec$minphase), ]
      
      # Remove series with more than 10% of zero returns
      idx <- apply(data$X[rownames(data$X) <= today, ], 2, function(x) { sum(x == 0) / length(x) <  0.1 } )
      X_train <- data$X[rownames(data$X) <= today, idx]
      bbstates <- bbstates[rownames(X_train), ]
      
      
      # Compute regime specific stats
      
      # All
      covmat <- cov(X_train)
      
      # Bear
      p_bear <- (1 - bbstates) / 2
      stats_bear <- cov.wt(x = X_train, 
                           wt = as.numeric( p_bear / sum(p_bear) ),
                           cor = FALSE, 
                           center = TRUE, 
                           method = "ML")
      mu_bear <- stats_bear$center
      covmat_bear <- make.positive.definite( stats_bear$cov )
      
      # Bull
      p_bull <- (1 + bbstates) / 2
      stats_bull <- cov.wt(x = X_train, 
                           wt = as.numeric( p_bull / sum(p_bull) ),
                           cor = FALSE, 
                           center = TRUE, 
                           method = "ML")
      mu_bull <- stats_bull$center
      covmat_bull <- make.positive.definite( stats_bull$cov )
      
      # Minimum torsion series
      S_all <- getter( mtca(X_train), "sources" )
      S_bear <- timeSeries( scale(X_train, mu_bear, FALSE) %*% 
                              torsion(covmat_bear, method = "mt"), time(X_train) )
      S_bull <- timeSeries( scale(X_train, mu_bull, FALSE) %*% 
                              torsion(covmat_bull, method = "mt"), time(X_train) )
      
      # Weighted minimum torsion series
      if ( !is.null(data$wmat) ) {
        wmat <- data$wmat[rownames(X_train), ]
        S_all_wght <- timeSeries( (scale(X_train, TRUE, FALSE) * wmat) %*% 
                                    getter(MT, "torsion"), time(X_train) )
        S_bear_wght <- timeSeries( (scale(X_train, mu_bear, FALSE) * wmat) %*% 
                                     torsion(covmat_bear, method = "mt"), time(X_train) )
        S_bull_wght <- timeSeries( (scale(X_train, mu_bull, FALSE) * wmat) %*% 
                                     torsion(covmat_bull, method = "mt"), time(X_train) )
      } else {
        S_all_wght <- S_all
        S_bear_wght <- S_bear
        S_bull_wght <- S_bull
      }
      
      # GARCH MD in MT space
      cvar_bear <- cvar_bull <- cvar_all <- S_all * NA
      cvar_fcast_bear <- matrix( NA, nrow = spec$garch_ctrl$steps, ncol = ncol(S_all) )
      cvar_fcast_bull <- cvar_fcast_all <- cvar_fcast_bear
      for ( j in 1:ncol(S_all) ) {
        fit_bear <- try( garch( S_bear[ ,j], spec$garch_ctrl ) )
        fit_bull <- try( garch( S_bull[ ,j], spec$garch_ctrl ) )
        fit_all <- try( garch( S_all[ ,j], spec$garch_ctrl ) )
        # if ( !inherits(fit_bear, "try-error") ) {
        if (  fit_bear@fit$convergence == 0 & fit_bull@fit$convergence == 0 ) {
          cvar_bear[ ,j] <- getCondVar( fit_bear )
          fcast_bear <- uGarchForecast( fitORspec = fit_bear )
          # cvar_fcast_bear[ ,j] <- getCondVar(fcast_bear)
          cvar_fcast_bear[ ,j] <- fcast_bear@forecast$sigmaFor^2
          cvar_bull[ ,j] <- getCondVar( fit_bull ) 
          fcast_bull <- uGarchForecast( fitORspec = fit_bull )
          # cvar_fcast_bull[ ,j] <- getCondVar(fcast_bull)
          cvar_fcast_bull[ ,j] <- fcast_bull@forecast$sigmaFor^2
        } else {
          cvar_bear[ ,j] <- ema( S_bear[ ,j], 0.1 )
          cvar_fcast_bear[ ,j] <- rep( as.numeric(cvar_bear[nrow(cvar_bear), j]), nrow(cvar_fcast_bear) )  
          cvar_bull[ ,j] <- ema( S_bull[ ,j], 0.1 )
          cvar_fcast_bull[ ,j] <- rep( as.numeric(cvar_bull[nrow(cvar_bull), j]), nrow(cvar_fcast_bull) )
        }
        # if ( !inherits(fit_all, "try-error" ) ) {
        if (  fit_all@fit$convergence == 0 ) {
          cvar_all[ ,j] <- getCondVar( fit_all )
          fcast_all <- uGarchForecast( fitORspec = fit_all )
          # cvar_fcast_all[ ,j] <- getCondVar(fcast_all)
          cvar_fcast_all[ ,j] <- fcast_all@forecast$sigmaFor^2
        } else {
          cvar_all[ ,j] <- ema( S_all[ ,j], 0.1 )
          cvar_fcast_all[ ,j] <- rep( as.numeric(cvar_all[nrow(cvar_all), j]), nrow(cvar_fcast_all) )
        }
        
      }
      
      # Scale insample mdistances
      md_bear_tmp <- apply( scale(cvar_bear, FALSE, apply(S_bear, 2, var)), 1, sum)
      md_bear_tmp <- timeSeries( md_bear_tmp, time(cvar_bear) )
      md_bull_tmp <- apply( scale(cvar_bull, FALSE, apply(S_bull, 2, var)), 1, sum)
      md_bull_tmp <- timeSeries( md_bull_tmp, time(cvar_bull) )
      md_all_tmp <- apply( scale(cvar_all, FALSE, apply(S_all, 2, var)), 1, sum)
      md_all_tmp <- timeSeries( md_all_tmp, time(cvar_all) )
      
      MD <- cbind( bear = md_bear_tmp, 
                   bull = md_bull_tmp,
                   all = md_all_tmp )
      MD_scl <- scale( MD, FALSE, TRUE )
      md_delta_tmp <- MD_scl[ ,"bear"] - MD_scl[ ,"bull"]
      
      # Scale forecasts
      md_bear_fcast_tmp <- apply( scale(cvar_fcast_bear, FALSE, apply(S_bear, 2, var)), 1, sum)
      md_bear_fcast_tmp <- as.timeSeries( md_bear_fcast_tmp )
      md_bull_fcast_tmp <- apply( scale(cvar_fcast_bull, FALSE, apply(S_bull, 2, var)), 1, sum)
      md_bull_fcast_tmp <- as.timeSeries( md_bull_fcast_tmp )
      md_all_fcast_tmp <- apply( scale(cvar_fcast_all, FALSE, apply(S_all, 2, var)), 1, sum)
      md_all_fcast_tmp <- as.timeSeries( md_all_fcast_tmp )
      
      MD_fcast <- cbind( bear = md_bear_fcast_tmp, 
                         bull = md_bull_fcast_tmp,
                         all = md_all_fcast_tmp )
      MD_fcast_scl <- scale( MD_fcast, FALSE, attr(MD_scl, "scaled:scale") )
      md_delta_fcast_tmp <- MD_fcast_scl[ ,"bear"] - MD_fcast_scl[ ,"bull"]
      
      
      md_bear[today, 1] <- md_bear_tmp[today, ]
      md_bull[today, 1] <- md_bull_tmp[today, ]
      md_all[today, 1] <- md_all_tmp[today, ]
      md_bear_scl[today, 1] <- MD_scl[today, "bear"]
      md_bull_scl[today, 1] <- MD_scl[today, "bull"]
      md_all_scl[today, 1] <- MD_scl[today, "all"]
      md_delta[today, 1] <- md_delta_tmp[today, ]
      md_bear_scl[today, -1] <- as.numeric(MD_fcast_scl[ ,"bear"])
      md_bull_scl[today, -1] <- as.numeric(MD_fcast_scl[ ,"bull"])
      md_all_scl[today, -1] <- as.numeric(MD_fcast_scl[ ,"all"])
      md_delta[today, -1] <- as.numeric(md_delta_fcast_tmp)
      states[today, ] <- bbstates[today, ]
    }
    
    return( TRUE ) 
  }
  
  
  
  
  
  # --------------------------------------------------------------------------
  #' @title Reference Class BBTGARCH
  #' @description setRefClass call for object of class BBTGARCH.
  #' @import methods
  #' @include class_sensor.R
  #' @export BBTGARCH
  # --------------------------------------------------------------------------
  BBTGARCH <- setRefClass( Class = "BBTGARCH", 
                        contains = "BBTurbulence",
                        methods = list( setCtrl = BBTGARCH.setCtrl,
                                        computeSignal = BBTGARCH.computeSignal,
                                        computeSignalOOS = BBTGARCH.computeSignalOOS,
                                        computeSignalPrediction = BBTGARCH.computeSignalPrediction,
                                        computeSignalInsample = BBTGARCH.computeSignalInsample) )
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  #' @export
  # --------------------------------------------------------------------------
  BBTEWMA.setCtrl <- function( universe, method, verbose = FALSE )
  {
    BBT <- BBTurbulence$new()
    BBT$setCtrl( universe = universe,
                 method = method,
                 verbose = verbose )
    spec <<- BBT$spec
    spec$name <<- paste0( "bbtewma_", spec$method, "_", spec$universe )
    return( TRUE )
  }
  
  
  # --------------------------------------------------------------------------
  #' @export
  # --------------------------------------------------------------------------
  BBTEWMA.computeSignalInsample <- function( )
  {
    
    if ( isTRUE(spec$verbose) ) {
      cat( "computing", spec$name, "signal insample...\n")
    }
    
    # Bears n' bulls
    BBSObj <- spec$BBS$copy()
    BBSObj$data$X_level <- cumulated(data$X_bm, "discrete")
    if ( isTRUE(spec$robust) ) {
      BBSObj$runRobust( minphase_vec = spec$minphase_vec,
                        mincycle_vec = spec$mincycle_vec )
    } else {
      BBSObj$run()
    }
    data$BBS <<- BBSObj
    bbstates <- na.omit( BBSObj$output$states )
    # # Omit last phase
    # bbstates <- bbstates[1:(nrow(bbstates)-spec$minphase), ]
    
    dates <- intersect(rownames(data$X), rownames(bbstates))
    X_train <- data$X[dates, ]
    bbstates <- bbstates[dates, ]
    
    # Compute regime specific stats
    
    # All
    covmat <- cov(X_train)
    
    # Bear
    p_bear <- (1 - bbstates) / 2
    stats_bear <- cov.wt(x = X_train, 
                         wt = as.numeric( p_bear / sum(p_bear) ),
                         cor = FALSE, 
                         center = TRUE, 
                         method = "ML")
    mu_bear <- stats_bear$center
    covmat_bear <- make.positive.definite( stats_bear$cov )
    
    # Bull
    p_bull <- (1 + bbstates) / 2
    stats_bull <- cov.wt(x = X_train, 
                         wt = as.numeric( p_bull / sum(p_bull) ),
                         cor = FALSE, 
                         center = TRUE, 
                         method = "ML")
    mu_bull <- stats_bull$center
    covmat_bull <- make.positive.definite( stats_bull$cov )
    
    # Minimum torsion series
    MT <-  mtca(X_train)
    S_all <- getter( MT, "sources" )
    S_bear <- timeSeries( scale(X_train, mu_bear, FALSE) %*% 
                            torsion(covmat_bear, method = "mt"), time(X_train) )
    S_bull <- timeSeries( scale(X_train, mu_bull, FALSE) %*% 
                            torsion(covmat_bull, method = "mt"), time(X_train) )
    
    # Weighted minimum torsion series
    if ( !is.null(data$wmat) ) {
      wmat <- data$wmat[rownames(X_train), ]
      S_all_wght <- timeSeries( (scale(X_train, TRUE, FALSE) * wmat) %*% 
                                  getter(MT, "torsion"), time(X_train) )
      S_bear_wght <- timeSeries( (scale(X_train, mu_bear, FALSE) * wmat) %*% 
                                   torsion(covmat_bear, method = "mt"), time(X_train) )
      S_bull_wght <- timeSeries( (scale(X_train, mu_bull, FALSE) * wmat) %*% 
                                   torsion(covmat_bull, method = "mt"), time(X_train) )
    } else {
      S_all_wght <- S_all
      S_bear_wght <- S_bear
      S_bull_wght <- S_bull
    }
    
    # EWMA MD in MT space
    ewma_alpha <- spec$ewma_alpha
    S_sq_bear_ewma <- ema( S_bear_wght^2, alpha = ewma_alpha )
    S_sq_bull_ewma <- ema( S_bull_wght^2, alpha = ewma_alpha )
    S_sq_all_ewma <- ema( S_all_wght^2, alpha = ewma_alpha )
    md_bear <- apply( scale(S_sq_bear_ewma, FALSE, apply(S_bear, 2, var)), 1, sum )
    md_bear <- timeSeries( md_bear, time(S_sq_bear_ewma) )
    md_bull <- apply( scale(S_sq_bull_ewma, FALSE, apply(S_bull, 2, var)), 1, sum )
    md_bull <- timeSeries( md_bull, time(S_sq_bull_ewma) )
    md_all <- apply( scale(S_sq_all_ewma, FALSE, apply(S_all, 2, var)), 1, sum )
    md_all <- timeSeries( md_all, time(S_sq_all_ewma) )
    
    MD <- cbind( bear = md_bear, 
                 bull = md_bull,
                 all = md_all )
    MD_scl <- scale( MD, FALSE, TRUE )
    md_delta <- MD_scl[ ,"bear"] - MD_scl[ ,"bull"]
    
    ans <- cbind( states = bbstates,
                  bear = md_bear,
                  bull = md_bull,
                  all = md_all,
                  bear_scl = MD_scl[ ,"bear"],
                  bull_scl = MD_scl[ ,"bull"],
                  all_scl = MD_scl[ ,"all"],
                  signal = md_delta )
    colnames(ans) <-  c("states", "bear", "bull", "all", "bear_scl", "bull_scl", "all_scl", "delta")
    signal$insample <<- ans
    
    if ( isTRUE(spec$verbose) ) {
      cat( " done.\n")
    }
    
    return( TRUE )   
  }
  
  
  # --------------------------------------------------------------------------
  #' @export
  # --------------------------------------------------------------------------
  BBTEWMA.computeSignal <- function( dates = NULL )
  {
    
    stopifnot( is.character(dates) )
    
    # Initialize outputs
    md_all <- data$X_bm[dates, ] * NA
    md_bear <- md_bull <- md_delta <- md_all
    md_all_scl <- md_bear_scl <- md_bull_scl <- md_all
    states <- md_all
    
    # Append signals on exit
    on.exit( { 
      lSignal <- list( na.omit( cbind( md_bear = md_bear,
                                       md_bull = md_bull,
                                       md_bear_scl = md_bear_scl,
                                       md_bull_scl = md_bull_scl,
                                       signal = md_delta,
                                       states = states ) ) )
      names(lSignal) <- spec$method
      .self$appendSignal( value = lSignal )
    } )
    
    # Loop over dates
    for ( today in dates ) {
      
      if ( isTRUE(spec$verbose) ) {
        cat( "Running BBTEWMA.computeSignal: date =", today, "\n" )
      }
      
      # Bears n' bulls
      BBSObj <- spec$BBS$copy()
      BBSObj$data$X_level <- cumulated(data$X_bm[(rownames(data$X_bm) <= today), ], "discrete")
      if ( isTRUE(spec$robust) ) {
        BBSObj$runRobust( minphase_vec = spec$minphase_vec,
                          mincycle_vec = spec$mincycle_vec )
      } else {
        BBSObj$run()
      }
      data$BBS <<- BBSObj
      bbstates <- na.omit( BBSObj$output$states )
      # # Omit last phase
      # bbstates <- bbstates[1:(nrow(bbstates)-spec$minphase), ]
      
      X_train <- data$X[rownames(data$X) <= today, ]
      bbstates <- bbstates[rownames(X_train), ]
      
      
      # Compute regime specific stats
      
      # All
      covmat <- cov(X_train)
      
      # Bear
      p_bear <- (1 - bbstates) / 2
      stats_bear <- cov.wt(x = X_train, 
                           wt = as.numeric( p_bear / sum(p_bear) ),
                           cor = FALSE, 
                           center = TRUE, 
                           method = "ML")
      mu_bear <- stats_bear$center
      covmat_bear <- make.positive.definite( stats_bear$cov )
      
      # Bull
      p_bull <- (1 + bbstates) / 2
      stats_bull <- cov.wt(x = X_train, 
                           wt = as.numeric( p_bull / sum(p_bull) ),
                           cor = FALSE, 
                           center = TRUE, 
                           method = "ML")
      mu_bull <- stats_bull$center
      covmat_bull <- make.positive.definite( stats_bull$cov )
      
      # Minimum torsion series
      MT <-  mtca(X_train)
      S_all <- getter( MT, "sources" )
      S_bear <- timeSeries( scale(X_train, mu_bear, FALSE) %*% 
                              torsion(covmat_bear, method = "mt"), time(X_train) )
      S_bull <- timeSeries( scale(X_train, mu_bull, FALSE) %*% 
                              torsion(covmat_bull, method = "mt"), time(X_train) )
      
      # Weighted minimum torsion series
      if ( !is.null(data$wmat) ) {
        wmat <- data$wmat[rownames(X_train), ]
        S_all_wght <- timeSeries( (scale(X_train, TRUE, FALSE) * wmat) %*% 
                                    getter(MT, "torsion"), time(X_train) )
        S_bear_wght <- timeSeries( (scale(X_train, mu_bear, FALSE) * wmat) %*% 
                                     torsion(covmat_bear, method = "mt"), time(X_train) )
        S_bull_wght <- timeSeries( (scale(X_train, mu_bull, FALSE) * wmat) %*% 
                                     torsion(covmat_bull, method = "mt"), time(X_train) )
      } else {
        S_all_wght <- S_all
        S_bear_wght <- S_bear
        S_bull_wght <- S_bull
      }
      
      
      # EWMA MD in MT space
      ewma_alpha <- spec$ewma_alpha
      S_sq_bear_ewma <- ema( S_bear_wght^2, alpha = ewma_alpha )
      S_sq_bull_ewma <- ema( S_bull_wght^2, alpha = ewma_alpha )
      S_sq_all_ewma <- ema( S_all_wght^2, alpha = ewma_alpha )
      md_bear_tmp <- apply( scale(S_sq_bear_ewma, FALSE, apply(S_bear, 2, var)), 1, sum )
      md_bear_tmp <- timeSeries( md_bear_tmp, time(S_sq_bear_ewma) )
      md_bull_tmp <- apply( scale(S_sq_bull_ewma, FALSE, apply(S_bull, 2, var)), 1, sum )
      md_bull_tmp <- timeSeries( md_bull_tmp, time(S_sq_bull_ewma) )
      md_all_tmp <- apply( scale(S_sq_all_ewma, FALSE, apply(S_all, 2, var)), 1, sum )
      md_all_tmp <- timeSeries( md_all_tmp, time(S_sq_all_ewma) )
      
      MD <- cbind( bear = md_bear_tmp, 
                   bull = md_bull_tmp,
                   all = md_all_tmp )
      MD_scl <- scale( MD, FALSE, TRUE )
      md_delta_tmp <- MD_scl[ ,"bear"] - MD_scl[ ,"bull"]
      
      md_bear[today, ] <- md_bear_tmp[today, ]
      md_bull[today, ] <- md_bull_tmp[today, ]
      md_all[today, ] <- md_all_tmp[today, ]
      md_bear_scl[today, ] <- MD_scl[today, "bear"]
      md_bull_scl[today, ] <- MD_scl[today, "bull"]
      md_all_scl[today, ] <- MD_scl[today, "all"]
      md_delta[today, ] <- md_delta_tmp[today, ]
      states[today, ] <- bbstates[today, ]
      
    }
    
    return( TRUE ) 
  }
  
  # --------------------------------------------------------------------------
  #' @title Reference Class BBTEWMA
  #' @description setRefClass call for object of class BBTEWMA.
  #' @import methods
  #' @include class_sensor.R
  #' @export BBTEWMA
  # --------------------------------------------------------------------------
  BBTEWMA <- setRefClass( Class = "BBTEWMA", 
                           contains = "BBTurbulence",
                           methods = list( setCtrl = BBTEWMA.setCtrl,
                                           computeSignal = BBTEWMA.computeSignal,
                                           computeSignalInsample = BBTEWMA.computeSignalInsample) )
  
  
  