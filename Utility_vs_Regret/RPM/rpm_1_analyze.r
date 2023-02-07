  
  
  ############################################################################
  ### BacktetDAA - RPM1 - ANALYZE
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     29.06.2020
  # First version:    28.05.2020
  # --------------------------------------------------------------------------
  
  
  # This file is to analyze DAA backtests using the RPM methodology.

  
  
  require(RP)
  require(slolz)
  require(DAARC)
  wd <- "R:/Asset_Management/Research_Projects/Multi_Asset/DAA/DAARC/BacktestDAA/"
  wd_warehouse <- paste0(wd, "RPM/waRehouse/")
  source( paste0(wd, "Source/custom_functions.R") )
  
  
  BacktestDAAWIP <- setRefClass( Class = "BacktestDAAWIP", 
                                 contains = "BacktestDAA",
                                 methods = list() )
  
  
  
  # Simulate
  file_names <- list.files( path = wd_warehouse, 
                            pattern = ".rds" )
  portfolio_names <- substr(file_names, 1, nchar(file_names) - 4)
  
  
  lWmat <- list()
  lLambda <- list()
  lES <- lEU <- lER <- list()
  for ( i in seq(along = portfolio_names) ) {
    
    BT <- loadBacktest( wd = wd_warehouse,
                        name = portfolio_names[i],
                        b_new = TRUE )
    lWmat[[i]] <- BT$output$weights
    lLambda[[i]] <- BT$output$lambda
    lES[[i]] <- BT$output$exp_state
    lEU[[i]] <- BT$output$exp_util
    lER[[i]] <- BT$output$exp_regret
  }
  names(lWmat) <- names(lLambda) <- 
    names(lES) <- names(lEU) <- names(lER) <- portfolio_names
  
  lWmat <- c( list(buynhold = lWmat[[1]][ ,1] * 0 + 1), lWmat )
  
  
  # Simulations
  FUN <- function(x, X_bm, tc)
  {
    tmp <- signalTesting.byTrading( X = X_bm,
                                    sig = x,
                                    n_lag = 1,
                                    tc = tc )
    return( tmp )
  }
  lSim <- lapply( lWmat, FUN = FUN, X_bm = BT$data$X, tc = 0.0004 )
  sim <- do.call( cbind, lSim )
  lSim_nc <- lapply( lWmat, FUN = FUN, X_bm = BT$data$X, tc = 0 )
  sim_nc <- do.call( cbind, lSim_nc )
  
  # Statistics
  stats_fields <- c("cumret", "sds", "sharpe", "maxDD")
  lStats <- descStats( na.omit(sim) )
  stats <- t( lStats$stats[stats_fields, ] )
  lStats_nc <- descStats( na.omit(sim_nc) )
  stats_nc <- t( lStats_nc$stats[stats_fields, ] )
  
  barplot( stats[ ,"maxDD"], beside = TRUE )
  
  
  plot.stats( lStats, sortby = NULL )  
  
  
  cbind(stats, stats_nc)
  
  
  
  plot( as.simTS(sim_nc[ ,1:5]) )
  
  
  
  
  

  i = which( names(lWmat) == "rpm_mvu1_bbturbulence.md_delta" )
  BT <- loadBacktest( wd = wd_warehouse,
                      name = names(lWmat)[i],
                      b_new = TRUE )
  plot(BT$data$signal)
  plot( lWmat[[i]], main = names(lWmat)[i] )
  
  
  
  plot( lEU[[2]], plot.type = "single" )
  EU <- do.call( cbind, lapply( lEU, FUN = function(x) { x[ ,ncol(x)] } ) )
  plot( EU, plot.type = "single" )
  abline( h = 0 )
  
  wmat <- do.call( cbind, lapply( lWmat, FUN = function(x) { x[ ,"w_star"] } ) )
  wmat_nr <- do.call( cbind, lapply( lWmat, FUN = function(x) { x[ ,"w_star_nr"] } ) )
  colnames(wmat) <- colnames(wmat_nr) <- names(lWmat)
  
  plot( wmat, plot.type = "single" )
  
  plot( apply(wmat, 1, mean) )
  lines( apply(wmat_nr, 1, mean), col = 2 )
  
  
  
  
  
  sig <- wmat
  y <- BT$data$X[ ,"BM"]
  test <- signalTesting.byTrading(X = y,
                                  sig = sig,
                                  n_lag = 1,
                                  tc = 0)
  colnames(test) <- colnames(sig)
  X_tmp <- na.omit(cbind(bm = y, test))
  plot( log(cumulated(X_tmp, "discrete")), plot.type = "single" )
  
  stats_fields <- c("cumret", "sds", "maxDD")
  lStats <- descStats(X_tmp) 
  t(lStats$stats[stats_fields, ])
  
  
  
  
  
  
  
  
  
    
    
  