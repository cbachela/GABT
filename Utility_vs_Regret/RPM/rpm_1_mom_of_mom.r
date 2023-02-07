  
  
  ############################################################################
  ### BacktetDAA - RPM1 - ANALYZE
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     29.06.2020
  # First version:    28.05.2020
  # --------------------------------------------------------------------------
  
  
  # This file is to analyze DAA backtests using the RPM methodology.
  
  require(SIDRC)
  
  wd <- "R:/Asset_Management/Research_Projects/Multi_Asset/SID/BacktestSIDRC/"
  
  source("R:/Asset_Management/Research_Projects/Multi_Asset/SID/BacktestSID/Source/BacktestSID_returns.R" )
  source("R:/Asset_Management/Research_Projects/Multi_Asset/SID/BacktestSID/Source/BacktestSID_portfolio.R" )
  source("R:/Asset_Management/Research_Projects/Multi_Asset/SID/BacktestSID/Source/BacktestSID_zzz.R" )
  source( paste0(wd, "Source/BacktestSIDRC_zzz.R") )
  
  
  
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
  
  
  colnames(sim) <- paste0("sim_", 1:ncol(sim))
  
  
  
  
  
  # Instantiate Backtest class
  # BT <- BacktestDAAWIP$new()
  # BT$setCtrl()
  # BT$updateData()
  
  
  Constraints <- constraints( selection = colnames(sim) )
  addConstraint(Constraints) <- budgetConstraint()
  addConstraint(Constraints) <- boxConstraint()
  
  
  BT2 <- BacktestSIDRC$new()
  # debugonce(aggWeekly)
  # BT$data$X <- aggWeekly( sim, compounding = "discrete", day = "Wed" )
  BT2$data$X_est <- sim
  BT2$data$X_sim <- sim
  rebdates <- rownames(BT2$data$X_est)[ seq(from = 252, to = nrow(BT2$data$X_est), by = 1) ]
  BT2$setCtrl( width = 252,
              alpha = 0.01,
              tau = 39 * 5,
              riskaversion = 1,
              Constraints = Constraints,
              NAMES_ADAPTIVE_EST = colnames(sim),
              # wd_warehouse = paste0(wd, "Momentum/MOFM/waRehouse/"),
              rebdates = rebdates )

  
  # BacktestDAAWIP$methods( mom = BacktestDAA.momvariancePortfolio )
  # BacktestDAAWIP$methods( runDAA = BacktestDAA.runDAA )
  # BT$runDAA( rebdates = BT$spec$rebdates, method = "mom" )
  
  # debugonce(BT$momPortfolio)
  # debugonce(BT$trainingData)
  BT2$runSID( method = "mom" )
  
  
  
  
  
  # Extract weights of single asset investment strategy
  
  lWmat_strat <- lWmat
  theta_mat <- BT2$output$weights
  # strat_array <- list2array( lWmat_strat )
  strat_mat <- timeSeries( do.call( cbind, lWmat_strat ), 
                           rownames(lWmat_strat[[1]]) )
  
  wmat <- theta_mat[ ,1] * NA
  for ( today in rownames(theta_mat) ) {
    
    coeff <- as.numeric( theta_mat[today, ] )
    wmat_tmp <- strat_mat[today, ]
    wmat[today, ] <- coeff %*% t(wmat_tmp)
  }
  
   
  plot( wmat )
  
  
  sim_mom <- signalTesting.byTrading( X = BT$data$X,
                                      sig = wmat,
                                      n_lag = 1,
                                      tc = 0.004 )
  sim_mom_nc <- signalTesting.byTrading( X = BT$data$X,
                                         sig = wmat,
                                         n_lag = 1,
                                         tc = 0 )
  
  
  sim_all <- na.omit( cbind( bm = BT$data$X, 
                             mom = sim_mom, 
                             mom_nc = sim_mom_nc, 
                             sim ) )
  
  
  
  # Statistics
  stats_fields <- c("cumret", "sds", "sharpe", "maxDD")
  lStats <- descStats( sim_all )
  stats <- t( lStats$stats[stats_fields, ] )
  lStats_nc <- descStats( sim_nc )
  stats_nc <- t( lStats_nc$stats[stats_fields, ] )
  
  barplot( stats[ ,"maxDD"], beside = TRUE )
  
  
  
  plot( as.simTS(sim_all[ ,1:5]) )
  
  stats
  stats_nc
  
  

  
  
  
  
  
  
