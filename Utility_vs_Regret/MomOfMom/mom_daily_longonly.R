  
  
  ############################################################################
  ### DAARC BACKTESTS - MOMENTUM
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     28.05.2020
  # First version:    28.05.2020
  # --------------------------------------------------------------------------
  
  
  require(RP)
  # require(AttM)
  require(DAARC)
  wd <- "R:/Asset_Management/Research_Projects/Multi_Asset/DAA/DAARC/"
  wd_warehouse <- paste0(wd, "Backtests/Momentum/waRehouse/")
  source( paste0(wd, "Source/custom_functions.R") )
  
  
  
  
  # Define rebalancing dates
  X <- DAA.updateData()
  # rebdates <- rownames(X)[ seq(from = 252, to = nrow(X), by = 100) ]
  rebdates <- rownames(X)[ seq(from = 252 * 2, to = nrow(X), by = 1) ]
  
  
  DAA <- DAARC$new()
  DAA$setCtrl()
  DAA$spec
  DAA$data <- X
  
  
  X_bm <- DAA$data[ ,"BM"]
  
  
  
  # --------------------------------------------------------------------------
  # BACKTEST
  # --------------------------------------------------------------------------
  
  
  source( paste0(wd, "Source/custom_functions.R") )
  source( "H:/R/AttM/R/class_Backtest.R" )
  source( "H:/R/AttM/R/class_BacktestBase.R" )
  
  
  # BT <- BacktestBase$new( rebdates = rebdates )
  BT <- BacktestBase$new( )
  BT$data <- list( X = X )
  BT$spec$rebdates <- rebdates
  BT$spec$mom_ctrl_width <- 252 * 2
  
  tau_vec <- 6 * (1:10) * 5
  BacktestBase$methods( mom = BacktestDAA.mom_lo )
 
  for ( i in seq(along = tau_vec) ) {
    
    BT$spec$mom_ctrl_tau <- tau_vec[i]
    # debugonce( BT$mom )
    # BT$mom( rebdate = rebdates[1] )
    try( lapply( rebdates, BT$mom ) )
    name <- paste0( "mom_lo_d_", tau_vec[i] )
    BT$save( wd = wd_warehouse,
             name = name,
             output_only = FALSE )
  }
  
  
  
  # --------------------------------------------------------------------------
  # SIMULATIONS
  # --------------------------------------------------------------------------
  
  
  # Load backtests and compute simulations

  filenames <- list.files(path = wd_warehouse,
                          pattern = ".rds")
  filenames <- unlist( lapply(filenames, FUN = function(x) { substring(x, 1, nchar(x)-4) }) )
  Names <- filenames[ grepl(pattern = "mom_lo_d", x = filenames) ]

  sim_bt <- NULL
  mumat <- NULL
  for ( i in seq(along = Names) ) {
    Obj <- AttM::loadBacktest(wd = wd_warehouse, name = Names[i])
    mu <- Obj$output$mu
    mumat <- cbind(mumat, mu)
    sig <- (1 + sign(mu-1)) / 2
    sim_tmp <- signalTesting.byTrading(X = X_bm,
                                      sig = sig,
                                      n_lag = 1,
                                      tc = 0)
    sim_bt <- cbind(sim_bt, sim_tmp)
  }
  colnames(sim_bt) <- colnames(mumat) <- Names
  
  sim <- na.omit( cbind(X_bm, sim_bt) )
  
  lStats <- descStats(sim)
  stats_fields <- c("cumret", "sds", "sharpe", "maxDD")
  t(lStats$stats[stats_fields, ])
  
  plot( as.simTS(sim) )
  
  colors <- fBasics::divPalette(n = ncol(mumat), "RdYlGn")
  plot(mumat, plot.type = "single", col = colors)
  abline(h = 1)
  
  
  
  
  
  
 
    
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  