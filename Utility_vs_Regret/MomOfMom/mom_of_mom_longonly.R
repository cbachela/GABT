  
  
  ############################################################################
  ### DAARC BACKTESTS - MOMENTUM OF STRATEGY PORTFOLIOS
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     30.05.2020
  # First version:    28.05.2020
  # --------------------------------------------------------------------------
  
  
  require(RP)
  # require(AttM)
  require(DAARC)
  wd <- "R:/Asset_Management/Research_Projects/Multi_Asset/DAA/DAARC/"
  wd_warehouse <- paste0(wd, "Backtests/Momentum/waRehouse/")

  
  source( paste0(wd, "Source/custom_functions.R") )
  source( "H:/R/AttM/R/class_Backtest.R" )
  source( "H:/R/AttM/R/class_BacktestBase.R" )
  
  
  
  
  
  # --------------------------------------------------------------------------
  # LOAD BACKTESTS AND COMPUTE SIMULATIONS
  # --------------------------------------------------------------------------
  
  X <- DAA.updateData()
  X_bm <- X[ ,"BM"]
  
  
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
  
  
  
  
  
  # --------------------------------------------------------------------------
  # MOM-VARIANCE OPTIMIZATION OF MOM BACKTESTS
  # --------------------------------------------------------------------------
  
  Constraints <- constraints( selection = colnames(sim) )
  addConstraint(Constraints) <- budgetConstraint()
  addConstraint(Constraints) <- boxConstraint()
  
  BT <- BacktestBase$new( )
  BacktestBase$methods( momvariancePortfolio = BacktestDAA.momvariancePortfolio )
  BT$data <- list( X = sim,
                   X_sim = sim )
  rebdates <- rownames(sim)[seq(from = 252, to = nrow(sim), by = 1)]
  BT$spec$Constraints <- Constraints
  BT$spec$rebdates <- rebdates
  BT$spec$mom_ctrl_width <- 252
  BT$spec$mom_ctrl_tau <- 13 * 5

  # debugonce( BT$momvariancePortfolio )
  # debugonce( BT$initialWeights )
  # debugonce( BT$floatingWeights )
  
  try( lapply( rebdates, BT$momvariancePortfolio ) )
  name <- "mom_of_mom_lo"
  BT$save( wd = wd_warehouse,
           name = name,
           output_only = FALSE )
  
  
  
  source( paste0(wd, "Source/custom_functions.R") )
  source( "H:/R/AttM/R/class_Backtest.R" )
  source( "H:/R/AttM/R/class_BacktestBase.R" )
  
  
  # Simulate
  
  sim_mofm <- simPortfolio( X = BT$data$X_sim,
                            wghts = BT$output$weights,
                            fc = 0,
                            vc = 0 )
  colnames(sim_mofm) <- "mofm"
  sim <- na.omit( cbind(BT$data$X_sim, sim_mofm) )
  
  colnames(sim)
  colors <- c("black", fBasics::divPalette(n = ncol(sim)-1, "RdYlGn"))
  plot( log(cumulated(sim, "discrete")), col = colors, plot.type = "single" )
  lines( log(cumulated(sim[ ,"mofm"], "discrete")), lwd = 3 )
  
  
  lStats <- descStats(sim)
  t(lStats$stats)
  
  
  
  