  

  ############################################################################
  ### BBTurbulence - GARCH SIMULATION
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     16.10.2020
  # First version:    10.10.2020
  # --------------------------------------------------------------------------
  

  # Idea:
  # Simulate garch data and compute turbulence measures oos and check 
  # if it would generate a good trading strategy.
  # Can only work with asymmetric garch, otherwise there is no reason
  # to expect that returns in turbulent periods coincide with negative 
  # market phases.
  # Problem: no straight forward function to generate multivariate asymmetric garch
  # data.
  
  
  # --------------------------------------------------------------------------
  # REQUIRE
  # --------------------------------------------------------------------------
  
  require(rugarch)
  require(DAARC)
  
  wd <- "R:/Asset_Management/Research_Projects/Multi_Asset/DAA/DAARC/Sensor_Testing/Univariate/Good_And_Bad_Turbulence/"
  source( paste0(wd, "Source/source_functions.R") )
  
  
 
  
  # --------------------------------------------------------------------------
  # DATA
  # --------------------------------------------------------------------------
  
  BBT <- loadSensor( sensor_name = "bbturbulence_base_dm" )
  # BBT <- loadSensor( sensor_name = "bbturbulence_scl2_dm" )
  BBT$spec$b_scl <- FALSE
  BBT$computeSignalInsample()
  Y <- log(1 + BBT$data$X)
  X_bm <- BBT$data$X_bm
  
  Turb <- Turbulence$new()
  Turb$setCtrl( method = "base", universe = "dm" )
  # Turb$setCtrl( method = "scl2", universe = "dm" )
  Turb$spec$b_scl <- FALSE
  Turb$updateData()
  Turb$computeSignalInsample()
  
  
  # --------------------------------------------------------------------------
  # PARAMETERS
  # --------------------------------------------------------------------------
  
  # BBQ parameters
  minphase_m <- 2
  mincycle_m <- 5
  minphase_w <- minphase_m * 4
  mincycle_w <- mincycle_m * 4
  minphase_d <- minphase_m * 4 * 5
  mincycle_d <- mincycle_m * 4 * 5
  theta <- 0.15
  
  # GARCH parameters
  M1 <- apply(Y, 2, mean)
  M2 <- cov(Y)
  n_sample = 10^3 * 10
  asTS = TRUE
  start_date = NULL
  frqncy = "day"
  seed = 1234
  variance_model = list(model = "gjrGARCH",
                        garchOrder = c(1, 1),
                        submodel = NULL, 
                        external.regressors = NULL, 
                        variance.targeting = FALSE)
  mean_model = list(armaOrder = c(0, 0), 
                    include.mean = FALSE, 
                    archm = FALSE,
                    archpow = 1, 
                    arfima = FALSE, 
                    external.regressors = NULL, 
                    archex = FALSE)
  distribution_model = "sstd"
  # start_pars = list(alpha1 = 0.05,
  #                   beta1 = 0.99,
  #                   gamma1 = 0.02)
  start_pars = list()
  # fixed_pars = list(alpha1 = 0.05,
  #                   beta1 = 0.93,
  #                   gamma1 = 0.02)
  fixed_pars = list()
  
  
  
  spec <- ugarchspec( variance.model = variance_model,
                      mean.model = mean_model,
                      distribution.model = "sstd",
                      start.pars = start_pars,
                      fixed.pars = fixed_pars )
  
  
  
  # --------------------------------------------------------------------------
  # GARCH SIMULATION
  # --------------------------------------------------------------------------
  
  # S <- mtca( Y )
  S <- pca( Y )
  lfit <- apply( log(1 + S), 2, ugarchfit, spec = spec )
  sim <- rmvgarch2( l_fit_or_spec = lfit,
                    M1 = M1, 
                    M2 = M2,
                    n_sample = n_sample,
                    asTS = TRUE,
                    start_date = NULL,
                    frqncy = frqncy,
                    seed = 999 )
  attr(sim, "seed")
  attr(sim, "seed_vec")
  
  
  
  descStats(sim, descStatsSpec(what = c("basics", "portmanteauSQ")) )
  plot.timeSeries( log(cumulated(sim, "discrete")) )
  
  
  # For simplicity, take eqw portfolio as benchmark so that we 
  # can work with unweighted Mahalanobis distances.
  X_synt <- timeSeries( sim[1:nrow(Y), ], time(Y) )
  X_bm_synt <- apply(X_synt, 1, mean)
  
  
  
  
  # --------------------------------------------------------------------------
  # COMPUTE DIFFERENT TURBULENCE MEASURES ON SYNTHETIC DATA - INSAMPLE
  # --------------------------------------------------------------------------
  
  # Raw turbulence
  MDSynt <- Turbulence$new()
  MDSynt$setCtrl( method = "base", universe = "dm" )
  MDSynt$spec$iso <- colnames(X_synt)
  MDSynt$spec$b_scl <- FALSE
  MDSynt$data <- list( X = X_synt,
                       S = X_synt,
                       X_bm = X_bm_synt )
  MDSynt$computeSignalInsample()
  
  # Relative turbulence
  BBTSynt <- BBT$copy()
  BBTSynt$spec$universe <- "synthetic_garch"
  BBTSynt$data <- list( X = X_synt,
                        S = X_synt,
                        X_bm = X_bm_synt )
  BBTSynt$computeSignalInsample()
  
  # Relative time series turbulence
  BBTSyntEwma <- BBTEWMA$new()
  BBTSyntEwma$setCtrl( method = "base", universe = "dm" )
  BBTSyntEwma$spec$iso <- colnames(X_synt)
  BBTSyntEwma$spec$b_scl <- FALSE
  BBTSyntEwma$data <- list( X = X_synt,
                            X_bm = X_bm_synt )
  BBTSyntEwma$computeSignalInsample()
  
  
  
  
  MDSynt$signal$insample
  
  
  
  
  
  # --------------------------------------------------------------------------
  # DESCRIPTIVE STATISTICS
  # --------------------------------------------------------------------------
  
  
  DS <- dirichletSampling( Y_train = X_bm_synt,
                           # X_train = BBTSyntEwma$signal$insample[ ,"all"],
                           X_train = MDSynt$signal$insample[ ,"sig_md"],
                           weights_fun = "l1",
                           sclfct = NULL,
                           scl_by_entropy = TRUE )
  ldens <- lapply( DS, density )
  plot.ldensity( ldens, fillin = FALSE )
  
  
  
  DS <- dirichletSampling( Y_train = X_bm_synt,
                           # X_train = BBTSynt$signal$insample[ ,"all"],
                           X_train = BBTSynt$signal$insample[ ,"delta"],
                           weights_fun = "l1",
                           sclfct = NULL,
                           scl_by_entropy = TRUE )
  ldens <- lapply( DS, density )
  plot.ldensity( ldens, fillin = FALSE )
  
  
  
  DS <- dirichletSampling( Y_train = X_bm_synt,
                           # X_train = BBTSyntEwma$signal$insample[ ,"all"],
                           X_train = BBTSyntEwma$signal$insample[ ,"delta"],
                           weights_fun = "l1",
                           sclfct = NULL,
                           scl_by_entropy = TRUE )
  ldens <- lapply( DS, density )
  plot.ldensity( ldens, fillin = FALSE )
  
  
  # On real data
  DS <- dirichletSampling( Y_train = BBT$data$X_bm,
                           X_train = BBT$signal$insample[ ,"delta"],
                           weights_fun = "l1",
                           sclfct = NULL,
                           scl_by_entropy = TRUE )
  ldens <- lapply( DS, density )
  plot.ldensity( ldens, fillin = FALSE )
  
  
  
  
  # --------------------------------------------------------------------------
  # BACKTESTS
  # --------------------------------------------------------------------------
  
  
  sig1 <- round( MDSynt$signal$insample$sig_md_ema0.089 )
  sig2 <- (1 + sign(BBTSynt$signal$insample[ ,"delta"])) / 2
  sig3 <- (1 + sign(BBTSyntEwma$signal$insample[ ,"delta"])) / 2
  signals <- cbind( mdsynt = sig1,
                    bbtsynt = sig2,
                    bbtsyntewma = sig3 )
  
  y <- X_bm_synt
  colnames(y) <- "bm"
  test <- signalTesting.byTrading( X = y,
                                   sig = signals,
                                   n_lag = 1,
                                   tc = 0.004 )
  colnames(test) <- colnames(signals)
  test_nc <- signalTesting.byTrading( X = y,
                                      sig = signals,
                                      n_lag = 1,
                                      tc = 0 )
  colnames(test_nc) <- colnames(signals)
  X_tmp <- na.omit(cbind(y, test))
  X_tmp_nc <- na.omit(cbind(y, test_nc))
  
  colors <- c(1, fBasics::rainbowPalette(ncol(test)))
  stats_fields <- c("cumret", "means", "sds", "sharpe", "maxDD")
  
  plot( as.simTS(X_tmp), col = colors )
  t( descStats( X_tmp )$stats[stats_fields, ] )
  drawDownStats(X_tmp)
  
  plot( as.simTS(X_tmp_nc), col = colors )
  t( descStats( X_tmp_nc )$stats[stats_fields, ] )
  drawDownStats(X_tmp_nc)
  
  
  
  
  
  # --------------------------------------------------------------------------
  # COMPARE SYNTHETIC GARCH TURBULENCE TO REAL TURBULENCE
  # --------------------------------------------------------------------------
  
  Obj <- BBT
  ObjSynt <- Obj$copy()
  ObjSynt$spec$universe <- "garch"
  ObjSynt$data <- list( X = X_synt,
                        S = X_synt,
                        X_bm = X_bm_synt )
  # debugonce(ObjSynt$computeSignalInsample )
  ObjSynt$computeSignalInsample()
  
  
  md <- BBT$signal$insample[ ,"all"]
  md_synt <- timeSeries( ObjSynt$signal$insample[1:nrow(md), "all"], time(md) )
  MD <- cbind(md, md_synt)
  ldens <- apply( MD, 2, density )
  
  plot.ldensity( ldens, fillin = FALSE )
  plot( MD, plot.type = "single" )
  
  
  DS <- dirichletSampling( Y_train = Obj$data$X_bm,
                           X_train = Obj$signal$insample[ ,"all"],
                           weights_fun = "l1",
                           sclfct = NULL,
                           scl_by_entropy = TRUE )
  ldens <- lapply( DS, density )
  plot.ldensity( ldens, fillin = FALSE )
  
  
  
  DS <- dirichletSampling( Y_train = ObjSynt$data$X_bm,
                           X_train = ObjSynt$signal$insample[ ,"all"],
                           weights_fun = "l1",
                           sclfct = NULL,
                           scl_by_entropy = TRUE )
  ldens <- lapply( DS, density )
  plot.ldensity( ldens, fillin = FALSE )
  
  
  
  plot( x = as.numeric(Obj$data$X_bm[rownames(Obj$signal$insample), ]),
        y = as.numeric(Obj$signal$insample[ ,"all"]) )
  plot( x = as.numeric(ObjSynt$data$X_bm),
        y = as.numeric(ObjSynt$signal$insample[ ,"all"]) )
  
  
  
  
  sig <- (1 + sign(ObjSynt$signal$insample[ ,"delta"])) / 2
  X_bm <- ObjSynt$data$X_bm
  test <- signalTesting.byTrading( X = X_bm,
                                   sig = sig,
                                   n_lag = 2,
                                   tc = 0.004 )
  colnames(test) <- colnames(sig)
  test_nc <- signalTesting.byTrading( X = X_bm,
                                      sig = sig,
                                      n_lag = 2,
                                      tc = 0 )
  colnames(test_nc) <- colnames(sig)
  X_tmp <- na.omit(cbind(X_bm, test))
  X_tmp_nc <- na.omit(cbind(X_bm, test_nc))
  
  colors <- c(1, fBasics::rainbowPalette(ncol(test)))
  stats_fields <- c("cumret", "means", "sds", "sharpe", "maxDD")
  
  plot( as.simTS(X_tmp), col = colors )
  t( descStats( X_tmp )$stats[stats_fields, ] )
  drawDownStats(X_tmp)
  
  plot( as.simTS(X_tmp_nc), col = colors )
  t( descStats( X_tmp_nc )$stats[stats_fields, ] )
  drawDownStats(X_tmp_nc)
  
  
  
  
  
  
  
  ############
  
  MDSYNT <- Turbulence$new()
  MDSYNT$setCtrl( method = "base", universe = "dm" )
  MDSYNT$spec$iso <- colnames(X_synt)
  MDSYNT$spec$b_scl <- FALSE
  MDSYNT$data <- list( X = X_synt,
                       S = X_synt,
                       X_bm = X_bm_synt )
  MDSYNT$computeSignalInsample()
  
  
  plot( Turb$signal$insample[ ,c("md", "sig_md", "sig_md_ema0.04")] )
  plot( MDSYNT$signal$insample[ ,c("md", "sig_md", "sig_md_ema0.04")] )
  
  
  lMD <- list( md = Turb$signal$insample[ ,"md"],
               md_synt = MDSYNT$signal$insample[ ,"md"] )
  ldens <- lapply( lMD,
                   FUN = density )
  plot.ldensity( ldens, fillin = FALSE ) 
  
  
  
  # --------------------------------------------------------------------------
  # BEAR BULL EWMA MAHALANOBIS ANALYSIS ON SIMULATED GARCH DATA
  # --------------------------------------------------------------------------
  
  BBTSYNT <- BBTEWMA$new()
  BBTSYNT$setCtrl( method = "base", universe = "dm" )
  BBTSYNT$spec$iso <- colnames(X_synt)
  BBTSYNT$spec$b_scl <- FALSE
  BBTSYNT$data <- list( X = X_synt,
                        X_bm = X_bm_synt )
  BBTSYNT$computeSignalInsample()
  
  plot( BBTSYNT$signal$insample )
  
  
  
  sig <- (1 + sign(BBTSYNT$signal$insample[ ,"delta"])) / 2
  X_bm <- BBTSYNT$data$X_bm
  test <- signalTesting.byTrading( X = X_bm,
                                   sig = sig,
                                   n_lag = 2,
                                   tc = 0.004 )
  colnames(test) <- colnames(sig)
  test_nc <- signalTesting.byTrading( X = X_bm,
                                      sig = sig,
                                      n_lag = 2,
                                      tc = 0 )
  colnames(test_nc) <- colnames(sig)
  X_tmp <- na.omit(cbind(X_bm, test))
  X_tmp_nc <- na.omit(cbind(X_bm, test_nc))
  
  colors <- c(1, fBasics::rainbowPalette(ncol(test)))
  stats_fields <- c("cumret", "means", "sds", "sharpe", "maxDD")
  
  plot( as.simTS(X_tmp), col = colors )
  t( descStats( X_tmp )$stats[stats_fields, ] )
  drawDownStats(X_tmp)
  
  plot( as.simTS(X_tmp_nc), col = colors )
  t( descStats( X_tmp_nc )$stats[stats_fields, ] )
  drawDownStats(X_tmp_nc)
  
  
  
  
  
  
  

  
  
  
  