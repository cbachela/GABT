  
  
  ############################################################################
  ### GOOD AND BAD TURBULENCE - FAMA FRENCH INDUSTRY PORTFOLIOS
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     09.03.2021
  # First version:    09.03.2021
  # --------------------------------------------------------------------------
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # REQUIRE
  # --------------------------------------------------------------------------
  
  require(RP)
  require(BSS)
  require(simolz)
  require(DAARC)
  
  wd <- "R:/Asset_Management/Research_Projects/Multi_Asset/DAA/DAARC/Sensor_Testing/Univariate/Good_And_Bad_Turbulence/"
  source( paste0(wd, "Source/source_functions.R") )
  

  
  
  # --------------------------------------------------------------------------
  # LOAD DATA
  # --------------------------------------------------------------------------
  
  # debugonce(get.file)
  X <- get.file( wd = paste0(wd, "Data/"), filename = "FF_Industry_10", sep = ";" )
  X <- X[isWeekday(time(X)), ]
  
  head(X)
  plot(X)
  
  
  # Benchmark 
  X_bm <- apply( X, 1, mean )
  colnames(X_bm) <- "bm"
  
  
  
  descStats(cbind(X, X_bm))
  
  
 
  
  
  
  # Object <- loadSensor( sensor_name = "bbtewma_base_dm" )
  # Obj <- DAARC::BBTEWMA$new()
  # Obj$setCtrl( universe = "dm",
  #              method = "base" )
  # Obj$spec$universe <- "ff10"
  # Obj$spec$width <- 1000
  # Obj$spec$name <- "bbtewma_base_ff10"
  # Obj$data <- list( X = X,
  #                   BBS = Object$data$BBS,
  #                   X_bm = X_bm,
  #                   capw = NULL )
  # Obj$computeSignalInsample()
  
  Obj <- loadSensor( sensor_name = "bbtewma_base_ff10" )
  # debugonce( Obj$computeSignal )
  Obj$updateSignal()
  Obj$save()
  
  
  sig_tmp <- Obj$signal$insample[ ,"delta"]
  sig_tmp <- Obj$getSignal()
  sig_tmp <- window(sig_tmp, "1935-01-01", end(sig_tmp))
  
  test <- signalTesting.byTrading( X = Obj$data$X_bm,
                                   sig = (1 + sign(sig_tmp)) / 2,
                                   n_lag = 1, 
                                   tc = 0 )
  X_tmp <- na.omit(cbind(Obj$data$X_bm, test))
  plot( as.simTS(X_tmp) )
  
  stats_fields <- c("cumret", "sds", "sharpe", "maxDD")
  descStats( X_tmp )$stats[stats_fields, ]
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # BBQ ALGO SETTINGS
  # --------------------------------------------------------------------------
  
  BBSF <- BBSFuzzy$new()
  BBSF$setCtrl()
  BBSF$spec$universe <- "ff10"
  BBSF$spec$round <- 0
  BBSF$data <- list( X_bm = X_bm,
                     X = X_bm )
  
  BBSF$computeSignalInsample()
  head(BBSF$signal$insample)
  
  BBSF$getSignal
  BBSF$spec$round
  
  BBSF$plotStates( type = "is" )
  
  mu_bear <- meanGeo( X = BBSF$data$X_bm[which(BBSF$getSignal(method = "insample") == -1), ] )
  mu_bull <- meanGeo( X = BBSF$data$X_bm[which(BBSF$getSignal(method = "insample") == 1), ] )

  mu_bear; mu_bull
  
  
  
  
  
  # --------------------------------------------------------------------------
  # TURBULENCE
  # --------------------------------------------------------------------------
  
  MDBase <- Turbulence$new()
  MDBase$setCtrl()
  MDBase$spec$b_scl <- FALSE
  MDBase$spec$name <- "turbulence_base_ff10"
  MDBase$spec$universe <- "ff10"
  MDBase$data <- list(X = X, 
                      S = X,
                      X_bm = X_bm)
  
  
  # MDBase <- loadSensor( sensor_name = "turbulence_base_ff10" )
  # MDBase$update()

  MDBase$computeSignalInsample()
  # MDBase$updateSignal()
  # MDBase$save()
  
  
  signal <- MDBase$signal$insample[ ,"md"]
  signal <- MDBase$getSignal()
  plot(signal)
  
  
  plot( cbind(log(cumulated(X_bm, "discrete")), 
              MDBase$signal$insample[ ,1:5]) )
  
  
  plot( MDBase$signal$insample[ ,"md"])
  plot( MDBase$signal$base[ ,"md"])
  
  
  
  tailleft( MDBase$signal$base )
  
  
  
  
  signal <- ema( MDBase$signal$insample[ ,"sig_md"], 1)
  signal <- ema( MDBase$getSignal(), 0.1)
  sig <- round(signal)
  test <- signalTesting.byTrading( X = X_bm,
                                   sig = sig,
                                   n_lag = 1,
                                   tc = 0 )
  X_tmp <- na.omit( cbind( bm = X_bm, bt = test) ) 
  plot( log(cumulated(X_tmp, "discrete")), plot.type = "single" )
  
  plot( log(cumulated(tail(X_tmp, 10000), "discrete")), plot.type = "single" )
  
  descStats(X_tmp)
  drawDownStats(X_tmp)
  
  
  
  
  
  # --------------------------------------------------------------------------
  # TURBULENCE ON MINIMUM TORSION SV
  # --------------------------------------------------------------------------
  
  MDMTSV <- Turbulence$new()
  MDMTSV$setCtrl( method = "mtsv" )
  MDMTSV$spec$universe <- "ff10"
  MDMTSV$spec$b_scl <- FALSE
  MDMTSV$spec$name <- "turbulence_mtsv_ff10"
  MDMTSV$data <- list(X = X, 
                      S = X,
                      X_bm = X_bm)
  
  
  # MDMTSV <- loadSensor( sensor_name = "turbulence_mtsv_ff10" )
  
  
  MDMTSV$computeSignalInsample()
  MDMTSV$updateSignal()
  # MDMTSV$save()
  
  
  signal <- MDMTSV$signal$insample[ ,"md"]
  # signal <- MDMTSV$getSignal()
  plot(signal)
  
  
  plot( cbind(log(cumulated(X_bm, "discrete")), 
              MDMTSV$signal$insample[ ,1:5]) )
  
  
  plot( MDMTSV$signal$insample[ ,"md"])
  plot( MDMTSV$signal$base[ ,"md"])
  
  
  
  
  
  signal <- ema( MDMTSV$signal$insample[ ,"sig_md"], 1)
  # signal <- ema( MDMTSV$getSignal(), 0.1)
  sig <- round(signal)
  test <- signalTesting.byTrading( X = X_bm,
                                   sig = sig,
                                   n_lag = 1,
                                   tc = 0 )
  X_tmp <- na.omit( cbind( bm = X_bm, bt = test) ) 
  plot( log(cumulated(X_tmp, "discrete")), plot.type = "single" )
  
  plot( log(cumulated(tail(X_tmp, 10000), "discrete")), plot.type = "single" )
  
  descStats(X_tmp)
  drawDownStats(X_tmp)
  
  
  
  
  
  # --------------------------------------------------------------------------
  # BEAR / BULL TURBULENCE
  # --------------------------------------------------------------------------
  
  # Base
  BBTBase <- BBTurbulence$new()
  BBTBase$setCtrl( method = "base" )
  BBTBase$spec$universe <- "ff10"
  BBTBase$spec$name <- "bbturbulence_base_ff10"
  BBTBase$spec$iso <- colnames(X)
  BBTBase$spec$BBS <- BBSF
  BBTBase$data <- list( X_bm = X_bm,
                        X = X,
                        BBS = BBSF )
  # BBTBase <- loadSensor( sensor_name = "bbturbulence_base_ff10" )
  # BBTBase$computeSignalInsample()
  # BBTBase$update()
  
  # debugonce(BBTBase$computeSignalInsample)
  BBTBase$computeSignalInsample()
  BBTBase$updateSignal()
  
  
  edit(BBTBase$computeSignalBase)
  
  
  
  plot( ema(BBTBase$signal$insample[ ,"delta"], 0.1) )
  abline(h = 0)
  
  signal <- ema( BBTBase$signal$insample[ ,"delta"], 0.1)
  sig <- (1 + sign(signal)) / 2
  test <- signalTesting.byTrading( X = X_bm,
                                   sig = sig,
                                   n_lag = 1,
                                   tc = 0.004 )
  X_tmp <- na.omit( cbind( bm = X_bm, bt = test) ) 
  plot( log(cumulated(X_tmp, "discrete")), plot.type = "single" )
  
  plot( log(cumulated(tail(X_tmp, 10000), "discrete")), plot.type = "single" )
  
  descStats(X_tmp)
  drawDownStats(X_tmp, top = 10)
  
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # MINIMUM TORSION EWMA BBT
  # --------------------------------------------------------------------------
  
  # Ewma
  BBTEwma <- BBTEWMA$new()
  BBTEwma$setCtrl( method = "base", universe = "dm" )
  BBTEwma$spec$universe <- "ff10"
  BBTEwma$spec$name <- "bbtewma_base_ff10"
  BBTEwma$spec$ewma_alpha <- 0.1
  BBTEwma$spec$iso <- colnames(X)
  BBTEwma$spec$name <- paste0(BBTEwma$spec$name, "_a", BBTEwma$spec$ewma_alpha)
  BBTEwma$spec$BBS <- BBSF
  BBTEwma$spec$width <- 252 * 5
  BBTEwma$data <- list( X_bm = X_bm,
                        X = X )
  # BBTEwma <- loadSensor( sensor_name = "bbtewma_base_ff10_a0.1" )
  # BBTEwma$computeSignalInsample()
  # BBTEwma$update()
  
  BBTEwma$computeSignalInsample()
  # debugonce(BBTEwma$computeSignal)
  # debugonce(BBTEwma$updateSignal)
  # BBTEwma$updateSignal()
  
  
  signal <- BBTEwma$getSignal()
  plot(signal)
  
  
  signal <- ema( BBTEwma$signal$insample[ ,"delta"], 1)
  signal <- ema( BBTEwma$getSignal(), 0.1 )
  sig <- (1 + sign(signal)) / 2
  test <- signalTesting.byTrading( X = X_bm,
                                   sig = sig,
                                   n_lag = 1,
                                   tc = 0 )
  X_tmp <- na.omit( cbind( bm = X_bm, bt = test) ) 
  plot( log(cumulated(X_tmp, "discrete")), plot.type = "single" )
  
  plot( log(cumulated(tail(X_tmp, 10000), "discrete")), plot.type = "single" )
  
  descStats(X_tmp)
  drawDownStats(X_tmp)
  
  
  plot(signal)
  abline(h = 0)    
  
  plot( BBTEwma$data$BBS$getSignal(method = "insample") )
  
  
  
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # MINIMUM TORSION GARCH BBT
  # --------------------------------------------------------------------------

  BBTgarch <- DAARC::BBTGARCH$new()
  BBTgarch$setCtrl( method = "base", universe = "dm" )
  BBTgarch$spec$universe <- "ff10"
  BBTgarch$spec$name <- "BBTgarch_base_ff10"
  BBTgarch$spec$iso <- colnames(X)
  BBTgarch$spec$BBS <- BBSF
  BBTgarch$data <- list( X_bm = X_bm,
                        X = X )
  # BBTgarch <- loadSensor( sensor_name = "bbtgarch_base_ff10_a0.1" )
  # BBTgarch$computeSignalInsample()
  # BBTgarch$update()
  
  BBTgarch$computeSignalInsample()
  
  
  
  signal <- ema( BBTgarch$signal$insample[ ,"delta"], 1)
  sig <- (1 + sign(signal)) / 2
  test <- signalTesting.byTrading( X = X_bm,
                                   sig = sig,
                                   n_lag = 1,
                                   tc = 0 )
  X_tmp <- na.omit( cbind( bm = X_bm, bt = test) ) 
  plot( log(cumulated(X_tmp, "discrete")), plot.type = "single" )
  
  plot( log(cumulated(tail(X_tmp, 10000), "discrete")), plot.type = "single" )
  
  descStats(X_tmp)
  drawDownStats(X_tmp)
  
  
  plot(signal)
  abline(h = 0)    
  
  
  
  
  
  # --------------------------------------------------------------------------
  # MINIMUM TORSION SV BBT
  # --------------------------------------------------------------------------
  
  require(stochvol)
  
  BBTsv <- DAARC::BBTSV$new()
  BBTsv$setCtrl( method = "base", universe = "dm" )
  BBTsv$spec$universe <- "ff10"
  BBTsv$spec$name <- "bbtsv_base_ff10"
  BBTsv$spec$iso <- colnames(X)
  BBTsv$spec$BBS <- BBSF
  BBTsv$data <- list( X_bm = X_bm,
                         X = X )
  # BBTsv <- loadSensor( sensor_name = "bbtsv_base_ff10" )
  # BBTsv$computeSignalInsample()
  # BBTsv$update()
  
  BBTsv$computeSignalInsample()
  
  
  
  signal <- ema( BBTsv$signal$insample[ ,"delta"], 1)
  sig <- (1 + sign(signal)) / 2
  test <- signalTesting.byTrading( X = X_bm,
                                   sig = sig,
                                   n_lag = 1,
                                   tc = 0 )
  X_tmp <- na.omit( cbind( bm = X_bm, bt = test) ) 
  plot( log(cumulated(X_tmp, "discrete")), plot.type = "single" )
  
  plot( log(cumulated(tail(X_tmp, 10000), "discrete")), plot.type = "single" )
  
  descStats(X_tmp)
  drawDownStats(X_tmp)
  
  
  plot(signal)
  abline(h = 0)    
  
  
  
  
  # Bootstrap
  
  require(RP)
  require(slolz)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  