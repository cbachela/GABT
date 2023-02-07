    
  
  
  ############################################################################
  ### GOOD AND BAD TURBULENCE - IN-SAMPLE - SCALE VS. BASE
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     22.10.2020
  # First version:    10.10.2020
  # --------------------------------------------------------------------------
    
  
  
  
  # --------------------------------------------------------------------------
  # BBTurbulence
  # --------------------------------------------------------------------------
  
  # Findings:
  # scl and scl2 behave similar whereas base performs badly.
  
  
  require(DAARC)
  
  
  Obj1 <- Turbulence$new()
  Obj1$setCtrl(method = "scl2")
  # Obj1$updateData()
  # lData <- Obj1$data
  Obj1$data <- lData
  Obj1$computeSignalInsample()
  
  Obj2 <- Obj1$copy()
  Obj2$spec$method <- "base"
  Obj2$data$wmat <- NULL
  Obj2$computeSignalInsample()
  
  Obj3 <- Obj1$copy()
  Obj3$spec$method <- "scl"
  Obj3$computeSignalInsample()
  
  
  
  sig <- cbind( scl2 = (1 + Obj1$signal$insample[ ,"states"]) / 2,
                base = (1 + Obj2$signal$insample[ ,"states"]) / 2 )
  test <- signalTesting.byTrading( X = Obj1$data$X_bm,
                                   sig = sig,
                                   n_lag = 1,
                                   tc = 0.004 )
  colnames(test) <- colnames(sig)
  test_nc <- signalTesting.byTrading( X = Obj1$data$X_bm,
                                      sig = sig,
                                      n_lag = 1,
                                      tc = 0 )
  colnames(test_nc) <- colnames(sig)
  X_tmp <- na.omit(cbind(BBT$data$X_bm, test, test_nc))
  plot( as.simTS(X_tmp) )
  
  
  
  
  signals <- cbind( scl2 = Obj1$signal$insample[ ,"md_sig"],
                    base = Obj2$signal$insample[ ,"delta"],
                    scl = Obj3$signal$insample[ ,"delta"] )
  
  test <- signalTesting.byTrading( X = Obj1$data$X_bm,
                                   sig = sig,
                                   n_lag = 1,
                                   tc = 0.004 )
  colnames(test) <- colnames(sig)
  test_nc <- signalTesting.byTrading( X = Obj1$data$X_bm,
                                      sig = sig,
                                      n_lag = 1,
                                      tc = 0 )
  colnames(test_nc) <- colnames(sig)
  X_tmp <- na.omit(cbind(BBT$data$X_bm, test, test_nc))
  plot( as.simTS(X_tmp) )
  descStats( X_tmp )
  
  
  
  
  TD <- trainingData( Y_train = Obj1$data$X_bm,
                      X_train = signals )
  lDS <- list()
  for ( j in 1:ncol(TD$X_train) ) {
    lDS[[j]] <- dirichletSampling( Y_train = TD$Y_train,
                                   X_train = TD$X_train[ ,j],
                                   weights_fun = "l1",
                                   sclfct = NULL,
                                   scl_by_entropy = TRUE )
  }
  
  lldens <- lapply( lDS, FUN = function(x) { lapply( x, density ) } )
  names(lldens) <- colnames(signals)
  
  plot.ldensity( lldens[[1]], main = names(lldens)[1])
  plot.ldensity( lldens[[2]], main = names(lldens)[2])
  plot.ldensity( lldens[[3]], main = names(lldens)[3])
  
  
  
  
  # --------------------------------------------------------------------------
  # BBTurbulence
  # --------------------------------------------------------------------------
  
  # Findings:
  # scl and scl2 behave similar whereas base performs badly.
  
  
  require(DAARC)
  
  
  Obj1 <- BBTurbulence$new()
  Obj1$setCtrl(method = "scl2")
  # Obj1$updateData()
  # lData <- Obj1$data
  Obj1$data <- lData
  Obj1$computeSignalInsample()
  
  Obj2 <- Obj1$copy()
  Obj2$spec$method <- "base"
  Obj2$data$wmat <- NULL
  Obj2$computeSignalInsample()
  
  Obj3 <- Obj1$copy()
  Obj3$spec$method <- "scl"
  Obj3$computeSignalInsample()
  
  
  
  sig <- cbind( scl2 = (1 + Obj1$signal$insample[ ,"states"]) / 2,
                base = (1 + Obj2$signal$insample[ ,"states"]) / 2 )
  test <- signalTesting.byTrading( X = Obj1$data$X_bm,
                                   sig = sig,
                                   n_lag = 1,
                                   tc = 0.004 )
  colnames(test) <- colnames(sig)
  test_nc <- signalTesting.byTrading( X = Obj1$data$X_bm,
                                      sig = sig,
                                      n_lag = 1,
                                      tc = 0 )
  colnames(test_nc) <- colnames(sig)
  X_tmp <- na.omit(cbind(BBT$data$X_bm, test, test_nc))
  plot( as.simTS(X_tmp) )
  
  
  
  
  signals <- cbind( scl2 = Obj1$signal$insample[ ,"delta"],
                    base = Obj2$signal$insample[ ,"delta"],
                    scl = Obj3$signal$insample[ ,"delta"] )
  sig <- (1 + sign(sig)) / 2
  
  test <- signalTesting.byTrading( X = Obj1$data$X_bm,
                                   sig = sig,
                                   n_lag = 1,
                                   tc = 0.004 )
  colnames(test) <- colnames(sig)
  test_nc <- signalTesting.byTrading( X = Obj1$data$X_bm,
                                      sig = sig,
                                      n_lag = 1,
                                      tc = 0 )
  colnames(test_nc) <- colnames(sig)
  X_tmp <- na.omit(cbind(BBT$data$X_bm, test, test_nc))
  plot( as.simTS(X_tmp) )
  descStats( X_tmp )

  
  
  
  TD <- trainingData( Y_train = Obj1$data$X_bm,
                      X_train = signals )
  lDS <- list()
  for ( j in 1:ncol(TD$X_train) ) {
    lDS[[j]] <- dirichletSampling( Y_train = TD$Y_train,
                                   X_train = TD$X_train[ ,j],
                                   weights_fun = "l1",
                                   sclfct = NULL,
                                   scl_by_entropy = TRUE )
  }
  
  lldens <- lapply( lDS, FUN = function(x) { lapply( x, density ) } )
  names(lldens) <- colnames(signals)
  
  plot.ldensity( lldens[[1]], main = names(lldens)[1])
  plot.ldensity( lldens[[2]], main = names(lldens)[2])
  plot.ldensity( lldens[[3]], main = names(lldens)[3])
  
  
  
  
  
  