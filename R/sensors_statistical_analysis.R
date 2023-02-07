  
  
  ############################################################################
  ### GOOD AND BAD TURBULENCE - SENSORS
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     31.01.2021
  # First version:    10.10.2020
  # --------------------------------------------------------------------------
  
  
  # --------------------------------------------------------------------------
  # REQUIRE
  # --------------------------------------------------------------------------
  
  require(RP)
  require(garcholz)
  require(DAARC)
  
  wd <- "R:/Asset_Management/Research_Projects/Multi_Asset/DAA/DAARC/Sensor_Testing/Univariate/Good_And_Bad_Turbulence/"
  source( paste0(wd, "Source/source_functions.R") )
  
  
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # DEFINE SIGNALS
  # --------------------------------------------------------------------------
  

  BBT <- loadSensor( sensor_name = "bbturbulence_base_dm" )
  BBTScl2 <- loadSensor( sensor_name = "bbturbulence_scl2_dm" )
  BBTEwma <- loadSensor( sensor_name = "bbtewma_base_dm" )
  BBTEwma05 <- loadSensor( sensor_name = "bbtewma05_base_dm" )
  BBTEwma01 <- loadSensor( sensor_name = "bbtewma01_base_dm" )
  BBTEwmaScl2 <- loadSensor( sensor_name = "bbtewma_scl2_dm" )
  BBTEwmaScl205 <- loadSensor( sensor_name = "bbtewma05_scl2_dm" )
  BBTEwmaScl201 <- loadSensor( sensor_name = "bbtewma01_scl2_dm" )
  BBTGarch <- loadSensor( sensor_name = "bbtgarch_base_dm" )
  BBTGarchScl2 <- loadSensor( sensor_name = "bbtgarch_scl2_dm" )
  
  NBCObj <- loadSensor( sensor_name = "nbc" )
  # NBCObj$update()
  # signals <- NBCObj$getSignal()
  # X_bm <- NBCObj$data$X_bm
  
  
  # BBT$update()
  # BBTScl2$update()
  # BBTEwma$update()
  # BBTEwma05$update()
  # BBTEwma01$update()
  # BBTEwmaScl2$update()
  # BBTEwmaScl205$update()
  # BBTEwmaScl201$update()
  # BBTGarch$update()
  # BBTGarchScl2$update()
  
  
  signals_nbc <- NBCObj$getSignal()
  
  signals_bbt <- cbind( base = BBT$getSignal(),
                        scl2 = BBTScl2$getSignal() )
  
  signals_bbtewma <- cbind( base = BBTEwma$getSignal(),
                            base05 = BBTEwma05$getSignal(),
                            base01 = BBTEwma01$getSignal(),
                            scl2 = BBTEwmaScl2$getSignal(),
                            scl205 = BBTEwmaScl205$getSignal(),
                            scl201 = BBTEwmaScl201$getSignal() )
  
  signals_bbtgarch <- cbind( base = BBTGarch$getSignal(),
                             scl2 = BBTGarchScl2$getSignal() )
  
  
  signals <- signals_bbtewma
  
  
  colors <- 1:ncol(signals)
  plot(tail(signals, 100), plot.type = "single", col = colors)
  legend("topleft", colnames(signals), lwd = 2, col = colors, text.col = colors, bty = "n" )
  abline(h = 0)
  
  
  
  
  # --------------------------------------------------------------------------
  # ACCURACY
  # --------------------------------------------------------------------------
  
  X_bm <- BBTEwma$data$X_bm
  sig <- (1 + sign(signals)) / 2
  
  FUN <- function(i)
  {
    TD <- trainingData( Y_train = (1 + sign(X_bm + 1e-16)) / 2, 
                        X_train = sig[ ,i],
                        # X_train = sig[ ,1] * 0 + runif(nrow(sig), 0, 1),
                        n_lag = 50 )
    brierScores( y.hat = TD$X_train, 
                 y.true = TD$Y_train,
                 method = "sq" )
  }
  setNames( unlist( lapply( 1:ncol(sig), FUN = FUN ) ), colnames(sig) )
  
  
  plot( x = TD$X_train, y = TD$Y_train )
  
  
  
  
  
  require(caret)
  
  #Creates vectors having data points
  expected_value <- factor(c(1,0,1,0,1,1,1,0,0,1))
  predicted_value <- factor(c(1,0,0,1,1,1,0,0,0,1))
  
  #Creating confusion matrix
  Cmat <- confusionMatrix( data = predicted_value, 
                           reference = expected_value )
  
  #Display results 
  Cmat
  
  
  
  TD <- trainingData( Y_train = (1 + sign(X_bm + 1e-16)) / 2, 
                      X_train = sig,
                      n_lag = 1 )
  
  Cmat <- confusionMatrix( data = factor(as.numeric(TD$X_train)),
                           reference = factor(as.numeric(TD$Y_train)) )
  
  Cmat
  
  
  
  
  
  