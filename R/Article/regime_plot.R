  
  ############################################################################
  ### GOOD AND BAD TURBULENCE - REGIME PLOT
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     03.06.2022
  # First version:    03.06.2022
  # --------------------------------------------------------------------------
  
  
  
  
  require(DAARC)

  
  
  BBT <- loadSensor( sensor_name = "bbturbulence_base_dm" )
  # BBT$updateData()
 
  BBSObj <- BBT$data$BBS$copy()
  BBSObj$data$X_level <- cumulated(BBT$data$X_bm, "discrete")
  BBSObj$runRobust()
  
  BBSObj$plotStates()
  BBSObj$phaseStats()
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  