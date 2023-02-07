  
  
  ############################################################################
  ### GOOD AND BAD TURBULENCE - WEIGHTED VS. UNWEIGHTED
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     22.01.2021
  # First version:    10.10.2020
  # --------------------------------------------------------------------------

  
  # --------------------------------------------------------------------------
  # REQUIRE
  # --------------------------------------------------------------------------
  
  require(RP)
  require(DAARC)
  
  wd <- "R:/Asset_Management/Research_Projects/Multi_Asset/DAA/DAARC/Sensor_Testing/Univariate/Good_And_Bad_Turbulence/"
  source( paste0(wd, "Source/source_functions.R") )
  
  
  
  # --------------------------------------------------------------------------
  # RAW TURBULENCE
  # --------------------------------------------------------------------------
  
  # Base
  Turb <- Turbulence$new()
  Turb$setCtrl( method = "base", universe = "dm" )
  # Turb$updateData()
  # lData <- Turb$data
  Turb$data <- lData
  Turb$computeSignalInsample()
  
  # Scl2
  TurbScl2 <- Turbulence$new()
  TurbScl2$setCtrl( method = "scl2", universe = "dm" )
  # TurbScl2$updateData()
  # lData_scl2 <- TurbScl2$data
  TurbScl2$data <- lData_scl2
  # debugonce( TurbScl2$computeSignalInsample )
  # debugonce( TurbScl2$prepareData )
  TurbScl2$computeSignalInsample()
  
  # Ewma
  TurbEwma <- Turbulence$new()
  TurbEwma$setCtrl( method = "mtewma", universe = "dm" )
  TurbEwma$data <- lData
  # debugonce( TurbEwma$computeSignalInsample )
  TurbEwma$computeSignalInsample()
  
  
  # Scl2 ewma
  TurbScl2Ewma <- Turbulence$new()
  TurbScl2Ewma$setCtrl( method = "mtewma", universe = "dm" )
  TurbScl2Ewma$data <- lData
  TurbScl2Ewma$data$wmat <- TurbScl2$data$wmat
  TurbScl2Ewma$computeSignalInsample()
  
  
  
  signals <- cbind( base = Turb$signal$insample[ ,"md"],
                    scl2 = TurbScl2$signal$insample[ ,"md"],
                    ewma = TurbEwma$signal$insample[ ,"md"],
                    scl2ewma = TurbScl2Ewma$signal$insample[ ,"md"])
  
  plot( signals )
  plot( signals, plot.type = "single" )
  plot( signals[ ,1:2], plot.type = "single" )
  
  
  boxplot( as.data.frame(signals) )
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # RELATIVE TURBULENCE
  # --------------------------------------------------------------------------
  
  # Base
  BBTBase <- BBTurbulence$new()
  BBTBase$setCtrl( method = "base", universe = "dm" )
  # BBTBase$updateData()
  # lD <- BBTBase$data
  BBTBase$data <- lD
  BBTBase$computeSignalInsample()  
  
  # Scl2
  BBTScl2 <- BBTurbulence$new()
  BBTScl2$setCtrl( method = "base", universe = "dm" )
  # BBTScl2$updateData()
  # lD_scl2 <- BBTScl2$data
  BBTScl2$data <- lD_scl2
  BBTScl2$computeSignalInsample()  
  
  # Ewma
  BBTEwma <- BBTEWMA$new()
  BBTEwma$setCtrl( method = "base", universe = "dm" )
  BBTEwma$spec$ewma_alpha <- 0.5
  BBTEwma$data <- lD
  BBTEwma$computeSignalInsample()
  
  # Scl2 Ewma
  BBTScl2Ewma <- BBTEWMA$new()
  BBTScl2Ewma$setCtrl( method = "base", universe = "dm" )
  BBTScl2Ewma$spec$ewma_alpha <- 0.5
  BBTScl2Ewma$data <- lD_scl2
  BBTScl2Ewma$computeSignalInsample()

  
  
  Name <- "delta"
  signals <- cbind( base = BBTBase$signal$insample[ ,Name],
                    scl2 = BBTScl2$signal$insample[ ,Name],
                    ewma = BBTEwma$signal$insample[ ,Name],
                    scl2ewma = BBTScl2Ewma$signal$insample[ ,Name])
  
  plot( signals )
  plot( signals, plot.type = "single" )
  plot( signals[ ,1:2], plot.type = "single" )
  
  
  boxplot( as.data.frame(signals) )
  
  
  plot( BBTScl2Ewma$signal$insample[ ,c("bear_scl", "bull_scl")], plot.type = "single" )
  plot( BBTEwma$signal$insample[ ,c("bear_scl", "bull_scl")], plot.type = "single" )
  
  
    
  