  
  
  ############################################################################
  ### SHINY - BEAR AND BULL STATES - PREPARE DATA
  ############################################################################
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version      03.07.2020
  # First version     03.07.2020
  # --------------------------------------------------------------------------
  
  
  require(BBSolz)
  require(DAARC)
  
  wd <- "R:/Asset_Management/R/Shiny/Signal_Testing/Bear_Bull_States/"
  
  # MSCI World Index
  X_bm <- rodbcGetOLZDBReturns( assetName = "NDDUWI" )
  
  # VIX
  vix_ret <- rodbcGetOLZDBReturns(assetName = "VIX",
                                  refCcy = "USD",
                                  frqncy = "daily",
                                  startDate = "1990-03-01",
                                  compounding = "continuous",
                                  na.rm = "r")
  vix_base <- 21.8999999999999  # This is the level of the vix as of 01.03.1990
  vix <- cumulated(vix_ret, compounding = "continuous") * vix_base
 
  
  
  # BBS
  BBS <- loadSensor(sensor_name = "bbs")
  BBSFuzzy <- loadSensor(sensor_name = "bbsfuzzy")

  
  
  
  
  # Lunde Timmermann filter
  setpar_filtering_alg( tr_bear = -30, 
                        tr_bull = 30 )
  bull <- run_filtering_alg( index = X_level * 100)
  LTF <- as.timeSeries( matrix(bull * 2 - 1, ncol = 1,
                               dimnames = list(rownames(X_level), 
                                               "LTF")) )
  
  
  
  
  
  ###
  # PAT BVG Backtests
  universe <- "dm_ex_ch"
  wd_cons <- "R:/Asset_Management/Research_Projects/Equity/Constraints/MSCI_Min_Vol_Constraints/"
  wd_trerr <- "R:/Asset_Management/Research_Projects/Equity/Optimization/Multi_Objective_Dynamic/Variance_TrackingError/"
  env_cons <- readRDS( file = paste0(wd_cons, "waRehouse/data4report.rds") )
  env_trerr <- readRDS( file = paste0(wd_trerr, "waRehouse/data4report.rds") )
  X_olz <- rodbcGetOLZDBReturns( assetName = "OLZEWEI SW",
                                 refCcy = "CHF",
                                 frqncy = "daily" )
  sim <- cbind( env_cons$lBM[[ universe ]],
                OLZ_fund = X_olz,
                env_cons$lSim_noESG[[ universe ]],
                env_trerr$lSim_noESG[[ universe ]] )
  
  
  
  
  
  # Save environment
  global_data <- new.env()
  global_data$sim <- sim
  global_data$X_bm <- X_bm
  global_data$vix <- vix
  # global_data$BBS <- BBS$signal$insample
  global_data$BBS <- LTF
  global_data$BBSFuzzy <- BBSFuzzy$signal$insample
  global_data$LTF <- LTF
  
  
  
  saveRDS( global_data, file = paste0(wd, "Data/global_data.rds") )
  
  
  
 
  
  
  
  
  
  
  