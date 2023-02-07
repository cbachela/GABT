  
  
  #############################################################
  ### SHINY FRONT END -  PORTFOLIO GEOMETRY - GLOBAL FILE   ###
  #############################################################
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version      20.04.2018
  # First version     20.04.2018
  # --------------------------------------------------------------------------


    
  # --------------------------------------------------------------------------
  # REQUIRE
  # --------------------------------------------------------------------------
  
  # require(GPO)
  require(RP)
  require(slolz)
  require(rgl)
  # require(shinyalert)

  source("utility_functions.R")
  

  
    

  # env_3d <- readRDS(file = "data/data_3d.Rds")
  env_3d <- readRDS(file = "R:/Asset_Management/R/Shiny/Portfolio_Geometry/data_3d.rds")
  env_4d <- readRDS(file = "data/data_4d.Rds")
  
  env <- env_3d
  
  X <- env_3d$X
  sim <- env_3d$sim
  cdf_mat <- env_3d$cdf_mat
  S_ls <- env_3d$S_ls
  S_lo <- env_3d$S_lo
  fp_ls <- env_3d$fp_ls
  fp_lo <- env_3d$fp_lo
  wmat_ls <- env_3d$wmat_ls
  wmat_lo <- env_3d$wmat_lo
  lCovmat <- env_3d$lCovmat
  lY <- env_3d$lY
  
 
  
  
  
  
  