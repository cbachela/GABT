  
  
  ############################################################################
  ### SHINY - BEAR AND BULL STATES - GLOBAL
  ############################################################################
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version      03.07.2020
  # First version     03.07.2020
  # --------------------------------------------------------------------------
  

  
  # --------------------------------------------------------------------------
  # REQUIRE
  # --------------------------------------------------------------------------
  
  require(shiny)
  require(shinyWidgets)
  require(DT)
  require(BBSolz)
  require(DAARC)
 
  
  wd <- "R:/Asset_Management/R/Shiny/Signal_Testing/Bear_Bull_States/"
  global_data <- readRDS( file = paste0(wd, "Data/global_data.rds") )
  
  
  
  
  
  