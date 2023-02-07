

#############################################################
### SHINY FRONT END - PORTFOLIO GEOMETRY - USER INTERFACE ###
#############################################################

# --------------------------------------------------------------------------
# Cyril Bachelard
# This version      20.04.2018
# First version     20.04.2018
# --------------------------------------------------------------------------


  
  
  
  # --------------------------------------------------------------------------
  # Define User Interface (UI)
  # --------------------------------------------------------------------------
  shinyUI(pageWithSidebar(
    
    # Application title
    headerPanel("Portfolio Geometry"),
    
    # Sidebar with controls
    sidebarPanel(
      
      wellPanel(
        
        # Show
        selectInput(inputId = 'what'
                    , label = 'Show'
                    , choices = list('Return Series',
                                     'Utility',
                                     'Constraints',
                                     'Covariance',
                                     'Efficient Frontier')
                    # , selected = 'Return Series')
                    , selected = "Utility")
        
        , conditionalPanel(
            condition = "input.what == 'Utility'"
            
            
          # Assets
          , selectInput(inputId = 'assets'
                        , label = 'Asset'
                        , choices = as.list(colnames(X))
                        , selected = colnames(X)[1]
                        , multiple = FALSE)
          
          # Utility function
          , selectInput(inputId = 'utility_function'
                        , label = 'Utility Function'
                        , choices = list('Mean-Variance'
                                         , 'Minimum-Variance'
                                         , 'Kahneman-Tversky'
                                         , 'Omega')
                                         # , 'Entropy')
                        , selected = 'Mean-Variance')
          
          # Parameters of utility funcitons
          , conditionalPanel(
            condition = "input.utility_function == 'Mean-Variance'"
            
            , sliderInput(inputId = "lambda"
                          , label = "Risk Aversion"
                          , min = 0
                          , max = 10
                          , value = 5
                          , step = 0.2)
          )
          
          , conditionalPanel(
            condition = "input.utility_function == 'Kahneman-Tversky'"
            
            , sliderInput(inputId = "RP"
                          , label = "Reference Point"
                          , min = -0.1
                          , max = 0.1
                          , value = 0
                          , step = 0.01)
            
            , sliderInput(inputId = "a"
                          , label = "Parameter a"
                          , min = 0
                          , max = 1
                          , value = 0.5
                          , step = 0.01)
            
            , sliderInput(inputId = "b"
                          , label = "Parameter b"
                          , min = 1
                          , max = 10
                          , value = 2.25
                          , step = 0.25)
          )
          
          , conditionalPanel(
            condition = "input.utility_function == 'Omega'"
            
            , sliderInput(inputId = "threshold"
                          , label = "Threshold"
                          , min = -0.1
                          , max = 0.1
                          , value = 0
                          , step = 0.001)
          )
        )
        
        , conditionalPanel(
            condition = "input.what == 'Covariance'"
          
          , selectInput(inputId = "shrinkage_target"
                        , label = "Shrinkage Target"
                        , choices = c("dcv", "ccor")
                        , selected = "dcv"
                        , multiple = FALSE)
          
          , sliderInput(inputId = "shrinkage_intensity"
                          , label = "Shrinkage Intensity"
                          , min = 0
                          , max = 1
                          , value = 0
                          , step = 0.1)
        )
        
        , conditionalPanel(
            condition = "input.what == 'Constraints'"
          
            , checkboxInput(inputId = "hide_S_ls"
                            #, label = "Show Long Only"
                            , label = "Constrained Set Only"
                            , value = FALSE)
            
            , radioButtons(inputId = "dimension"
                           , label = "Dimensionality"
                           , choices = c(3, 4)
                           , selected = 3)
            
            , checkboxInput(inputId = "upper"
                            , label = "Upper Bounds"
                            , value = FALSE)
            
            , conditionalPanel(
              condition = "input.upper == '1'"

              , sliderInput(inputId = "upper_x1"
                            , label = "Upper Bound on Asset 1"
                            , min = 0.34
                            , max = 1
                            , value = 1
                            , step = 0.1)
              , sliderInput(inputId = "upper_x2"
                            , label = "Upper Bound on Asset 2"
                            , min = 0.34
                            , max = 1
                            , value = 1
                            , step = 0.1)
              , sliderInput(inputId = "upper_x3"
                            , label = "Upper Bound on Asset 3"
                            , min = 0.34
                            , max = 1
                            , value = 1
                            , step = 0.1)
              , sliderInput(inputId = "upper_x23"
                            , label = "Upper Bound on Assets 2 and 3"
                            , min = 0.5
                            , max = 1
                            , value = 1
                            , step = 0.1)

            )
         
            , radioButtons(inputId = "cons"
                           , label = "Constraints"
                           , choices = c("Budget"
                                         ,"Long Only"
                                         , "Turnover"
                                         , "Tracking Error"
                                         , "VaR"
                                         , "Expected Shortfall")
                           , selected = "Budget"
                           , inline = FALSE)
            
            , conditionalPanel(
              condition = "input.cons == 'Turnover'"
              
              , checkboxInput(inputId = "add_diamond"
                              , label = "Add Diamond"
                              , value = FALSE)
            )
            , conditionalPanel(
              condition = "input.cons == 'Tracking Error'"
              
              , checkboxInput(inputId = "add_ellipsoid"
                              , label = "Add Ellipsoid"
                              , value = FALSE)
            )
        )
        
        
        # # Start date
        # , selectInput(inputId = 'startdate'
        #               , label = 'Start Date'
        #               , choices = rownames(X)
        #               , selected = rownames(X)[1])
        # 
        # # End date
        # , selectInput(inputId = 'enddate'
        #               , label = 'End Date'
        #               , choices = rownames(X)
        #               , selected = rownames(X)[nrow(X)])
        
      )
      
      # , wellPanel(
      #   
      #   # estimate
      #   #useShinyalert(),
      #   actionButton(inputId = 'goButton'
      #                , label = 'Estimate')
      #   
      #   # # download weights
      #   # , downloadButton('downloadData', 'Download Weights')
      #   
      #   )
      # 
      ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Plot"
                 , h3( textOutput(outputId = "title.plot_2d") )
                 , plotOutput( outputId = "plot_2d" )
                 # , tableOutput(outputId="tmpTable")
                 , h3( textOutput( outputId = "title.plot_3d") )
                 , plotOutput( outputId = "plot_3d" )
        ),
        tabPanel("Table"
                 , h3( textOutput( outputId = "caption") )
                 , verbatimTextOutput( outputId = "tmpPrint" )
        )
      ))  
    
  ))
  
  
  
  
  
  
  
  
  
  
  
