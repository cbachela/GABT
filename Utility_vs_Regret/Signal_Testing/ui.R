  
  
  ############################################################################
  ### SHINY  -  SIGNAL TESTING  -  USER INTERFACE FILE  
  ############################################################################
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version      20.06.2019
  # First version     20.06.2019
  # --------------------------------------------------------------------------
  
  
  # Conditional Density
  # Trading
  
  
  
  # --------------------------------------------------------------------------
  # Define User Interface (UI)
  # --------------------------------------------------------------------------
  shinyUI(pageWithSidebar(
    
    # Application title
    headerPanel("Signal Processing"),
    
    # Sidebar with controls
    sidebarPanel(
      
      wellPanel(
        
        # What 
        selectInput(inputId = 'what'
                    , label = 'Show'
                    , choices = list('Signal',
                                     'In/Out Backtest',
                                     'Bear/Bull Classification',
                                     'Conditional Density',
                                     'FCMDM',
                                     'Logistic Regression',
                                     'Performance After Signal',
                                     'Importance Sampling')
                    , selected = "Signal")
        
        # Dependent variable
        , selectInput(inputId = 'dependent_variable'
                      , label = 'Dependent Variable'
                      , choices = c('BM Returns', 
                                    'Bear / Bull States',
                                    'Fuzzy Bear / Bull States')
                      , selected = 'BM Returns')
        
        
        # Sensor
        , selectInput(inputId = 'sensor'
                      , label = 'Sensor'
                      , choices = names(DAA$sensors)
                      , selected = names(DAA$sensors)[1])
        
        # Signal
        , htmlOutput(outputId = "signal")
        
        
        , conditionalPanel(
          condition = "input.what == 'Signal'"
          , htmlOutput(outputId = "x_eval")
        )
        
        , conditionalPanel(
          condition = "input.what == 'FCMDM'"
          , sliderInput(inputId = "n_centers"
                        , label = "Number of cluster"
                        , min = 2
                        , max = 10
                        , value = 2
                        , step = 1)
          , radioButtons(inputId = "fcmdm_which_alpha"
                         , label = "Dirichlet parameter"
                         , choices = c("memb", "pmat", "pmat_times_nrow")
                         , selected = "pmat_times_nrow")
          , radioButtons(inputId = "fcmdm_which_plot"
                         , label = "Choose plot"
                         , choices = c("Signal and cluster center",
                                       "Membership function")
                         , selected = "Signal and cluster center")
        )
        
        , conditionalPanel(
          condition = "input.what == 'Logistic Regression'"
          , sliderInput(inputId = "x_eval_logit"
                        , label = "Signal level for evaluation"
                        , min = 0
                        , max = 1
                        , value = 0
                        , step = 0.01)
        )
        
        , conditionalPanel(
          condition = "input.what == 'In/Out Backtest'"
          , checkboxInput(inputId = "tc"
                        , label = "Transaction Costs"
                        , value = TRUE)
        )
        
        , conditionalPanel(
          condition = "input.what == 'Bear/Bull Classification' || 
                       input.what == 'Logistic Regression'"
          , sliderInput(inputId = "minphase"
                        , label = "Minimum phase length"
                        , min = 5
                        , max = 5 * 4 * 5
                        , value = 5
                        , step = 5)
          , sliderInput(inputId = "mincycle"
                        , label = "Minimum cycle length"
                        , min = 5
                        , max = 5 * 4 * 10
                        , value = 5 * 4 * 5
                        , step = 5)
        )
       
        , conditionalPanel(
          condition = "input.what == 'Performance After Signal'"
          , sliderInput(inputId = "n_lag_pas"
                        , label = "Lag"
                        , min = 0
                        , max = 60
                        , value = 1
                        , step = 1)
        )
        # # Lag
        # , conditionalPanel(
        #   condition = "input.what != 'Performance After Signal'"
        #   , sliderInput(inputId = "n_lag"
        #                 , label = "Lag"
        #                 , min = 0
        #                 , max = 10
        #                 , value = 1
        #                 , step = 1)
        # )
        
        # Aggregate returns
        , conditionalPanel(
          condition = "input.what == 'Importance Sampling'"
          , sliderInput(inputId = "n_agg"
                        , label = "Aggregate Returns over x days"
                        , min = 0
                        , max = 42
                        , value = 21
                        , step = 1)
        )
        

      )

    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Plot"
                 , h3( textOutput(outputId = "title.plot_1") )
                 , plotOutput( outputId = "plot_1" )
                 # , tableOutput(outputId="tmpTable")
                 , h3( textOutput( outputId = "title.plot_2") )
                 , plotOutput( outputId = "plot_2" )
        ),
        tabPanel("Table"
                 , h3( textOutput( outputId = "caption") )
                 , verbatimTextOutput( outputId = "tmpPrint" )
        )
      ))  
    
  ))
  
  
  
  
  
  
  
  
  
  
  
  
  