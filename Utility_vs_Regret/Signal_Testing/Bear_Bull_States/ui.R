  
  
  ############################################################################
  ### SHINY - BEAR AND BULL STATES - USER INTERFACE
  ############################################################################
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version      03.07.2020
  # First version     03.07.2020
  # --------------------------------------------------------------------------
  
  
  
  
  # --------------------------------------------------------------------------
  # Define User Interface (UI)
  # --------------------------------------------------------------------------
  shinyUI(pageWithSidebar(
    
    # Application title
    headerPanel("Portfolio Characteristics"),
    
    # Sidebar with controls
    sidebarPanel(
      wellPanel(
        
        
        # Series
        selectInput( inputId = 'series'
                     , label = 'Series'
                     , choices = list('MSCI World', 
                                      "VIX")
                     , selected = 'MSCI World' )
        
        # Strategies
        # , selectInput( inputId = 'strategies'
        #                , label = 'Strategies'
        #                , choices = colnames(global_data$sim)
        #                , multiple = TRUE
        #                , selected = c("bm", "OLZ_fund") )
        , pickerInput( inputId = 'strategies'
                       , label = 'Strategies'
                       , multiple = TRUE
                       # , options = pickerOptions(maxOptions = ncol(global_data$sim))
                       , options = list( 'actions-box' = TRUE )
                       , choices = colnames(global_data$sim)
                       , selected = c("bm", "OLZ_fund") )
        
        # Pattern recognition algorithm
        , selectInput( inputId = 'algo'
                       , label = 'Pattern recognition algorithm'
                       , choices = list('Pagan',
                                        'Lunde-Timmermann')
                       , selected = 'Pagan' )
        
        , conditionalPanel(
          condition = "input.algo == 'Pagan'"
          
            # Parameters of BBS Algo
            , sliderInput( inputId = "minphase"
                           , label = "Minimum phase length (in weeks)"
                           , min = 1
                           , max = 1 * 4 * 5
                           , value = 1
                           , step = 1 )
          
            , sliderInput( inputId = "mincycle"
                            , label = "Minimum cycle length (in weeks)"
                            , min = 1
                            , max = 1 * 4 * 10
                            , value = 1 * 4 * 5
                            , step = 1 )
            
            # , sliderInput(inputId = "drop_last_phase"
            #               , label = "Drop last phase"
            #               , min = 1
            #               , max = 1 * 4 * 5
            #               , value = 1
            #               , step = 1)
        )
        
        , conditionalPanel(
          condition = "input.algo == 'Lunde-Timmermann'"
            , sliderInput( inputId = "th_down"
                           , label = "Percentage decrease"
                           , min = -30
                           , max = -10
                           , value = -20
                           , step = 1 ) 
          , sliderInput( inputId = "th_up"
                         , label = "Percentage increase"
                         , min = 10
                         , max = 30
                         , value = 20
                         , step = 1 ) 
        )
          
      )
      , wellPanel(
        
        # # estimate
        # actionButton(inputId = 'goButton'
        #              , label = 'Estimate')
        
        # # download weights
        # , downloadButton('downloadData', 'Download Weights')
        
      )),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Plot"
                 , h3( textOutput(outputId = "title.plot1") )
                 , plotOutput( outputId = "plot1" )
                 , h3( textOutput( outputId = "title.plot2") )
                 , plotOutput( outputId = "plot2" )
        ),
        tabPanel("Table"
                 # , h3( textOutput( outputId = "caption") )
                 # , verbatimTextOutput( outputId = "tmpPrint" )
                 , h2("Desriptive Statistcs")
                 , DT::dataTableOutput( outputId = "tmpTable" )
        )
      ))
    
  ))
  
  
  
  
  
  
  
  
  
  
  
