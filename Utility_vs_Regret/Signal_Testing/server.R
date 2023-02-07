  
  
  ############################################################################
  ### SHINY  -  SIGNAL TESTING  -  SERVER FILE  
  ############################################################################
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version      20.06.2019
  # First version     20.06.2019
  # --------------------------------------------------------------------------
  
  
  
  # --------------------------------------------------------------------------
  # SERVER
  # --------------------------------------------------------------------------
  shinyServer(function(input, output) {
    
    # DATA
    reactiveData <- reactive({
      customData( input = input )
    })
    
    # # build a reactive expression that only invalidates 
    # # when the value of input$goButton becomes out of date 
    # # (i.e., when the button is pressed)
    # rebalance <- eventReactive(input$goButton, {
    #   customRebalance(input = input)
    # })
    
    
    # Update input
    output$signal <- renderUI({ 
      sensor_name <- ifelse( is.null(input$sensor), "momentum", input$sensor )
      choices <- "signal"
      if ( sensor_name == "bcpolz" ) {
        choices <- c( choices, colnames(DAA$sensors[[sensor_name]]$signal$base$momat) )
      } else {
        choices <- c( choices, colnames(DAA$sensors[[sensor_name]]$signal$base) )
      }
      choices <- unique(choices)
      selected <- choices[1]
      selectInput(inputId = "signal"
                  , label = "Signal"
                  , choices = choices
                  , selected = selected)
    })
    
    output$x_eval <- renderUI({ 
      # data_object <- reactiveData()
      sensor_name <- tolower(input$sensor[1])
      Sensor <- DAA$sensors[[sensor_name]]
      signal_name <- input$signal[1]
      signal <- as.numeric(Sensor$signal$base[ ,signal_name])
      data_object = list(signal = signal)
      sliderInput(inputId = "x_eval"
                  , label = "Signal level for evaluation"
                  , min = round(min(data_object$signal, na.rm = TRUE), 2)
                  , max = round(max(data_object$signal, na.rm = TRUE), 2)
                  , value = round(min(data_object$signal, na.rm = TRUE), 2)
                  , step = round((max(data_object$signal, na.rm = TRUE) - 
                                    min(data_object$signal, na.rm = TRUE)) / 10, 2))
    })
    
    
    # First chart
    output$plot_1 <- renderPlot({
      data_object <- reactiveData()
      customPlot( input = input,
                  data = data_object )
    })
    output$title.plot_1 <- renderText({
      txt <- switch(input$what,
                    "vrp_raw" = "VRP RAW",
                    "Weights" = "Weights")
    })
    
    
    # Second chart
    output$plot_2 <- renderPlot({
      data_object <- reactiveData()
      customPlot2( input = input,
                   data = data_object )
    })
    output$title.plot_2 <- renderText({
      txt <- switch(input$what,
                    "vrp_raw" = "VRP RAW",
                    "Weights" = "Weights")
    })
    
    # TABLES
    # # Return a table
    # output$tmpTable <- renderTable({
    #   # input$goButton
    #   customTable(input = input)
    # }, digits = 4)
    # observeEvent(input$goButton, {
    #   # Show a modal when the button is pressed
    #   shinyalert("Apero!"
    #              , "Thank you for your attention."
    #              , type = "success"
    #              , animation = FALSE
    #              
    #              )
    # })
    
    
    # Return a table
    output$tmpPrint <- renderPrint({
      data_object <- reactiveData()
      tbl <- customTable(input = input, 
                         data = data_object)
      print(tbl)
    })
    
    # Title table
    output$caption <- renderText({
      txt <- paste0("Statistics")
      return( txt )
    })
    
    # # DOWNLOAD
    # output$downloadData <- downloadHandler(
    #   filename = function() {
    #     paste(Sys.Date(), '_Weights_', input$universe, '.csv', sep="")
    #   },
    #   content = function(file) {
    #     write.csv(customDownload(input = input
    #                              , object = reactiveData())
    #                              , file)
    #   }
    # )
    
    
  })
