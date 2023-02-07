  
  
  ############################################################################
  ### SHINY - BEAR AND BULL STATES - SERVER SIDE
  ############################################################################
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version      03.07.2020
  # First version     03.07.2020
  # --------------------------------------------------------------------------
  


  
  # --------------------------------------------------------------------------
  # CUSTOM FUNCTIONS
  # --------------------------------------------------------------------------
  customData <- function( input )
  {
    
    # Series
    if ( input$series == "MSCI World" ) {
      X_level <- cumulated(global_data$X_bm, "discrete")
    } else if ( input$series == "VIX" ) {
      X_level <- global_data$vix
    }
    
    # States
    if ( input$algo == "Pagan" ) {
      BBS <- bbs( X = X_level,
                  minphase = input$minphase * 5,
                  mincycle = input$mincycle * 5,
                  e = 0 )
    } else {
      setpar_filtering_alg( tr_bear = input$th_down, 
                            tr_bull = input$th_up )
      bull <- run_filtering_alg( index = X_level * 100)
      BBS <- as.timeSeries( matrix(bull * 2 - 1, ncol = 1,
                                   dimnames = list(rownames(X_level), 
                                                   "LTF")) )
    }
    
    global_data$BBS <<- BBS
    global_data$X_level <<- X_level
    
    ans <- list( X_level = X_level,
                 BBS = BBS )
   
    return( ans )
  }
  
  
  # --------------------------------------------------------------------------
  customPlot <- function( input )
  {
    bbsPlot( input = input )
  }
  
  # --------------------------------------------------------------------------
  bbsPlot <- function( input )
  {
    plot(global_data$X_level)
    abline(v = time(global_data$BBS)[which(global_data$BBS == 1)], col = "lightgreen")
    abline(v = time(global_data$BBS)[which(global_data$BBS == -1)], col = "tomato")
    lines(global_data$X_level)
  }
  
  # --------------------------------------------------------------------------
  bbsfuzzyPlot <- function( input, data )
  {
    states_fuzzy <- global_data$BBSFuzzy
    states_fuzzy_unique <- sort(unique(states_fuzzy))
    # barplot( states_fuzzy_unique )
    colors <- fBasics::divPalette(n = length(states_fuzzy_unique), "RdYlGn")
    plot(data$X_level)
    for ( i in seq(along = states_fuzzy_unique) ) {
      idx <- which(states_fuzzy == states_fuzzy_unique[i])
      abline(v = time(states_fuzzy)[ idx ], col = colors[i], lwd = 4)
    }
    lines(data$X_level)
  }
  
 
  
  # --------------------------------------------------------------------------
  customTable <- function( input )
  {
    
    Y <- na.omit( global_data$sim[ ,unique(c("bm", input$strategies))] )
    BBS <- global_data$BBS
    dates <- intersect( rownames(Y), rownames(BBS) )
    Y <- Y[dates, ]
    BBS <- BBS[dates, ]
    dates_down <- rownames(BBS)[which(BBS == -1)]
    dates_up <- rownames(BBS)[which(BBS == 1)]
    
    stats_fields <- c("means", "sds", "sharpe", "maxDD", "te")
    spec <- descStatsSpec(bmName = "bm")
    lStats_down <- descStats(Y[intersect(rownames(Y), dates_down), ], spec)
    lStats_up <- descStats(Y[intersect(rownames(Y), dates_up), ], spec)
    stats_down <- lStats_down$stats
    stats_up <- lStats_up$stats
    
    # lStats_bear$desc; lStats_bull$desc
    # t(stats_bear[stats_fields, ])    
    # t(stats_bull[stats_fields, ]) 
    
    
    
    field <- "means"
    stats <- cbind( Downphase = stats_down[field, ],
                    Upphase = stats_up[field, ] )
    stats <- round(stats, 3)
    # stats
    # 
    # stats_rnk <- apply( stats[-c(3:4), ], 2, rank )
    # cbind( stats_rnk, apply(stats_rnk, 1, sum) )
    
    return( stats )
  }
  
  
  
  # --------------------------------------------------------------------------
  # SERVER
  # --------------------------------------------------------------------------
  shinyServer( function( input, output ) {
    
    # DATA
    reactiveData <- reactive({
      customData( input = input )
    })
    
    
    # CHARTS
    output$plot1 <- renderPlot({
      # input$goButton
      reactiveData()
      customPlot( input = input )
    })
    output$title.plot1 <- renderText({
      txt <- "Bear and Bull States"
    })
    
    # CHARTS
    output$plot2 <- renderPlot({
      # input$goButton
      reactiveData()
      # bbsfuzzyPlot( input = input )
      bbsPlot( input = input )
     
    })
    output$title.plot2 <- renderText({
      txt <- "Fuzzy Bear and Bull States"
    })
    
    # TABLES
    # Return a table
    # output$tmpTable <- renderTable({
    #   # input$goButton
    #   customTable( input = input )
    # }, digits = 4)
    output$tmpTable <- DT::renderDataTable({
      reactiveData()
      customTable(input = input)
    })

    
    # Title table
    output$caption <- renderText({
      txt <- paste0("xxx")
      return( txt )
    })
   
    
    
  })
  
  
  
  