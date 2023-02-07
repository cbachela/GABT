  
  
  ############################################################################
  ### SHINY  -  SIGNAL TESTING  -  GLOBAL FILE  
  ############################################################################
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version      20.06.2019
  # First version     20.06.2019
  # --------------------------------------------------------------------------
  
  
  
  # --------------------------------------------------------------------------
  # REQUIRE
  # --------------------------------------------------------------------------
  
  require(ggridges)
  require(ggplot2)
  require(viridis)
  require(hrbrthemes)
  
  require(boot)
  require(BBSolz)
  require(visolz)
  require(slolz)
  require(DAARC)
  # require(shinyalert)
  # source("utility_functions.R")
  
  wd <- "R:/Asset_Management/R/Shiny/Signal_Testing/Data/"
  global_data <- readRDS( file = paste0(wd, "data.rds") )
  lMuSigma <- global_data$lMuSigma
  DAA <- global_data$DAA
 

  
  
  
  # --------------------------------------------------------------------------
  # CUSTOM FUNCTIONS
  # --------------------------------------------------------------------------

  # customData
  # customPlot
  # customPlot2
  # conditionalDensityPlot
  # signalPlot
  # bbsPlot
  # bbsPlot2
  # importanceSamplingPlot
  # importanceSamplingPlot2
  # fcmdmPlot
  # fcmdmPlot2
  # logitRegressionPlot
  # logisticRegressionPlot2
  # backtestTradingPlot
  # backtestTradingPlot2
  # pasPlot
  # customTable
  

  
  
  
  
  
  # --------------------------------------------------------------------------
  customData <- function( input )
  {

    # Dependend variable
    y <- switch( input$dependent_variable,
                 "BM Returns" = DAA$data[isWeekday(rownames(DAA$data)), ],
                 "Bear / Bull States" = (1 + global_data$BBS) / 2,
                 "Fuzzy Bear / Bull States" = global_data$states_fuzzy )
    
    # Independent variable
    sensor_name <- tolower(input$sensor[1])
    Sensor <- DAA$sensors[[sensor_name]]
    signal_name <- input$signal[1]
    if ( signal_name == "signal" & sensor_name != "asyturb" ) {   # ~~~~~~~~~~~~~ 
      signal <- Sensor$getSignal()
    } else {
      signal <- Sensor$signal$base[ ,signal_name]
    }
    
    TD <- trainingData(Y_train = y,
                       X_train = signal)
    ans <- list(y = TD$Y_train,
                signal = TD$X_train)
   
    if ( input$what == "In/Out Backtest" ) {
        
      n_lag <- input$n_lag
      tc <- ifelse( isTRUE(input$tc), 0.0045, 0 )
      y <- ans$y
      signal <- ans$signal
      sig <- lag( (1 + sign(signal)) / 2, 0, trim = TRUE )
      test <- signalTesting.byTrading(X = y,
                                      sig = sig,
                                      n_lag = n_lag,
                                      tc = tc)
      X_bt <- na.omit(cbind(BM = y, Backtest = test))
      ans$sig <- sig
      ans$X_bt <- X_bt
    }
    
    if ( input$what %in% c("Importance Sampling",
                           "Bear/Bull Classification",
                           "Logistic Regression") ) {
      n_lag <- input$n_lag
      apply_roll_y <- list(fun = NULL, width = NULL, by = NULL)
      if ( input$what == "Importance Sampling" & input$n_agg > 0 ) {
        # FUN <- function(X) { exp(apply(log(1 + X), 2, sum)) - 1 }
        FUN <- function(X) { apply(X, 2, sum) }
        apply_roll_y = list(fun = FUN,
                            width = input$n_agg,
                            by = input$n_agg)
      }
      if ( input$what %in% c("Bear/Bull Classification",
                             "Logistic Regression") ) {
        
        y_tmp <- DAA$data[isWeekday(rownames(DAA$data)), ]
        y_level <- log(cumulated(exp(y_tmp)-1, "discrete"))
        BBS <- bbs(y_level, 
                   minphase = input$minphase, 
                   mincycle = input$mincycle, 
                   e = 0)
        global_data$BBS <<- BBS
        
        if ( input$what == "Logistic Regression") {
          reg <- regression( Y_train = (1 + global_data$BBS) / 2,
                             X_train = ans$signal,
                             n_lag = n_lag,
                             type = "logit" )
          ans$reg <- reg
          ans$signal <- reg$y
          ans$y <- DAA$data[isWeekday(rownames(DAA$data)), ]
        }
      }
      lBoot <- importanceSampling(Y_train = ans$y,
                                  X_train = ans$signal,
                                  n_lag = n_lag,
                                  apply_roll_y = apply_roll_y)
      ans$lBoot <- lBoot
    }
    
    if ( input$what == "FCMDM" ) {
      
      cm <- cluster.cmeans(x = ans$signal,  
                           centers = input$n_centers,
                           method = "cmeans", 
                           m = 2)
      memb <- cm$membership
      pmat <- apply( memb, 2, function(x) { x / sum(x) } )
      
      ans$cm <- cm
      ans$memb <- memb
      ans$pmat <- pmat

    }
    if ( input$what == "Performance After Signal" ) {
      
      # sigAggFUN <- parse( text = "apply(X, 2, sum)" )
      sigAggFUN <- NULL
      lPAS <- pas(X = DAA$data[isWeekday(rownames(DAA$data)), ], 
                  sig = ans$signal, 
                  n_lags = c(input$n_lag_pas, 1), 
                  compounding = "discrete",
                  sigAggFUN = sigAggFUN )
      ans$lPAS <- lPAS
    }
    
    return( ans )
  }
  
  
  # --------------------------------------------------------------------------
  customPlot <- function( input, data )  
  {

    if ( input$what == "Signal" ) {
      signalPlot( input, data )
    }
    if ( input$what == "In/Out Backtest" ) {
      backtestTradingPlot( input, data )
    }
    if ( input$what == "Bear/Bull Classification" ) {
      bbsPlot( input, data )
    }
    if ( input$what == "Performance After Signal" ) {
      pasPlot(input, data )
    }
    if ( input$what == "Importance Sampling" ) {
      importanceSamplingPlot( input, data )
    }
    if ( input$what == "FCMDM" ) {
      fcmdmPlot( input, data )
    }
    if ( input$what == "Logistic Regression" ) {
      logisticRegressionPlot( input, data )
    }
    if ( input$what == "Conditional Density" ) {
      conditionalDensityPlot( input, data )
    }
    
  }
  
  # --------------------------------------------------------------------------
  customPlot2 <- function( input, data )  
  {
    if ( input$what == "Signal" ) {
      signalPlot2( input, data )
    }
    if ( input$what == "In/Out Backtest" ) {
      backtestTradingPlot2( input, data )
    }
    if ( input$what == "Bear/Bull Classification" ) {
      bbsPlot2( input, data )
    }
    if ( input$what == "Performance After Signal" ) {
      pasPlot(input, data )
    }
    if ( input$what == "Importance Sampling" ) {
      importanceSamplingPlot2( input, data )
    }
    if ( input$what == "FCMDM" ) {
      fcmdmPlot2( input, data )
    }
    if ( input$what == "Logistic Regression" ) {
      logisticRegressionPlot2( input, data )
    }
    
  }
  
  # --------------------------------------------------------------------------
  conditionalDensityPlot <- function( input, data )
  {
    TD <- trainingData(Y_train = data$y,
                       X_train = data$signal)
    y <- TD$Y_train
    x <- TD$X_train
    x_margin <- seq(from = min(x), to = max(x), length.out = 10)
    fit <- hdrcde::cde(y = y, 
                       x = x, 
                       deg = 1,
                       x.margin = x_margin, 
                       nymargin = 100)
    
    plot(fit)
  }
  
  # --------------------------------------------------------------------------
  signalPlot <- function( input, data )
  {
   
    # tmp <- cbind(data$y, log(cumulated(data$y, "discrete")), data$signal)
    # colnames(tmp) <- c("BM", "BM_cumulated", paste0(input$sensor, "_", input$signal)) 
    # HC <- hc_lineChart(tmp,
    #                    type = "stock",
    #                    compare = "none",
    #                    title = "Signal",
    #                    ylab = "Level",
    #                    connectNulls = TRUE)
    # return( HC )
    
    IS <- importanceSampling(Y_train = data$y,
                             X_train = data$signal,
                             X_eval = input$x_eval)
    w <- timeSeries( matrix(IS[[1]]$weights, ncol = 1), rownames(IS[[1]]$data) )
    Boot <- boot(data$y, statistic = meanBoot, R = 10^3)
    
    # tmp <- cbind(data$signal,
    #              w,
    #              data$y, 
    #              log(cumulated(data$y, "discrete")))
    # colnames(tmp) <- c(paste0(input$sensor, "_", input$signal), "CPMF", "BM", "BM_cumulated")
    # plot(tmp)
    
    
    shades_of_grey <- grey.colors(n = length(w))
    idx <- 1
    color <- "tomato"
    par(mfrow = c(2, 2))
    plot(data$signal, main = "Signal")
    abline(h = IS[[idx]]$X_eval, col = color)
    lines(data$signal)
    plot(w, ylab = "CPMF", main = "Conditional Probability Mass Function")
    abline(v = time(w), col = rev(shades_of_grey)[rank(as.numeric(w))])
    lines(w)
    plot(IS[[idx]]$data, col = "white", main = "Dependent Variable")
    points(IS[[idx]]$data, col = rev(shades_of_grey)[rank(as.numeric(w))], pch = 19)
    abline(h = 0, col = "lightgrey")
    lDens <- apply( cbind(Boot$t, IS[[1]]$t), 2, density )
    plot( lDens[[1]], main = "Conditional (red) and Unconditional (black)\n Distribution of Mean Estimator",
          ylim = range(unlist(lapply(lDens, FUN = function(x) { x$y }))),
          xlim = range(unlist(lapply(lDens, FUN = function(x) { x$x }))) )
    lines( lDens[[2]], col = color )
    abline(v = mean(data$y))
    abline(v = mean(IS[[1]]$t), col = 2)
    
  }
  # --------------------------------------------------------------------------
  signalPlot2 <- function( input, data )
  {
    # tmp <- cbind(data$y, log(cumulated(data$y, "discrete")), data$signal)
    # colnames(tmp) <- c("BM", "BM_cumulated", paste0(input$sensor, "_", input$signal))
    # plot( tmp, plot.type = "multiple" )
  }
  
  
  # --------------------------------------------------------------------------
  bbsPlot <- function( input, data )
  {
    plot(data$y_level)
    abline(v = time(global_data$BBS)[which(global_data$BBS == 1)], col = "lightgreen")
    abline(v = time(global_data$BBS)[which(global_data$BBS == -1)], col = "tomato")
    lines(data$y_level)
  }
  
  # --------------------------------------------------------------------------
  bbsPlot2 <- function( input, data )
  {
    importanceSamplingPlot( input = input,
                            data = data )
  }
  
  
  # --------------------------------------------------------------------------
  importanceSamplingPlot <- function( input, data )
  {

    Boot <- boot(data = data$lBoot[[1]]$data, 
                 statistic = meanBoot, 
                 R = 10^3)  
    
    lDens <- lapply( data$lBoot, function(x) { density(x$t) } )
    colors <- fBasics::divPalette(n = length(lDens), "Spectral")
    plot( lDens[[1]],
          xlim = range(unlist(lapply(lDens, function(x) { x$x }))),
          ylim = range(unlist(lapply(lDens, function(x) { x$y }))),
          main = "Conditional estimator distributions")
    for ( i in 1:length(lDens) ) {
      lines(lDens[[i]], col = colors[i])
    }    
    lines(density(Boot$t))
    legend("topleft", c(names(data$lBoot), "ordinary bootstrap"), lwd = 2, 
           col = c(colors, "black"), text.col = c(colors, "black"), bty = "n")
    
  }
  
  
  # --------------------------------------------------------------------------
  importanceSamplingPlot2 <- function( input, data )
  {
    mu <- unlist(lapply(data$lBoot, FUN = function(x) { mean(x$t) }))
    colors <- fBasics::divPalette(n = length(mu), "RdYlGn")
    barplot(mu, col = colors)
  }
  
  # --------------------------------------------------------------------------
  fcmdmPlot <- function( input, data )
  {
    colors <- fBasics::divPalette(n = ncol(data$pmat), "RdYlGn")
    if ( input$fcmdm_which_plot == "Signal and cluster center" ) {
      centers <- sort(data$cm$cm$centers)
      plot( data$signal )
      for ( k in 1:length(centers) ) {
        abline(h = centers[k], col = colors[k])
      }
    } else if ( input$fcmdm_which_plot == "Membership function" ) {
      plot( data$memb, col = colors )
    }
    
  }
  # --------------------------------------------------------------------------
  fcmdmPlot2 <- function( input, data )
  {

    # easyGgplot2::ggplot2.density(data = data$DF, 
    #                              xName = "mean", 
    #                              groupName = 'method', 
    #                              fillGroupDensity = TRUE,
    #                              addMeanLine = TRUE,
    #                              legendPosition = "top" )

    # lDens <- lapply( data$lDF, FUN = density )
    
    n_sim <- 10^3
    if ( input$fcmdm_which_alpha == "pmat" ) {
      mat <- data$pmat
      mat <- cbind(rep(1/nrow(mat), nrow(mat)), mat)
    } else if ( input$fcmdm_which_alpha == "pmat_times_nrow" ) {
      mat <- data$pmat * nrow(data$pmat)
      mat <- cbind(rep(1, nrow(mat)), mat)
    } else {
      mat <- data$memb
      mat <- cbind(rep(1, nrow(mat)), mat)
    }
    
    
    
    # Posterior distributions
    y <- as.numeric(data$y[rownames(data$pmat), ])
    m <- matrix(NA, nrow = n_sim, ncol = ncol(mat),
                dimnames = list(NULL, c("const", paste0("cluster_", 1:(ncol(mat)-1)))))
    for ( k in 1:ncol(mat) ) {
      samples_tmp <- MCMCpack::rdirichlet(n = n_sim, alpha = mat[ ,k])
      m[ ,k] <- apply(samples_tmp, 1, function(x) { y %*% x } )
    }
  
    lDens <- apply( m, 2, density )
    colors <- c("black", fBasics::divPalette(n = length(lDens)-1, "RdYlGn"))
    plot( lDens[[1]],
          xlim = range(unlist(lapply(lDens, function(x) { x$x }))),
          ylim = range(unlist(lapply(lDens, function(x) { x$y }))),
          main = "Posterior distributions")
    for ( i in 1:length(lDens) ) {
      lines(lDens[[i]], col = colors[i])
    }    
    legend("topleft", names(lDens), lwd = 2, 
           col = colors, text.col = colors, bty = "n")
    
    
  }
  
  
  # --------------------------------------------------------------------------
  logisticRegressionPlot <- function( input, data )
  {
    plot( data$reg$z, data$reg$y )
  }
  # --------------------------------------------------------------------------
  logisticRegressionPlot2 <- function( input, data )
  {
    
    IS <- importanceSampling(Y_train = data$y,
                             X_train = data$signal,
                             X_eval = input$x_eval_logit)
    w <- timeSeries( matrix(IS[[1]]$weights, ncol = 1), rownames(IS[[1]]$data) )
    Boot <- boot(data$y, statistic = meanBoot, R = 10^3)

    shades_of_grey <- grey.colors(n = length(w))
    idx <- 1
    color <- "tomato"
    par(mfrow = c(2, 2))
    plot(data$signal, main = "Signal")
    abline(h = IS[[idx]]$X_eval, col = color)
    lines(data$signal)
    plot(w, ylab = "CPMF", main = "Conditional Probability Mass Function")
    abline(v = time(w), col = rev(shades_of_grey)[rank(as.numeric(w))])
    lines(w)
    plot(IS[[idx]]$data, col = "white", main = "Returns")
    points(IS[[idx]]$data, col = rev(shades_of_grey)[rank(as.numeric(w))], pch = 19)
    abline(h = 0, col = "lightgrey")
    lDens <- apply( cbind(Boot$t, IS[[1]]$t) * 252, 2, density )
    plot( lDens[[1]], main = "Conditional (red) and Unconditional (black)\n Distribution of Mean Estimator",
          ylim = range(unlist(lapply(lDens, FUN = function(x) { x$y }))),
          xlim = range(unlist(lapply(lDens, FUN = function(x) { x$x }))) )
    lines( lDens[[2]], col = color )
    abline(v = mean(data$y) * 252)
    abline(v = mean(IS[[1]]$t) * 252, col = 2)
  
  }
  
  # --------------------------------------------------------------------------
  backtestTradingPlot <- function( input, data )
  {
    
    plot( log(cumulated(data$X_bt, "discrete")), plot.type = "single" )
    abline(v = time(data$sig)[which(data$sig == 1)], col = 3)
    abline(v = time(data$sig)[which(data$sig == 0)], col = 2)
    lines(log(cumulated(data$X_bt[ ,1], "discrete")), col = 1)
    lines(log(cumulated(data$X_bt[ ,2], "discrete")), col = 4)
    legend("topleft", lwd = 2, colnames(data$X_bt),
           col = c(1, 4), text.col = c(1, 4), bty = "n")

    # HC <- hc_lineChart(cumulated(data$X_bt, "discrete"),
    #                    type = "stock",
    #                    compare = "percent",
    #                    title = "Performance",
    #                    ylab = "Cumulative Return",
    #                    connectNulls = TRUE)
    # return( HC )

  }
  # --------------------------------------------------------------------------
  backtestTradingPlot2 <- function( input, data )
  {
    plot.simTS( as.simTS(data$X_bt), col = c(1, 4) )
  }
  
  
  
  # --------------------------------------------------------------------------
  pasPlot <- function( input, data )  
  {
    
    lPAS <- data$lPAS
    
    idx <- 1
    Y_train <- lPAS[[idx]]$Y_train
    X_train <- lPAS[[idx]]$X_train
    reg <- regression( Y_train = Y_train,
                       X_train = X_train,
                       type = "ols")
    summary(reg$reg)
    
    plot( x = as.numeric(X_train), y = as.numeric(Y_train), pch = 19 )
    abline(h = 0, col = "grey")
    abline(v = 0, col = "grey")
    abline(reg$reg, col = 2)
    
    
  }
  
  
  
  # --------------------------------------------------------------------------
  customTable <- function( input, data )
  {
    
    if ( input$what == "Bear/Bull Classification" ) {
      
      
    }
    
  }
  
  
  
  
  
  
