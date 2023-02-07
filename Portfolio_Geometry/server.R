  
  
  ##########################################################
  ### SHINY FRONT END - PORTFOLIO GEOMETRY - SERVER SIDE ###
  ##########################################################
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version      20.04.2018
  # First version     20.04.2018
  # --------------------------------------------------------------------------


  # CUSTOM FUNCTIONS:
  #
  # cutomData
  # customTable
  # customPlotUtility
  # customPlotPortfolio
  # efficientFrontierPlot
  # efficientFrontierPlot_3AssetSpace
  # colorCoding
  # custom2DPlot
  # custom3DPlot
  
  
  # --------------------------------------------------------------------------
  # CUSTOM FUNCTIONS
  # --------------------------------------------------------------------------
  
  # --------------------------------------------------------------------------
  customData <- function( input )
  {
    asset_name <- input$assets[1]
    # x_tmp <- window(X[ ,asset_name], input$startdate, input$enddate)
    x_tmp <- X[ ,asset_name]
    
    # FUN <- match.fun( paste0(input$utility_function, "Utility"))
    
    if ( input$utility_function == "Minimum-Variance" ) {
      y <- minvarianceUtility(x_tmp) * (-1)
    } else if ( input$utility_function == "Mean-Variance" ) {
      y <- mvUtility(x_tmp, 
                     lambda = input$lambda) 
    } else if ( input$utility_function == "Kahneman-Tversky" ) {
      y <- ktUtility(x_tmp, 
                     RP = input$RP, 
                     a = input$a, 
                     b = input$b)
    } else if ( input$utility_function == "Omega" ) {
      y <- omegaRatio(x = x_tmp, cdf = density(x)$y, tauvec = x_tmp)
      y[is.na(y)] <- 0
      
    } else if ( input$utility_function == "Entropy" ) {
      y <- entropyUtility(x_tmp, lambda = input$lambda)
      x_tmp <- mean(x_tmp)
    } else {
      stop("@ customData: Invalid choice of utility function")
    }
    
    ans <- list(x = x_tmp,
                y = y)
    return( ans )
  }
  
  # --------------------------------------------------------------------------
  customTable <- function(input)
  {
    # X_tmp <- window(X, input$startdate, input$enddate)
    # X_tmp <- X
    # lStats <- descStats(X_tmp, 
    #                     descStatsSpec(what = c("basics", "maxDD")))
    lStats <- env$lStats
    ans <- lStats
    return( ans )
  }
  
  # --------------------------------------------------------------------------
  customPlotUtility <- function(input, object = NULL)
  {
    what <- tolower(input$what)
    if ( is.null(object) ) {
      return( NULL )
    } else {
      if ( what == 'utility' ) {
     
        plot( hist(object$x, breaks = round(nrow(X) / 6)), main = "", xaxt = "n", yaxt = "n", 
              ylab = "", xlab = "", col = "lightblue")
        # lines(density(object$x), lwd = 2, col = "steelblue")
        axis(side = 4)
        mtext("Frequency", side = 4, line = 3, cex = 1.1) 
        par(new = TRUE)
        plot( as.numeric(object$x), 13 * as.numeric(object$y),
              col = "red",
              main = "",
              ylab = "Utility",
              xlab = "Returns")
        # axis(side = 1)
        # mtext(side = 1, line = 3, 'Utility')
        abline(h = 0, col = "darkgrey")
        abline(v = 0, col = "darkgrey")
        # legend("topleft", col = color, text.col = color, lty=1, lwd=2, cex=1.1, 
        #        legend=c("MSCI World Index", "200 Day Moving Average"), bty="n")

      } else {
        return( NULL )
      }
    }
  }


  
  # --------------------------------------------------------------------------
  customPlotPortfolio <- function(input)
  {

    # sim_tmp <- window(sim, input$startdate, input$enddate)
    sim_tmp <- sim

    if ( input$utility_function == "Minimum-Variance" ) {
      U <- apply( 13 * sim_tmp, 2, minvarianceUtility )
    } else if ( input$utility_function == "Mean-Variance" ) {
      U <- apply( 13 * sim_tmp, 2, mvUtility, lambda = input$lambda )
    } else if ( input$utility_function == "Kahneman-Tversky" ) {
      U <- apply( 13 * sim_tmp, 2, ktUtility, RP = input$RP, a = input$a, b = input$b )
    } else {
      stop("@ customPlotPortfolio: Invalid choice of utility function")
    }

    tmp <- apply(U, 2, mean)
    idx <- which(tmp == max(tmp))

    ans <- setNames( S[idx, ], colnames(S_lo) )
    barplot( ans )

  }
  
  
  
  # --------------------------------------------------------------------------
  efficientFrontierPlot <- function(input)
  {
    # Efficient frontier in risk-return space
    
    wmat_ls <- attr(fp_ls, "Weights")
    w_ls <- getWeights(attr(fp_ls, "GMV"))
    tmp <- apply(wmat_ls, 1, function(x) { sum(abs(x - w_ls))})
    idx_ls <- which(tmp == 0)
    
    wmat_lo <- attr(fp_lo, "Weights")
    w_lo <- getWeights(attr(fp_lo, "GMV"))
    tmp <- apply(wmat_lo, 1, function(x) { sum(abs(x - w_lo))})
    idx_lo <- which(tmp == 0)
    
    plot( fp_ls[ ,"targetRisk"], 
          fp_ls[ ,"targetReturn"], pch = 19,
          xlab = "Risk", ylab = "Return", col = "grey35", type = "o")
    points( fp_lo[ ,"targetRisk"], 
            fp_lo[ ,"targetReturn"], pch = 19, col = "tomato", type = "o")
    points(fp_ls[-c(1:idx_ls), "targetRisk"], 
           fp_ls[-c(1:idx_ls), "targetReturn"], pch = 19)
    points(fp_lo[-c(1:idx_lo), "targetRisk"], 
           fp_lo[-c(1:idx_lo), "targetReturn"], pch = 19, col = 2)
    points(fp_ls[idx_ls, "targetRisk"], 
           fp_ls[idx_ls, "targetReturn"], pch = 19, cex = 2)
    points(fp_lo[idx_lo, "targetRisk"], 
           fp_lo[idx_lo, "targetReturn"], pch = 19, cex = 2, col = 2)
    
    
    }
  

  
  # --------------------------------------------------------------------------
  efficientFrontierPlot_3AssetSpace <- function(input)
  {
    # Efficient frontier plot in asset space
    wmat_ls <- attr(fp_ls, "Weights")
    w_ls <- getWeights(attr(fp_ls, "GMV"))
    tmp <- apply(wmat_ls, 1, function(x) { sum(abs(x - w_ls))})
    idx_ls <- which(tmp == 0)
    
    wmat_lo <- attr(fp_lo, "Weights")
    w_lo <- getWeights(attr(fp_lo, "GMV"))
    tmp <- apply(wmat_lo, 1, function(x) { sum(abs(x - w_lo))})
    idx_lo <- which(tmp == 0)
    
    covmat <- cov(X)
    Names <- colnames(X)
    FUN <- function(x) { as.numeric( t(x) %*% covmat %*% x ) * (-1) }
    obj_values <- apply(S_lo, 1, FUN = FUN)
    color <- colorCoding(Set = S_lo, 
                         obj_values = obj_values)
    plot3d(S_ls[ ,1], S_ls[ ,2], S_ls[ ,3], 
           xlab = Names[1], ylab = Names[2], zlab = Names[3], 
           col = "grey", pch = 19, cex = 2, size = 1)
    points3d(S_lo[ ,1], S_lo[ ,2], S_lo[ ,3], col = color, size = 3)
    points3d(wmat_ls[c(1:idx_ls), 1], wmat_ls[c(1:idx_ls), 2], wmat_ls[c(1:idx_ls), 3], size = 5, col = "grey35")
    points3d(wmat_lo[c(1:idx_lo), 1], wmat_lo[c(1:idx_lo), 2], wmat_lo[c(1:idx_lo), 3], size = 7, col = "tomato")
    points3d(wmat_ls[-c(1:idx_ls), 1], wmat_ls[-c(1:idx_ls), 2], wmat_ls[-c(1:idx_ls), 3], size = 5, col = 1)
    points3d(wmat_lo[-c(1:idx_lo), 1], wmat_lo[-c(1:idx_lo), 2], wmat_lo[-c(1:idx_lo), 3], size = 7, col = 2)
    points3d(wmat_ls[idx_ls, 1], wmat_ls[idx_ls, 2], wmat_ls[idx_ls, 3], size = 12, col = 1)
    points3d(wmat_lo[idx_lo, 1], wmat_lo[idx_lo, 2], wmat_lo[idx_lo, 3], size = 15, col = 2)
  
  }
  

  
  # --------------------------------------------------------------------------
  colorCoding <- function(Set, 
                          obj_values = NULL, 
                          FUN, 
                          col = c("lightblue", "darkblue"),
                          n_steps = 5)
  {
    txt <- colorRampPalette(colors = col)(n_steps)
    n <- nrow(Set) / length(txt)
    colors <- rep(NA, nrow(Set))
    for ( i in seq(along = txt) ) {
      colors[((n*i)-n+1):(n*i)] <- rep(txt[i], n)
    }
    color <- rep(NA, nrow(Set))
    if ( is.null(obj_values) ) {
      obj_values <- apply(Set, 1, FUN)
    }
    color[order(obj_values)] <- colors
    
    return( color )
  }


  
  # --------------------------------------------------------------------------
  custom2DPlot <- function(input, object)
  {
     
    X_tmp <- X
    
    if ( input$what == "Utility" ) {
      customPlotUtility(input = input, 
                        object = object)
    
    } else if ( input$what == "Covariance" ) {
      # covmat <- lCovmat[[ input$shrinkage_intensity * 10 + 1 ]]
      covmat <- cov.shrinkage(X = X_tmp, 
                              lambda = input$shrinkage_intensity,
                              target_method = input$shrinkage_target)
      varcorImagePlot(covmat = covmat )
      
    } else if ( input$what == "Return Series" ) {
      plot( X_tmp )

    } else if (input$what == "Efficient Frontier" ) {
      efficientFrontierPlot(input = input)
      
    }
    
  }
  
  

  
  # --------------------------------------------------------------------------
  custom3DPlot <- function(input, object)
  {
    
    Names <- colnames(X)
    # X_tmp <- window(X[ ,1:3], input$startdate, input$enddate)
    X_tmp <- X
    # covmat <- lCovmat[[ input$shrinkage_intensity * 10 + 1 ]]
    covmat <- cov.shrinkage(X = X_tmp, 
                            lambda = input$shrinkage_intensity,
                            target_method = input$shrinkage_target)
    mu <- meanGeo(X_tmp)
    
    
    # Plot utility levels
    if ( input$what == "Utility" ) {
     
      clear3d()
      plot3d(S_ls[ ,1], S_ls[ ,2], S_ls[ ,3], 
             xlab = Names[1], ylab = Names[2], zlab = Names[3], 
             col = "grey", pch = 19, cex = 2, size = 1)

      # Colorcoding
      if ( input$utility_function == "Mean-Variance" ) {
        FUN <- function(x, lambda) 
        { 
          as.numeric( t(x) %*% mu ) - lambda * as.numeric( t(x) %*% covmat %*% x )
        }
        obj_values <- apply(S_lo, 1, FUN = FUN, lambda = input$lambda)
      
      } else if ( input$utility_function == "Minimum-Variance" ) {
        FUN <- function(x) { as.numeric( t(x) %*% covmat %*% x ) * (-1) }
        obj_values <- apply(S_lo, 1, FUN = FUN)
    
      } else if ( input$utility_function == "Kahneman-Tversky" ) {
        FUN <- function(x, RP, a, b) { ktUtility(x, RP = RP, a = a, b = b) }
        tmp <- apply(sim, 2, FUN = FUN, RP = input$RP, a = input$a, b = input$b)
        obj_values <- apply(tmp, 2, mean) 
        
      } else if ( input$utility_function == "Entropy" ) {
        FUN <- function(x, lambda) { entropyUtility(x, lambda) }
        obj_values <- apply(sim, 2, FUN = FUN, lambda = input$lambda)
      
      } else if ( input$utility_function == "Omega" ) {
        # FUN <- function(x, threshold) { as.numeric( omegaRatio(x, cdf = density(x, n = length(x))$y, tauvec = threshold)) } 
        # # obj_values <- apply(sim[ ,1:3], 2, FUN = FUN, threshold = input$threshold)
        # obj_values <- numeric(ncol(sim))
        # for (i in 1:ncol(sim)) {
        #   obj_values[i] <- FUN(sim[ ,i], threshold = input$threshold)
        # }
        obj_values <- numeric(ncol(sim))
        for ( i in 1:ncol(sim) ) {
          obj_values[i] <- omegaRatio(x = sim[ ,i], 
                                      cdf = cdf_mat[ ,i], 
                                      tauvec = input$threshold)
        }
        
      }
       
      
      color <- colorCoding(Set = S_lo, 
                           obj_values = obj_values)
      
      # Add colored long-only portfolios
      points3d(S_lo[ ,1], S_lo[ ,2], S_lo[ ,3], col = color, size = 4)
      
      
    # Plot covariance ellipsoid
    } else if ( input$what == "Covariance" ) {

      
      # Intersecting ellipse
      N <- ncol(X_tmp)
      y <- c(1, rep(0, N-1))
      v <- rep(1, N)
      g <- 1
      E <- ellipsoid(q = rep(0, N),
                     Q = covmat)
      H <- hyperplane(v = v, 
                      b = g)
      tmat <- alignVector(y = y, x = v)
      tmat_inv <- solve(tmat)
      z <- (g * tmat %*% v) / as.numeric( t(v) %*% v )     
      E_prime <- ePlus( E = eMult(E = E, Amat = tmat), 
                       bvec = -z )
      
      q_prime <- getCentre(E_prime)
      Q_prime <- getShape(E_prime)
      
      
      r_tilde <- q_prime[-1] + q_prime[1] * (solve(Q_prime[-1, -1]) %*% Q_prime[-1, 1])
      r_prime <- c(0, r_tilde)
      
      
      ones <- rep(1, N)
      covmat_inv <- solve(covmat)
      w_ls <- covmat_inv %*% ones / as.numeric( t(ones) %*% covmat_inv %*% ones)
      w_th <- twoAssetPortfolio(w1 = w_ls, 
                                w2 = rep(1/N, N), 
                                covmat = covmat, 
                                alpha = 0.05)$w
      
      w_th_prime <- (tmat %*% w_th - z)[-1]
      h <- t(w_th_prime - r_tilde) %*% Q_prime[-1, -1] %*% (w_th_prime - r_tilde)
                       
      
      R_tilde <- Q_prime[-1, -1] / as.numeric( h )
      R_prime <- cbind( rep(0, N), rbind(rep(0, N-1), R_tilde) )
      
      Y_prime <- runifE(S = R_tilde,
                        z_hat = r_tilde,
                        gamma_threshold = 1,
                        n_points = ncol(sim), 
                        b_pushy = FALSE)
      Y_prime <- cbind(0, Y_prime)
      Y <- t(apply(Y_prime, 1, function(x) { tmat_inv %*% (x + z) }))
      
      # Covariance ellipsoid
      Y_covmat <- runifE(S = covmat,
                         z_hat = rep(0, N),
                         gamma_threshold = as.numeric(t(w_th) %*% covmat %*% w_th),
                         n_points = ncol(sim), 
                         b_pushy = FALSE)
      
      
      # Colorcoding
      color <- colorCoding(Set = S_lo, 
                           FUN = function(x) { as.numeric( t(x) %*% covmat %*% x ) * (-1) } )
      
      # Plot
      clear3d()
      plot3d(S_ls[ ,1], S_ls[ ,2], S_ls[ ,3], 
             xlab = Names[1], ylab = Names[2], zlab = Names[3], 
             col = "lightgrey", pch = 19, cex = 2, size = 1.5)
      points3d(S_lo[ ,1], S_lo[ ,2], S_lo[ ,3], col = color, size = 3.5)
      # i <- input$shrinkage_intensity * 10 + 1
      # points3d(lY[[i]][ ,1], lY[[i]][ ,2], lY[[i]][ ,3], 
      #            col = "darkorange", size = 1)
      points3d(Y_covmat[ ,1], Y_covmat[ ,2], Y_covmat[ ,3], 
               col = "darkorange", size = ifelse(N == 3, 1, 2))
      points3d(Y[ ,1], Y[ ,2], Y[ ,3], 
               col = "darkorange", size = 4)
      points3d(w_ls[1], w_ls[2], w_ls[3], 
               col = "darkorange", size = 8)
      
      
    
    # Plot feasible set
    } else if ( input$what == "Constraints" ) {
      
  
      if ( input$dimension == 4 ) {
        env <- env_4d
        
      } else {
        env <- env_3d
      }
     
      S_lo <- env$S_lo
      S_ls <- env$S_ls
      
      if ( isTRUE(input$hide_S_ls) ) {
        S_ls <- NULL
      }
      
      if ( isTRUE(input$upper) ) {

        # Check upper bounds 
        idx <- cbind( (S_lo[ ,1] <= input$upper_x1),
                      (S_lo[ ,2] <= input$upper_x2),
                      (S_lo[ ,3] <= input$upper_x3))
      
        # Check group constraint
        idx <- cbind( idx, (apply(S_lo[ ,2:3], 1, sum) <= input$upper_x23) )
        idx <- apply(idx, 1, all)

      } else {
        idx <- 1:nrow(S_lo)
      }
      S_lo <- S_lo[idx, ]
      VaR <- env$VaR[idx]
      ES <- env$ES[idx]
      to <- env$to[idx]
      w_init <- env$w_init
      
      if ( input$cons == "Budget" ) {
        
        color <- "grey"
        S_lo <- NULL
        S_ls <- env$S_ls
        
      } else if ( input$cons == "Long Only" ) {
      
        color <- "orange"
        
      } else if ( input$cons == "Turnover" ) {
        
        color <- colorCoding(Set = S_lo, 
                             obj_values = to, 
                             col = c("darkred", "orange"),
                             n_steps = 5)
        
      } else if ( input$cons == "Tracking Error" ) {
      
        covmat <- cov(env$X)
        tr_err <- apply(S_lo, 1, function(x) { t(x - env$w_init) %*% covmat %*% (x - env$w_init) })
        
        color <- colorCoding(Set = S_lo, 
                             obj_values = tr_err, 
                             col = c("darkred", "orange"),
                             n_steps = 5)
        
      } else if ( input$cons == "VaR" ) {
        color <- colorCoding(Set = S_lo, 
                             obj_values = VaR, 
                             col = c("orange", "darkred"),
                             n_steps = 5)
        
      } else if ( input$cons == "Expected Shortfall" ) {
        color <- colorCoding(Set = S_lo, 
                             obj_values = ES, 
                             col = c("orange", "darkred"),
                             n_steps = 5)
      }
      
      # Plot
      clear3d()
      if ( !is.null(S_ls) ) {
        plot3d(S_ls[ ,1], S_ls[ ,2], S_ls[ ,3], 
               xlab = Names[1], ylab = Names[2], zlab = Names[3], 
               col = "lightgrey", pch = 19, cex = 2, size = 1.5)
        points3d(S_lo[ ,1], S_lo[ ,2], S_lo[ ,3], col = color, size = 3)
      } else {
        plot3d(S_lo[ ,1], S_lo[ ,2], S_lo[ ,3], 
               xlab = Names[1], ylab = Names[2], zlab = Names[3], 
               col = color, size = 5,
               ylim = c(-1, 1), xlim = c(-1, 1), zlim = c(-1, 1))
      }
      if ( input$cons %in% c("Turnover", "Tracking Error") ) {
        points3d(env$w_init[1], w_init[2], w_init[3], size = 7)
      }
      if ( isTRUE(input$add_diamond) ) {
        points3d(env$S_unc[env_3d$idx_to, 1],
                 env$S_unc[env_3d$idx_to, 2],
                 env$S_unc[env_3d$idx_to, 3], col = "red", size = 3)
      }
      if ( isTRUE(input$add_ellipsoid) ) {
        points3d(env$S_unc[env_3d$idx_te, 1],
                 env$S_unc[env_3d$idx_te, 2],
                 env$S_unc[env_3d$idx_te, 3], col = "red", size = 3)
      }
      

      
    } else if ( input$what == "Return Series" ) {
     
      
    } else if (input$what == "Efficient Frontier" ) {
      efficientFrontierPlot_3AssetSpace(input = input)
      
    }
    
  }
  
  
  
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
    
    
   
    # 2-D CHARTS
    output$plot_2d <- renderPlot({
      object <- reactiveData()
      custom2DPlot( input = input,
                    object = object )
    })
    output$title.plot_2d <- renderText({
      txt <- switch(input$what,
                    "Utility" = paste0(input$utility_function, " Utility Function"),
                    "Weights" = "Weights",
                    "Covariance" = "Correlation Image Plot",
                    "Return Series" = "Return Series",
                    "Efficient Frontier" = "Efficient Frontier in Risk-Return Space")
    })
    
    
    # 3-D CHARTS
    output$plot_3d <- renderPlot({
      if ( input$what == "Return Series" ) {
        plot(log(cumulated(X, "discrete")), plot.type = "single",
             ylab = "log(value)")
        legend("topleft", lwd = 2, colnames(X), 
               col = 1:ncol(X), text.col = 1:ncol(X), bty = "n")
      } else {
        custom3DPlot( input = input )
      }
    })
    output$title.plot_3d <- renderText({
      txt <- ifelse( input$what == "Return Series", "Cumulative Returns", "")
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

    
    # Return a summary print
    output$tmpPrint <- renderPrint({
      tbl <- customTable(input = input)
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
  