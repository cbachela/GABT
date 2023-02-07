  
  
  ############################################################################
  ### PREPARE DATA FOR SHINY SIGNAL TESTING
  ############################################################################
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version      20.06.2019
  # First version     20.06.2019
  # --------------------------------------------------------------------------
  
  
  require(DAARC)
  
  
  DAA <- loadSensor( sensor_name = "daa" )
  DAA$updateData()
  
  WhiteNoise <- loadSensor(sensor_name = "whitenoise")
  DAA$sensors[["whitenoise"]] <- WhiteNoise
  
  names(DAA$sensors)
  

  ls(DAA$sensors[[1]]$posterior)
  musigma <- DAA$sensors[[1]]$evaluatePosterior()
  
  
  lMuSigma <- lapply(DAA$sensors, FUN = function(x) { x$evaluatePosterior() })
  names(lMuSigma) <- names(DAA$sensors)
  
  
  
  
  sensor_names <- c("momentum", 
                    "vrp", 
                    # "stlfsi", 
                    "bcpolz", 
                    "turbulence",
                    "eigen")
  signals <- NULL
  for ( j in seq(along = sensor_names) ) {
    
    sensor_name <- sensor_names[j]
    Object <- loadSensor(sensor_name = sensor_name)
    Object$update()
    DAA$sensors[[sensor_name]] <- Object
    signals <- cbind(signals, 
                     Object$getSignal())
  }  
  
  DAA$updateData()
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # Bear / Bull Labeling (BBS Algo)
  # --------------------------------------------------------------------------
  
  mincycle <- 5 * 4 * 7 
  minphase <- 5 * 4 * 2
  k_peak <- k_trough <- l_peak <- l_trough <- minphase * 1
  e <- 0
  
  X_bm <- DAA$data[isWeekday(rownames(DAA$data)), ]
  X_level <- log(cumulated(X_bm, "discrete"))
  BBS <- bbs(X = X_level,
             mincycle = mincycle, 
             minphase = minphase,
             k.peak = k_peak,
             l.peak = l_peak,
             k.trough = k_trough,
             l.trough = l_trough,
             e = e)
  states <- (1 + BBS) / 2
  
  
  plot(X_level)
  abline(v = time(BBS)[which(BBS == 1)], col = "green")
  abline(v = time(BBS)[which(BBS == -1)], col = "red")
  lines(X_level)
  
  
  
  # --------------------------------------------------------------------------
  # ASYMETRIC TURBULENCE SIGNAL ON SENSOR SERIES
  # --------------------------------------------------------------------------
  
  
  
  # Bear markets
  p_bear <- (1 - BBS) / 2
  TD <- trainingData(Y_train = p_bear, 
                     X_train = signals)
  mu_bear <- t(TD$X_train) %*% (TD$Y_train / sum(TD$Y_train))
  scnd_mom = t(TD$X_train) %*% (TD$X_train * ((TD$Y_train / sum(TD$Y_train)) %*% matrix(1, 1, ncol(TD$X_train))) )
  scnd_mom = (scnd_mom + t(scnd_mom)) / 2
  sigma_bear = scnd_mom - mu_bear %*% t(mu_bear)
  
  
  # Bull markets
  p_bull <- (1 + BBS) / 2
  TD <- trainingData(Y_train = p_bull, 
                     X_train = signals)
  mu_bull <- t(TD$X_train) %*% (TD$Y_train / sum(TD$Y_train))
  scnd_mom = t(TD$X_train) %*% (TD$X_train * ((TD$Y_train / sum(TD$Y_train)) %*% matrix(1, 1, ncol(TD$X_train))) )
  scnd_mom = (scnd_mom + t(scnd_mom)) / 2
  sigma_bull = scnd_mom - mu_bull %*% t(mu_bull)
  
  md_all <- mahalanobis.ts(X = TD$X_train)
  md_bear <- mahalanobis.ts(X = TD$X_train, center = mu_bear, covmat = sigma_bear)
  md_bull <- mahalanobis.ts(X = TD$X_train, center = mu_bull, covmat = sigma_bull)
  MD <- cbind(all = md_all,
              bear = md_bear,
              bull = md_bull)
  MD_scl <- scale(MD, FALSE, TRUE)
  md_delta <- MD_scl[ ,"bear"] - MD_scl[ ,"bull"]
  md_delta_ema <- ema(MD_scl[ ,"bear"], 0.05) - ema(MD_scl[ ,"bull"], 0.05)
  
  sig_tmp <- cbind(md_delta, md_delta_ema)
  colnames(sig_tmp) <- c("signal", "signal_ema")
  
  DAA$sensors$asyturb <- list(data = list(),
                              signal = list(base = sig_tmp))
  
  
  
  # --------------------------------------------------------------------------
  # Fuzzy Bear / Bull Labeling (BBS Algo)
  # --------------------------------------------------------------------------
  
  # Run algo on daily data, then aggregate fuzzy signal
  
  X_level <- log(cumulated(X_bm, "discrete"))
  minphase <- 5 * 4 * 2
  k_peak <- k_trough <- l_peak <- l_trough <- minphase * 2
  e <- 0
  
  mincycle_vec <- round(seq(from = 5*4, to = 5*4*7, length.out = 10))
  minphase_vec <- round(seq(from = 5, to = 5*4, length.out = 10))
  
  states_mat <- X_level[ ,rep(1, length(minphase_vec) * length(mincycle_vec))] * NA
  k <- 1
  for ( j in seq(along = mincycle_vec) ) {
    for ( i in seq(along = minphase_vec) ) {
      mincycle <- mincycle_vec[j]
      minphase <- minphase_vec[i]
      BBS_tmp <- bbs(X = X_level,
                     mincycle = mincycle, 
                     minphase = minphase,
                     # k.peak = k_peak,
                     # l.peak = l_peak,
                     # k.trough = k_trough,
                     # l.trough = l_trough,
                     e = e)
      states_mat[ ,k] <- (1 + BBS_tmp) / 2
      k <- k + 1
    }
  }
  
  states_fuzzy_raw <- apply(states_mat, 1, mean)
  # states_fuzzy_ema <- ema(states_fuzzy_raw, alpha = 0.05)
  # states_fuzzy_raw_m <- applyRoll(Data = states_fuzzy_raw,
  #                                 Width = 21,
  #                                 By = 1,
  #                                 FUN = mean)

  states_fuzzy <- states_fuzzy_raw
  
  
  
  par(mfrow = c(2, 1))
  
  states_fuzzy_unique <- sort(unique(states_fuzzy))
  colors <- fBasics::divPalette(n = length(states_fuzzy_unique), "RdYlGn")
  X_level <- log(cumulated(X_bm, "discrete"))
  
  plot(X_level, main = "Crisp Bear / Bull States")
  abline(v = time(BBS)[which(BBS == 1)], col = colors[length(colors)])
  abline(v = time(BBS)[which(BBS == -1)], col = colors[1])
  lines(X_level, col = "blue")
  
  plot(X_level, main = "Fuzzy Bear / Bull States")
  for ( i in seq(along = states_fuzzy_unique) ) {
    abline( v = time(states_fuzzy)[which(states_fuzzy == states_fuzzy_unique[i])],
            col = colors[i] )
  }
  lines(X_level, col = "blue")
  
  
  dev.off()
  plot( states_fuzzy )
  
  
  
 
  
  
  # --------------------------------------------------------------------------
  # SAVE
  # --------------------------------------------------------------------------
  
  env <- new.env()
  env$DAA <- DAA
  env$lMuSigma <- lMuSigma
  env$BBS <- BBS
  env$states_fuzzy <- states_fuzzy
  saveRDS(object = env, file = "R:/Asset_Management/R/Shiny/Signal_Testing/Data/data.rds")
  
  

  
  
  source("R:/Asset_Management/R/Shiny/Signal_Testing/global.R")
  

  input <- list(what = "In/Out Backtest",
                sensor = "VRP",
                dependent_variable = "BM Returns",
                signal = "signal",
                n_lag = 1,
                n_lag_pas = 40)
  
  # debugonce(customData)
  data <- customData( input )
  plot( data$signal )
  
  input$minphase <- 5 * 4 * 2
  input$mincycle <- 5 * 4 * 5
  bbsPlot( input, data )
  
  # debugonce(customPlot)
  # debugonce(pasPlot)
  # debugonce(logisticRegressionPlot)
  # debugonce(regression)
  customPlot( input = input, data = data )  
  
  
  
  test <- function( input, data ) 
  {
    
    # Parameters
    n_lag <- input$n_lag
    # width <- by <- input$n_agg
    width <- by <- 1
    
    # debugonce(trainingData)
    FUN <- function(X) { exp(apply(log(1 + X), 2, sum)) - 1 }
    lTD <- trainingData(Y_train = data$y,
                        X_train = data$signal,
                        n_lag = n_lag,
                        apply_roll_y = list(fun = FUN, 
                                            width = width, 
                                            by = by))
    
    X_eval_vec <- seq(from = min(lTD$X_train), to = max(lTD$X_train), length.out = 10)
    Dist <- timeSeries( t( apply( lTD$X_train, 1, function(x) { abs(x - X_eval_vec) }) ), time(lTD$X_train) )
    FUN <- function(x) 
    {
      wghts = NULL
      if ( length(unique(x)) > 1 ) {
        x <- x / max(x)
        wghts <- as.numeric( (1 - x) / sum(1 - x) )
      } 
      return( wghts )
    }
    wmat <- apply(Dist, 2, FUN)


    idx <- 10
    dates <- intersect( rownames(X), rownames(wmat))
    y <- X[dates, ] * wmat[dates, idx] * nrow(wmat)

    plot(y)
    dd <- drawdowns(y)

    plot(dd)

    plot(density(abs(dd)))
    
 
    
    
    
  }
  
  
  
  
  
  
  
