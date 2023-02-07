  
  
  ##################################################################
  ### SHINY FRONT END -  PORTFOLIO GEOMETRY - DATA PREPERATION   ###
  ##################################################################
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version      20.04.2018
  # First version     20.04.2018
  # --------------------------------------------------------------------------

  
  
  
  # --------------------------------------------------------------------------
  # REQUIRE
  # --------------------------------------------------------------------------
  
  require(RP)
  require(GPO)
  require(rgl)
  
  source("R:/Asset_Management/R/Shiny/Utility_Maximization/utility_functions.R")
  source("R:/Asset_Management/Research Projects/Random_Portfolios/random_portfolios/ellipsoidal_toolbox.R")
  
  


  # --------------------------------------------------------------------------
  # PARAMETERS
  # --------------------------------------------------------------------------
  
  # Which assets to choose
  idx <- c(4, 5, 6) # , 7)   #// Add or omit last one if you want 4D or 3D ~~~~~~~~~~~~~~~~~~~~~~~~

  
  # --------------------------------------------------------------------------
  # LOAD DATA
  # --------------------------------------------------------------------------

  Data <- GPO::getMSCIData(universe = "dm",
                           frqncy = "w",
                           day = "Wed",
                           compounding = "discrete",
                           heading = "ISO3",
                           ccy = "Local")

  # Number of random portfolios to be generated
  if ( length(idx) == 3 ) {
    nsim <- 10^4 * 4    
  } else if (length(idx) == 4 ) {
    nsim <- 10^4 * 10
  }
  
  X <- Data[ ,idx]
  N <- ncol(X)
  # colnames(X) <- c("Asset 1", "Asset 2", "Asset 3")
  alpha <- 0.05   # variance upper bound, i.e. sigma(w) <= sigma(w_gmv) * (1 + alpha)
  
  
  lStats <- descStats(X, descStatsSpec(what = c("basics", "maxDD")))
  lStats
  
  
  Names <- colnames(X)
  covmat <- cov(X)
  mu <- meanGeo(X)
  sigmaFUN <- function(x) { as.numeric( t(x) %*% covmat %*% x ) }
  
  # Generate constraints object  
  cons_ls <- cons_lo <- constraints(selection = Names[1:N])
  addConstraint(cons_ls) <- budgetConstraint()
  addConstraint(cons_ls) <- boxConstraint(name = "LongShort")
  addConstraint(cons_lo) <- budgetConstraint()
  addConstraint(cons_lo) <- boxConstraint(name = "LongOnly", 
                                          lower = setNames(rep(0, N), Names[1:N]), 
                                          upper = setNames(rep(1, N), Names[1:N]) )
  
  # Minimum variance portfolio
  GPS <- gps(Data = X)
  GPO <- gpo(GPS = GPS)
  w_gmv <- getWeights(GPO)
  
  # QEQW portfolio
  GPS <- gps(Data = X,
             Constraints = cons_lo,
             Covariance = covCtrl(method = "dcv"))
  GPO <- gpo(GPS = GPS)
  x0 = getWeights(GPO)
  
  # Random weights
  S_ls <- rpolytope(Constraints = cons_ls,
                    nsim = nsim,
                    algo = "hitandrun",
                    thin = 10^2,
                    x0 = x0)
  S_lo <- rpolytope(Constraints = cons_lo,
                    nsim = nsim * 1,
                    algo = "hitandrun",
                    thin = 10^2,
                    x0 = x0)
  # S_unc <- matrix( runif(nsim * N * 100, -0.5, 1), nsim * 100, N)
  
  
  
  # plot3d(S_unc[ ,1], S_unc[ ,2], S_unc[ ,3], col = color,
  #        xlab = Names[1], ylab = Names[2], zlab = Names[3], size = 3)
  
  
  
  
  # S_ls <- rpolytope(Constraints = cons_ls, 
  #                   nsim = nsim * 4, 
  #                   algo = "hitandrun",
  #                   thin = 10^2,
  #                   x0 = x0)
  # idx_lo <- apply(S_ls, 1, function(x) { all(x >= 0) })
  # S_ls <- S_ls[-idx_lo, ]
  
  
  
  
  # Add turnover constraint (linearized)
  w_init <- setNames( c( (1:ncol(X)) / sum(1:ncol(X)) ) , colnames(X) )
  
  S_unc <- cbind( runif(nsim * 30, w_init[1] - 0.5, w_init[1] + 0.5),
                  runif(nsim * 30, w_init[2] - 0.5, w_init[2] + 0.5),
                  runif(nsim * 30, w_init[3] - 0.5, w_init[3] + 0.5) )
  colnames(S_unc) <- colnames(X)
  # Turnover
  to <- apply(S_unc, 1, function(x) { sum(abs(x - w_init))})
  idx_to <- which(to < 0.36)
  # Tracking error
  te <- apply(S_unc, 1, function(x) { t(x - w_init) %*% covmat %*% (x - w_init) })
  idx_te <- which(te < quantile(te, 0.1))

  
  
  
  
  GPS_tmp <- GPS
  # addConstraint(GPS_tmp) <- budgetConstraint(rhs = 2, sense = "<=")
  # addConstraint(GPS_tmp) <- boxConstraint(name = "LongShort")
  addConstraint(GPS_tmp) <- turnoverConstraint(rhs = 0.2, 
                                               w_init = w_init, 
                                               linearize = TRUE)
  GPP <- gpp(GPS_tmp)
  GPP_lin <- linearize(GPP, method = "double")
  GPO_lin <- gpsolve(GPP_lin)
  ( w_lin <- getWeights(GPO_lin) )
  ( x0 <- GPO_lin@Result$x )
 
  cons <- getConstraints(GPP_lin)
  eps <- 1e-05
  cons@linear$rhs["turnover_budget"] <- cons@linear$rhs["turnover_budget"] + eps
  
  S_lin <- rpolytope(Constraints = cons,
                     nsim = nsim * 5, 
                     algo = "hitandrun",
                     x0 = x0)
  S_tob <- S_lin[ ,1:N]
  

  toFUN <- function(x, x0) { sum(abs(x - x0)) }
  to <- apply(S_lo, 1, toFUN, x0 = w_init)
  tob <- apply(S_tob, 1, toFUN, x0 = w_init)
  max(to); max(tob)
  
  
  
  
  
  # Colorcode
  txt <- colorRampPalette(colors = c("darkred", "orange"))(5)
  n <- nrow(S_lo) / length(txt)
  colors <- rep(NA, nrow(S_lo))
  for ( i in seq(along = txt)) {
    colors[((n*i)-n+1):(n*i)] <- rep(txt[i], n)
  }
  # colors <- fBasics::divPalette(n = nrow(S_lo), "Spectral")
  color <- rep(NA, nrow(S_lo))
  color[order(to)] <- colors
  
  
  
  # Chart
  plot3d(S_lo[ ,1], S_lo[ ,2], S_lo[ ,3], col = color,
         xlab = Names[1], ylab = Names[2], zlab = Names[3], size = 3)
  points3d(w_init[1], w_init[2], w_init[3], size = 10)
  # points3d(S_lin[ ,1], S_lin[ ,2], S_lin[ ,3])
  
  
  
  
  sds <- apply(S_lo, 1, sigmaFUN)
  q <- quantile(sds, seq(from = 0, to = 1, length.out = 10)) 
  col_txt <- rev(colorRampPalette(colors = c("lightblue", "darkblue"))(10))
  eps <- 1e-06
  lIdx <- list()
  for (i in seq(along = q)) {
    lIdx[[i]] <- which( abs(sds - q[i]) <= (eps * i) )
  }
  
  for (i in 1:length(lIdx)) {
    points3d(S_lo[lIdx[[i]], 1], S_lo[lIdx[[i]], 2], S_lo[lIdx[[i]], 3], col = col_txt[i], size = 4)
  }
  # points3d(S_lo[ ,1], S_lo[ ,2], S_lo[ ,3], col = color, size = 5)
  
  
  
  
  
  # Tracking Error
  tr_err <- apply(S_lo, 1, function(x) { t(x - w_init) %*% covmat %*% (x - w_init) })
  
  
  
  # Simulations
  sim <- timeSeries( apply(S_lo, 1, function(x) { X %*% x }), time(X) )
  colnames(sim) <- paste0("sim", 1:ncol(sim))

  
  
  # Value at Risk
  VaR <- ValueAtRisk(sim, method = "hist", p = 0.99 )[[2]][1, ]
  
  
  # Expected Shortfall
  ES <- ExpectedShortfall(Data = sim, method = "hist", p = 0.99 )[[2]][1, ]
  
  
  # CDF on sim for Omega ratio
  cdf_mat <- apply( sim, 2, function(x) { density(x, n = length(x))$y } )
  
  
  
  
  
  # Utility
  
  # Parameters
  nsim <- 10^4
  # a <- 0.88
  # b <- 2.25
  a <- 0.5
  b <- 1.9
  
  
  X_m <- aggMonthly(X, compounding = "continuous")
  # sim <- ( exp(X_m) - 1) %*% t(S_lo)
  sim <- X_m %*% t(S_lo)
  
  RP <- mean(X_m[ ,3]) * (-1)
  # RP <- -0.05 
  RP <- 0
  
  # U <- apply(sim, 2, ktUtility, 
  #            RP = RP, 
  #            a = a, 
  #            b = b)
  
  U <- apply(sim , 2, acUtility,
             b1 = lambda2b(0.5)$b1,
             b2 = lambda2b(0.5)$b2)
  
  avg_util <- apply(U, 2, mean)
  
  range(avg_util)
  idx_maxU <- which(avg_util == max(avg_util))
  S_lo[idx_maxU, ]
  avg_util[idx_maxU]
  
  
  
  
  
  # Plot utility function
  
  x <- as.numeric(sim[ ,1])
  y <- ktUtility(mu = x, RP = RP, a = a, b = b)
  
  plot(sort(x), sort(y), type = "l", lwd = 5)
  abline(h = 0, col = "darkgrey")
  abline(v = 0, col = "darkgrey")
  
  
  
  
  
  # Colorcoding
  # txt <- fBasics::divPalette(n = 50, name = "Spectral")
  txt <- colorRampPalette(colors = c("lightblue", "darkblue"))(5)
  n <- nrow(S_lo) / length(txt)
  colors <- rep(NA, nrow(S_lo))
  for ( i in seq(along = txt)) {
    colors[((n*i)-n+1):(n*i)] <- rep(txt[i], n)
  }
  # colors <- fBasics::divPalette(n = nrow(S_lo), "Spectral")
  color <- rep(NA, nrow(S_lo))
  # color[order(avg_util)] <- colors
  sds <- apply(S_lo, 1, sigmaFUN)
  color[order(sds)] <- colors
  
  
  
  
  # # Using solnp
  # 
  # # edit(gpp.nlp.solnp)
  # # edit(gpp.nlp.donlp2)
  # 
  # 
  # objFun <- function(x) 
  # {
  #   sim <- X_m %*% x
  #   U <- ktUtility(mu = sim, RP = RP, a = a, b = b)
  #   # U <- mean(sim)
  #   return( -mean(U) )
  # }
  # eqFun <- function(x) { sum(x) }
  # ineqFun <- sigmaFUN
  # 
  # 
  # lower <- getConstraints(cons_lo, "bounds")$lower
  # upper <- getConstraints(cons_lo, "bounds")$upper
  # idx_budget <- which(rownames(getConstraints(cons_lo, "linear")$Amat) == "budget")
  # budget = NULL
  # if ( length(idx_budget) > 0 ) {
  #   budget <- getConstraints(cons_lo, "linear")$rhs[idx_budget]
  # }
  # sigma_th <- ineqFun(x = w_gmv) * (1 + alpha)
  # 
  # 
  # x_init <- x0
  # # x_init <- S_lo[1, ]
  # opt <- solnp( pars = x_init, 
  #               fun = objFun, 
  #               eqfun = eqFun,
  #               eqB = budget, 
  #               ineqfun = ineqFun, 
  #               ineqLB = ineqFun(x = w_gmv), 
  #               ineqUB = sigma_th, 
  #               LB = lower, 
  #               UB = upper, 
  #               control = list() )
  # 
  # w_solnp <- opt$pars
  # opt$values
  
  
  
  
  
  
  
  
  
  
  E <- ellipsoid(q = rep(0, ncol(X)),
                 Q = solve(covmat))
  H <- hyperplane(v = rep(1, ncol(X)),
                  g = 1)
  I <- hpintersection(E = E, H = H)
  I_nm1 <- hpintersection.Nminus1(E = E, H = H)
  
  
  
  
  # Uniform samples from ellipoid
  
  sigma_th <- sigmaFUN(x = w_gmv) * (1 + alpha)
  Y_covmat <- runifE( S = covmat, 
                      z_hat = rep(0, ncol(covmat)),
                      gamma_threshold = sigma_th, 
                      n_points = nsim * 5, 
                      b_pushy = FALSE )
  
  
  
  # Uniform samples from n-1 ellipsoid (intersection of E and H)
  
  tmat <- I_nm1$tmat
  tmat_inv <- solve(tmat)
  f <- I_nm1$f
  w_prime <- I_nm1$w_prime
  W_prime <- I_nm1$W_prime
  
  
  # Find portfolio with variance sigma_th  
  # // 
  # // have to write procedure. Here we give a hardcoded exampple
  # // 
  w2 <- w_gmv * 0
  w2[1] <- 1
  sigmaFUN(w2)
  tmp <- twoAssetPortfolio(w1 = w_gmv, 
                           w2 = w2, 
                           covmat = covmat, 
                           alpha = alpha)
  w_th <- tmp$w
  sigmaFUN(w_th)
  
  v_th <- as.numeric( tmat %*% w_th - f )  
  th <- t(v_th - w_prime)[-1] %*% solve(W_prime[-1, -1]) %*% (v_th - w_prime)[-1]
  
  Y_prime <- runifE( S = solve(W_prime[-1, -1]), 
                     z_hat = w_prime[-1],
                     gamma_threshold = as.numeric( th ), 
                     n_points = nsim, 
                     b_pushy = FALSE )
  Y_prime <- cbind(0, Y_prime)
  colnames(Y_prime) <- Names
  Y <- t( apply(Y_prime, 1, function(x) { tmat_inv %*% (x + f) } ) )
  colnames(Y) <- Names
  
  Y_prime2 <- runifE( S = solve(W_prime[-1, -1]), 
                      z_hat = w_prime[-1],
                      gamma_threshold = as.numeric( th ), 
                      n_points = nsim * 1, # ~~~~~~~~~~~~~ 
                      b_pushy = TRUE )    #!
  Y_prime2 <- cbind(0, Y_prime2)
  colnames(Y_prime2) <- Names
  Y2 <- t( apply(Y_prime2, 1, function(x) { tmat_inv %*% (x + f) } ) )
  colnames(Y2) <- Names
  
  
  # U <- apply(Y_prime2, 1, objFun)
  # idx <- which(U == min(U))
  # ( w_rp <- Y2[idx, ] )
  
  # objFun(x = w_rp)
  # objFun(x = opt$pars)
  
  # w_rp
  # opt$pars
  
  
  
  
  eps <- 1e-02
  plot3d(S_ls[ ,1], S_ls[ ,2], S_ls[ ,3], col = "grey",
         xlab = Names[1], ylab = Names[2], zlab = Names[3], pch = 19, cex = 2, size = 1)
  points3d(S_lo[ ,1], S_lo[ ,2], S_lo[ ,3], col = color, size = 3)
  # points3d(w_solnp[1], w_solnp[2], w_solnp[3], col = 2, size = 10, lwd = 3)
  # points3d(w_rp[1], w_rp[2], w_rp[3], col = 3, size = 10, lwd = 3)
  # points3d(w_gmv[1], w_gmv[2], w_gmv[3], col = "orange", size = 10)
  points3d(Y_covmat[ ,1], Y_covmat[ ,2], Y_covmat[ ,3], col = "orange", size = 3)
  points3d(Y[ ,1] + eps, Y[ ,2] + eps, Y[ ,3] + eps, col = "darkorange", size = 3)
  # points3d(Y2[ ,1] + eps, Y2[ ,2] + eps, Y2[ ,3] + eps, col = "darkorange", size = 3)
  
  
  
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # EFFICIENT FRONTIER PLOT
  # --------------------------------------------------------------------------
  
  frqncy <- xts::periodicity(X)
  scalefactor <- switch(frqncy$scale,
                        daily   = 252,
                        weekly  = 52,
                        monthly = 12, 
                        yearly  = 1)
  
  GPS_lo <- gps(Data = X, 
                Solver = solverCtrl(utility = list(n_frontierpoints = 200,
                                                   scalefactor = scalefactor)))
  GPO_lo <- gpo(GPS = GPS_lo)
  w_lo <- getWeights(GPO_lo)
  GPS_ls <- GPS_lo
  addConstraint(GPS_ls) <- boxConstraint(name = "LongShort")
  fp_lo <- frontierPoints(GPS = GPS_lo)
  fp_ls <- frontierPoints(GPS = GPS_ls)
  wmat_lo <- attr(fp_lo, "Weights")
  wmat_ls <- attr(fp_ls, "Weights")
  
  
  
  
  
  
  
  
  
  
  # --------------------------------------------------------------------------
  # SHRINKING THE COVARIANCE MATRIX
  # --------------------------------------------------------------------------
  
  
  lambda_vec <- seq(from = 0, to = 1, length.out = 10)
  
  
  lCovmat <- list()
  for ( i in seq(along = lambda_vec) ) {
    
    covmat_shrink <- covariance(Data = X, 
                                covCtrl(method = "shrinkage",
                                        ellipsis = list(target = "dcv", 
                                                        lambda = lambda_vec[i])))
    lCovmat[[i]] <- covmat_shrink
  }
  names(lCovmat) <- paste0("lambda=", lambda_vec)
  
  
  plot(lCovmat[[1]])
  plot(lCovmat[[2]])
  plot(lCovmat[[3]])
  plot(lCovmat[[4]])
  plot(lCovmat[[5]])
  plot(lCovmat[[10]])
  
  
  alpha <- 0.05
  H <- hyperplane(v = rep(1, ncol(X)),
                  g = 1)
  lY <- lYI <- list()
  for ( i in seq(along = lambda_vec) ) {
    
    # Analytic solution to long-short minimum variance
    ones <- rep(1, ncol(lCovmat[[i]]))
    covmat_inv <- solve(lCovmat[[i]])
    w_ls <- covmat_inv %*% ones / as.numeric( t(ones) %*% covmat_inv %*% ones)
    sigma_ls <- as.numeric( t(w_ls) %*% lCovmat[[i]] %*% w_ls )
    sigma_th <- sigma_ls * (1 + alpha)
    
    lY[[i]] <- runifE( S = lCovmat[[i]], 
                       z_hat = rep(0, ncol(covmat)),
                       gamma_threshold = sigma_th, 
                       n_points = nsim * 5, 
                       b_pushy = FALSE )
    
    E <- ellipsoid(q = rep(0, ncol(X)),
                   Q = covmat_inv)
    I <- hpintersection(E = E, H = H)
    
    lYI[[i]] <- runifE(S = solve(make.positive.definite(getShape(I))), 
                       z_hat = getCentre(I),
                       gamma_threshold = sigma_th,
                       n_points = nsim,
                       b_pushy = FALSE)
    
    
  }  
  
  
  
  # Efficient frontier plot in asset space
  plot3d(S_ls[ ,1], S_ls[ ,2], S_ls[ ,3], col = "grey",
         xlab = Names[1], ylab = Names[2], zlab = Names[3], pch = 19, cex = 2, size = 1)
  points3d(S_lo[ ,1], S_lo[ ,2], S_lo[ ,3], col = color, size = 3)
  points3d(Y_covmat[ ,1], Y_covmat[ ,2], Y_covmat[ ,3], col = "orange", size = 3)
  points3d(lY[[5]][ ,1], lY[[5]][ ,2], lY[[5]][ ,3], col = "darkorange", size = 3)
  points3d(lY[[10]][ ,1], lY[[10]][ ,2], lY[[10]][ ,3], col = "red", size = 3)
  points3d(lYI[[1]][ ,1], lYI[[1]][ ,2], lYI[[1]][ ,3], col = "green", size = 3)
  points3d(lYI[[5]][ ,1], lYI[[5]][ ,2], lYI[[5]][ ,3], col = "green", size = 3)
  points3d(lYI[[10]][ ,1], lYI[[10]][ ,2], lYI[[10]][ ,3], col = "darkgreen", size = 3)
  


  
  
  
  # --------------------------------------------------------------------------
  # SAVE
  # --------------------------------------------------------------------------
  
  env <- new.env()
  
  # Data 
  env$X <- X
  env$sim <- sim
  
  # Specs
  env$cons_lo <- cons_lo
  
  # Random portfolios
  env$S_ls <- S_ls
  env$S_lo <- S_lo
  env$S_tob <- S_tob
  env$S_unc <- S_unc
  env$idx_to <- idx_to
  env$idx_te <- idx_te
  
  # Statistics
  env$lStats <- lStats
  env$sds <- sds
  env$to <- to
  env$tr_err <- tr_err
  env$VaR <- VaR
  env$ES <- ES
  env$w_init <- w_init
  env$nsim <- nsim
  
  # Efficient Frontier
  env$fp_lo <- fp_lo
  env$fp_ls <- fp_ls
  env$wmat_lo <- wmat_lo
  env$wmat_ls <- wmat_ls
  
  # Covariances
  env$lCovmat <- lCovmat
  
  # Ellipsoids
  env$lY <- lY
  
  
  if ( ncol(X) == 3 ) {
    saveRDS(object = env, file = "R:/Asset_Management/R/Shiny/Portfolio_Geometry/data_3d.rds")
  } else if ( ncol(X) == 4 ) {
    saveRDS(object = env, file = "R:/Asset_Management/R/Shiny/Portfolio_Geometry/data_4d.rds")
  }  
  # --------------------------------------------------------------------------
  
  
  
  
  
  
  
