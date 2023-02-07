  
  
  ############################################################################
  ### GOOD AND BAD TURBULENCE - GAS
  ############################################################################
  
  
  # --------------------------------------------------------------------------
  # Cyril Bachelard
  # This version:     30.06.2021
  # First version:    30.06.2021
  # --------------------------------------------------------------------------
  
  
  require(betategarch)
  require(GAS)
  require(DAARC)
  
  
  
  
  # --------------------------------------------------------------------------
  # Data
  # --------------------------------------------------------------------------
  
  X <- getMSCIData( universe = "dm", 
                    frqncy = "d", 
                    ccy = "Local" )
  X <- log(1 + X[isWeekday(time(X)), ])
  X <- X[ ,1:5]
  
  
  
  # --------------------------------------------------------------------------
  # GAS package
  # --------------------------------------------------------------------------
  
  
  DistLabels()
  DistInfo( DistLabel = "norm" )
  
  gas_spec <- UniGASSpec( Dist = "norm",
                          # Dist = "std", 
                          ScalingType = "Identity",
                          # ScalingType = "Inv",
                          GASPar = list(location = TRUE, 
                                        scale = TRUE,
                                        skewness = FALSE, 
                                        shape = FALSE, 
                                        shape2 = FALSE) )
  cvar <- X * NA
  lFit <- list()
  for ( j in 1:ncol(X) ) {
    fit <- UniGASFit( GASSpec = gas_spec, 
                          data = as.numeric(X[ ,j]), 
                          fn.optimizer = fn.optim,
                          # fn.optimizer = fn.solnp,
                          Compute.SE = FALSE )
    lFit[[j]] <- fit
    cvar[ ,j] <- getMoments(fit)[1:nrow(X), "M2"]
    
  }
  
  # for ( j in 1:length(lFit) ) {
  #   cvar[ ,j] <- getFilteredParameters(lFit[[j]])[1:nrow(X), "scale"]
  # }
  
  
  plot( cvar )
  
  
  
  # --------------------------------------------------------------------------
  # betategarch package
  # --------------------------------------------------------------------------
  
  fit <- tegarch( y = X[ ,1], 
                  asym = TRUE, 
                  skew = TRUE, 
                  components = 1, 
                  initial.values = NULL,
                  lower = NULL, 
                  upper = NULL,
                  hessian = TRUE, 
                  lambda.initial = NULL,
                  c.code = TRUE, 
                  logl.penalty = NULL, 
                  aux = NULL )
  
  fit
  
  # Graph of fitted volatility (conditional standard deviation):
  plot( fitted(fit) )
  
  
  
  cvar <- X * NA
  lFit <- list()
  for ( j in 1:ncol(X) ) {
    fit <- tegarch( y = X[ ,j], 
                    asym = TRUE, 
                    skew = TRUE, 
                    components = 1, 
                    initial.values = NULL,
                    lower = NULL, 
                    upper = NULL,
                    hessian = TRUE, 
                    lambda.initial = NULL,
                    c.code = TRUE, 
                    logl.penalty = NULL, 
                    aux = NULL )
    lFit[[j]] <- fit
    cvar[ ,j] <- fitted(fit)
    
  }
  
  plot( cvar )
  
  
  
  