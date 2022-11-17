library(nimble)
# function to compute the analysis distribution from the forecasted and observed states 
GEF <- function(Forecast, Observed, H, verbose=1) {
  library(nimble)
  
  GEF.nimble <-  nimbleCode({
    
    # Sorting out qs
    qq ~ dunif(0.0001, 5)
    
    # Do we have more than one forecast variable? 
    if (N > 1){
      
      # Use Q to inflate the diagonals
      pf_theta[1:N, 1:N] <-  pf[1:N, 1:N] + ((qq-1) * diag(diag(pf[1:N, 1:N])))
      
      # X model
      X.mod[1:N] ~ dmnorm(muf[1:N], cov = pf_theta[1:N, 1:N])
      
    }else{
      
      # Use Q to inflate the diagonals
      pf_theta <-  pf * qq
      
      # X model
      X.mod ~ dnorm(muf, var = pf_theta)
    }
    
    # Now apply H 
    if (YN > 1){
      
      # If we have more than one observation, but only one forecast var, we repeat it 
      if (N < 2){
        tmpX  <- X.mod
        for (i in 1:YN) {
          Xs[i] <- tmpX
        }
      }else{
        
        # Otherwise, we just map it to the correct forecast using H
        for (i in 1:YN) {
          tmpX[i]  <- X.mod[H[i]]
          Xs[i] <- tmpX[i]
        }
      }
    }else{
      
      # If we have one observation and one forecast var, then it's straightforward
      if (N == 1){
        tmpX  <- X.mod
        Xs <- tmpX
      }else{
        
        # Otherwise, we have PDA and will handle those variables later, but we still need to map since X.mod is a vector
        tmpX  <- X.mod[H]
        Xs <- tmpX
      }
    }
    
    # Now consider data model 
    # Do we have more than one observation?
    if (YN > 1){
      
      # Likelihood
      y.obs[1:YN] ~ dmnorm(Xs[1:YN], prec = r[1:YN, 1:YN])
      
    }else{
      
      # Likelihood
      y.obs ~ dnorm(Xs, var = r)
      
    }
    
    # Putting the forecast variables without obs in Xall 
    # Otherwise we just use Xs
    if (nNotH > 0){
      
      if (nNotH > 1){
        for (i in 1:nNotH) {
          tmpXmod[i]  <- X.mod[NotH[i]]
          Xall[MapXall[i]] <- tmpXmod[i]
        }
      }else{
        tmpXmod <- X.mod[NotH]
        Xall[MapXall] <- tmpXmod
      }
      
      # Add back values from Xs
      if (YN > 1){
        for (i in 1:YN) {
          tmpXH[i]  <- Xs[i]
          Xall[i] <- tmpXH[i]
        }
      }else{
        tmpXH <- Xs
        Xall[1] <- tmpXH
      }
    }
    
  })
  # Fix H so it just maps the Obs to Forecast variables (will be a vector)
  H = apply(H, 1, function(row){
    which(row == 1)
  })
  if (is.list(H)) H <- H[[1]]
  
  # Forecast inputs
  Q <- Forecast$Q # process error
  X <- Forecast$X # states from all ensembles 
  R <- NULL
  # Observed inputs
  if (is.null(R)) R <- Observed$R
  Y <- Observed$Y
  # forecast prior info 
  mu.f <- as.numeric(apply(X, 2, mean, na.rm = TRUE))
  if (ncol(X)==1){Pf <- cov(X)[1,1]} else{Pf <- cov(X)}
  if (length(Y) ==1){r <- R[1,1]} else{r <- solve(R)}

  if(is.null(Q)) Q <-1
  
  # Figure out H mapping 
  # Do we have unobserved variables in Forecast?
  NotH <- which(!(1:ncol(X) %in% H))
  if (length(NotH)==0){
    Xall = mu.f[H]
    Xs = Xall
    NotH = 0
    nNotH = 0
    MapXall = 0
    DEmap = c(H)
  }else{
    Xall = c(mu.f[H], mu.f[NotH])
    Xs = mu.f[H]
    nNotH = length(NotH)
    MapXall = (length(H)+1):length(Xall)
    DEmap = c(H, NotH)
  }
  
  ## First step will create the C++ NIMBLE libraries that is needed and takes time
  ## on the else of this if which will be t>1, only the variables inside the model will
  ## be updated and the model will be run using the same library

  # In the first step we just build the model and in other steps we look to see if it exists and 
  # if yes we only update the variables in the model 
  if(!exists("Cmcmc")) {
    
    # Initial values
    inits.pred <-
      list(
        X.mod = as.vector(mu.f),
        qq = Q,
        Xall = Xall,
        Xs = Xs
      ) 
    
    #dimensions.tobit <- list(X.mod = ncol(X))
    
    # Constants defined in the model
    constants.tobit <-list( N = ncol(X),# number of state variables
                            YN = length(Y), # number of obs
                            H = H,
                            NotH = NotH,
                            nNotH = nNotH, 
                            MapXall = MapXall
    )
    
    # Data used for setting the likelihood and other stuff
    data.tobit <- list( muf = as.vector(mu.f),
                        pf = Pf,
                        y.obs = Y,
                        r = r
    )
    
    #------ This list could be appended by some info from obs.list for this day to accommodate other nimble models
    # Both the nimble model (here called GEF.nimble) and data.tobit object can be updated/appended to accommodate more complex models. 
    #  To make this step more generalization, we know that everything can go into the SDA (different data/model types) but 
    # It always spits out one important thing and that is Xa and Pa. These four lists and the nimble model could be 
    # appended or be brought in bythe obs.list at different time steps to model all sorts of complex bayesian models and data types while not 
    # hard coding any of that here. 

    #-------------------------------------------------------- NIMBLE MODEL ----------------------------
    # This is the first step in making the nimble model - Nimble does some preliminary checks on the code    
    model_pred <- nimbleModel(GEF.nimble,
                              data = data.tobit,
                              #dimensions = dimensions.tobit,
                              constants = constants.tobit,
                              inits = inits.pred,
                              name = 'base')
    
    
    model_pred$initializeInfo()
    
    ## Adding X.mod,q,r as data for building model.
    conf <- configureMCMC(model_pred, print=TRUE)
    
    # This will need to be changed in the case that we have varying Hs in the obs.list. 
    # For example, if we have 10 cm value, but no 20 cm value. 
    if(constants.tobit$nNotH > 0){
      conf$addMonitors(c("Xs","qq","Xall"))
    }else{
      conf$addMonitors(c("Xs","qq"))
    }
    
    # We will put this in the global env
    Rmcmc <<- buildMCMC(conf)
    Cmodel <<- compileNimble(model_pred, showCompilerOutput=TRUE)
    Cmcmc <<- compileNimble(Rmcmc, project = model_pred, showCompilerOutput=TRUE)
    
  }else {
    Cmodel$y.obs <- Y
    Cmodel$qq <- Q # update from the prior q 
    Cmodel$muf <- mu.f
    Cmodel$pf <- Pf
    Cmodel$r <- solve(R) #precision
  }
  
  # Run the MCMC
  DEsamples <- runMCMC(Cmcmc, niter = 100000, nburnin = 30000, nchains = 1, samplesAsCodaMCMC = TRUE)

  
  #---------------- Write down some GEF_Diagnostics
  if(verbose==1){
    if(!dir.exists(file.path(getwd(),"GEF_Diagnostics"))){
      dir.create(file.path(getwd(),"GEF_Diagnostics"))
    }
    save(DEsamples, file=file.path(getwd(),"GEF_Diagnostics",paste0("Chains_", Observed$Date, ".RData")))
  }

  ## -------------------- update parameters and estimate mua and pa
  if (nNotH == 0){
    iX <- grep("Xs", colnames(DEsamples), fixed = TRUE)
  }else{
    iX <- grep("Xall", colnames(DEsamples), fixed = TRUE)
  }
  vals = iX[sapply(1:ncol(X), function(ind){which(DEmap==ind)[1]})]
  
  if (length(vals) > 1){
    mu.a <- colMeans(DEsamples[,vals])
    Pa   <- cov(DEsamples[,vals])
    Pa[is.na(Pa)] <- 0
  }else{
    mu.a <- mean(DEsamples[,vals])
    Pa   <- var(DEsamples[,vals])
    Pa[is.na(Pa)] <- 0
  }
  
  Q <- mean(DEsamples[, "qq"])
  if (ncol(X) > 1){
    Pf_theta <- Pf + (Q-1) * diag(diag(Pf))
  }else{
    Pf_theta <- Pf * Q
  }
  
  return(list(
    mu.f = mu.f,
    Pf = Pf,
    mu.a = as.numeric(mu.a),
    Pa = Pa,
    Pf_theta = Pf_theta,
    H=H,
    Q=Q
  ))
}

# function to compute the analysis distribution from the forecasted and observed states 
EnKF <- function(Forecast, Observed, R=NULL, H, theta = 1, extraArg = NULL, ...) {
  
  # Forecast inputs
  Q <- NULL
  X <- Forecast$X # states from all ensembles 
  
  # Observed inputs
  if (is.null(R)) R <- Observed$R
  Y <- Observed$Y
  
  # EnKF ---------------------------------------------
  
  # forecast prior info 
  mu.f <- as.numeric(apply(X, 2, mean, na.rm = TRUE))
  Pf <- cov(X)
  
  # MK: We might want to rethink this hack in the context of SDA with soil moisture
  # 0.1 is a pretty high variance value for the forecast
  # diag(Pf)[which(diag(Pf) == 0)] <- 0.1 # hack for zero variance
  
  # remove covariance values from Pf is ensemble size is too small to allow for good estimation 
  if (nrow(X) < 30) Pf <- diag(diag(Pf))
  
  # process error
  if (!is.null(Q)) {
    Pf <- Pf + Q
  }
  
  # inflate forecast variance (just the diagonals) 
  Pf_theta <- Pf
  diag(Pf_theta) <- diag(Pf) * diag(theta)
  
  # Kalman gain calculation with inflation factor
  K <- Pf_theta %*% t(H) %*% solve((R + H %*% Pf_theta %*% t(H)))
  
  # compute analysis distribution
  mu.a <- mu.f + K %*% (Y - H %*% mu.f)
  Pa   <- (diag(ncol(X)) - K %*% H) %*% Pf_theta
  
  return(list(
    mu.f = mu.f,
    Pf = Pf,
    mu.a = as.numeric(mu.a),
    Pa = Pa,
    Pf_theta = Pf_theta,
    H=H,
    Q=NULL
  ))
}

# adjust analysis using original likelihood of forecast (not inflated) 
adj.ens<-function(Pf, X, mu.f, mu.a, Pa){
  
  vars = colnames(X)
  X = as.matrix(X)
  
  S_f  <- svd(Pf)
  L_f  <- S_f$d
  V_f  <- S_f$v
  
  # normalize
  Z <- X*0
  
  for(i in seq_len(nrow(X))){
    Z[i,] <- 1/sqrt(L_f) * t(V_f)%*%(X[i,]-mu.f)
  }
  Z[is.na(Z)]<-0
  Z[is.infinite(Z)] <- 0
  
  # analysis
  S_a  <- svd(Pa)
  
  L_a  <- S_a$d
  V_a  <- S_a$v
  
  # analysis ensemble 
  X_a <- X*0
  for(i in seq_len(nrow(X))){
    # she decomposed Pa - then it's putting it back together but with a different Z which comes from the likelihood of that ens    
    X_a[i,] <- V_a %*% diag(sqrt(L_a), ncol = ncol(X)) %*%Z[i,] + mu.a
  }
  
  #if(sum(mu.a - colMeans(X_a)) > 1 | sum(mu.a - colMeans(X_a)) < -1) logger.warn('Problem with ensemble adjustment (1)')
  #if(sum(diag(Pa) - diag(cov(X_a))) > 5 | sum(diag(Pa) - diag(cov(X_a))) < -5) logger.warn('Problem with ensemble adjustment (2)')
  
  analysis <- as.matrix(X_a)
  
  # The following is a hack to ensure that soil moisture does not go below zero
  # also need to remember that soil moisture cannot go above 1, though this doesn't typically happen 
  SWS = grep('sw',vars)
  if (length(SWS) > 0){
    if (length(SWS) > 1){
      analysis[,SWS] <- apply(analysis[,SWS], 2,
                              function(x) {
                                #x[x <= 0] <- median(x)
                                x[x <= 0] <- 0.1 # MK: the air-dry values for layers 2, 3, and 4 are 0.089, 0.087, 0.087 
                                x[x >= 0.8] <- 0.44 # MK: the saturated values for layers 2, 3, and 4 are 0.437, 0.435, 0.435
                                return(x)
                              })
    }else{
      x <- analysis[,SWS]
      x[x <= 0] <- 0.1 # MK: the air-dry values for layers 2, 3, and 4 are 0.089, 0.087, 0.087 
      x[x >= 0.8] <- 0.44 # MK: the saturated values for layers 2, 3, and 4 are 0.437, 0.435, 0.435
      analysis[,SWS] <-  x
    }
  }
  
  # Ks also cannot go below zero
  Ks = grep('ks',vars)
  if (length(Ks) > 0){
    analysis[,Ks] <- apply(analysis[,Ks], 2,
                           function(x) {
                             x[x <= 0] <- 0.1 # MK: from the "reasonable" range of values for this parameter in real soils
                             return(x)
                           })
  }
  
  return(analysis)
}

# function to apply method by Miyoshi et. al (2013) for inflation factor and R matrix 
inflate <- function(Observed, X, Xa, H, R, theta, p = 0.5){
  
  # determine innovations for analysis and forecast
  obs.mat = matrix(rep(Observed$Y, nrow(X)), nrow = nrow(X), byrow = TRUE)
  
  # check to see if we are adjusting variables outside of obs data 
  # i.e. if we have more forecast variables than observed 
  flag = FALSE
  if (ncol(X) > ncol(obs.mat)){
   
    # which forecast are in obs?
    inds = which(apply(H, 2, function(col){any(col == 1)}))
    
    fn = ncol(X)
    
    # adjust all data to the correct dimensionos
    X <- X[,inds]
    Xa <- Xa[,inds]
    H <- H[,inds]
    theta <- theta[inds,inds]
    flag = TRUE
  }
  
  dof = obs.mat - X 
  daf = obs.mat - Xa
  
  # we can run into matrix notation issues if the X has only one variable
  if (is.null(dim(X))) dim(X) <- c(length(X), 1)
  
  # update R using E(daf x dof) = cov(daf, dof) + mean(daf)*mean(dof)
  Ru = cov(daf, dof) + (apply(daf, 2, mean) %*% t(apply(dof, 2, mean)))
  
  # then make sure R has no covariance values
  diags = diag(Ru)
  Ru = Ru * 0 
  diag(Ru) = diags
  
  # estimate updated inflation factor
  exp = cov(dof) + (apply(dof, 2, mean) %*% t(apply(dof, 2, mean)))
  denom = solve(H %*% (cov(X) * theta) %*% t(H))
  theta_update = (exp - Ru) * denom
  
  # smooth inflation factor using past estimated and newly estimated
  theta = ((1-p) * theta) + (p * theta_update)
  Rnew = ((1-p) * R) + (p * Ru)
  
  # make sure the R is not negative
  if(Rnew > 1e-6) {
    R <- Rnew
  }
  
  if (flag){
    theta_bigger = diag(fn)
    theta_bigger[inds, inds] = theta
    theta = theta_bigger
  }
  
  # set minimum R variance values to be around accuracy ranges for sensors (at the very least)
  # vars = colnames(R)
  # diags = diag(R)
  # SMs = grep('sw',vars)
  # if (length(SMs) > 0) {
  #   toChange = which(diags[SMs] < 0.01^2)
  #   diags[SMs[toChange]] = 0.01^2
  # }
  # diag(R) = diags
  # TO DO: we probably should add a lower limit for other variables in R 
  
  # make sure theta isn't shrinking the forecast uncertainty and remove all inflation of covariance values
  tempDiag = diag(theta)
  tempDiag[which(tempDiag < 10)] = 10
  theta = diag(length(tempDiag))
  diag(theta) = tempDiag
  
  return(list(theta = theta, R = R, Rnew= Rnew, Ru= Ru))
}

# data assimilation workflow with switches
data_assimilation <- function(Forecast, Observed, R, theta, miyoshi,DA.method, p = 0.5, Rcov){

  
  if (!is.null(Forecast$X)){
    if (nrow(Forecast$X) > 0){
      
      # form H based on available data and forecast 
      H <- Observed$H #%>% as.matrix()
      
      Forecast$X <- Forecast$X %>% select(Observed$Forecast.name)

      if(is.numeric(R)) {
        # should R covariance values be included? if no, remove them 
        if (Rcov == 0){
          Rnames = colnames(Observed$R)
          Observed$R = diag(diag(Observed$R), ncol = length(Observed$Y))
          colnames(Observed$R) = Rnames
          rownames(Observed$R) = Rnames
        }
      }

      # Based on the method, call the appropriate SDA function
      # Method=1 EnKF (Miyoshi=1 includes Miyoshi algorithm)
      # Method=2 Generalized Ensemble Filter (GEF)
      
      if(DA.method==1) {
         if (is.null(theta)) theta <- matrix(10, ncol(Forecast$X), ncol(Forecast$X))
        # first data assimilation without adjustment or later with Miyoshi algorithm
        if (miyoshi == 0){
          
          DAR <- EnKF(Forecast, Observed, H = H, theta = theta)
          DAResult <- adj.ens(DAR$Pf, Forecast$X, DAR$mu.f, DAR$mu.a, DAR$Pa)
          
        }else{
          
          # in the case that this is the first day of assimilation, R and theta will be NULL
          if (is.null(R)) R <- Observed$R
         
          
          DAR <- EnKF(Forecast, Observed, R = R, H = H, theta = theta)
          
          DAResult <- adj.ens(DAR$Pf, Forecast$X, DAR$mu.f, DAR$mu.a, DAR$Pa)
          
          infl <- inflate(Observed = Observed, X = as.matrix(Forecast$X),
                          Xa = DAResult, H = H, R = R, theta = theta, p = p)
          
          Observed$R <- R 
          DAR$Theta <- theta
          DAR$Inflate <- infl
          # set up values for next time step
          R <- infl$R # This is reestimation of the R
          theta <- infl$theta # this inflates forecast
        }
        

        
      }else if (DA.method==2) {
        
        # in the case that this is the first day of assimilation, R and theta will be NULL
        #if (is.null(R)) R <- Observed$R
        #if (is.null(theta)) theta <- matrix(1,ncol(Forecast$X),ncol(Forecast$X))
        
        # Turned off use of Miyoshi inflation script for use with GEF
        #Observed$R <- R

        if(Observed$H == "RTM") {
          # Source the model
          source('RTM_GEF.R')
          #DAR <- callr::r(GEF.RTM, args=list(Forecast, Observed, H = H))
          DAR <- GEF.RTM (Forecast, Observed, H = H)
        }else {
          #DAR <- callr::r(GEF, args=list(Forecast, Observed, H = H))
          DAR <- GEF(Forecast, Observed, H = H)
        }

        
        colnames(DAR$Pa) <- colnames(DAR$Pf)
        rownames(DAR$Pa) <- rownames(DAR$Pf)
        
        DAResult <- adj.ens(DAR$Pf, Forecast$X, DAR$mu.f, DAR$mu.a, DAR$Pa)
        
        colnames(DAR$Pa)
        
        #infl <- inflate(Observed = Observed, X = as.matrix(Forecast$X),
        #                Xa = DAResult, H = H, R = R, theta = theta, p = p)
        
        #DAR$Theta_Myioshi <- theta
        
        # set up values for next time step
        #R <- infl$R
        #theta <- infl$theta
      }else if (DA.method==4) {
        library(lubridate)
        
        # Create dir for writing the 4DEnVar files
        if(!dir.exists("4DEnVar_Files")) dir.create("4DEnVar_Files")
        if(dir.exists(file.path("4DEnVar_Files", gsub("-","_", today)))) unlink(file.path("4DEnVar_Files", gsub("-","_", today)), recursive = TRUE)
        if(!dir.exists(file.path("4DEnVar_Files", gsub("-","_", today)))) dir.create(file.path("4DEnVar_Files", gsub("-","_", today)))
        
        #save.image(file.path("4DEnVar_Files", gsub("-","_", today), "Test4DVar.RData"))
        #-----------------------------------
        today.vec <- (strsplit(today, '-')[[1]])
        
        today_date_frmt <- Date_APSIM_to_R(today.vec[2] %>% as.numeric(),
                        today.vec[1] %>% as.numeric()) %>% as.Date()
        
        #------------ Finding obs for this time window
        window_dates <- lubridate::interval(lubridate::ymd(today_date_frmt-window_length),
                                            lubridate::ymd(today_date_frmt))
        
        
        window_obs<- obs.list %>%
          purrr::keep(~lubridate::ymd(.x$Date) %within% window_dates)
        
        if(length(window_obs) > window_length){
          diff <- length(window_obs) - window_length
          window_obs <- window_obs[(1+(diff)):length(window_obs)]
        }
        
        Ys <- names(window_obs) %>%
          purrr::map_dfr(~ data.frame(Y = window_obs[[.x]]$Y,
                                  R = window_obs[[.x]]$R %>% as.numeric(), 
                                  Date= .x
                                  )
                     )

        obs_days <- names(window_obs) %>%
          purrr::map_dbl(~ strsplit(.x,"-")[[1]][1] %>% as.numeric)
        
        obs_years <- names(window_obs) %>%
          purrr::map_dbl(~ strsplit(.x,"-")[[1]][2] %>% as.numeric)
        #------------ Finding Predict for this time window
        X <- data.all %>%
          dplyr::filter(day %in% obs_days & year %in% obs_years) %>%
          dplyr::select(day, year, window_obs[[1]]$Forecast.name, ensemble, sat1, dul1, ll151, bd1, swcon1) %>%
          mutate(Date=paste0(day,"-", year)) %>%
          left_join(Ys, by="Date")
        
        hx <- X %>%
          split(.$ensemble) %>%
          purrr::map_dfr(~ .x[["sw1"]]) %>%
          as.matrix() %>%
          `colnames<-`(NULL)
        nrow_states <- nrow(hx)
        #------ Add soil params to Xb
        #Background params
        Pb <- X %>%
          split(.$ensemble) %>%
          purrr::map_dfr(~ .x[['sat1']]) %>%
          apply(2, mean)%>%
          bind_rows(
            X %>%
              split(.$ensemble) %>%
              purrr::map_dfr(~ .x[['dul1']]) %>%
              apply(2, mean) 
          ) %>%
          bind_rows(
            X %>%
              split(.$ensemble) %>%
              purrr::map_dfr(~ .x[['bd1']]) %>%
              apply(2, mean)
          )%>%
          bind_rows(
            X %>%
              split(.$ensemble) %>%
              purrr::map_dfr(~ .x[['ll151']]) %>%
              apply(2, mean) 
          ) %>%
          as.matrix() %>%
          `colnames<-`(NULL)
        
        
        Xb <- rbind(hx, Pb)
        #--------------------------------------------------------------------------------------------------

        #xb     --- the background ensemble of initial state and/or parameters. Each ensmble is a column. rows are in time. The analysis is the posterior of this.
        write.table(Xb, file=file.path("4DEnVar_Files", gsub("-","_", today), "0xb.dat"), col.names = FALSE, row.names = FALSE)
        
        #hx the ensmble of model predicted observations. Each column is a ens of model simulations through the time window . Each row is a time step.
        write.table(hx, file=file.path("4DEnVar_Files", gsub("-","_", today), "0hx.dat"), col.names = FALSE, row.names = FALSE)
        
        #hx_bar --- the model predicted observations for the mean of xb (n_obs rows)  
        write.table(apply(hx, 1, mean), file=file.path("4DEnVar_Files", gsub("-","_", today), "0hxbar.dat"), col.names = FALSE, row.names = FALSE)
        
        #R      --- the observation uncertainty covariance matrix (n_obs rows; n_obs cols)  
        write.table(diag(Ys$R), file=file.path("4DEnVar_Files", gsub("-","_", today), "0R.dat"), col.names = FALSE, row.names = FALSE)
        
        #y      --- the observations (n_obs rows)  
        write.table(Ys$Y, file=file.path("4DEnVar_Files", gsub("-","_", today), "0y.dat"), col.names = FALSE, row.names = FALSE)
        #----------------------------------------- Run 4DEnVar
        #"/pysims/data/4DEnVar_engine-main/4DEnVar"
        output_4DVar <- system(paste0("/pysims/data/4DEnVar_engine-main/4DEnVar"," ",
                                      file.path("4DEnVar_Files", gsub("-","_", today), "0xb.dat")," ",
                                      file.path("4DEnVar_Files", gsub("-","_", today), "0hx.dat")," ",
                                      file.path("4DEnVar_Files", gsub("-","_", today), "0y.dat")," ",
                                      file.path("4DEnVar_Files", gsub("-","_", today), "0R.dat")," ",
                                      file.path("4DEnVar_Files", gsub("-","_", today), "0hxbar.dat")
                                      ), 
                     intern = TRUE)
        
        
        if(output_4DVar > 0) {
          # I'm extracting the last time step
          # 4 parameters + 1 dash line 
          # Towards the end of the growing season sometimes there is not enough obs for the whole time window. 
          # Use whatever that is available to extract the same states/params
          sel_ind <- window_length
          if(nrow_states < window_length) sel_ind <- nrow_states
          
          
          print("RAW 4DEnVar-------------")
          print(output_4DVar)
          Xa <- output_4DVar[((sel_ind+4)+1+sel_ind):length(output_4DVar) ] %>% 
            purrr::map_dfc(~.x %>%strsplit(" ") %>% unlist() %>% as.numeric) %>%
            as.matrix() %>%
            `colnames<-`(c("sw1","sat1","dul1","bd1","ll151"))
          
          mu.a <- mean(Xa)
          Pa <- var(Xa)
          DAResult <- Xa
          DAR <- list(
            mu.a = as.numeric(mu.a),
            Pa = Pa,
            Xa = Xa, 
            Xb = Xb, 
            hx = hx, 
            output_4DVar = output_4DVar, 
            X = X
          )
          
          # print("Xa-------------")
          # print(Xa)
        }
        
      }
      
      # gather info on SDA for this time step
      DAR$Adj <- DAResult
      DAR$Obs <- Observed
      DAR$Forecast <- Forecast
      
      return(list(Xa = DAResult, DAR = DAR, R = R, theta = theta))

    }
  }else{
    return(NULL)
  }
}

# function to get the obs.mean and obs.cov for the date when there is known data available
get_Observed <- function(today, obs.list, dateinfo){
  
  ind <- which(names(obs.list) == today)
  if(length(ind)==0) {
    Observed <- list()
  }else{
     Observed <- obs.list[[ind]]
  }
 
  return(Observed)
}

# function to get the forecast ensemble set up 
get_Forecast <- function(data, date, yr){
  

  
  X <-  data %>%
    filter(day == date, year == yr)%>%
    distinct(ensemble,.keep_all = TRUE) %>% 
    dplyr::select(-ensemble) 
  
 
  
  return(X)
}

# function to set up values for next time step 
nextStep <- function(Xa, Maintained){
  
  results <- colnames(Xa) %>%
    purrr::map(function(col.name) {
      #I remove the numbers from the name and find the state variables with the 
      #exact same name with no numbers - for example for sw1 this gives me all
      # sw1-sw8
      tmp <- Maintained [, grep(paste0(gsub('[[:digit:]]+', '', col.name),"$"),
                                gsub('[[:digit:]]+', '', colnames(Maintained))
      )]
      
      #original names
      orgi.name <- colnames(tmp)
      
      #those that will be replaced
      rep.col.xa <- grep(paste0(gsub('[[:digit:]]+', '', col.name),"$"),
                         gsub('[[:digit:]]+', '', colnames(Xa))
      )
      
      if (is.null(ncol(tmp))){
        tmp <- data.frame(x = tmp)
        colnames(tmp) <- col.name
      } 
      
      # Where would the rep.col.xa will go
      rep.col.main <- purrr::map_dbl(colnames(Xa)[rep.col.xa], ~grep(.x, colnames(tmp)))
        
      #replace the column from Xa with the the one from maintained
      tmp[, rep.col.main] <- Xa[, rep.col.xa]

      tmp[ tmp < 0 ] = 0
      
      tmp%>%
        `colnames<-`(c(orgi.name))
      
      
    }) %>%
    setNames(colnames(Xa))
  
  
  return(results)
}


submit_result_Sqlite <- function(){

  if(!exists('tmp_met')) {
  #read met to get lat/lon
  tmp_met <- readLines("met00000.met")
  # process string
  latlon.df <- data.frame(
    longitude = stringr::str_extract(tmp_met[4], "-[[:digit:]]+.[[:digit:]]+"),
    latitude =stringr::str_extract(tmp_met[3], "[[:digit:]]+.[[:digit:]]+")
  )
  }

# Write to the db
    tryCatch(
      expr = {
        if(!exists('con')) {
        # Start working with the DB ---------------------
        con <<-  DBI::dbConnect(RSQLite::SQLite(),
                               file.path("DA.sqlite3"))
        }
        #  write my model outputs
        DBI::dbWriteTable(con, "Outputs", cbind(latlon.df, data.all) %>%
                            mutate(Date=Date_APSIM_to_R(Time$Year, Time$Day)), append = TRUE)
        #Close the connection
        #DBI::dbDisconnect(con)
      },error = function(e) {
       print(e)
      }
    )

  
}

Date_APSIM_to_R <- function(year, day){
  return(as.character(as.Date(as.Date(day, origin = paste0(year -1, "-12-31")))))
}