###################################################
## Comparaison des WAIC des différents modèles
###################################################

library(R2jags)

# https://gist.github.com/oliviergimenez/70ad17910a62635ff6a062f8ec34292f

# Models
# ~ detection ~ abondance

# Null model 1 : ~ 1 ~ random effect
mnull1.jags <- function(){
  #Set priors for parameters to be estimated
  delta ~ dnorm(0, 0.001)
  mu.pa ~ dnorm(0, 0.001)
  mu.N ~ dnorm(0, 0.001)

  # Random effect of site on abundance
  for(k in 1:nsite){
    rand.N[k] ~ dnorm(0, tau.N)
  }
  sigma.N ~ dunif(0, 10)
  tau.N  <- pow(sigma.N, -2)
  
  ## Build design matrix for N and ps
  for (i in 1:n){
    ## Area effect on detection probability
    logit(pa[i]) <- mu.pa
    logit(pb[i]) <- mu.pa+delta
    
    ## try for dependant observers
    piMat[i, 1] <- ifelse(prim_ran[i]==1, pa[i],
                          pb[i]
    ) # observed by A or B
    
    piMat[i, 2] <- ifelse(sec_ran[i]==1, pa[i] * (1-pb[i]),
                          pb[i] * (1-pa[i])
    ) # observed by B and missed by A or observed by A missed by B
    
    ## Constrains on population size
    log(N[i]) <- mu.N + rand.N[site[i]]
    
    ## Likelihood for the double-observer survey
    for(j in 1:2){
      C[i, j]      ~ dpois(lambda[i, j])
      lambda[i, j] <- N[i] * piMat[i, j]
      fit[i, j]    <- exp(N[i] * piMat[i, j])
    }
  }
}

# Null model 2 : ~ 1 ~ spatial autocorrelation
mnull2.jags <- function(){
  #Set priors for parameters to be estimated
  delta ~ dnorm(0, 0.001)
  mu.pa ~ dnorm(0, 0.001)
  
  # priors for the spatial component
  global.mu_lanslevillard ~ dnorm(0, 0.01)
  #global.mu_lanslevillard <- 0
  global.tau_lanslevillard ~ dgamma(0.001, 0.001)
  for(i in 1:35){
    mu_lanslevillard[i] ~ dnorm(global.mu_lanslevillard, global.tau_lanslevillard)
  }     
  
  global.mu_beaufort ~ dnorm(0, 0.01)
  #global.mu_beaufort <- 0
  global.tau_beaufort ~ dgamma(0.001, 0.001)
  for(i in 1:35){
    mu_beaufort[i] ~ dnorm(global.mu_beaufort, global.tau_beaufort)
  }                 
  
  # Spatial component
  lambda_spat ~ dgamma(1, 0.1)
  sigmasq <- 1/sigmasq.inv
  sigmasq.inv ~ dgamma(2, 1)
  
  ## Build design matrix for N and ps
  for (i in 1:n){
    ## Area effect on detection probability
    logit(pa[i]) <- mu.pa
    logit(pb[i]) <- mu.pa+delta
    
    ## try for dependant observers
    piMat[i, 1] <- ifelse(prim_ran[i]==1, pa[i],
                          pb[i]
    ) # observed by A or B
    
    piMat[i, 2] <- ifelse(sec_ran[i]==1, pa[i] * (1-pb[i]),
                          pb[i] * (1-pa[i])
    ) # observed by B and missed by A or observed by A missed by B
    
    ## Constrains on population size
    log(N[i]) <- CT[site[i]]
    
    ## Likelihood for the double-observer survey
    for(j in 1:2){
      C[i, j]      ~ dpois(lambda[i, j])
      lambda[i, j] <- N[i] * piMat[i, j]
      fit[i, j]    <- exp(N[i] * piMat[i, j])
    }
  }
  
  # Spatial component Les beaufort
  CT[1:35] ~ dmnorm.vcov(mu_beaufort[1:35], D.covar_beaufort[1:35, 1:35])
  ## turning covariances into precisions
  for(i in 1:35){
    D.covar_beaufort[i, i] <- sigmasq # diagonale de la matrice de covariance
    for(j in 1:(i - 1)){
      # constrain covariance matrix by the distance matrix
      D.covar_beaufort[i, j] <- sigmasq * exp(-(lambda_spat * D_beaufort[i, j]))
      D.covar_beaufort[j, i] <- D.covar_beaufort[i, j] # symmetrical matrix
    }
  }
  
  # Spatial component lanslevillard
  CT[36:70] ~ dmnorm.vcov(mu_lanslevillard[1:35], D.covar_lanslevillard[1:35, 1:35])
  ## turning covariances into precisions
  for(i in 1:35){
    D.covar_lanslevillard[i, i] <- sigmasq # diagonale de la matrice de covariance
    for(j in 1:(i - 1)){
      # constrain covariance matrix by the distance matrix
      D.covar_lanslevillard[i, j] <- sigmasq * exp(-(lambda_spat * D_lanslevillard[i, j]))
      D.covar_lanslevillard[j, i] <- D.covar_lanslevillard[i, j] # symmetrical matrix
    }
  }
}

# Model 1 : ~ secteur ~ spatial autocorrelation
msecteurP.jags <- function(){
  #Set priors for parameters to be estimated
delta ~ dnorm(0, 0.001)

# Fixed parameters for detection
for(k in 1:npar.secteur){
  beta.secteurP[k] ~ dnorm(0, 0.0001)
}

# priors for the spatial component
global.mu_lanslevillard ~ dnorm(0, 0.01)
#global.mu_lanslevillard <- 0
global.tau_lanslevillard ~ dgamma(0.001, 0.001)
for(i in 1:35){
  mu_lanslevillard[i] ~ dnorm(global.mu_lanslevillard, global.tau_lanslevillard)
}     

global.mu_beaufort ~ dnorm(0, 0.01)
#global.mu_beaufort <- 0
global.tau_beaufort ~ dgamma(0.001, 0.001)
for(i in 1:35){
  mu_beaufort[i] ~ dnorm(global.mu_beaufort, global.tau_beaufort)
}                 

# Spatial component
lambda_spat ~ dgamma(1, 0.1)
sigmasq <- 1/sigmasq.inv
sigmasq.inv ~ dgamma(2, 1)

## Build design matrix for N and ps
for (i in 1:n){
  ## Area effect on detection probability
  logit(pa[i]) <- inprod(beta.secteurP[], Xsecteur[i,])
  logit(pb[i]) <- inprod(beta.secteurP[], Xsecteur[i,])+delta
  
  ## try for dependant observers
  piMat[i, 1] <- ifelse(prim_ran[i]==1, pa[i],
                        pb[i]
  ) # observed by A or B
  
  piMat[i, 2] <- ifelse(sec_ran[i]==1, pa[i] * (1-pb[i]),
                        pb[i] * (1-pa[i])
  ) # observed by B and missed by A or observed by A missed by B
  
  ## Constrains on population size
  log(N[i]) <- CT[site[i]]

  ## Likelihood for the double-observer survey
  for(j in 1:2){
    C[i, j]      ~ dpois(lambda[i, j])
    lambda[i, j] <- N[i] * piMat[i, j]
    fit[i, j]    <- exp(N[i] * piMat[i, j])
  }
}

# Spatial component Les beaufort
CT[1:35] ~ dmnorm.vcov(mu_beaufort[1:35], D.covar_beaufort[1:35, 1:35])
## turning covariances into precisions
for(i in 1:35){
  D.covar_beaufort[i, i] <- sigmasq # diagonale de la matrice de covariance
  for(j in 1:(i - 1)){
    # constrain covariance matrix by the distance matrix
    D.covar_beaufort[i, j] <- sigmasq * exp(-(lambda_spat * D_beaufort[i, j]))
    D.covar_beaufort[j, i] <- D.covar_beaufort[i, j] # symmetrical matrix
  }
}

# Spatial component lanslevillard
CT[36:70] ~ dmnorm.vcov(mu_lanslevillard[1:35], D.covar_lanslevillard[1:35, 1:35])
## turning covariances into precisions
for(i in 1:35){
  D.covar_lanslevillard[i, i] <- sigmasq # diagonale de la matrice de covariance
  for(j in 1:(i - 1)){
    # constrain covariance matrix by the distance matrix
    D.covar_lanslevillard[i, j] <- sigmasq * exp(-(lambda_spat * D_lanslevillard[i, j]))
    D.covar_lanslevillard[j, i] <- D.covar_lanslevillard[i, j] # symmetrical matrix
  }
}
}

# Model 2 : ~ secteur+derangement ~ spatial autocorrelation
mderangementP.jags <- function(){
  ## Set priors for parameters to be estimated
  delta ~ dnorm(0, 0.001)
  
  # Fixed parameters for detection
  for(k in 1:npar.secteur.derangement){
    beta.secteur.derangementP[k] ~ dnorm(0, 0.0001)
  }
  
  # priors for the spatial component
  global.mu_lanslevillard ~ dnorm(0, 0.01)
  global.tau_lanslevillard ~ dgamma(0.001, 0.001)
  for(i in 1:35){
    mu_lanslevillard[i] ~ dnorm(global.mu_lanslevillard, global.tau_lanslevillard)
  }     
  
  global.mu_beaufort ~ dnorm(0, 0.01)
  global.tau_beaufort ~ dgamma(0.001, 0.001)
  for(i in 1:35){
    mu_beaufort[i] ~ dnorm(global.mu_beaufort, global.tau_beaufort)
  }                 
  
  # Spatial component
  lambda_spat ~ dgamma(1, 0.1)
  sigmasq <- 1/sigmasq.inv
  sigmasq.inv ~ dgamma(2, 1)
  
  ## Build design matrix for N and ps
  for (i in 1:n){
    ## Area effect on detection probability
    logit(pa[i]) <- inprod(beta.secteur.derangementP[], Xsecteur.derangement[i,])
    logit(pb[i]) <- inprod(beta.secteur.derangementP[], Xsecteur.derangement[i,])+delta
    
    ## try for dependant observers
    piMat[i, 1] <- ifelse(prim_ran[i]==1, pa[i],
                          pb[i]
    ) # observed by A or B
    
    piMat[i, 2] <- ifelse(sec_ran[i]==1, pa[i] * (1-pb[i]),
                          pb[i] * (1-pa[i])
    ) # observed by B and missed by A or observed by A missed by B
    
    ## Constrains on population size
    log(N[i]) <- CT[site[i]] 
    
    ## Likelihood for the double-observer survey
    for(j in 1:2){
      C[i, j]      ~ dpois(lambda[i, j])
      lambda[i, j] <- N[i] * piMat[i, j]
      fit[i, j]    <- exp(N[i] * piMat[i, j])
    }
    
  }
  
  # Spatial component Les beaufort
  CT[1:35] ~ dmnorm.vcov(mu_beaufort[1:35], D.covar_beaufort[1:35, 1:35])
  ## turning covariances into precisions
  for(i in 1:35){
    D.covar_beaufort[i, i] <- sigmasq # diagonale de la matrice de covariance
    for(j in 1:(i - 1)){
      # constrain covariance matrix by the distance matrix
      D.covar_beaufort[i, j] <- sigmasq * exp(-(lambda_spat * D_beaufort[i, j]))
      D.covar_beaufort[j, i] <- D.covar_beaufort[i, j] # symmetrical matrix
    }
  }
  
  # Spatial component lanslevillard
  CT[36:70] ~ dmnorm.vcov(mu_lanslevillard[1:35], D.covar_lanslevillard[1:35, 1:35])
  ## turning covariances into precisions
  for(i in 1:35){
    D.covar_lanslevillard[i, i] <- sigmasq # diagonale de la matrice de covariance
    for(j in 1:(i - 1)){
      # constrain covariance matrix by the distance matrix
      D.covar_lanslevillard[i, j] <- sigmasq * exp(-(lambda_spat * D_lanslevillard[i, j]))
      D.covar_lanslevillard[j, i] <- D.covar_lanslevillard[i, j] # symmetrical matrix
    }
  }
}


# Model 3 : ~ secteur+orientation ~ spatial autocorrelation
morientationP.jags <- function(){
  ## Set priors for parameters to be estimated
  delta ~ dnorm(0, 0.001)
  
  # Fixed parameters for detection
  for(k in 1:npar.secteur.orientation){
    beta.secteur.orientationP[k] ~ dnorm(0, 0.0001)
  }
  
  # priors for the spatial component
  global.mu_lanslevillard ~ dnorm(0, 0.01)
  global.tau_lanslevillard ~ dgamma(0.001, 0.001)
  for(i in 1:35){
    mu_lanslevillard[i] ~ dnorm(global.mu_lanslevillard, global.tau_lanslevillard)
  }     
  
  global.mu_beaufort ~ dnorm(0, 0.01)
  global.tau_beaufort ~ dgamma(0.001, 0.001)
  for(i in 1:35){
    mu_beaufort[i] ~ dnorm(global.mu_beaufort, global.tau_beaufort)
  }                 
  
  # Spatial component
  lambda_spat ~ dgamma(1, 0.1)
  sigmasq <- 1/sigmasq.inv
  sigmasq.inv ~ dgamma(2, 1)
  
  ## Build design matrix for N and ps
  for (i in 1:n){
    ## Area effect on detection probability
    logit(pa[i]) <- inprod(beta.secteur.orientationP[], Xsecteur.orientation[i,])
    logit(pb[i]) <- inprod(beta.secteur.orientationP[], Xsecteur.orientation[i,])+delta
    
    ## try for dependant observers
    piMat[i, 1] <- ifelse(prim_ran[i]==1, pa[i],
                          pb[i]
    ) # observed by A or B
    
    piMat[i, 2] <- ifelse(sec_ran[i]==1, pa[i] * (1-pb[i]),
                          pb[i] * (1-pa[i])
    ) # observed by B and missed by A or observed by A missed by B
    
    ## Constrains on population size
    log(N[i]) <- CT[site[i]] 
    
    ## Likelihood for the double-observer survey
    for(j in 1:2){
      C[i, j]      ~ dpois(lambda[i, j])
      lambda[i, j] <- N[i] * piMat[i, j]
      fit[i, j]    <- exp(N[i] * piMat[i, j])
    }
    
  }
  
  # Spatial component Les beaufort
  CT[1:35] ~ dmnorm.vcov(mu_beaufort[1:35], D.covar_beaufort[1:35, 1:35])
  ## turning covariances into precisions
  for(i in 1:35){
    D.covar_beaufort[i, i] <- sigmasq # diagonale de la matrice de covariance
    for(j in 1:(i - 1)){
      # constrain covariance matrix by the distance matrix
      D.covar_beaufort[i, j] <- sigmasq * exp(-(lambda_spat * D_beaufort[i, j]))
      D.covar_beaufort[j, i] <- D.covar_beaufort[i, j] # symmetrical matrix
    }
  }
  
  # Spatial component lanslevillard
  CT[36:70] ~ dmnorm.vcov(mu_lanslevillard[1:35], D.covar_lanslevillard[1:35, 1:35])
  ## turning covariances into precisions
  for(i in 1:35){
    D.covar_lanslevillard[i, i] <- sigmasq # diagonale de la matrice de covariance
    for(j in 1:(i - 1)){
      # constrain covariance matrix by the distance matrix
      D.covar_lanslevillard[i, j] <- sigmasq * exp(-(lambda_spat * D_lanslevillard[i, j]))
      D.covar_lanslevillard[j, i] <- D.covar_lanslevillard[i, j] # symmetrical matrix
    }
  }
}

# Model 4 : ~ secteur+meteo ~ spatial autocorrelation
mmeteoP.jags <- function(){
  ## Set priors for parameters to be estimated
  delta ~ dnorm(0, 0.001)
  
  # Fixed parameters for detection
  for(k in 1:npar.secteur.meteo){
    beta.secteur.meteoP[k] ~ dnorm(0, 0.0001)
  }
  
  # priors for the spatial component
  global.mu_lanslevillard ~ dnorm(0, 0.01)
  global.tau_lanslevillard ~ dgamma(0.001, 0.001)
  for(i in 1:35){
    mu_lanslevillard[i] ~ dnorm(global.mu_lanslevillard, global.tau_lanslevillard)
  }     
  
  global.mu_beaufort ~ dnorm(0, 0.01)
  global.tau_beaufort ~ dgamma(0.001, 0.001)
  for(i in 1:35){
    mu_beaufort[i] ~ dnorm(global.mu_beaufort, global.tau_beaufort)
  }                 
  
  # Spatial component
  lambda_spat ~ dgamma(1, 0.1)
  sigmasq <- 1/sigmasq.inv
  sigmasq.inv ~ dgamma(2, 1)
  
  ## Build design matrix for N and ps
  for (i in 1:n){
    ## Area effect on detection probability
    logit(pa[i]) <- inprod(beta.secteur.meteoP[], Xsecteur.meteo[i,])
    logit(pb[i]) <- inprod(beta.secteur.meteoP[], Xsecteur.meteo[i,])+delta
    
    ## try for dependant observers
    piMat[i, 1] <- ifelse(prim_ran[i]==1, pa[i],
                          pb[i]
    ) # observed by A or B
    
    piMat[i, 2] <- ifelse(sec_ran[i]==1, pa[i] * (1-pb[i]),
                          pb[i] * (1-pa[i])
    ) # observed by B and missed by A or observed by A missed by B
    
    ## Constrains on population size
    log(N[i]) <- CT[site[i]] 
    
    ## Likelihood for the double-observer survey
    for(j in 1:2){
      C[i, j]      ~ dpois(lambda[i, j])
      lambda[i, j] <- N[i] * piMat[i, j]
      fit[i, j]    <- exp(N[i] * piMat[i, j])
    }
    
  }
  
  # Spatial component Les beaufort
  CT[1:35] ~ dmnorm.vcov(mu_beaufort[1:35], D.covar_beaufort[1:35, 1:35])
  ## turning covariances into precisions
  for(i in 1:35){
    D.covar_beaufort[i, i] <- sigmasq # diagonale de la matrice de covariance
    for(j in 1:(i - 1)){
      # constrain covariance matrix by the distance matrix
      D.covar_beaufort[i, j] <- sigmasq * exp(-(lambda_spat * D_beaufort[i, j]))
      D.covar_beaufort[j, i] <- D.covar_beaufort[i, j] # symmetrical matrix
    }
  }
  
  # Spatial component lanslevillard
  CT[36:70] ~ dmnorm.vcov(mu_lanslevillard[1:35], D.covar_lanslevillard[1:35, 1:35])
  ## turning covariances into precisions
  for(i in 1:35){
    D.covar_lanslevillard[i, i] <- sigmasq # diagonale de la matrice de covariance
    for(j in 1:(i - 1)){
      # constrain covariance matrix by the distance matrix
      D.covar_lanslevillard[i, j] <- sigmasq * exp(-(lambda_spat * D_lanslevillard[i, j]))
      D.covar_lanslevillard[j, i] <- D.covar_lanslevillard[i, j] # symmetrical matrix
    }
  }
}

# Model 5 : ~ secteur+session ~ spatial autocorrelation
msessionP.jags <- function(){
  ## Set priors for parameters to be estimated
  delta ~ dnorm(0, 0.001)
  
  # Fixed parameters for detection
  for(k in 1:npar.secteur.session){
    beta.secteur.sessionP[k] ~ dnorm(0, 0.0001)
  }
  
  # priors for the spatial component
  global.mu_lanslevillard ~ dnorm(0, 0.01)
  global.tau_lanslevillard ~ dgamma(0.001, 0.001)
  for(i in 1:35){
    mu_lanslevillard[i] ~ dnorm(global.mu_lanslevillard, global.tau_lanslevillard)
  }     
  
  global.mu_beaufort ~ dnorm(0, 0.01)
  global.tau_beaufort ~ dgamma(0.001, 0.001)
  for(i in 1:35){
    mu_beaufort[i] ~ dnorm(global.mu_beaufort, global.tau_beaufort)
  }                 
  
  # Spatial component
  lambda_spat ~ dgamma(1, 0.1)
  sigmasq <- 1/sigmasq.inv
  sigmasq.inv ~ dgamma(2, 1)
  
  ## Build design matrix for N and ps
  for (i in 1:n){
    ## Area effect on detection probability
    logit(pa[i]) <- inprod(beta.secteur.sessionP[], Xsecteur.session[i,])
    logit(pb[i]) <- inprod(beta.secteur.sessionP[], Xsecteur.session[i,])+delta

    ## try for dependant observers
    piMat[i, 1] <- ifelse(prim_ran[i]==1, pa[i],
                          pb[i]
    ) # observed by A or B
    
    piMat[i, 2] <- ifelse(sec_ran[i]==1, pa[i] * (1-pb[i]),
                          pb[i] * (1-pa[i])
    ) # observed by B and missed by A or observed by A missed by B
    
    ## Constrains on population size
    log(N[i]) <- CT[site[i]] 
    
    ## Likelihood for the double-observer survey
    for(j in 1:2){
      C[i, j]      ~ dpois(lambda[i, j])
      lambda[i, j] <- N[i] * piMat[i, j]
      fit[i, j]    <- exp(N[i] * piMat[i, j])
    }
    
  }
  
  # Spatial component Les beaufort
  CT[1:35] ~ dmnorm.vcov(mu_beaufort[1:35], D.covar_beaufort[1:35, 1:35])
  ## turning covariances into precisions
  for(i in 1:35){
    D.covar_beaufort[i, i] <- sigmasq # diagonale de la matrice de covariance
    for(j in 1:(i - 1)){
      # constrain covariance matrix by the distance matrix
      D.covar_beaufort[i, j] <- sigmasq * exp(-(lambda_spat * D_beaufort[i, j]))
      D.covar_beaufort[j, i] <- D.covar_beaufort[i, j] # symmetrical matrix
    }
  }
  
  # Spatial component lanslevillard
  CT[36:70] ~ dmnorm.vcov(mu_lanslevillard[1:35], D.covar_lanslevillard[1:35, 1:35])
  ## turning covariances into precisions
  for(i in 1:35){
    D.covar_lanslevillard[i, i] <- sigmasq # diagonale de la matrice de covariance
    for(j in 1:(i - 1)){
      # constrain covariance matrix by the distance matrix
      D.covar_lanslevillard[i, j] <- sigmasq * exp(-(lambda_spat * D_lanslevillard[i, j]))
      D.covar_lanslevillard[j, i] <- D.covar_lanslevillard[i, j] # symmetrical matrix
    }
  }
}

# Model 6 : ~ secteur+derangement+meteo ~ spatial autocorrelation
mderangement.meteoP.jags <- function(){
  ## Set priors for parameters to be estimated
  delta ~ dnorm(0, 0.001)
  
  # Fixed parameters for detection
  for(k in 1:npar.secteur.derangement.meteo){
    beta.secteur.derangement.meteoP[k] ~ dnorm(0, 0.0001)
  }
  
  # priors for the spatial component
  global.mu_lanslevillard ~ dnorm(0, 0.01)
  global.tau_lanslevillard ~ dgamma(0.001, 0.001)
  for(i in 1:35){
    mu_lanslevillard[i] ~ dnorm(global.mu_lanslevillard, global.tau_lanslevillard)
  }     
  
  global.mu_beaufort ~ dnorm(0, 0.01)
  global.tau_beaufort ~ dgamma(0.001, 0.001)
  for(i in 1:35){
    mu_beaufort[i] ~ dnorm(global.mu_beaufort, global.tau_beaufort)
  }                 
  
  # Spatial component
  lambda_spat ~ dgamma(1, 0.1)
  sigmasq <- 1/sigmasq.inv
  sigmasq.inv ~ dgamma(2, 1)
  
  ## Build design matrix for N and ps
  for (i in 1:n){
    ## Area effect on detection probability
    logit(pa[i]) <- inprod(beta.secteur.derangement.meteoP[], Xsecteur.derangement.meteo[i,])
    logit(pb[i]) <- inprod(beta.secteur.derangement.meteoP[], Xsecteur.derangement.meteo[i,])+delta
    
    ## try for dependant observers
    piMat[i, 1] <- ifelse(prim_ran[i]==1, pa[i],
                          pb[i]
    ) # observed by A or B
    
    piMat[i, 2] <- ifelse(sec_ran[i]==1, pa[i] * (1-pb[i]),
                          pb[i] * (1-pa[i])
    ) # observed by B and missed by A or observed by A missed by B
    
    ## Constrains on population size
    log(N[i]) <- CT[site[i]] 
    
    ## Likelihood for the double-observer survey
    for(j in 1:2){
      C[i, j]      ~ dpois(lambda[i, j])
      lambda[i, j] <- N[i] * piMat[i, j]
      fit[i, j]    <- exp(N[i] * piMat[i, j])
    }
    
  }
  
  # Spatial component Les beaufort
  CT[1:35] ~ dmnorm.vcov(mu_beaufort[1:35], D.covar_beaufort[1:35, 1:35])
  ## turning covariances into precisions
  for(i in 1:35){
    D.covar_beaufort[i, i] <- sigmasq # diagonale de la matrice de covariance
    for(j in 1:(i - 1)){
      # constrain covariance matrix by the distance matrix
      D.covar_beaufort[i, j] <- sigmasq * exp(-(lambda_spat * D_beaufort[i, j]))
      D.covar_beaufort[j, i] <- D.covar_beaufort[i, j] # symmetrical matrix
    }
  }
  
  # Spatial component lanslevillard
  CT[36:70] ~ dmnorm.vcov(mu_lanslevillard[1:35], D.covar_lanslevillard[1:35, 1:35])
  ## turning covariances into precisions
  for(i in 1:35){
    D.covar_lanslevillard[i, i] <- sigmasq # diagonale de la matrice de covariance
    for(j in 1:(i - 1)){
      # constrain covariance matrix by the distance matrix
      D.covar_lanslevillard[i, j] <- sigmasq * exp(-(lambda_spat * D_lanslevillard[i, j]))
      D.covar_lanslevillard[j, i] <- D.covar_lanslevillard[i, j] # symmetrical matrix
    }
  }
}


# Model 7 : ~ secteur ~ altitude+spatial autocorrelation
maltitudeN.jags <- function(){
  #Set priors for parameters to be estimated
  delta ~ dnorm(0, 0.001)
  beta.altN ~ dnorm(0, 0.0001)
  
  # Fixed parameters for detection
  for(k in 1:npar.secteur){
    beta.secteurP[k] ~ dnorm(0, 0.0001)
  }
  
  # priors for the spatial component
  global.mu_lanslevillard ~ dnorm(0, 0.01)
  #global.mu_lanslevillard <- 0
  global.tau_lanslevillard ~ dgamma(0.001, 0.001)
  for(i in 1:35){
    mu_lanslevillard[i] ~ dnorm(global.mu_lanslevillard, global.tau_lanslevillard)
  }     
  
  global.mu_beaufort ~ dnorm(0, 0.01)
  #global.mu_beaufort <- 0
  global.tau_beaufort ~ dgamma(0.001, 0.001)
  for(i in 1:35){
    mu_beaufort[i] ~ dnorm(global.mu_beaufort, global.tau_beaufort)
  }                 
  
  # Spatial component
  lambda_spat ~ dgamma(1, 0.1)
  sigmasq <- 1/sigmasq.inv
  sigmasq.inv ~ dgamma(2, 1)
  
  ## Build design matrix for N and ps
  for (i in 1:n){
    ## Area effect on detection probability
    logit(pa[i]) <- inprod(beta.secteurP[], Xsecteur[i,])
    logit(pb[i]) <- inprod(beta.secteurP[], Xsecteur[i,])+delta
    
    ## try for dependant observers
    piMat[i, 1] <- ifelse(prim_ran[i]==1, pa[i],
                          pb[i]
    ) # observed by A or B
    
    piMat[i, 2] <- ifelse(sec_ran[i]==1, pa[i] * (1-pb[i]),
                          pb[i] * (1-pa[i])
    ) # observed by B and missed by A or observed by A missed by B
    
    ## Constrains on population size
    log(N[i]) <- beta.altN*altitude[i] + CT[site[i]] 
    
    ## Likelihood for the double-observer survey
    for(j in 1:2){
      C[i, j]      ~ dpois(lambda[i, j])
      lambda[i, j] <- N[i] * piMat[i, j]
      fit[i, j]    <- exp(N[i] * piMat[i, j])
    }
  }
  
  # Spatial component Les beaufort
  CT[1:35] ~ dmnorm.vcov(mu_beaufort[1:35], D.covar_beaufort[1:35, 1:35])
  ## turning covariances into precisions
  for(i in 1:35){
    D.covar_beaufort[i, i] <- sigmasq # diagonale de la matrice de covariance
    for(j in 1:(i - 1)){
      # constrain covariance matrix by the distance matrix
      D.covar_beaufort[i, j] <- sigmasq * exp(-(lambda_spat * D_beaufort[i, j]))
      D.covar_beaufort[j, i] <- D.covar_beaufort[i, j] # symmetrical matrix
    }
  }
  
  # Spatial component lanslevillard
  CT[36:70] ~ dmnorm.vcov(mu_lanslevillard[1:35], D.covar_lanslevillard[1:35, 1:35])
  ## turning covariances into precisions
  for(i in 1:35){
    D.covar_lanslevillard[i, i] <- sigmasq # diagonale de la matrice de covariance
    for(j in 1:(i - 1)){
      # constrain covariance matrix by the distance matrix
      D.covar_lanslevillard[i, j] <- sigmasq * exp(-(lambda_spat * D_lanslevillard[i, j]))
      D.covar_lanslevillard[j, i] <- D.covar_lanslevillard[i, j] # symmetrical matrix
    }
  }
}

# Model 7 : ~ secteur ~ habitat+spatial autocorrelation
mhabitatN.jags <- function(){
  ## Set priors for parameters to be estimated
  delta ~ dnorm(0, 0.001)
  
  # Fixed parameters for detection
  for(k in 1:npar.secteur){
    beta.secteurP[k] ~ dnorm(0, 0.0001)
  }
  
  for(k in 1:npar.habitat){
    beta.habitatN[k] ~ dnorm(0, 0.0001)
  }
  
  # priors for the spatial component
  global.mu_lanslevillard ~ dnorm(0, 0.01)
  global.tau_lanslevillard ~ dgamma(0.001, 0.001)
  for(i in 1:35){
    mu_lanslevillard[i] ~ dnorm(global.mu_lanslevillard, global.tau_lanslevillard)
  }     
  
  global.mu_beaufort ~ dnorm(0, 0.01)
  global.tau_beaufort ~ dgamma(0.001, 0.001)
  for(i in 1:35){
    mu_beaufort[i] ~ dnorm(global.mu_beaufort, global.tau_beaufort)
  }                 
  
  # Spatial component
  lambda_spat ~ dgamma(1, 0.1)
  sigmasq <- 1/sigmasq.inv
  sigmasq.inv ~ dgamma(2, 1)
  
  ## Build design matrix for N and ps
  for (i in 1:n){
    ## Area effect on detection probability
    logit(pa[i]) <- inprod(beta.secteurP[], Xsecteur[i,])
    logit(pb[i]) <- inprod(beta.secteurP[], Xsecteur[i,])+delta
    
    ## try for dependant observers
    piMat[i, 1] <- ifelse(prim_ran[i]==1, pa[i],
                          pb[i]
    ) # observed by A or B
    
    piMat[i, 2] <- ifelse(sec_ran[i]==1, pa[i] * (1-pb[i]),
                          pb[i] * (1-pa[i])
    ) # observed by B and missed by A or observed by A missed by B
    
    ## Constrains on population size
    log(N[i]) <- inprod(beta.habitatN[], Xhabitat[i, ]) + CT[site[i]] 
    
    ## Likelihood for the double-observer survey
    for(j in 1:2){
      C[i, j]      ~ dpois(lambda[i, j])
      lambda[i, j] <- N[i] * piMat[i, j]
      fit[i, j]    <- exp(N[i] * piMat[i, j])
    }
    
  }
  
  # Spatial component Les beaufort
  CT[1:35] ~ dmnorm.vcov(mu_beaufort[1:35], D.covar_beaufort[1:35, 1:35])
  ## turning covariances into precisions
  for(i in 1:35){
    D.covar_beaufort[i, i] <- sigmasq # diagonale de la matrice de covariance
    for(j in 1:(i - 1)){
      # constrain covariance matrix by the distance matrix
      D.covar_beaufort[i, j] <- sigmasq * exp(-(lambda_spat * D_beaufort[i, j]))
      D.covar_beaufort[j, i] <- D.covar_beaufort[i, j] # symmetrical matrix
    }
  }
  
  # Spatial component lanslevillard
  CT[36:70] ~ dmnorm.vcov(mu_lanslevillard[1:35], D.covar_lanslevillard[1:35, 1:35])
  ## turning covariances into precisions
  for(i in 1:35){
    D.covar_lanslevillard[i, i] <- sigmasq # diagonale de la matrice de covariance
    for(j in 1:(i - 1)){
      # constrain covariance matrix by the distance matrix
      D.covar_lanslevillard[i, j] <- sigmasq * exp(-(lambda_spat * D_lanslevillard[i, j]))
      D.covar_lanslevillard[j, i] <- D.covar_lanslevillard[i, j] # symmetrical matrix
    }
  }
}

# Model 8 : ~ secteur+derangement ~ altitude+spatial autocorrelation
mderangementP.altitudeN.jags <- function(){
  ## Set priors for parameters to be estimated
  delta ~ dnorm(0, 0.001)
  beta.altN ~ dnorm(0, 0.001)
  
  # Fixed parameters for detection
  for(k in 1:npar.secteur.derangement){
    beta.secteur.derangementP[k] ~ dnorm(0, 0.0001)
  }
  
  # priors for the spatial component
  global.mu_lanslevillard ~ dnorm(0, 0.01)
  global.tau_lanslevillard ~ dgamma(0.001, 0.001)
  for(i in 1:35){
    mu_lanslevillard[i] ~ dnorm(global.mu_lanslevillard, global.tau_lanslevillard)
  }     
  
  global.mu_beaufort ~ dnorm(0, 0.01)
  global.tau_beaufort ~ dgamma(0.001, 0.001)
  for(i in 1:35){
    mu_beaufort[i] ~ dnorm(global.mu_beaufort, global.tau_beaufort)
  }                 
  
  # Spatial component
  lambda_spat ~ dgamma(1, 0.1)
  sigmasq <- 1/sigmasq.inv
  sigmasq.inv ~ dgamma(2, 1)
  
  ## Build design matrix for N and ps
  for (i in 1:n){
    ## Area effect on detection probability
    logit(pa[i]) <- inprod(beta.secteur.derangementP[], Xsecteur.derangement[i,])
    logit(pb[i]) <- inprod(beta.secteur.derangementP[], Xsecteur.derangement[i,])+delta
    
    ## try for dependant observers
    piMat[i, 1] <- ifelse(prim_ran[i]==1, pa[i],
                          pb[i]
    ) # observed by A or B
    
    piMat[i, 2] <- ifelse(sec_ran[i]==1, pa[i] * (1-pb[i]),
                          pb[i] * (1-pa[i])
    ) # observed by B and missed by A or observed by A missed by B
    
    ## Constrains on population size
    log(N[i]) <- beta.altN*altitude[i] + CT[site[i]] 
    
    ## Likelihood for the double-observer survey
    for(j in 1:2){
      C[i, j]      ~ dpois(lambda[i, j])
      lambda[i, j] <- N[i] * piMat[i, j]
      fit[i, j]    <- exp(N[i] * piMat[i, j])
    }
    
  }
  
  # Spatial component Les beaufort
  CT[1:35] ~ dmnorm.vcov(mu_beaufort[1:35], D.covar_beaufort[1:35, 1:35])
  ## turning covariances into precisions
  for(i in 1:35){
    D.covar_beaufort[i, i] <- sigmasq # diagonale de la matrice de covariance
    for(j in 1:(i - 1)){
      # constrain covariance matrix by the distance matrix
      D.covar_beaufort[i, j] <- sigmasq * exp(-(lambda_spat * D_beaufort[i, j]))
      D.covar_beaufort[j, i] <- D.covar_beaufort[i, j] # symmetrical matrix
    }
  }
  
  # Spatial component lanslevillard
  CT[36:70] ~ dmnorm.vcov(mu_lanslevillard[1:35], D.covar_lanslevillard[1:35, 1:35])
  ## turning covariances into precisions
  for(i in 1:35){
    D.covar_lanslevillard[i, i] <- sigmasq # diagonale de la matrice de covariance
    for(j in 1:(i - 1)){
      # constrain covariance matrix by the distance matrix
      D.covar_lanslevillard[i, j] <- sigmasq * exp(-(lambda_spat * D_lanslevillard[i, j]))
      D.covar_lanslevillard[j, i] <- D.covar_lanslevillard[i, j] # symmetrical matrix
    }
  }
}

# Model 9 : ~ secteur+derangement ~ habitat+spatial autocorrelation
mderangementP.habitatN.jags <- function(){
  ## Set priors for parameters to be estimated
  delta ~ dnorm(0, 0.001)

  # Fixed parameters for detection
  for(k in 1:npar.secteur.derangement){
    beta.secteur.derangementP[k] ~ dnorm(0, 0.0001)
  }
  
  for(k in 1:npar.habitat){
    beta.habitatN[k] ~ dnorm(0, 0.0001)
  }
  
  # priors for the spatial component
  global.mu_lanslevillard ~ dnorm(0, 0.01)
  global.tau_lanslevillard ~ dgamma(0.001, 0.001)
  for(i in 1:35){
    mu_lanslevillard[i] ~ dnorm(global.mu_lanslevillard, global.tau_lanslevillard)
  }     
  
  global.mu_beaufort ~ dnorm(0, 0.01)
  global.tau_beaufort ~ dgamma(0.001, 0.001)
  for(i in 1:35){
    mu_beaufort[i] ~ dnorm(global.mu_beaufort, global.tau_beaufort)
  }                 
  
  # Spatial component
  lambda_spat ~ dgamma(1, 0.1)
  sigmasq <- 1/sigmasq.inv
  sigmasq.inv ~ dgamma(2, 1)
  
  ## Build design matrix for N and ps
  for (i in 1:n){
    ## Area effect on detection probability
    logit(pa[i]) <- inprod(beta.secteur.derangementP[], Xsecteur.derangement[i,])
    logit(pb[i]) <- inprod(beta.secteur.derangementP[], Xsecteur.derangement[i,])+delta
    
    ## try for dependant observers
    piMat[i, 1] <- ifelse(prim_ran[i]==1, pa[i],
                          pb[i]
    ) # observed by A or B
    
    piMat[i, 2] <- ifelse(sec_ran[i]==1, pa[i] * (1-pb[i]),
                          pb[i] * (1-pa[i])
    ) # observed by B and missed by A or observed by A missed by B
    
    ## Constrains on population size
    log(N[i]) <- inprod(beta.habitatN[], Xhabitat[i, ]) + CT[site[i]] 
    
    ## Likelihood for the double-observer survey
    for(j in 1:2){
      C[i, j]      ~ dpois(lambda[i, j])
      lambda[i, j] <- N[i] * piMat[i, j]
      fit[i, j]    <- exp(N[i] * piMat[i, j])
    }
    
  }
  
  # Spatial component Les beaufort
  CT[1:35] ~ dmnorm.vcov(mu_beaufort[1:35], D.covar_beaufort[1:35, 1:35])
  ## turning covariances into precisions
  for(i in 1:35){
    D.covar_beaufort[i, i] <- sigmasq # diagonale de la matrice de covariance
    for(j in 1:(i - 1)){
      # constrain covariance matrix by the distance matrix
      D.covar_beaufort[i, j] <- sigmasq * exp(-(lambda_spat * D_beaufort[i, j]))
      D.covar_beaufort[j, i] <- D.covar_beaufort[i, j] # symmetrical matrix
    }
  }
  
  # Spatial component lanslevillard
  CT[36:70] ~ dmnorm.vcov(mu_lanslevillard[1:35], D.covar_lanslevillard[1:35, 1:35])
  ## turning covariances into precisions
  for(i in 1:35){
    D.covar_lanslevillard[i, i] <- sigmasq # diagonale de la matrice de covariance
    for(j in 1:(i - 1)){
      # constrain covariance matrix by the distance matrix
      D.covar_lanslevillard[i, j] <- sigmasq * exp(-(lambda_spat * D_lanslevillard[i, j]))
      D.covar_lanslevillard[j, i] <- D.covar_lanslevillard[i, j] # symmetrical matrix
    }
  }
}

# Model 10 : ~ secteur+derangement ~ altitude*secteur+spatial autocorrelation
mderangementP.secteur.altitudeN.jags <- function(){
  ## Set priors for parameters to be estimated
  delta ~ dnorm(0, 0.001)

  # Fixed parameters for detection
  for(k in 1:npar.secteur.derangement){
    beta.secteur.derangementP[k] ~ dnorm(0, 0.0001)
  }
  
  # Fixed parameters for abundance
  for(k in 1:2){
    beta.altN[k] ~ dnorm(0, 0.0001)
  }
  
  # priors for the spatial component
  global.mu_lanslevillard ~ dnorm(0, 0.01)
  global.tau_lanslevillard ~ dgamma(0.001, 0.001)
  for(i in 1:35){
    mu_lanslevillard[i] ~ dnorm(global.mu_lanslevillard, global.tau_lanslevillard)
  }     
  
  global.mu_beaufort ~ dnorm(0, 0.01)
  global.tau_beaufort ~ dgamma(0.001, 0.001)
  for(i in 1:35){
    mu_beaufort[i] ~ dnorm(global.mu_beaufort, global.tau_beaufort)
  }                 
  
  # Spatial component
  lambda_spat ~ dgamma(1, 0.1)
  sigmasq <- 1/sigmasq.inv
  sigmasq.inv ~ dgamma(2, 1)
  
  ## Build design matrix for N and ps
  for (i in 1:n){
    ## Area effect on detection probability
    logit(pa[i]) <- inprod(beta.secteur.derangementP[], Xsecteur.derangement[i,])
    logit(pb[i]) <- inprod(beta.secteur.derangementP[], Xsecteur.derangement[i,])+delta
    
    ## try for dependant observers
    piMat[i, 1] <- ifelse(prim_ran[i]==1, pa[i],
                          pb[i]
    ) # observed by A or B
    
    piMat[i, 2] <- ifelse(sec_ran[i]==1, pa[i] * (1-pb[i]),
                          pb[i] * (1-pa[i])
    ) # observed by B and missed by A or observed by A missed by B
    
    ## Constrains on population size
    log(N[i]) <- beta.altN[secteur[i]]*altitude[i] + CT[site[i]] 
    
    ## Likelihood for the double-observer survey
    for(j in 1:2){
      C[i, j]      ~ dpois(lambda[i, j])
      lambda[i, j] <- N[i] * piMat[i, j]
      fit[i, j]    <- exp(N[i] * piMat[i, j])
    }
    
  }
  
  # Spatial component Les beaufort
  CT[1:35] ~ dmnorm.vcov(mu_beaufort[1:35], D.covar_beaufort[1:35, 1:35])
  ## turning covariances into precisions
  for(i in 1:35){
    D.covar_beaufort[i, i] <- sigmasq # diagonale de la matrice de covariance
    for(j in 1:(i - 1)){
      # constrain covariance matrix by the distance matrix
      D.covar_beaufort[i, j] <- sigmasq * exp(-(lambda_spat * D_beaufort[i, j]))
      D.covar_beaufort[j, i] <- D.covar_beaufort[i, j] # symmetrical matrix
    }
  }
  
  # Spatial component lanslevillard
  CT[36:70] ~ dmnorm.vcov(mu_lanslevillard[1:35], D.covar_lanslevillard[1:35, 1:35])
  ## turning covariances into precisions
  for(i in 1:35){
    D.covar_lanslevillard[i, i] <- sigmasq # diagonale de la matrice de covariance
    for(j in 1:(i - 1)){
      # constrain covariance matrix by the distance matrix
      D.covar_lanslevillard[i, j] <- sigmasq * exp(-(lambda_spat * D_lanslevillard[i, j]))
      D.covar_lanslevillard[j, i] <- D.covar_lanslevillard[i, j] # symmetrical matrix
    }
  }
}

# Model 11 : ~ secteur+derangement ~ altitude*secteur
mderangementP.secteur.altitudeN.ssCT.jags <- function(){
  ## Set priors for parameters to be estimated
  delta ~ dnorm(0, 0.001)
  
  # Fixed parameters for detection
  for(k in 1:npar.secteur.derangement){
    beta.secteur.derangementP[k] ~ dnorm(0, 0.0001)
  }
  
  # Fixed parameters for abundance
  for(k in 1:2){
    beta.altN[k] ~ dnorm(0, 0.0001)
  }
  
  ## Build design matrix for N and ps
  for (i in 1:n){
    ## Area effect on detection probability
    logit(pa[i]) <- inprod(beta.secteur.derangementP[], Xsecteur.derangement[i,])
    logit(pb[i]) <- inprod(beta.secteur.derangementP[], Xsecteur.derangement[i,])+delta
    
    ## try for dependant observers
    piMat[i, 1] <- ifelse(prim_ran[i]==1, pa[i],
                          pb[i]
    ) # observed by A or B
    
    piMat[i, 2] <- ifelse(sec_ran[i]==1, pa[i] * (1-pb[i]),
                          pb[i] * (1-pa[i])
    ) # observed by B and missed by A or observed by A missed by B
    
    ## Constrains on population size
    log(N[i]) <- beta.altN[secteur[i]]*altitude[i] 
    
    ## Likelihood for the double-observer survey
    for(j in 1:2){
      C[i, j]      ~ dpois(lambda[i, j])
      lambda[i, j] <- N[i] * piMat[i, j]
      fit[i, j]    <- exp(N[i] * piMat[i, j])
    }
    
  }
}

# Priors
npl_marmottes <- model.matrix(data = site.covs_tot, ~ secteur)
npl_marmottes_secteur_orientation <- model.matrix(data = site.covs_tot, ~ secteur+orientation)
npl_marmottes_secteur_derangement <- model.matrix(data = site.covs_tot, ~ secteur+derangement_new)
npl_marmottes_secteur_meteo <- model.matrix(data = site.covs_tot, ~ secteur+meteo)
npl_marmottes_secteur_session <- model.matrix(data = site.covs_tot, ~ secteur+session)
npl_marmottes_habitat   <- model.matrix(data = site.covs_tot, ~ typologie_habitat_new)
npl_marmottes_secteur_derangement_meteo <- model.matrix(data = site.covs_tot, ~ secteur+derangement_new+meteo)
stops_marmottes <- length(observations1[,1])

sectf2 <- site.covs_tot$secteur
sectf2 <- ifelse(sectf2=="Beaufort", 1, 2)

#prim_ran <- sample(c(1:2), dim(observations1)[1], replace=T)
prim_ran <- ifelse(observations1$obs_prim=="Manolo"|observations1$obs_prim=="Léane", 1, 2)
sec_ran <- ifelse(prim_ran==1, 2, 1)

inits.mnull1 <- function(){
  list(
    mu.pa = rnorm(1),
    delta=rnorm(1),
    mu.N = rnorm(1)
  )
}

inits.mnull2 <- function(){
  list(
    lambda_spat = 0.05,
    sigmasq.inv =  1/0.3,
    global.tau_lanslevillard = 0.05,
    global.tau_beaufort = 0.05
  )
}

inits.msecteurP <- function(){
  list(
    beta.secteurP = rnorm(dim(npl_marmottes)[2]),
    lambda_spat = 0.05,
    sigmasq.inv =  1/0.3,
    global.tau_lanslevillard = 0.05,
    global.tau_beaufort = 0.05
  )
}

inits.mderangementP <- function(){
  list(
    beta.secteur.derangementP = rnorm(dim(npl_marmottes_secteur_derangement)[2]),
    lambda_spat = 0.05,
    sigmasq.inv =  1/0.3,
    global.tau_lanslevillard = 0.05,
    global.tau_beaufort = 0.05
  )
}

inits.morientationP <- function(){
  list(
    beta.secteur.orientationP = rnorm(dim(npl_marmottes_secteur_orientation)[2]),
    lambda_spat = 0.05,
    sigmasq.inv =  1/0.3,
    global.tau_lanslevillard = 0.05,
    global.tau_beaufort = 0.05
  )
}

inits.mmeteoP <- function(){
  list(
    beta.secteur.meteoP = rnorm(dim(npl_marmottes_secteur_meteo)[2]),
    lambda_spat = 0.05,
    sigmasq.inv =  1/0.3,
    global.tau_lanslevillard = 0.05,
    global.tau_beaufort = 0.05
  )
}

inits.msessionP <- function(){
  list(
    beta.secteur.sessionP = rnorm(dim(npl_marmottes_secteur_session)[2]),
    lambda_spat = 0.05,
    sigmasq.inv =  1/0.3,
    global.tau_lanslevillard = 0.05,
    global.tau_beaufort = 0.05
  )
}

inits.mderangement.meteoP <- function(){
  list(
    beta.secteur.derangement.meteoP = rnorm(dim(npl_marmottes_secteur_derangement_meteo)[2]),
    lambda_spat = 0.05,
    sigmasq.inv =  1/0.3,
    global.tau_lanslevillard = 0.05,
    global.tau_beaufort = 0.05
  )
}

inits.maltitudeN <- function(){
  list(
    beta.secteurP = rnorm(dim(npl_marmottes)[2]),
    lambda_spat = 0.05,
    sigmasq.inv =  1/0.3,
    global.tau_lanslevillard = 0.05,
    global.tau_beaufort = 0.05
  )
}

inits.mhabitatN <- function(){
  list(
    beta.secteurP = rnorm(dim(npl_marmottes)[2]),
    lambda_spat = 0.05,
    sigmasq.inv =  1/0.3,
    global.tau_lanslevillard = 0.05,
    global.tau_beaufort = 0.05
  )
}

inits.mderangementP.altitudeN <- function(){
  list(
    beta.secteur.derangementP = rnorm(dim(npl_marmottes_secteur_derangement)[2]),
    lambda_spat = 0.05,
    sigmasq.inv =  1/0.3,
    global.tau_lanslevillard = 0.05,
    global.tau_beaufort = 0.05
  )
}

inits.mderangementP.habitatN <- function(){
  list(
    beta.secteur.derangementP = rnorm(dim(npl_marmottes_secteur_derangement)[2]),
    lambda_spat = 0.05,
    sigmasq.inv =  1/0.3,
    global.tau_lanslevillard = 0.05,
    global.tau_beaufort = 0.05
  )
}

inits.mderangementP.secteur.altitudeN <- function(){
  list(
    beta.secteur.derangementP = rnorm(dim(npl_marmottes_secteur_derangement)[2]),
    lambda_spat = 0.05,
    sigmasq.inv =  1/0.3,
    global.tau_lanslevillard = 0.05,
    global.tau_beaufort = 0.05
  )
}

inits.mderangementP.secteur.altitudeN.ssCT <- function(){
  list(
    beta.secteur.derangementP = rnorm(dim(npl_marmottes_secteur_derangement)[2])
  )
}

# Parameters to save
params.mnull1 <- c("lambda", "N", "pa", "pb", "sigma.N",
                   "mu.pa", "rand.N", "mu.N")

params.mnull2 <- c("lambda", "N", "pa", "pb", "lambda_spat", "sigmasq",
                   "D.covar_lanslevillard", "D.covar_beaufort", "CT", "mu_lanslevillard",
                   "mu_beaufort")

params.msecteurP <- c("lambda", "N", "pa", "pb", "lambda_spat", "sigmasq",
              "D.covar_lanslevillard", "D.covar_beaufort", "CT", "mu_lanslevillard",
              "mu_beaufort")

params.morientationP <- c("lambda", "N", "pa", "pb", "lambda_spat", "sigmasq",
                      "D.covar_lanslevillard", "D.covar_beaufort", "CT", "mu_lanslevillard",
                      "mu_beaufort")

params.mderangementP <- c("lambda", "N", "pa", "pb", "lambda_spat", "sigmasq",
                      "D.covar_lanslevillard", "D.covar_beaufort", "CT", "mu_lanslevillard",
                      "mu_beaufort")

params.mmeteoP <- c("lambda", "N", "pa", "pb", "lambda_spat", "sigmasq",
                          "D.covar_lanslevillard", "D.covar_beaufort", "CT", "mu_lanslevillard",
                          "mu_beaufort")

params.msessionP <- c("lambda", "N", "pa", "pb", "lambda_spat", "sigmasq",
                          "D.covar_lanslevillard", "D.covar_beaufort", "CT", "mu_lanslevillard",
                          "mu_beaufort")

params.mderangement.meteoP <- c("lambda", "N", "pa", "pb", "lambda_spat", "sigmasq",
                          "D.covar_lanslevillard", "D.covar_beaufort", "CT", "mu_lanslevillard",
                          "mu_beaufort")

params.maltitudeN <- c("lambda", "N", "pa", "pb", "lambda_spat", "sigmasq",
               "D.covar_lanslevillard", "D.covar_beaufort", "CT", "mu_lanslevillard",
               "mu_beaufort", "beta.altN")

params.mhabitatN <- c("lambda", "N", "pa", "pb", "lambda_spat", "sigmasq",
                      "D.covar_lanslevillard", "D.covar_beaufort", "CT", "mu_lanslevillard",
                      "mu_beaufort", "beta.habitatN")

params.mderangementP.altitudeN <- c("lambda", "N", "pa", "pb", "lambda_spat", "sigmasq",
                          "D.covar_lanslevillard", "D.covar_beaufort", "CT", "mu_lanslevillard",
                          "mu_beaufort", "beta.altN")

params.mderangementP.habitatN <- c("lambda", "N", "pa", "pb", "lambda_spat", "sigmasq",
                                    "D.covar_lanslevillard", "D.covar_beaufort", "CT", "mu_lanslevillard",
                                    "mu_beaufort", "beta.habitatN")

params.mderangementP.secteur.altitudeN <- c("lambda", "N", "pa", "pb", "lambda_spat", "sigmasq",
                                    "D.covar_lanslevillard", "D.covar_beaufort", "CT", "mu_lanslevillard",
                                    "mu_beaufort", "beta.altN")

params.mderangementP.secteur.altitudeN.ssCT <- c("lambda", "N", "pa", "pb", "beta.altN")

# Data
jags.data.mnull1 <- list(
  C = y_tot,
  n = stops_marmottes,
  prim_ran=prim_ran,
  sec_ran=sec_ran,
  #altitude = site.covs_tot$saltitude[1:544,],
  nsite = length(unique(site.covs_tot$site)),                  
  site = site.covs_tot$site
)

jags.data.mnull2 <- list(
  C = y_tot,
  n = stops_marmottes,
  prim_ran=prim_ran,
  sec_ran=sec_ran,
  site = site.covs_tot$site,
  D_beaufort = distance_mat_geo_beaufort / 1000,
  D_lanslevillard = distance_mat_geo_lanslevillard / 1000
)

jags.data.msecteurP <- list(
  C = y_tot,
  n = stops_marmottes,
  Xsecteur = npl_marmottes,
  npar.secteur = dim(npl_marmottes)[2],
  prim_ran=prim_ran,
  sec_ran=sec_ran,
  site = site.covs_tot$site,
  D_beaufort = distance_mat_geo_beaufort / 1000,
  D_lanslevillard = distance_mat_geo_lanslevillard / 1000
)

jags.data.mderangementP <- list(
  C = y_tot,
  n = stops_marmottes,
  Xsecteur.derangement = npl_marmottes_secteur_derangement,
  npar.secteur.derangement = dim(npl_marmottes_secteur_derangement)[2],
  prim_ran=prim_ran,
  sec_ran=sec_ran,
  site = site.covs_tot$site,
  D_beaufort = distance_mat_geo_beaufort / 1000,
  D_lanslevillard = distance_mat_geo_lanslevillard / 1000
)

jags.data.morientationP <- list(
  C = y_tot,
  n = stops_marmottes,
  Xsecteur.orientation = npl_marmottes_secteur_orientation,
  npar.secteur.orientation = dim(npl_marmottes_secteur_orientation)[2],
  prim_ran=prim_ran,
  sec_ran=sec_ran,
  site = site.covs_tot$site,
  D_beaufort = distance_mat_geo_beaufort / 1000,
  D_lanslevillard = distance_mat_geo_lanslevillard / 1000
)

jags.data.mmeteoP <- list(
  C = y_tot,
  n = stops_marmottes,
  Xsecteur.meteo = npl_marmottes_secteur_meteo,
  npar.secteur.meteo = dim(npl_marmottes_secteur_meteo)[2],
  prim_ran=prim_ran,
  sec_ran=sec_ran,
  site = site.covs_tot$site,
  D_beaufort = distance_mat_geo_beaufort / 1000,
  D_lanslevillard = distance_mat_geo_lanslevillard / 1000
)

jags.data.msessionP <- list(
  C = y_tot,
  n = stops_marmottes,
  Xsecteur.session = npl_marmottes_secteur_session,
  npar.secteur.session = dim(npl_marmottes_secteur_session)[2],
  prim_ran=prim_ran,
  sec_ran=sec_ran,
  site = site.covs_tot$site,
  D_beaufort = distance_mat_geo_beaufort / 1000,
  D_lanslevillard = distance_mat_geo_lanslevillard / 1000
)

jags.data.mderangement.meteoP <- list(
  C = y_tot,
  n = stops_marmottes,
  Xsecteur.derangement.meteo = npl_marmottes_secteur_derangement_meteo,
  npar.secteur.derangement.meteo = dim(npl_marmottes_secteur_derangement_meteo)[2],
  prim_ran=prim_ran,
  sec_ran=sec_ran,
  site = site.covs_tot$site,
  D_beaufort = distance_mat_geo_beaufort / 1000,
  D_lanslevillard = distance_mat_geo_lanslevillard / 1000
)

jags.data.maltitudeN <- list(
  C = y_tot,
  n = stops_marmottes,
  Xsecteur = npl_marmottes,
  npar.secteur = dim(npl_marmottes)[2],
  prim_ran=prim_ran,
  sec_ran=sec_ran,
  site = site.covs_tot$site,
  altitude = site.covs_tot$altitude,
  D_beaufort = distance_mat_geo_beaufort / 1000,
  D_lanslevillard = distance_mat_geo_lanslevillard / 1000
)

jags.data.mhabitatN <- list(
  C = y_tot,
  n = stops_marmottes,
  Xsecteur = npl_marmottes,
  Xhabitat = npl_marmottes_habitat,
  npar.secteur = dim(npl_marmottes)[2],
  npar.habitat = dim(npl_marmottes_habitat)[2],
  prim_ran=prim_ran,
  sec_ran=sec_ran,
  site = site.covs_tot$site,
  D_beaufort = distance_mat_geo_beaufort / 1000,
  D_lanslevillard = distance_mat_geo_lanslevillard / 1000
)

jags.data.mderangementP.altitudeN <- list(
  C = y_tot,
  n = stops_marmottes,
  Xsecteur.derangement = npl_marmottes_secteur_derangement,
  npar.secteur.derangement = dim(npl_marmottes_secteur_derangement)[2],
  prim_ran=prim_ran,
  sec_ran=sec_ran,
  site = site.covs_tot$site,
  altitude = site.covs_tot$altitude,
  D_beaufort = distance_mat_geo_beaufort / 1000,
  D_lanslevillard = distance_mat_geo_lanslevillard / 1000
)

jags.data.mderangementP.habitatN <- list(
  C = y_tot,
  n = stops_marmottes,
  Xsecteur.derangement = npl_marmottes_secteur_derangement,
  npar.secteur.derangement = dim(npl_marmottes_secteur_derangement)[2],
  Xhabitat = npl_marmottes_habitat,
  npar.habitat = dim(npl_marmottes_habitat)[2],
  prim_ran=prim_ran,
  sec_ran=sec_ran,
  site = site.covs_tot$site,
  D_beaufort = distance_mat_geo_beaufort / 1000,
  D_lanslevillard = distance_mat_geo_lanslevillard / 1000
)

jags.data.mderangementP.secteur.altitudeN <- list(
  C = y_tot,
  n = stops_marmottes,
  Xsecteur.derangement = npl_marmottes_secteur_derangement,
  npar.secteur.derangement = dim(npl_marmottes_secteur_derangement)[2],
  prim_ran=prim_ran,
  sec_ran=sec_ran,
  site = site.covs_tot$site,
  secteur = sectf2,
  altitude = site.covs_tot$altitude,
  D_beaufort = distance_mat_geo_beaufort / 1000,
  D_lanslevillard = distance_mat_geo_lanslevillard / 1000
)

jags.data.mderangementP.secteur.altitudeN.ssCT <- list(
  C = y_tot,
  n = stops_marmottes,
  Xsecteur.derangement = npl_marmottes_secteur_derangement,
  npar.secteur.derangement = dim(npl_marmottes_secteur_derangement)[2],
  prim_ran=prim_ran,
  sec_ran=sec_ran,
  site = site.covs_tot$site,
  secteur = sectf2,
  altitude = site.covs_tot$altitude
)

# JAGS
fit.mnull1.jags <- jags(data = jags.data.mnull1, 
                        inits = inits.mnull1, 
                        parameters.to.save = params.mnull1, 
                        model.file = mnull1.jags,
                        n.chains = 4, 
                        n.iter = 500, 
                        n.burnin = 100, 
                        n.thin = 1)
save(fit.mnull1.jags, file="Output/fit.mnull1.jags.Rdata")


fit.mnull2.jags <- jags(data = jags.data.mnull2, 
                    inits = inits.mnull2, 
                    parameters.to.save = params.mnull2, 
                    model.file = mnull2.jags,
                    n.chains = 4, 
                    n.iter = 500, 
                    n.burnin = 100, 
                    n.thin = 1)
save(fit.mnull2.jags, file="Output/fit.mnull2.jags.Rdata")


fit.msecteurP.jags <- jags(data = jags.data.msecteurP, 
                inits = inits.msecteurP, 
                parameters.to.save = params.msecteurP, 
                model.file = msecteurP.jags,
                n.chains = 4, 
                n.iter = 500, 
                n.burnin = 100, 
                n.thin = 1)
save(fit.msecteurP.jags, file="Output/fit.msecteurP.jags.Rdata")


fit.mderangementP.jags <- jags(data = jags.data.mderangementP, 
                    inits = inits.mderangementP, 
                    parameters.to.save = params.mderangementP, 
                    model.file = mderangementP.jags,
                    n.chains = 4, 
                    n.iter = 500, 
                    n.burnin = 100, 
                    n.thin = 1)
save(fit.mderangementP.jags, file="Output/fit.mderangementP.jags.Rdata")


fit.morientationP.jags <- jags(data = jags.data.morientationP, 
                              inits = inits.morientationP, 
                              parameters.to.save = params.morientationP, 
                              model.file = morientationP.jags,
                              n.chains = 4, 
                              n.iter = 500, 
                              n.burnin = 100, 
                              n.thin = 1)
save(fit.morientationP.jags, file="Output/fit.morientationP.jags.Rdata")


fit.mmeteoP.jags <- jags(data = jags.data.mmeteoP, 
                               inits = inits.mmeteoP, 
                               parameters.to.save = params.mmeteoP, 
                               model.file = mmeteoP.jags,
                               n.chains = 4, 
                               n.iter = 500, 
                               n.burnin = 100, 
                               n.thin = 1)
save(fit.mmeteoP.jags, file="Output/fit.mmeteoP.jags.Rdata")


fit.msessionP.jags <- jags(data = jags.data.msessionP, 
                               inits = inits.msessionP, 
                               parameters.to.save = params.msessionP, 
                               model.file = msessionP.jags,
                               n.chains = 4, 
                               n.iter = 500, 
                               n.burnin = 100, 
                               n.thin = 1)
save(fit.msessionP.jags, file="Output/fit.msessionP.jags.Rdata")

fit.mderangement.meteoP.jags <- jags(data = jags.data.mderangement.meteoP, 
                               inits = inits.mderangement.meteoP, 
                               parameters.to.save = params.mderangement.meteoP, 
                               model.file = mderangement.meteoP.jags,
                               n.chains = 4, 
                               n.iter = 500, 
                               n.burnin = 100, 
                               n.thin = 1)
save(fit.mderangement.meteoP.jags, file="Output/fit.mderangement.meteoP.jags.Rdata")

fit.maltitudeN.jags <- jags(data = jags.data.maltitudeN, 
                              inits = inits.maltitudeN, 
                              parameters.to.save = params.maltitudeN, 
                              model.file = maltitudeN.jags,
                              n.chains = 4, 
                              n.iter = 500, 
                              n.burnin = 100, 
                              n.thin = 1)
save(fit.maltitudeN.jags, file="Output/fit.maltitudeN.jags.Rdata")


fit.mhabitatN.jags <- jags(data = jags.data.mhabitatN, 
                           inits = inits.mhabitatN, 
                           parameters.to.save = params.mhabitatN, 
                           model.file = mhabitatN.jags,
                           n.chains = 4, 
                           n.iter = 500, 
                           n.burnin = 100, 
                           n.thin = 1)
save(fit.mhabitatN.jags, file="Output/fit.mhabitatN.jags.Rdata")

#caterplot(fit.mhabitatN.jags, "beta.habitatN", labels.loc="axis")
#fit.mhabitatN.mcmc <- as.mcmc(fit.mhabitatN.jags)
#plot(fit.msecteurP.mcmc)

fit.mderangementP.altitudeN.jags <- jags(data = jags.data.mderangementP.altitudeN, 
                               inits = inits.mderangementP.altitudeN, 
                               parameters.to.save = params.mderangementP.altitudeN, 
                               model.file = mderangementP.altitudeN.jags,
                               n.chains = 4, 
                               n.iter = 500, 
                               n.burnin = 100, 
                               n.thin = 1)

save(fit.mderangementP.altitudeN.jags, file="Output/fit.mderangementP.altitudeN.jags.Rdata")

fit.mderangementP.habitatN.jags <- jags(data = jags.data.mderangementP.habitatN, 
                                         inits = inits.mderangementP.habitatN, 
                                         parameters.to.save = params.mderangementP.habitatN, 
                                         model.file = mderangementP.habitatN.jags,
                                         n.chains = 4, 
                                         n.iter = 500, 
                                         n.burnin = 100, 
                                         n.thin = 1)
save(fit.mderangementP.habitatN.jags, file="Output/fit.mderangementP.habitatN.jags.Rdata")

fit.mderangementP.secteur.altitudeN.jags <- jags(data = jags.data.mderangementP.secteur.altitudeN, 
                                         inits = inits.mderangementP.secteur.altitudeN, 
                                         parameters.to.save = params.mderangementP.secteur.altitudeN, 
                                         model.file = mderangementP.secteur.altitudeN.jags,
                                         n.chains = 4, 
                                         n.iter = 500, 
                                         n.burnin = 100, 
                                         n.thin = 1)
save(fit.mderangementP.secteur.altitudeN.jags, file="Output/fit.mderangementP.secteur.altitudeN.jags.Rdata")

fit.mderangementP.secteur.altitudeN.ssCT.jags <- jags(data = jags.data.mderangementP.secteur.altitudeN.ssCT, 
                                                 inits = inits.mderangementP.secteur.altitudeN.ssCT, 
                                                 parameters.to.save = params.mderangementP.secteur.altitudeN.ssCT, 
                                                 model.file = mderangementP.secteur.altitudeN.ssCT.jags,
                                                 n.chains = 4, 
                                                 n.iter = 500, 
                                                 n.burnin = 100, 
                                                 n.thin = 1)
save(fit.mderangementP.secteur.altitudeN.jags, file="Output/fit.mderangementP.secteur.altitudeN.jags.Rdata")

# WAIC
samples.mnull1 <- jags.samples(fit.mnull1.jags$model, 
                               c("WAIC","deviance","N"), 
                               type = "mean", 
                               n.iter = 500,
                               n.burnin = 100,
                               n.thin = 1)
save(samples.mnull1, file="Output/samples.mnull1.Rdata")


samples.mnull2 <- jags.samples(fit.mnull2.jags$model, 
                           c("WAIC","deviance","N"), 
                           type = "mean", 
                           n.iter = 500,
                           n.burnin = 100,
                           n.thin = 1)
save(samples.mnull2, file="Output/samples.mnull2.Rdata")


samples.msecteurP <- jags.samples(fit.msecteurP.jags$model, 
                           c("WAIC","deviance","N"), 
                           type = "mean", 
                           n.iter = 500,
                           n.burnin = 100,
                           n.thin = 1)
save(samples.msecteurP, file="Output/samples.msecteurP.Rdata")


samples.mderangementP <- jags.samples(fit.mderangementP.jags$model, 
                           c("WAIC","deviance","N"), 
                           type = "mean", 
                           n.iter = 500,
                           n.burnin = 100,
                           n.thin = 1)
save(samples.mderangementP, file="Output/samples.mderangementP.Rdata")


samples.morientationP <- jags.samples(fit.morientationP.jags$model, 
                           c("WAIC","deviance","N"), 
                           type = "mean", 
                           n.iter = 500,
                           n.burnin = 100,
                           n.thin = 1)
save(samples.morientationP, file="Output/samples.morientationP.Rdata")


samples.mmeteoP <- jags.samples(fit.mmeteoP.jags$model, 
                                      c("WAIC","deviance","N"), 
                                      type = "mean", 
                                      n.iter = 500,
                                      n.burnin = 100,
                                      n.thin = 1)
save(samples.mmeteoP, file="Output/samples.mmeteoP.Rdata")


samples.msessionP <- jags.samples(fit.msessionP.jags$model, 
                                      c("WAIC","deviance","N"), 
                                      type = "mean", 
                                      n.iter = 500,
                                      n.burnin = 100,
                                      n.thin = 1)
save(samples.msessionP, file="Output/samples.msessionP.Rdata")

samples.mderangement.meteoP <- jags.samples(fit.mderangement.meteoP.jags$model, 
                                      c("WAIC","deviance","N"), 
                                      type = "mean", 
                                      n.iter = 500,
                                      n.burnin = 100,
                                      n.thin = 1)
save(samples.mderangement.meteoP, file="Output/samples.mderangement.meteoP.Rdata")


samples.maltitudeN <- jags.samples(fit.maltitudeN.jags$model, 
                                     c("WAIC","deviance","N"), 
                                     type = "mean", 
                                     n.iter = 500,
                                     n.burnin = 100,
                                     n.thin = 1)
save(samples.maltitudeN, file="Output/samples.maltitudeN.Rdata")


samples.mhabitatN <- jags.samples(fit.mhabitatN.jags$model, 
                                  c("WAIC","deviance","N"), 
                                  type = "mean", 
                                  n.iter = 500,
                                  n.burnin = 100,
                                  n.thin = 1)
save(samples.mhabitatN, file="Output/samples.mhabitatN.Rdata")

samples.mderangementP.altitudeN <- jags.samples(fit.mderangementP.altitudeN.jags$model, 
                                      c("WAIC","deviance","N"), 
                                      type = "mean", 
                                      n.iter = 500,
                                      n.burnin = 100,
                                      n.thin = 1)
save(samples.mderangementP.altitudeN, file="Output/samples.mderangementP.altitudeN.Rdata")

samples.mderangementP.habitatN <- jags.samples(fit.mderangementP.habitatN.jags$model, 
                                                c("WAIC","deviance","N"), 
                                                type = "mean", 
                                                n.iter = 500,
                                                n.burnin = 100,
                                                n.thin = 1)
save(samples.mderangementP.habitatN, file="Output/samples.mderangementP.habitatN.Rdata")

samples.mderangementP.secteur.altitudeN <- jags.samples(fit.mderangementP.secteur.altitudeN.jags$model, 
                                                c("WAIC","deviance","N"), 
                                                type = "mean", 
                                                n.iter = 500,
                                                n.burnin = 100,
                                                n.thin = 1)
save(samples.mderangementP.secteur.altitudeN, file="Output/samples.mderangementP.secteur.altitudeN.Rdata")

samples.mderangementP.secteur.altitudeN.ssCT <- jags.samples(fit.mderangementP.secteur.altitudeN.ssCT.jags$model, 
                                                        c("WAIC","deviance","N"), 
                                                        type = "mean", 
                                                        n.iter = 500,
                                                        n.burnin = 100,
                                                        n.thin = 1)
save(samples.mderangementP.secteur.altitudeN.ssCT, file="Output/samples.mderangementP.secteur.altitudeN.ssCT.Rdata")


# Extract lambda spat for each model fit
load("Output/fit.mnull1.jags.Rdata")
load("Output/fit.mnull2.jags.Rdata")
load("Output/fit.msecteurP.jags.Rdata")
load("Output/fit.mderangementP.jags.Rdata")
load("Output/fit.morientationP.jags.Rdata")
load("Output/fit.mmeteoP.jags.Rdata")
load("Output/fit.msessionP.jags.Rdata")
load("Output/fit.mderangement.meteoP.jags.Rdata")
load("Output/fit.maltitudeN.jags.Rdata")
load("Output/fit.mhabitatN.jags.Rdata")
load("Output/fit.mderangementP.altitudeN.jags.Rdata")
load("Output/fit.mderangementP.habitatN.jags.Rdata")

fit.mnull1.jags$BUGSoutput$summary[c("lambda_spat","sigmasq"),]
fit.mnull2.jags$BUGSoutput$summary[c("lambda_spat","sigmasq"),]
fit.msecteurP.jags$BUGSoutput$summary[c("lambda_spat","sigmasq"),]
fit.mderangementP.jags$BUGSoutput$summary[c("lambda_spat","sigmasq"),]
fit.morientationP.jags$BUGSoutput$summary[c("lambda_spat","sigmasq"),]
fit.mmeteoP.jags$BUGSoutput$summary[c("lambda_spat","sigmasq"),]
fit.msessionP.jags$BUGSoutput$summary[c("lambda_spat","sigmasq"),]
fit.mderangement.meteoP.jags$BUGSoutput$summary[c("lambda_spat","sigmasq"),]
fit.maltitudeN.jags$BUGSoutput$summary[c("lambda_spat","sigmasq"),]
fit.mhabitatN.jags$BUGSoutput$summary[c("lambda_spat","sigmasq"),]
fit.mderangementP.altitudeN.jags$BUGSoutput$summary[c("lambda_spat","sigmasq"),]
fit.mderangementP.habitatN.jags$BUGSoutput$summary[c("lambda_spat","sigmasq"),]

# Extract Rhat for each model fit
mean(fit.mderangementP.habitatN.jags$BUGSoutput$summary[,"Rhat"])
mean(fit.mnull1.jags$BUGSoutput$summary[,"Rhat"])
mean(fit.mnull2.jags$BUGSoutput$summary[,"Rhat"])
mean(fit.msecteurP.jags$BUGSoutput$summary[,"Rhat"])
mean(fit.mderangementP.jags$BUGSoutput$summary[,"Rhat"])
mean(fit.morientationP.jags$BUGSoutput$summary[,"Rhat"])
mean(fit.mmeteoP.jags$BUGSoutput$summary[,"Rhat"])
mean(fit.msessionP.jags$BUGSoutput$summary[,"Rhat"])
mean(fit.mderangement.meteoP.jags$BUGSoutput$summary[,"Rhat"])
mean(fit.maltitudeN.jags$BUGSoutput$summary[,"Rhat"])
mean(fit.mhabitatN.jags$BUGSoutput$summary[,"Rhat"])
mean(fit.mderangementP.altitudeN.jags$BUGSoutput$summary[,"Rhat"])
mean(fit.mderangementP.habitatN.jags$BUGSoutput$summary[,"Rhat"])

