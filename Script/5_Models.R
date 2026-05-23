###################################################
## Comparaison des WAIC des différents modèles
###################################################

library(R2jags)

# https://gist.github.com/oliviergimenez/68ad17910a62635ff6a062f8ec34292f

# Models

# Model 1 
m1.altitude <- function(){
    ## Set priors for parameters to be estimated
    delta ~ dnorm(0, 0.001)
    
    # Fixed parameters for detection
    for(k in 1:npar.secteur){
      beta.secteurP[k] ~ dnorm(0, 0.0001)
    }

    for(k in 1:2){
      beta.altN[k] ~ dnorm(0, 0.0001)
    }
    
    for(k in 1:npar.habitatpsi){
      beta.habitatPsi[k] ~ dnorm(0, 0.0001)
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
      # Ajout composante zero-inflated :
      # https://biometry.github.io/APES//LectureNotes/2016-JAGS/ZeroInflation/ZeroInflation_JAGS.html
      # https://jbds.isdsa.org/public/journals/1/html/v2n2/shao/
      
      ## ZERO-INFLATION
      logit(psi[i]) <- inprod(beta.habitatPsi[], Xhabitatpsi[i, ])
      z[i] ~ dbern(psi[i])
      
      ## ABONDANCE 
      log(N[i]) <- beta.altN[secteur[i]]*altitude[i] + CT[site[i]] 
      
      N_eff[i] <- N[i] * z[i] + 0.00001
      
      ## Likelihood for the double-observer survey
      for(j in 1:2){
        C[i, j]      ~ dpois(lambda[i, j])
        lambda[i, j] <- N_eff[i] * piMat[i, j]
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


m2.altitude <- function(){
  ## Set priors for parameters to be estimated
  delta ~ dnorm(0, 0.001)
  
  # Fixed parameters for detection
  for(k in 1:npar.secteur){
    beta.secteurP[k] ~ dnorm(0, 0.0001)
  }
  
  for(k in 1:2){
    beta.altN[k] ~ dnorm(0, 0.0001)
  }
  
  for(k in 1:2){
    beta.altN2[k] ~ dnorm(0, 0.0001)
  }
  
  for(k in 1:npar.habitatpsi){
    beta.habitatPsi[k] ~ dnorm(0, 0.0001)
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
    # Ajout composante zero-inflated :
    # https://biometry.github.io/APES//LectureNotes/2016-JAGS/ZeroInflation/ZeroInflation_JAGS.html
    # https://jbds.isdsa.org/public/journals/1/html/v2n2/shao/
    
    ## ZERO-INFLATION
    logit(psi[i]) <- inprod(beta.habitatPsi[], Xhabitatpsi[i, ])
    z[i] ~ dbern(psi[i])
    
    ## ABONDANCE 
    log(N[i]) <- beta.altN[secteur[i]] * altitude[i] +
                 beta.altN2[secteur[i]] * altitude[i]^2 +
                 CT[site[i]] 
    
    N_eff[i] <- N[i] * z[i] + 0.00001
    
    ## Likelihood for the double-observer survey
    for(j in 1:2){
      C[i, j]      ~ dpois(lambda[i, j])
      lambda[i, j] <- N_eff[i] * piMat[i, j]
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


# Modèle 2
m3.habitat <- function(){
  ## Set priors for parameters to be estimated
  delta ~ dnorm(0, 0.001)
  
  # Fixed parameters for detection
  for(k in 1:npar.secteur){
    beta.secteurP[k] ~ dnorm(0, 0.0001)
  }
  
  for(k in 1:npar.habitat){
   beta.habitatN[k] ~ dnorm(0, 0.0001)
  }
  
  for(k in 1:npar.habitatpsi){
    beta.habitatPsi[k] ~ dnorm(0, 0.0001)
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
    # Ajout composante zero-inflated :
    # https://biometry.github.io/APES//LectureNotes/2016-JAGS/ZeroInflation/ZeroInflation_JAGS.html
    # https://jbds.isdsa.org/public/journals/1/html/v2n2/shao/
    
    ## ZERO-INFLATION
    logit(psi[i]) <- inprod(beta.habitatPsi[], Xhabitatpsi[i, ])
    z[i] ~ dbern(psi[i])
    
    ## ABONDANCE 
    log(N[i]) <- inprod(beta.habitatN[], Xhabitat[i, ]) + CT[site[i]]

    N_eff[i] <- N[i] * z[i] + 0.00001
    
    ## Likelihood for the double-observer survey
    for(j in 1:2){
      C[i, j]      ~ dpois(lambda[i, j])
      lambda[i, j] <- N_eff[i] * piMat[i, j]
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

# Priors
npl_marmottes <- model.matrix(data = site.covs_tot, ~ secteur)
npl_marmottes_habitat   <- model.matrix(data = site.covs_tot, ~ typologie_habitat_new)
npl_marmottes_habitat_secteur   <- model.matrix(data = site.covs_tot, ~ secteur+typologie_habitat_new)
stops_marmottes <- length(observations1[,1])

prim_ran <- ifelse(observations1$obs_prim=="Manolo"|observations1$obs_prim=="Léane", 1, 2)
sec_ran <- ifelse(prim_ran==1, 2, 1)

sectf2 <- site.covs_tot$secteur
sectf2 <- ifelse(sectf2=="Beaufort", 1, 2)

inits.m1 <-  function(){
  list  (
    beta.secteurP = rnorm(dim(npl_marmottes)[2]),
    lambda_spat = 0.05,
    sigmasq.inv =  1/0.3,
    global.tau_lanslevillard = 0.05,
    global.tau_beaufort = 0.05
  )
}

# Parameters to save
params.m1 <- c("lambda", "N_eff", "pa", "pb", "lambda_spat",
               "sigmasq", "D.covar_beaufort", "D.covar_lanslevillard",
               "CT", "mu_beaufort", "mu_lanslevillard",
               "beta.altN", "beta.secteurP", "beta.habitatPsi")

params.m2 <- c("lambda", "N_eff", "pa", "pb", "lambda_spat",
               "sigmasq", "D.covar_beaufort", "D.covar_lanslevillard",
               "CT", "mu_beaufort", "mu_lanslevillard",
               "beta.altN", "beta.altN2", "beta.secteurP", "beta.habitatPsi")

params.m3 <- c("lambda", "N_eff", "pa", "pb", "lambda_spat",
              "sigmasq", "D.covar_beaufort", "D.covar_lanslevillard",
              "CT", "mu_beaufort", "mu_lanslevillard",
              "beta.habitatN", "beta.secteurP", "beta.habitatPsi")

# Data
jags.data.m1 <- list(
  C = y_tot,
  n = stops_marmottes,
  z = ifelse(rowSums(y_tot) > 0, 1, 0),
  Xsecteur = npl_marmottes,
  Xhabitatpsi = npl_marmottes_habitat,
  npar.secteur = dim(npl_marmottes)[2],
  npar.habitatpsi = dim(npl_marmottes_habitat)[2],
  prim_ran=prim_ran,
  sec_ran=sec_ran,
  site = site.covs_tot$site,
  secteur = sectf2,
  altitude = site.covs_tot$saltitude[1:552,], # altitude centrée-réduite
  D_beaufort = distance_mat_geo_beaufort / 1000,
  D_lanslevillard = distance_mat_geo_lanslevillard / 1000
)

jags.data.m2 <- list(
  C = y_tot,
  n = stops_marmottes,
  z = ifelse(rowSums(y_tot) > 0, 1, 0),
  Xsecteur = npl_marmottes,
  Xhabitatpsi = npl_marmottes_habitat,
  npar.secteur = dim(npl_marmottes)[2],
  npar.habitatpsi = dim(npl_marmottes_habitat)[2],
  prim_ran=prim_ran,
  sec_ran=sec_ran,
  site = site.covs_tot$site,
  secteur = sectf2,
  altitude = site.covs_tot$saltitude[1:552,], # altitude centrée-réduite
  D_beaufort = distance_mat_geo_beaufort / 1000,
  D_lanslevillard = distance_mat_geo_lanslevillard / 1000
)

jags.data.m3 <- list(
  C = y_tot,
  n = stops_marmottes,
  z = ifelse(rowSums(y_tot) > 0, 1, 0),
  Xsecteur = npl_marmottes,
  Xhabitat = npl_marmottes_habitat_secteur,
  Xhabitatpsi = npl_marmottes_habitat,
  npar.secteur = dim(npl_marmottes)[2],
  npar.habitat = dim(npl_marmottes_habitat_secteur)[2],
  npar.habitatpsi = dim(npl_marmottes_habitat)[2],
  prim_ran=prim_ran,
  sec_ran=sec_ran,
  site = site.covs_tot$site,
  secteur = sectf2,
  D_beaufort = distance_mat_geo_beaufort / 1000,
  D_lanslevillard = distance_mat_geo_lanslevillard / 1000
)

# JAGS
fit.m1.jags <- jags(data = jags.data.m1, 
                        inits = inits.m1, 
                        parameters.to.save = params.m1, 
                        model.file = m1.altitude,
                        n.chains = 4, 
                        n.iter = 5000, 
                        n.burnin = 1000, 
                        n.thin = 1)
setwd("/Volumes/Boulot/FDC_Savoie/2025")
save(fit.m1.jags, file="fit.m1.jags.Rdata")

fit.m2.jags <- jags(data = jags.data.m2, 
                    inits = inits.m1, 
                    parameters.to.save = params.m2, 
                    model.file = m2.altitude,
                    n.chains = 4, 
                    n.iter = 5000, 
                    n.burnin = 1000, 
                    n.thin = 1)
setwd("/Volumes/Boulot/FDC_Savoie/2025")
save(fit.m2.jags, file="fit.m2.jags.Rdata")

fit.m3.jags <- jags(data = jags.data.m3, 
                    inits = inits.m1, 
                    parameters.to.save = params.m3, 
                    model.file = m3.habitat,
                    n.chains = 4, 
                    n.iter = 5000, 
                    n.burnin = 1000, 
                    n.thin = 1)
setwd("/Volumes/Boulot/FDC_Savoie/2025")
save(fit.m3.jags, file="fit.m3.jags.Rdata")

# WAIC
samples.m1 <- jags.samples(fit.m1.jags$model, 
                               c("WAIC","deviance","N"), 
                               type = "mean", 
                               n.iter = 5000,
                               n.burnin = 1000,
                               n.thin = 1)
setwd("/Volumes/Boulot/FDC_Savoie/2025")
save(samples.m1, file="samples.m1.Rdata")

samples.m2 <- jags.samples(fit.m2.jags$model, 
                           c("WAIC","deviance","N"), 
                           type = "mean", 
                           n.iter = 5000,
                           n.burnin = 1000,
                           n.thin = 1)
setwd("/Volumes/Boulot/FDC_Savoie/2025")
save(samples.m2, file="samples.m2.Rdata")

samples.m3 <- jags.samples(fit.m3.jags$model, 
                           c("WAIC","deviance","N"), 
                           type = "mean", 
                           n.iter = 5000,
                           n.burnin = 1000,
                           n.thin = 1)
setwd("/Volumes/Boulot/FDC_Savoie/2025")
save(samples.m3, file="samples.m3.Rdata")


# Extract lambda spat for each model fit
setwd("/Volumes/Boulot/FDC_Savoie/2025")
load("fit.m1.jags.Rdata")
load("fit.m2.jags.Rdata")
load("fit.m3.jags.Rdata")

fit.m1.jags$BUGSoutput$summary[c("lambda_spat","sigmasq"),]
fit.m2.jags$BUGSoutput$summary[c("lambda_spat","sigmasq"),]
fit.m3.jags$BUGSoutput$summary[c("lambda_spat","sigmasq"),]

# Extract Rhat for each model fit
mean(fit.m1.jags$BUGSoutput$summary[,"Rhat"])
mean(fit.m2.jags$BUGSoutput$summary[,"Rhat"])
mean(fit.m3.jags$BUGSoutput$summary[,"Rhat"])
