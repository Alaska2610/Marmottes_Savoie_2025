## Set priors for parameters to be estimated
delta ~ dnorm(0, 0.001)

# # Fixed parameters for detection
# for(k in 1:npar.secteur.derangement){
#     beta.secteur.derangementP[k] ~ dnorm(0, 0.0001)
# }

# Fixed parameters for detection
for(k in 1:npar.secteur){
  beta.secteurP[k] ~ dnorm(0, 0.0001)
}

# Fixed parameters for abundance
for(k in 1:2){
  beta.altN[k] ~ dnorm(0, 0.0001)
}

# # priors for the spatial component
# global.mu_beaufort ~ dnorm(0, 0.01)
# global.tau_beaufort ~ dgamma(0.001, 0.001)
# for(i in 1:35){
#     mu_beaufort[i] ~ dnorm(global.mu_beaufort, global.tau_beaufort)
# }     
# 
# global.mu_lanslevillard ~ dnorm(0, 0.01)
# global.tau_lanslevillard ~ dgamma(0.001, 0.001)
# for(i in 1:35){
#     mu_lanslevillard[i] ~ dnorm(global.mu_lanslevillard, global.tau_lanslevillard)
# }                 

# Spatial component
# lambda_spat ~ dgamma(1, 0.1)
# sigmasq <- 1/sigmasq.inv
# sigmasq.inv ~ dgamma(2, 1)

## Build design matrix for N and ps
for (i in 1:n){
  ## Area effect on detection probability
  # logit(pa[i]) <- inprod(beta.secteur.derangementP[], Xsecteur.derangement[i,])
  # logit(pb[i]) <- inprod(beta.secteur.derangementP[], Xsecteur.derangement[i,])+delta
  
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
  log(N[i]) <- beta.altN[secteur[i]]*altitude[i]
  #log(N[i]) <- beta.altN[secteur[i]]*altitude[i] + CT[site[i]] 
  
  ## Likelihood for the double-observer survey
  for(j in 1:2){
    C[i, j]      ~ dpois(lambda[i, j])
    lambda[i, j] <- N[i] * piMat[i, j]
    fit[i, j]    <- exp(N[i] * piMat[i, j])
  }
}

# # Spatial component Lanslevillard
#   CT[1:35] ~ dmnorm.vcov(mu_lanslevillard[1:35], D.covar_lanslevillard[1:35, 1:35])
#   ## turning covariances into precisions
#   for(i in 1:35){
#     D.covar_lanslevillard[i, i] <- sigmasq # diagonale de la matrice de covariance
#         for(j in 1:(i - 1)){
#             # constrain covariance matrix by the distance matrix
#             D.covar_lanslevillard[i, j] <- sigmasq * exp(-(lambda_spat * D_lanslevillard[i, j]))
#             D.covar_lanslevillard[j, i] <- D.covar_lanslevillard[i, j] # symmetrical matrix
#           }
#     }
# 
#   # Spatial component beaufort
#   CT[36:70] ~ dmnorm.vcov(mu_beaufort[1:35], D.covar_beaufort[1:35, 1:35])
#   ## turning covariances into precisions
#   for(i in 1:35){
#     D.covar_beaufort[i, i] <- sigmasq # diagonale de la matrice de covariance
#         for(j in 1:(i - 1)){
#             # constrain covariance matrix by the distance matrix
#             D.covar_beaufort[i, j] <- sigmasq * exp(-(lambda_spat * D_beaufort[i, j]))
#             D.covar_beaufort[j, i] <- D.covar_beaufort[i, j] # symmetrical matrix
#           }
#     }




jags_marmottes7 <- jags.model(
  model7.spec_marmottes,
  data = list(
    C = y_tot,
    n = stops_marmottes,
    Xsecteur = npl_marmottes_secteur,
    npar.secteur = dim(npl_marmottes_secteur)[2],
    prim_ran=prim_ran,
    sec_ran=sec_ran,
    secteur = sectf2,
    altitude = site.covs_tot$saltitude[1:552,] # altitude centrée-réduite
  ),
  inits = list(
    beta.secteurP = rnorm(dim(npl_marmottes_secteur)[2]),
    lambda_spat = 0.01,
    sigmasq.inv =  1/0.3,
    global.tau_beaufort = 0.05,
    global.tau_lanslevillard = 0.05
  ),
  n.chains = 4,
  n.adapt = 500
)


"beta.secteurP",
"beta.altN",
"N",
"lambda"


plot(out_marmottes7_coda[,c(#"global.mu_beaufort","global.mu_lanslevillard",
  "beta.altN[1]","beta.altN[2]",
  "beta.secteurP[1]",
  "beta.secteurP[2]")],
  #"lambda_spat","sigmasq")],
  trace = FALSE)


plot(out_marmottes7_coda[,c(#"global.mu_beaufort","global.mu_lanslevillard",
  "beta.altN[1]","beta.altN[2]",
  "beta.secteurP[1]",
  "beta.secteurP[2]")],
  #"lambda_spat","sigmasq")], 
  density = FALSE)


summary(out_marmottes7_coda[,c("beta.altN[1]","beta.altN[2]",
                               "beta.secteurP[1]",
                               "beta.secteurP[2]")])

gelman.diag(out_marmottes7_coda[,c("beta.altN[1]","beta.altN[2]",
                                   "beta.secteurP[1]",
                                   "beta.secteurP[2]")])