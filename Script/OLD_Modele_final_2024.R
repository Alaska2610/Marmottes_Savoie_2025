#####################################
## Modèle secteur+dérangement sur P
## CT sur N
## Test sur 2025 avec exactement les mêmes variables
#####################################
library(coda)
library(rjags)
library(qgraph)
library(ggmap)
#library(mcmcplots)
library(ggplot2)
library(dplyr)
library(tidyr)
library(unmarked)
library(gmodels)

CB <- TRUE
if(CB == TRUE){setwd("./Scripts")}else
{setwd("G:/FDC_Savoie/Marmottes_Savoie/")}

load("Output/observations1.Rdata")

# Merging some habitats
head(observations1)
table(observations1$habitat)
observations1$typologie_habitat_new <- observations1$habitat
observations1[is.element(observations1$habitat, "Prairie de fauche"),]$typologie_habitat_new <- "alpage_herbace"
observations1[is.element(observations1$habitat, "Alpage"),]$typologie_habitat_new <- "alpage_herbace"
observations1[is.element(observations1$habitat, "Foret claire"),]$typologie_habitat_new <- "pre_bois"
observations1[is.element(observations1$habitat, "Lande"),]$typologie_habitat_new <- "pre_bois"
observations1[is.element(observations1$habitat, "Eboulis /substrat rocheux"),]$typologie_habitat_new <- "eboulis"
table(observations1$typologie_habitat_new)

# Disturbance 0/1
table(observations1$derangement)
observations1$derangement_new <- observations1$derangement
observations1[is.element(observations1$derangement, c("Agricole","Autre","Touristique")),]$derangement_new <- "Derangement"
table(observations1$derangement_new)

observations1$saltitude <- scale(observations1$altitude)

hist(observations1$nb_obs_prim)
hist(observations1$nb_obs_sec)

## Unmarked Double-Observer sur toutes les obs
y_tot <- matrix(NA, dim(observations1), 2, byrow=T)
y_tot[,1] <- observations1$nb_obs_prim
y_tot[,2] <- observations1$nb_obs_sec-observations1$nb_obs_prim
y_tot[,2] <- ifelse(y_tot[,2]<0,0,y_tot[,2])
site.covs_tot <- data.frame(observations1[,c("secteur","meteo",
                                             "derangement","derangement_new",
                                             "habitat","typologie_habitat_new",
                                             "orientation","latitude","longitude",
                                             "julian_day","session",
                                             "altitude","saltitude","maille")])
observer_m    <- as.matrix(
  observations1[,c("obs_prim","obs_sec")]
)
site.covs_tot$site <- site.covs_tot$maille
table(site.covs_tot$maille, site.covs_tot$secteur)
# 1 à 291 : Beaufort, 308 à 517 : Lanslevillard
# Renommage des mailles pour qu'elles aillent de 1 à 70
order_maille <- unique(site.covs_tot[order(site.covs_tot$maille),]$maille)
for(i in order_maille){
  site.covs_tot[is.element(site.covs_tot$maille, i),]$site <- match(i,order_maille)
}

## Distance matrix between counting points (distance entre mailles)
site.covs_tot$latitude <- as.numeric(as.character(site.covs_tot$latitude))
site.covs_tot$longitude <- as.numeric(as.character(site.covs_tot$longitude))
## Average of the coordinates by site to have a single counting point per site
coord_lanslevillard <- site.covs_tot[is.element(site.covs_tot$secteur, "Lanslevillard"),c("latitude","longitude","site")]
coord_lanslevillard_mean <- aggregate(list(x=coord_lanslevillard$latitude, y=coord_lanslevillard$longitude), by=list(site=coord_lanslevillard$site), mean)
rownames(coord_lanslevillard_mean) <- coord_lanslevillard_mean$site
distance_mat_geo1_lanslevillard <- dist(coord_lanslevillard_mean[,c("x","y")]) # meters
distance_mat_geo_lanslevillard <- as.matrix(distance_mat_geo1_lanslevillard)

coord_beaufort <- site.covs_tot[is.element(site.covs_tot$secteur, "Beaufort"),c("latitude","longitude","site")]
coord_beaufort_mean <- aggregate(list(x=coord_beaufort$latitude, y=coord_beaufort$longitude), by=list(site=coord_beaufort$site), mean)
rownames(coord_beaufort_mean) <- coord_beaufort_mean$site
distance_mat_geo1_beaufort <- dist(coord_beaufort_mean[,c("x","y")]) # meters
distance_mat_geo_beaufort <- as.matrix(distance_mat_geo1_beaufort, labels=T)

# Test modèle sans effet de dérangement sur P et pas d'autocorrélation spatiale
model7.string_marmottes <-"
          model {

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
"

## Define parameters
npl_marmottes_secteur_derangement <- model.matrix(data = site.covs_tot, ~ secteur+derangement_new)
npl_marmottes_secteur <- model.matrix(data = site.covs_tot, ~ secteur)
stops_marmottes <- length(observations1[,1])

prim_ran <- ifelse(observations1$obs_prim=="Manolo"|observations1$obs_prim=="Léane", 1, 2)
sec_ran <- ifelse(prim_ran==1, 2, 1)

sectf2 <- site.covs_tot$secteur
sectf2 <- ifelse(sectf2=="Beaufort", 1, 2)

## Fit model
model7.spec_marmottes <- textConnection(model7.string_marmottes)
jags_marmottes7 <- jags.model(
  model7.spec_marmottes,
  data = list(
    C = y_tot,
    n = stops_marmottes,
    Xsecteur.derangement = npl_marmottes_secteur_derangement,
    npar.secteur.derangement = dim(npl_marmottes_secteur_derangement)[2],
    prim_ran=prim_ran,
    sec_ran=sec_ran,
    site = site.covs_tot$site,
    secteur = sectf2,
    altitude = site.covs_tot$saltitude[1:552,], # altitude centrée-réduite
    D_beaufort = distance_mat_geo_beaufort / 1000,
    D_lanslevillard = distance_mat_geo_lanslevillard / 1000
  ),
  inits = list(
    beta.secteur.derangementP = rnorm(dim(npl_marmottes_secteur_derangement)[2]),
    lambda_spat = 0.05,
    sigmasq.inv =  1/0.3,
    global.tau_lanslevillard = 0.05,
    global.tau_beaufort = 0.05
  ),
  n.chains = 4,
  n.adapt = 100
)
save(jags_marmottes7, file="Output/jags_marmottes7.Rdata")
# load("Output/jags_marmottes7.Rdata")

## Burn-in
update(jags_marmottes7, 500) #5000

## Save MCMC posteriors for estimated parameters
out_marmottes7_coda <- coda.samples(
  model = jags_marmottes7,
  variable.names =
    c(
      "beta.secteur.derangementP",
      "beta.altN",
      "lambda_spat",
      "sigmasq",
      "global.mu_beaufort",
      "global.mu_lanslevillard",
      "N",
      "lambda"
    ),
  n.iter = 500 #50000
)
save(out_marmottes7_coda, file="Output/out_marmottes7_coda.Rdata")
load("Output/out_marmottes7_coda.Rdata")

colnames(out_marmottes7_coda[[1]])

par(mfrow =  c(3, 3))
## Plot histograms
plot(out_marmottes7_coda[,c(#"global.mu_beaufort","global.mu_lanslevillard",
                            "beta.altN[1]","beta.altN[2]",
                            "beta.secteur.derangementP[1]",
                            "beta.secteur.derangementP[2]",
                            "beta.secteur.derangementP[3]")],
                            #"lambda_spat","sigmasq")],
     trace = FALSE)


## Plot trace
plot(out_marmottes7_coda[,c(#"global.mu_beaufort","global.mu_lanslevillard",
                            "beta.altN[1]","beta.altN[2]",
                            "beta.secteur.derangementP[1]",
                            "beta.secteur.derangementP[2]",
                            "beta.secteur.derangementP[3]")],
                            #"lambda_spat","sigmasq")],
     density = FALSE)


## Pb de convergence quand l'effet du dérangement sur la détection est intégré dans le modèle !
## Convergence ok sans dérangement

summary(out_marmottes7_coda[,c("beta.altN[1]","beta.altN[2]",
                               "beta.secteur.derangementP[1]",
                               "beta.secteur.derangementP[2]",
                               "beta.secteur.derangementP[3]")])


gelman.diag(out_marmottes7_coda[,c("beta.altN[1]","beta.altN[2]",
                                   "beta.secteur.derangementP[1]",
                                   "beta.secteur.derangementP[2]",
                                   "beta.secteur.derangementP[3]")])


par(mfrow=c(1,1))
par(mar=c(2,8,1,1))
caterplot(out_marmottes7_coda, "beta.secteur.derangementP", labels.loc="axis", 
          labels=c("Intercept (Beaufort/ \n Pas de dérangement)",
                   "Lanslevillard",
                   "Dérangement"),
          quantiles=list(outer=c(0.025,0.975),inner=c(0.05,0.95)))

mcmc_marmottes7[,c("beta.altN[1]","beta.altN[2]",
    "beta.secteur.derangementP[1]",
    "beta.secteur.derangementP[2]",
    "beta.secteur.derangementP[3]")]
mcmc_marmottes7 <- as.data.frame(as.matrix(out_marmottes7_coda))
mcmc_marmottes7$beaufort_0Derangement <- mcmc_marmottes7$`beta.secteur.derangementP[1]`
mcmc_marmottes7$lanslevillard_0Derangement <- mcmc_marmottes7$`beta.secteur.derangementP[1]`+mcmc_marmottes7$`beta.secteur.derangementP[2]`
mcmc_marmottes7$beaufort_1Derangement <- mcmc_marmottes7$`beta.secteur.derangementP[1]`+mcmc_marmottes7$`beta.secteur.derangementP[3]`
mcmc_marmottes7$lanslevillard_1Derangement <- mcmc_marmottes7$`beta.secteur.derangementP[1]`+mcmc_marmottes7$`beta.secteur.derangementP[2]`+mcmc_marmottes7$`beta.secteur.derangementP[3]`
mcmc_marmottes7_quant <- t(apply(mcmc_marmottes7, 2, quantile, p=c(0.025, 0.975)))
mcmc_marmottes7_mean <- data.frame(mean_mar=apply(mcmc_marmottes7, 2, mean))
mcmc_marmottes7_all <- cbind(mcmc_marmottes7_mean, mcmc_marmottes7_quant)
mcmc_marmottes7_all$variable <- rownames(mcmc_marmottes7_all)
colnames(mcmc_marmottes7_all)[2:3] <- c("CrI2.5","CrI97.5")
mcmc_marmottes7_all$pmean_mar <- plogis(mcmc_marmottes7_all$mean_mar)
mcmc_marmottes7_all$pCrI2.5 <- plogis(mcmc_marmottes7_all$CrI2.5)
mcmc_marmottes7_all$pCrI97.5 <- plogis(mcmc_marmottes7_all$CrI97.5)
mcmc_marmottes7_all[1667:1670,]
ggplot(mcmc_marmottes7_all[1667:1670,], aes(x=pmean_mar, y=variable))+
  geom_point()+
  geom_errorbar(aes(xmin=pCrI2.5, xmax=pCrI97.5), width=0.1)+
  xlim(c(0,1)) + ylab("") + xlab("Probabilité de détection, IC95%")


## Extraction des N et lambda par point de comptage + IC
mean_lambdaN <- summary(out_marmottes7_coda)$statistics
CI_lambdaN <- summary(out_marmottes7_coda)$quantiles

save(mean_lambdaN, file="Output/mean_lambdaN.Rdata")
save(CI_lambdaN, file="Output/CI_lambdaN.Rdata")
load("Output/mean_lambdaN.Rdata")
load("Output/CI_lambdaN.Rdata")

## Plot N ~ altitude
N_all <- cbind(mean_lambdaN[1:552,], CI_lambdaN[1:552,], site.covs_tot)
colnames(N_all)[c(5,9)] <- c("IC2.5","IC97.5")
N_sum_mailles <- aggregate(list(Mean=N_all$Mean, IC2.5=N_all$IC2.5, IC97.5=N_all$IC97.5),
                           by=list(secteur=N_all$secteur, numero_maille=N_all$maille),
                           function(x){sum(x)/2})
N_sum_mailles$densite <- N_sum_mailles$Mean/(pi*150*150/10000)
aggregate(densite~secteur, N_sum_mailles, ci)

Altitude_mean <- aggregate(altitude~secteur+maille, N_all, mean)
N_all_mailles <- cbind(N_sum_mailles, altitude=Altitude_mean[,"altitude"])
ggplot(N_all_mailles, aes(x=altitude, y=Mean, color=secteur))+
  geom_point()+
  geom_errorbar(aes(ymin=IC2.5, ymax=IC97.5))+
  labs(color="Sites d'étude")+
  ylab("Nombre de marmottes estimées par point de comptage")+
  xlab("Altitude (m)")+
  theme(legend.position="bottom")

## Plot lambda ~ C (observateur secondaire)
lambda_all_obsprim <- cbind(mean_lambdaN[561:1112,], CI_lambdaN[561:1112,],
                            nb_obs=y_tot[,1], site.covs_tot)
lambda_all_obsprim$observateur <- "primaire"
lambda_all_obssec <- cbind(mean_lambdaN[1113:1664,], CI_lambdaN[1113:1664,],
                           nb_obs=y_tot[,2], site.covs_tot)
lambda_all_obssec$observateur <- "secondaire"
lambda_all <- rbind(lambda_all_obsprim, lambda_all_obssec)
colnames(lambda_all)[c(5,9)] <- c("IC2.5","IC97.5")
lambda_sum_mailles <- aggregate(list(Mean=lambda_all$Mean, IC2.5=lambda_all$IC2.5,
                                     IC97.5=lambda_all$IC97.5, nb_obs=lambda_all$nb_obs),
                           by=list(secteur=lambda_all$secteur,
                                   numero_maille=lambda_all$maille,
                                   observateur=lambda_all$observateur),
                           function(x){sum(x)/2})
ggplot(lambda_sum_mailles, aes(x=nb_obs, y=Mean, color=secteur))+
  geom_abline(slope=1, linetype="dashed")+
  geom_point()+
  facet_grid(.~observateur)+
  geom_errorbar(aes(ymin=IC2.5, ymax=IC97.5))+
  labs(color="Sites d'étude")+
  ylab("Nombre de marmottes estimées par point de comptage \n multipliées par la probabilité de détection (lambda)")+
  xlab("Nombre de marmottes observées")+
  theme(legend.position="bottom")

cor(lambda_sum_mailles$Mean, lambda_sum_mailles$nb_obs)
summary(lm(Mean~nb_obs*observateur, lambda_sum_mailles))

