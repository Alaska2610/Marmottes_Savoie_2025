#####################################
## Modèle secteur+dérangement sur P
## CT sur N
## Test sur 2025 avec exactement les mêmes variables
#####################################
library(rjags)
library(geodist)
## La fonction dist() dans R base fonctionne également
library(qgraph)
library(ggmap)
library(mcmcplots)

CB <- TRUE
if(CB == TRUE){setwd("./Scripts")}else
{setwd("G:/FDC_Savoie/Marmottes_Savoie/")}

load("Output/observations1.Rdata")

# Disturbance 0/1
table(observations1$derangement)
observations1$derangement_new <- observations1$derangement
observations1[is.element(observations1$derangement, c("Agricole","Autre","Touristique")),]$derangement_new <- "Derangement"
table(observations1$derangement_new)

observations1$saltitude <- scale(observations1$altitude)

## Unmarked Double-Observer sur toutes les obs
y_tot <- matrix(NA, dim(observations1), 2, byrow=T)
y_tot[,1] <- observations1$nb_obs_prim
y_tot[,2] <- observations1$nb_obs_sec-observations1$nb_obs_prim
y_tot[,2] <- ifelse(y_tot[,2]<0,0,y_tot[,2])
site.covs_tot <- data.frame(observations1[,c("secteur","meteo",
                                             "derangement","derangement_new",
                                             "habitat",
                                             "orientation","latitude","longitude",
                                             "julian_day",
                                             "session","altitude","maille")])
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
                global.mu_beaufort ~ dnorm(0, 0.01)
                global.tau_beaufort ~ dgamma(0.001, 0.001)
                for(i in 1:35){
                    mu_beaufort[i] ~ dnorm(global.mu_beaufort, global.tau_beaufort)
                }     
                
                global.mu_lanslevillard ~ dnorm(0, 0.01)
                global.tau_lanslevillard ~ dgamma(0.001, 0.001)
                for(i in 1:35){
                    mu_lanslevillard[i] ~ dnorm(global.mu_lanslevillard, global.tau_lanslevillard)
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
                    #log(N[i]) <- beta.altN*altitude[i] + CT[site[i]] 
                    log(N[i]) <- beta.altN[secteur[i]]*altitude[i] + CT[site[i]] 

            ## Likelihood for the double-observer survey
                for(j in 1:2){
                    C[i, j]      ~ dpois(lambda[i, j])
                    lambda[i, j] <- N[i] * piMat[i, j]
                    fit[i, j]    <- exp(N[i] * piMat[i, j])
                }
                }

                # Spatial component Lanslevillard
                  CT[1:35] ~ dmnorm.vcov(mu_lanslevillard[1:35], D.covar_lanslevillard[1:35, 1:35])
                  ## turning covariances into precisions
                  for(i in 1:35){
                    D.covar_lanslevillard[i, i] <- sigmasq # diagonale de la matrice de covariance
                        for(j in 1:(i - 1)){
                            # constrain covariance matrix by the distance matrix
                            D.covar_lanslevillard[i, j] <- sigmasq * exp(-(lambda_spat * D_lanslevillard[i, j]))
                            D.covar_lanslevillard[j, i] <- D.covar_lanslevillard[i, j] # symmetrical matrix
                          }
                    }

                  # Spatial component beaufort
                  CT[36:70] ~ dmnorm.vcov(mu_beaufort[1:33], D.covar_beaufort[1:33, 1:33])
                  ## turning covariances into precisions
                  for(i in 1:35){
                    D.covar_beaufort[i, i] <- sigmasq # diagonale de la matrice de covariance
                        for(j in 1:(i - 1)){
                            # constrain covariance matrix by the distance matrix
                            D.covar_beaufort[i, j] <- sigmasq * exp(-(lambda_spat * D_beaufort[i, j]))
                            D.covar_beaufort[j, i] <- D.covar_beaufort[i, j] # symmetrical matrix
                          }
                    }

        }
"

## Define parameters
npl_marmottes_secteur_derangement <- model.matrix(data = site.covs_tot, ~ secteur+derangement_new)
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
    D_lanslevillard = distance_mat_geo_lanslevillard / 1000,
    D_beaufort = distance_mat_geo_beaufort / 1000
  ),
  inits = list(
    beta.secteur.derangementP = rnorm(dim(npl_marmottes_secteur_derangement)[2]),
    lambda_spat = 0.05,
    sigmasq.inv =  1/0.3,
    global.tau_beaufort = 0.05,
    global.tau_lanslevillard = 0.05
  ),
  n.chains = 4,
  n.adapt = 1000
)

setwd("G:/FDC_Savoie/Marmottes_Savoie/BDD/Rdata")
save(jags_marmottes7, file="jags_marmottes7.Rdata")
load("jags_marmottes7.Rdata")

## Once the model compile, run for a large number of iterations
update(jags_marmottes7, 5000)

## Save MCMC posteriors for estimated parameters
out_marmottes7_coda <- coda.samples(
  model = jags_marmottes7,
  variable.names =
    c(
      "beta.secteur.derangementP",
      "beta.altN",
      "lambda_spat",
      "sigmasq",
      "global.mu_lanslevillard",
      "global.mu_beaufort",
      "N",
      "lambda"
    ),
  n.iter = 50000
)

save(out_marmottes7_coda, file="Ouput/out_marmottes7_coda.Rdata")
load("Output/out_marmottes7_coda.Rdata")

par(mfrow =  c(3, 3))
## Plot histograms
plot(out_marmottes7_coda, trace = FALSE)
## Plot trace
plot(out_marmottes7_coda, density = FALSE)
# mixage ok 
summary(out_marmottes7_coda)
# Effet du dérangement 
gelman.diag(out_marmottes7_coda)

mcmc_marmottes7 <- as.data.frame(as.matrix(out_marmottes7_coda))
plot(density(mcmc_marmottes7$'global.mu_beaufort'))

# secteur + derangement_new sur P, altitude + CT sur N
# altitude centré-réduit
# exp(0.83)=2.3 => plus du doublement du nb de marmottes par augmentation d'1 unité
# 1 unité ~ 235m (sd(observations1$altitude))
# exp(0.83/235*100)=1.42 => 42% d'augmentation du nb de marmottes tous les 100m
# exp(0.18/235*100)=1.08 => IC inf
# exp(1.54/235*100)=1.93 => IC sup
# 1. Empirical mean and standard deviation for each variable,
# plus standard error of the mean:
#   
#                                 Mean     SD Naive SE Time-series SE
# beta.altN                     0.8309 0.3415 0.001707       0.053103
# beta.secteur.derangementP[1]  1.6318 0.3164 0.001582       0.004162
# beta.secteur.derangementP[2]  1.1529 0.5217 0.002609       0.005068
# beta.secteur.derangementP[3] -1.2608 0.5581 0.002790       0.012158
# global.mu_beaufort            -0.8240 0.9130 0.004565       0.189082
# global.mu_lanslevillard           -2.1365 0.7648 0.003824       0.084001
# lambda_spat                   1.9445 1.3420 0.006710       0.072650
# sigmasq                       1.9726 1.4385 0.007193       0.085435
# 
# 2. Quantiles for each variable:
#   
#                                 2.5%     25%     50%     75%   97.5%
# beta.altN                     0.1812  0.5995  0.8383  1.0388  1.5427
# beta.secteur.derangementP[1]  1.0204  1.4190  1.6279  1.8420  2.2632
# beta.secteur.derangementP[2]  0.1510  0.7971  1.1419  1.4960  2.2097
# beta.secteur.derangementP[3] -2.3360 -1.6391 -1.2620 -0.8826 -0.1596
# global.mu_beaufort            -2.8378 -1.2931 -0.8614 -0.3688  1.5828
# global.mu_lanslevillard           -3.6112 -2.5811 -2.1757 -1.7362 -0.3253
# lambda_spat                   0.3815  1.0439  1.6376  2.4623  5.4525
# sigmasq                       0.4749  1.1712  1.6675  2.3645  5.2835

par(mfrow=c(1,1))
par(mar=c(2,8,1,1))
caterplot(out_marmottes7_coda, "beta.secteur.derangementP", labels.loc="axis", 
          labels=c("Intercept (beaufort/ \n Pas de dérangement)",
                   "lanslevillard",
                   "Dérangement"),
          quantiles=list(outer=c(0.025,0.975),inner=c(0.05,0.95)))

# C'est correct de faire comme ça pour avoir les IC de pdet par secteur et dérangement ?
mcmc_marmottes7 <- as.data.frame(as.matrix(out_marmottes7_coda))
mcmc_marmottes7$beaufort_0Derangement <- mcmc_marmottes7$`beta.secteur.derangementP[1]`
mcmc_marmottes7$lanslevillard_0Derangement <- mcmc_marmottes7$`beta.secteur.derangementP[1]`+mcmc_marmottes7$`beta.secteur.derangementP[2]`
mcmc_marmottes7$beaufort_1Derangement <- mcmc_marmottes7$`beta.secteur.derangementP[1]`+mcmc_marmottes7$`beta.secteur.derangementP[3]`
mcmc_marmottes7$lanslevillard_1Derangement <- mcmc_marmottes7$`beta.secteur.derangementP[1]`+mcmc_marmottes7$`beta.secteur.derangementP[2]`++mcmc_marmottes7$`beta.secteur.derangementP[3]`
mcmc_marmottes7_quant <- t(apply(mcmc_marmottes7, 2, quantile, p=c(0.025, 0.975)))
mcmc_marmottes7_mean <- data.frame(mean_mar=apply(mcmc_marmottes7, 2, mean))
mcmc_marmottes7_all <- cbind(mcmc_marmottes7_mean, mcmc_marmottes7_quant)
mcmc_marmottes7_all$variable <- rownames(mcmc_marmottes7_all)
colnames(mcmc_marmottes7_all)[2:3] <- c("CrI2.5","CrI97.5")
mcmc_marmottes7_all$pmean_mar <- plogis(mcmc_marmottes7_all$mean_mar)
mcmc_marmottes7_all$pCrI2.5 <- plogis(mcmc_marmottes7_all$CrI2.5)
mcmc_marmottes7_all$pCrI97.5 <- plogis(mcmc_marmottes7_all$CrI97.5)
mcmc_marmottes7_all[1642:1645,]
ggplot(mcmc_marmottes7_all[1642:1645,], aes(x=pmean_mar, y=variable))+
  geom_point()+
  geom_errorbar(aes(xmin=pCrI2.5, xmax=pCrI97.5), width=0.1)+
  xlim(c(0,1)) + ylab("") + xlab("Probabilité de détection, IC95%")


## Extraction des N et lambda par point de comptage + IC
mean_lambdaN <- summary(out_marmottes7_coda)$statistics
CI_lambdaN <- summary(out_marmottes7_coda)$quantiles

setwd("G:/FDC_Savoie/Marmottes_Savoie/BDD/Rdata")
save(mean_lambdaN, file="mean_lambdaN.Rdata")
save(CI_lambdaN, file="CI_lambdaN.Rdata")
load("mean_lambdaN.Rdata")
load("CI_lambdaN.Rdata")

## Plot N ~ altitude
N_all <- cbind(mean_lambdaN[1:544,], CI_lambdaN[1:544,], site.covs_tot)
colnames(N_all)[c(5,9)] <- c("IC2.5","IC97.5")
N_sum_mailles <- aggregate(list(Mean=N_all$Mean, IC2.5=N_all$IC2.5, IC97.5=N_all$IC97.5), 
                           by=list(secteur=N_all$secteur, numero_maille=N_all$numero_maille),
                           function(x){sum(x)/2})
N_sum_mailles$densite <- N_sum_mailles$Mean/(pi*150*150/10000)
aggregate(densite~secteur, N_sum_mailles, ci)
Altitude_mean <- aggregate(altitude~secteur+numero_maille, N_all, mean)
N_all_mailles <- cbind(N_sum_mailles, altitude=Altitude_mean[,"altitude"])
ggplot(N_all_mailles, aes(x=altitude, y=Mean, color=secteur))+
  geom_point()+
  geom_errorbar(aes(ymin=IC2.5, ymax=IC97.5))+
  labs(color="Sites d'étude")+
  ylab("Nombre de marmottes estimées par point de comptage")+
  xlab("Altitude (m)")+
  theme(legend.position="bottom")

## Plot lambda ~ C (observateur secondaire)
lambda_all_obsprim <- cbind(mean_lambdaN[551:1094,], CI_lambdaN[551:1094,], nb_obs=y_tot[,"na"], site.covs_tot)
lambda_all_obsprim$observateur <- "primaire"
lambda_all_obssec <- cbind(mean_lambdaN[1095:1638,], CI_lambdaN[1095:1638,], nb_obs=y_tot[,"nb"], site.covs_tot)
lambda_all_obssec$observateur <- "secondaire"
lambda_all <- rbind(lambda_all_obsprim, lambda_all_obssec)
colnames(lambda_all)[c(5,9)] <- c("IC2.5","IC97.5")
lambda_sum_mailles <- aggregate(list(Mean=lambda_all$Mean, IC2.5=lambda_all$IC2.5, IC97.5=lambda_all$IC97.5, nb_obs=lambda_all$nb_obs), 
                           by=list(secteur=lambda_all$secteur, numero_maille=lambda_all$numero_maille, observateur=lambda_all$observateur),
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

## 545:551
mean_lambdaN[c(545:551,1640:1641),]
CI_lambdaN[c(545:551,1640:1641),]

## N corrected by Pavailability
setwd("G:/FDC_Savoie/Marmottes_Savoie/BDD/Rdata")
load("out_marmottes_param.Rdata")
load("out_marmottes7_coda.Rdata")
mcmc_marmottes7 <- as.data.frame(as.matrix(out_marmottes7_coda))
mcmc_marmottes_pdispo <- as.data.frame(as.matrix(out_marmottes_param))
for(i in 1:544){
  mcmc_marmottes7[,i] <- mcmc_marmottes7[,i]/mcmc_marmottes_pdispo
}

newNdispo <- cbind(t(apply(mcmc_marmottes7[,1:544], 2, quantile, p=c(0.025, 0.975))), apply(mcmc_marmottes7[,1:544], 2, mean))
newNdispo1 <- as.data.frame(cbind(newNdispo,
                                  site.covs_tot=site.covs_tot))
colnames(newNdispo1)[c(1:3)] <- c("IC2.5","IC97.5","MeanDispo")
newNdispo1_ag <- aggregate(list(MeanDispo=newNdispo1$MeanDispo, IC2.5=newNdispo1$IC2.5, IC97.5=newNdispo1$IC97.5), 
                                by=list(secteur=newNdispo1$site.covs_tot.secteur, numero_maille=newNdispo1$site.covs_tot.numero_maille),
                                function(x){sum(x)/2})
newNdispo1_ag$densite_meandispo <- newNdispo1_ag$MeanDispo/(pi*150*150/10000)
newNdispo1_ag$densite_icinf <- newNdispo1_ag[,"IC2.5"]/(pi*150*150/10000)
newNdispo1_ag$densite_icsup <- newNdispo1_ag[,"IC97.5"]/(pi*150*150/10000)
aggregate(densite_meandispo~secteur, newNdispo1_ag, ci)
aggregate(densite_icinf~secteur, newNdispo1_ag, ci)
aggregate(densite_icsup~secteur, newNdispo1_ag, ci)

alt_coord1 <- aggregate(list(altitude=N_marmottes7$site.covs_tot.altitude, x=N_marmottes7$site.covs_tot.x, y=N_marmottes7$site.covs_tot.y), 
                       by=list(numero_maille=N_marmottes7$site.covs_tot.numero_maille), mean)
newNdispo1_ag2 <- merge(newNdispo1_ag, alt_coord1, by="numero_maille")

N_marmottes7_lanslevillard_cordispo <- newNdispo1_ag2[which(newNdispo1_ag2$x<6.8),]
N_marmottes7_beaufort_cordispo <- newNdispo1_ag2[which(newNdispo1_ag2$x>6.8),]

ggmap::register_google(key="AIzaSyA53J3oEn4CEPw1xB7Grb2Ei_-AYYdcXes")
mapmarmot_lanslevillard_N <- get_map(location = c(lon = mean(N_marmottes7_lanslevillard_cordispo$x), lat = mean(N_marmottes7_lanslevillard_cordispo$y)), zoom = 13,
                                maptype = "terrain", scale = 2)

ggmap(mapmarmot_lanslevillard_N) +
  geom_point(data = N_marmottes7_lanslevillard_cordispo, aes(x = x, y = y, size=densite_meandispo))+
  ggtitle("lanslevillard")+ scale_size_area(limits=c(0,5.5))+
  labs(size="Densité")
#+scale_color_gradient(colours=rainbow(5))

mapmarmot_beaufort_N <- get_map(location = c(lon = mean(N_marmottes7_beaufort_cordispo$x), lat = mean(N_marmottes7_beaufort_cordispo$y)), zoom = 13,
                               maptype = "terrain", scale = 2)
ggmap(mapmarmot_beaufort_N) +
  geom_point(data = N_marmottes7_beaufort_cordispo, aes(x = x, y = y, size=densite_meandispo))+
  ggtitle("beaufort")+ scale_size_area(limits=c(0,5.5))+
  labs(size="Densité")


## Get more parameters
out_marmottes7_samples <- jags.samples(
  model = jags_marmottes7,
  variable.names =
    c(
      "lambda",
      "N",
      "pa",
      "pb",
      "beta.secteur.derangementP",
      "lambda_spat",
      "sigmasq",
      "D.covar_beaufort",
      "D.covar_lanslevillard",
      "CT",
      "mu_beaufort",
      "mu_lanslevillard",
      "beta.altN"
    ),
  n.iter = 10000
)

setwd("G:/FDC_Savoie/Marmottes_Savoie/BDD/Rdata")
save(out_marmottes7_samples, file="out_marmottes7_samples.Rdata")
load("out_marmottes7_samples.Rdata")

mean(apply(out_marmottes7_samples$mu_beaufort, 1, mean))
mean(apply(out_marmottes7_samples$mu_lanslevillard, 1, mean))
mean(apply(out_marmottes7_samples$lambda_spat, 1, mean))
quantile(out_marmottes7_samples$lambda_spat, c(0.025,0.975))
out_marmottes7_samples$N
mean(apply(out_marmottes7_samples$pa, 1, mean))
mean(apply(out_marmottes7_samples$pb, 1, mean))
quantile(out_marmottes7_samples$pa, c(0.025,0.975))
quantile(out_marmottes7_samples$pb, c(0.025,0.975))
apply(out_marmottes7_samples$beta.altN, 1, mean)
quantile(out_marmottes7_samples$beta.altN, c(0.025,0.975))
apply(out_marmottes7_samples$beta.secteur.derangementP, 1, mean)

N_marmottes7 <- as.data.frame(cbind(N=apply(out_marmottes7_samples$N,1,mean),
                                    pa=apply(out_marmottes7_samples$pa,1,mean),
                                    pb=apply(out_marmottes7_samples$pb,1,mean),
                                    site.covs_tot=site.covs_tot,
                                    y_tot=y_tot))
N_marmottes7_ag <- aggregate(N~site.covs_tot.numero_maille+site.covs_tot.periode+site.covs_tot.secteur, N_marmottes7, sum)
N_marmottes7_ag1 <- aggregate(N~site.covs_tot.numero_maille+site.covs_tot.secteur, N_marmottes7_ag, mean)
alt_coord <- aggregate(list(altitude=N_marmottes7$site.covs_tot.altitude, x=N_marmottes7$site.covs_tot.x, y=N_marmottes7$site.covs_tot.y), 
                       by=list(site.covs_tot.numero_maille=N_marmottes7$site.covs_tot.numero_maille), mean)
N_marmottes7_ag2 <- merge(N_marmottes7_ag1, alt_coord, by="site.covs_tot.numero_maille")
N_marmottes7_ag2$densite <- N_marmottes7_ag2$N/(pi*150*150/10000)
aggregate(densite~site.covs_tot.secteur, N_marmottes7_ag2, ci)
ggplot(N_marmottes7_ag2, aes(x=altitude, y=N, colour=site.covs_tot.secteur))+
  #geom_boxplot(aes(colour=site.covs_tot.secteur))+
  geom_point()+
  xlab("")+
  ylab("Nombre de marmottes estimé par point de comptage")+
  labs(color = "Sites d'étude")

N_marmottes7_lanslevillard <- N_marmottes7_ag2[which(N_marmottes7_ag2$x<6.8),]
N_marmottes7_beaufort <- N_marmottes7_ag2[which(N_marmottes7_ag2$x>6.8),]

ggmap::register_google(key="AIzaSyA53J3oEn4CEPw1xB7Grb2Ei_-AYYdcXes")
mapmarmot_lanslevillard_N <- get_map(location = c(lon = mean(N_marmottes7_lanslevillard$x), lat = mean(N_marmottes7_lanslevillard$y)), zoom = 13,
                                maptype = "terrain", scale = 2)

ggmap(mapmarmot_lanslevillard_N) +
  geom_point(data = N_marmottes7_lanslevillard, aes(x = x, y = y, size=densite))+
  ggtitle("lanslevillard")+ scale_size_area(limits=c(0,2.2))
#+scale_color_gradient(colours=rainbow(5))

mapmarmot_beaufort_N <- get_map(location = c(lon = mean(N_marmottes7_beaufort$x), lat = mean(N_marmottes7_beaufort$y)), zoom = 13,
                               maptype = "terrain", scale = 2)
ggmap(mapmarmot_beaufort_N) +
  geom_point(data = N_marmottes7_beaufort, aes(x = x, y = y, size=densite))+
  ggtitle("beaufort")+ scale_size_area(limits=c(0,2.2))
