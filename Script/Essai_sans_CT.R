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
model1.string_marmottes <-"
          model {
          
            ## Set priors for parameters to be estimated
  delta ~ dnorm(0, 0.001)
  
  # Fixed parameters for detection
  for(k in 1:npar.secteur){
    beta.secteurP[k] ~ dnorm(0, 0.0001)
  }
  
  for(k in 1:npar.habitat){
    beta.habitatN[k] ~ dnorm(0, 0.0001)
  }
  
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
    log(N[i]) <- inprod(beta.habitatN[], Xhabitat[i, ])
    
    ## Likelihood for the double-observer survey
    for(j in 1:2){
      C[i, j]      ~ dpois(lambda[i, j])
      lambda[i, j] <- N[i] * piMat[i, j]
      fit[i, j]    <- exp(N[i] * piMat[i, j])
    }
    
  }

          }
"

## Define parameters
npl_marmottes <- model.matrix(data = site.covs_tot, ~ secteur)
npl_marmottes_habitat   <- model.matrix(data = site.covs_tot, ~ typologie_habitat_new)
stops_marmottes <- length(observations1[,1])

prim_ran <- ifelse(observations1$obs_prim=="Manolo"|observations1$obs_prim=="Léane", 1, 2)
sec_ran <- ifelse(prim_ran==1, 2, 1)

sectf2 <- site.covs_tot$secteur
sectf2 <- ifelse(sectf2=="Beaufort", 1, 2)

## Fit model
model1.spec_marmottes <- textConnection(model1.string_marmottes)
jags_marmottes1 <- jags.model(
  model1.spec_marmottes,
  data = list(
    C = y_tot,
    n = stops_marmottes,
    Xsecteur = npl_marmottes,
    Xhabitat = npl_marmottes_habitat,
    npar.secteur = dim(npl_marmottes)[2],
    npar.habitat = dim(npl_marmottes_habitat)[2],
    prim_ran=prim_ran,
    sec_ran=sec_ran
  ),
  inits = list(
    beta.secteurP = rnorm(dim(npl_marmottes)[2])
  ),
  n.chains = 4,
  n.adapt = 100
)
save(jags_marmottes1, file="Output/jags_marmottes1.Rdata")

## Burn-in
update(jags_marmottes1, 500) #5000

## Save MCMC posteriors for estimated parameters
out_marmottes1_coda <- coda.samples(
  model = jags_marmottes1,
  variable.names =
    c("lambda", "N", 
      "beta.habitatN",
      "beta.secteurP"
    ),
  n.iter = 500 #50000
)
save(out_marmottes1_coda, file="Output/out_marmottes1_coda.Rdata")

colnames(out_marmottes1_coda[[1]])

par(mfrow =  c(3, 3))
## Plot histograms
plot(out_marmottes1_coda[,c(
  "beta.habitatN[1]",
  "beta.habitatN[2]",
  "beta.habitatN[3]",
  "beta.secteurP[1]",
  "beta.secteurP[2]")],
  trace = FALSE)


## Plot trace
plot(out_marmottes1_coda[,c(
  "beta.habitatN[1]",
  "beta.habitatN[2]",
  "beta.habitatN[3]",
  "beta.secteurP[1]",
  "beta.secteurP[2]")],
  density = FALSE)

summary(out_marmottes1_coda[,c("beta.habitatN[1]",
                               "beta.habitatN[2]",
                               "beta.habitatN[3]",
                               "beta.secteurP[1]",
                               "beta.secteurP[2]")])


gelman.diag(out_marmottes1_coda[,c("beta.habitatN[1]",
                                   "beta.habitatN[2]",
                                   "beta.habitatN[3]",
                                   "beta.secteurP[1]",
                                   "beta.secteurP[2]")])


par(mfrow=c(1,1))
par(mar=c(2,8,1,1))
caterplot(out_marmottes1_coda, "beta.habitatN", labels.loc="axis", 
          labels=c("Alpages",
                   "Eboulis",
                   "Pré/bois"),
          quantiles=list(outer=c(0.025,0.975),inner=c(0.05,0.95)))

caterplot(out_marmottes1_coda, "beta.secteurP", labels.loc="axis", 
          labels=c("Intercept (Beaufort)",
                   "Lanslevillard"),
          quantiles=list(outer=c(0.025,0.975),inner=c(0.05,0.95)))

mcmc_marmottes1 <- as.data.frame(as.matrix(out_marmottes1_coda))
mcmc_marmottes1[,c("beta.habitatN[1]",
                   "beta.habitatN[2]",
                   "beta.habitatN[3]")]
mcmc_marmottes1$Beaufort <- mcmc_marmottes1$`beta.secteurP[1]`
mcmc_marmottes1$Lanslevillard <- mcmc_marmottes1$`beta.secteurP[1]`+mcmc_marmottes1$`beta.secteurP[2]`
mcmc_marmottes1_quant <- t(apply(mcmc_marmottes1, 2, quantile, p=c(0.025, 0.975)))
mcmc_marmottes1_mean <- data.frame(mean_mar=apply(mcmc_marmottes1, 2, mean))
mcmc_marmottes1_all <- cbind(mcmc_marmottes1_mean, mcmc_marmottes1_quant)
mcmc_marmottes1_all$variable <- rownames(mcmc_marmottes1_all)
colnames(mcmc_marmottes1_all)[2:3] <- c("CrI2.5","CrI97.5")
mcmc_marmottes1_all$pmean_mar <- plogis(mcmc_marmottes1_all$mean_mar)
mcmc_marmottes1_all$pCrI2.5 <- plogis(mcmc_marmottes1_all$CrI2.5)
mcmc_marmottes1_all$pCrI97.5 <- plogis(mcmc_marmottes1_all$CrI97.5)
mcmc_marmottes1_all[1662:1663,]
ggplot(mcmc_marmottes1_all[1662:1663,], aes(x=pmean_mar, y=variable))+
  geom_point()+
  geom_errorbar(aes(xmin=pCrI2.5, xmax=pCrI97.5), width=0.1)+
  xlim(c(0,1)) + ylab("") + xlab("Probabilité de détection, IC95%")

## Extraction des N et lambda par point de comptage + IC
mean_lambdaN <- summary(out_marmottes1_coda)$statistics
CI_lambdaN <- summary(out_marmottes1_coda)$quantiles

## N ~ secteur et N ~ habitat
N_all <- cbind(mean_lambdaN[1:552,], CI_lambdaN[1:552,], site.covs_tot)
colnames(N_all)[c(5,9)] <- c("IC2.5","IC97.5")
N_sum_mailles <- aggregate(list(Mean=N_all$Mean, IC2.5=N_all$IC2.5, IC97.5=N_all$IC97.5),
                           by=list(secteur=N_all$secteur, numero_maille=N_all$maille, habitat=N_all$typologie_habitat_new),
                           function(x){sum(x)/2})
N_sum_mailles$densite <- N_sum_mailles$Mean/(pi*150*150/10000)
aggregate(densite~secteur, N_sum_mailles, ci)
aggregate(densite~habitat, N_sum_mailles, ci)

## Plot lambda ~ C (observateur secondaire)
lambda_all_obsprim <- cbind(mean_lambdaN[558:1109,], CI_lambdaN[558:1109,],
                            nb_obs=y_tot[,1], site.covs_tot)
lambda_all_obsprim$observateur <- "primaire"
lambda_all_obssec <- cbind(mean_lambdaN[1110:1661,], CI_lambdaN[1110:1661,],
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


## N corrected by Pavailability
load("Output/out_marmottes_parametrage.Rdata")
load("Output/out_marmottes1_coda.Rdata")
mcmc_marmottes1 <- as.data.frame(as.matrix(out_marmottes1_coda))
mcmc_marmottes_pdispo <- as.data.frame(as.matrix(out_marmottes_parametrage))
for(i in 1:552){
  mcmc_marmottes1[,i] <- mcmc_marmottes1[,i]/mcmc_marmottes_pdispo[1:2000,]
}

newNdispo <- cbind(t(apply(mcmc_marmottes1[,1:552], 2, quantile, p=c(0.025, 0.975))), apply(mcmc_marmottes1[,1:552], 2, mean))
newNdispo1 <- as.data.frame(cbind(newNdispo,
                                  site.covs_tot=site.covs_tot))
colnames(newNdispo1)[c(1:3)] <- c("IC2.5","IC97.5","MeanDispo")
newNdispo1_ag <- aggregate(list(MeanDispo=newNdispo1$MeanDispo, IC2.5=newNdispo1$IC2.5, IC97.5=newNdispo1$IC97.5),
                           by=list(secteur=newNdispo1$site.covs_tot.secteur, numero_maille=newNdispo1$site.covs_tot.site),
                           function(x){sum(x)/2})
newNdispo1_ag$densite_meandispo <- newNdispo1_ag$MeanDispo/(pi*150*150/10000)
newNdispo1_ag$densite_icinf <- newNdispo1_ag[,"IC2.5"]/(pi*150*150/10000)
newNdispo1_ag$densite_icsup <- newNdispo1_ag[,"IC97.5"]/(pi*150*150/10000)
aggregate(densite_meandispo~secteur, newNdispo1_ag, ci)
aggregate(densite_icinf~secteur, newNdispo1_ag, ci)
aggregate(densite_icsup~secteur, newNdispo1_ag, ci)



alt_coord1 <- aggregate(list(altitude=N_marmottes1$site.covs_tot.altitude, x=N_marmottes1$site.covs_tot.longitude, y=N_marmottes1$site.covs_tot.latitude), 
                        by=list(numero_maille=N_marmottes1$site.covs_tot.site), mean)
newNdispo1_ag2 <- merge(newNdispo1_ag, alt_coord1, by="numero_maille")

N_marmottes1_lanslevillard_cordispo <- newNdispo1_ag2[which(newNdispo1_ag2$x>6.8),]
N_marmottes1_beaufort_cordispo <- newNdispo1_ag2[which(newNdispo1_ag2$x<6.8),]

ggmap::register_google(key="AIzaSyA53J3oEn4CEPw1xB7Grb2Ei_-AYYdcXes")
mapmarmot_lanslevillard_N <- get_map(location = c(lon = mean(N_marmottes1_lanslevillard_cordispo$x), lat = mean(N_marmottes1_lanslevillard_cordispo$y)), zoom = 13,
                                     maptype = "terrain", scale = 2)

ggmap(mapmarmot_lanslevillard_N) +
  geom_point(data = N_marmottes1_lanslevillard_cordispo, aes(x = x, y = y, size=densite_meandispo))+
  ggtitle("Lanslevillard")+ scale_size_area(limits=c(0,5.5))+
  labs(size="Densité")
#+scale_color_gradient(colours=rainbow(5))

mapmarmot_beaufort_N <- get_map(location = c(lon = mean(N_marmottes1_beaufort_cordispo$x), lat = mean(N_marmottes1_beaufort_cordispo$y)), zoom = 13,
                                maptype = "terrain", scale = 2)
ggmap(mapmarmot_beaufort_N) +
  geom_point(data = N_marmottes1_beaufort_cordispo, aes(x = x, y = y, size=densite_meandispo))+
  ggtitle("Beaufort")+ scale_size_area(limits=c(0,5.5))+
  labs(size="Densité")


## Get more parameters
out_marmottes1_samples <- jags.samples(
  model = jags_marmottes1,
  variable.names =
    c(
      "lambda",
      "N",
      "pa",
      "pb",
      "beta.secteurP",
      "lambda_spat",
      "sigmasq",
      "D.covar_beaufort",
      "D.covar_lanslevillard",
      "CT",
      "mu_beaufort",
      "mu_lanslevillard",
      "beta.habitatN"
    ),
  n.iter = 500
)
save(out_marmottes1_samples, file="Output/out_marmottes1_samples.Rdata")
load("Output/out_marmottes1_samples.Rdata")

mean(apply(out_marmottes1_samples$mu_beaufort, 1, mean))
mean(apply(out_marmottes1_samples$mu_lanslevillard, 1, mean))
mean(apply(out_marmottes1_samples$lambda_spat, 1, mean))
quantile(out_marmottes1_samples$lambda_spat, c(0.025,0.975))
out_marmottes1_samples$N
mean(apply(out_marmottes1_samples$pa, 1, mean))
mean(apply(out_marmottes1_samples$pb, 1, mean))
quantile(out_marmottes1_samples$pa, c(0.025,0.975))
quantile(out_marmottes1_samples$pb, c(0.025,0.975))
apply(out_marmottes1_samples$beta.secteurP, 1, mean)

N_marmottes1 <- as.data.frame(cbind(N=apply(out_marmottes1_samples$N,1,mean),
                                    pa=apply(out_marmottes1_samples$pa,1,mean),
                                    pb=apply(out_marmottes1_samples$pb,1,mean),
                                    site.covs_tot=site.covs_tot,
                                    y_tot=y_tot))
N_marmottes1_ag <- aggregate(N~site.covs_tot.site+site.covs_tot.periode+site.covs_tot.secteur, N_marmottes1, sum)
N_marmottes1_ag1 <- aggregate(N~site.covs_tot.site+site.covs_tot.secteur, N_marmottes1_ag, mean)
alt_coord <- aggregate(list(altitude=N_marmottes1$site.covs_tot.altitude, x=N_marmottes1$site.covs_tot.x, y=N_marmottes1$site.covs_tot.y),
                       by=list(site.covs_tot.site=N_marmottes1$site.covs_tot.site), mean)
N_marmottes1_ag2 <- merge(N_marmottes1_ag1, alt_coord, by="site.covs_tot.site")
N_marmottes1_ag2$densite <- N_marmottes1_ag2$N/(pi*150*150/10000)
aggregate(densite~site.covs_tot.secteur, N_marmottes1_ag2, ci)
ggplot(N_marmottes1_ag2, aes(x=altitude, y=N, colour=site.covs_tot.secteur))+
  #geom_boxplot(aes(colour=site.covs_tot.secteur))+
  geom_point()+
  xlab("")+
  ylab("Nombre de marmottes estimé par point de comptage")+
  labs(color = "Sites d'étude")

N_marmottes1_lanslevillard <- N_marmottes1_ag2[which(N_marmottes1_ag2$x<6.8),]
N_marmottes1_beaufort <- N_marmottes1_ag2[which(N_marmottes1_ag2$x>6.8),]

ggmap::register_google(key="AIzaSyA53J3oEn4CEPw1xB7Grb2Ei_-AYYdcXes")
mapmarmot_lanslevillard_N <- get_map(location = c(lon = mean(N_marmottes1_lanslevillard$x), lat = mean(N_marmottes1_lanslevillard$y)), zoom = 13,
                                     maptype = "terrain", scale = 2)

ggmap(mapmarmot_lanslevillard_N) +
  geom_point(data = N_marmottes1_lanslevillard, aes(x = x, y = y, size=densite))+
  ggtitle("lanslevillard")+ scale_size_area(limits=c(0,2.2))
#+scale_color_gradient(colours=rainbow(5))

mapmarmot_beaufort_N <- get_map(location = c(lon = mean(N_marmottes1_beaufort$x), lat = mean(N_marmottes1_beaufort$y)), zoom = 13,
                                maptype = "terrain", scale = 2)
ggmap(mapmarmot_beaufort_N) +
  geom_point(data = N_marmottes1_beaufort, aes(x = x, y = y, size=densite))+
  ggtitle("beaufort")+ scale_size_area(limits=c(0,2.2))