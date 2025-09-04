setwd("G:/FDC_Savoie/Marmottes_Savoie/Scripts")
source("1_MiseEnForme.R")

################################################################################
# Visualisation des données -----------------------------------------------
################################################################################

  ## Paramétrage
ggplot(parametrage, aes(duree))+
  geom_histogram(bins=10)+
  facet_grid(.~secteur)+
  xlab("Durée d'observations de 30 marmottes")+
  ylab("Fréquence")

# shapiro.test(parametrage$duree)
# shapiro.test(log(parametrage$duree))

moy_param_log <- aggregate(log_duree~secteur, parametrage, ci)
# exp(0.67) ; exp(0.97)

summary(lm(log_duree~secteur, parametrage))
# summary(lm(log_duree~secteur, parametrage[which(parametrage$duree<30),]))
# Pas de différences significatives de durées d'observations des marmottes entre les 2 secteurs.


  ## Observations
# Création des tables pour les graphiques ci-après
tab_secteur <- table(observations1$secteur)
tab_maille_secteur <- aggregate(numero_maille~secteur, observations1, function(x){length(unique(x))})
# observations1[is.element(observations$numero_maille, c("455","456")),]
# 33 mailles à Bessans au lieu de 31. Ok pour tout garder ? Mailles 455 et 456 en plus. A confirmer avec la FDC.

tab_meteo_secteur <- as.data.frame(table(observations1$condition_meteo, observations1$secteur))
tab_habitat_secteur <- as.data.frame(table(observations1$typologie_habitat, observations1$secteur))
tab_habitat_derangement <- as.data.frame(table(observations1$derangement, observations1$secteur))

moy_obs_prim <- aggregate(obs_primaire_nb_observation~secteur, observations1, ci)
moy_obs_sec <- aggregate(obs_secondaire_nb_observation~secteur, observations1, ci)

# A-t-on une météo et habitat par maille et par période ?
# Non pour la météo
tab_meteo1 <- table(observations1$condition_meteo, observations1$numero_maille, observations1$periode)
# Non pour l'habitat
tab_habitat2 <- table(observations1$typologie_habitat, observations1$numero_maille, observations1$periode)

nb_obs_prim <- aggregate(obs_primaire_nb_observation~secteur, observations1, sum)
# 161 à Bessans, 87 aux Menuires
colnames(nb_obs_prim)[2] <- "nb_obs"
nb_obs_prim$observateur <- "primaire"
nb_obs_sec <- aggregate(obs_secondaire_nb_observation~secteur, observations1, sum)
# 201 à Bessans, 94 aux Menuires
colnames(nb_obs_sec)[2] <- "nb_obs"
nb_obs_sec$observateur <- "secondaire"
nb_obs_prim_sec <- rbind(nb_obs_prim, nb_obs_sec)

nb_obs_max_prim <- aggregate(obs_primaire_nb_observation~secteur, observations1, max)
nb_obs_min_prim <- aggregate(obs_primaire_nb_observation~secteur, observations1, min)

nb_periode_secteur <- aggregate(obs_secondaire_nb_observation~secteur*periode, observations1, sum)

nb_obs_mailles <- aggregate(obs_secondaire_nb_observation~secteur*periode*numero_maille, observations1, sum)
nb_obs_mailles1 <- aggregate(obs_secondaire_nb_observation~secteur*numero_maille, nb_obs_mailles, mean)
nb_obs_mailles1$densite <- nb_obs_mailles1$obs_secondaire_nb_observation/(pi*150*150/10000)
nb_obs_mailles2 <- aggregate(densite~secteur, nb_obs_mailles1, ci)
# Nombre moyen de marmottes observées par maille

# Surface prospectée par point (rayon d'observation de 150m)
surf_point <- pi*150*150/10000 # surface prospectée de ~7ha par maille
# surf_point*33
# surf_point*35
#####

# Graphs
ggplot(nb_obs_prim_sec, aes(y=nb_obs, x=secteur, fill=observateur))+
  geom_bar(stat="identity", position="dodge", width=0.5)+
  xlab("Secteur")+
  ylab("Nombre de marmottes observées")+
  labs(fill="Observateur")

ggplot(observations1, aes(obs_secondaire_nb_observation))+
  geom_histogram()+
  facet_grid(.~secteur)+
  xlab("Nombre de marmottes observées par l'observateur secondaire")+
  ylab("Fréquence")

ggplot(observations1, aes(diff_obs))+
  geom_histogram()+
  facet_grid(.~secteur)+
  xlab("Différences du nombre de marmottes entre \n l'observateur primaire et secondaire")+
  ylab("Fréquence")

# Répartition des quarts de mailles selon la météo
ggplot(tab_meteo_secteur, aes(x=Var1, y=Freq, fill=Var2))+
  geom_bar(stat="identity", position="dodge")+
  labs(fill = "Secteur")+
  xlab("Météo")+
  ylab("Fréquence")

# Répartition des quarts de mailles selon l'habitat
ggplot(tab_habitat_secteur, aes(x=Var1, y=Freq, fill=Var2))+
  geom_bar(stat="identity", position="dodge")+
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1, size=12))+
  labs(fill = "Secteur")+
  xlab("Habitat")+
  ylab("Fréquence")

# Répartition des quarts de mailles selon le dérangement
ggplot(tab_habitat_derangement, aes(x=Var1, y=Freq, fill=Var2))+
  geom_bar(stat="identity", position="dodge")+
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1, size=12))+
  labs(fill = "Sites d'étude")+
  xlab("Dérangement")+
  ylab("Fréquence")

# Nombre d'animaux comptés selon l'altitude
lm_alt <- lm(obs_secondaire_nb_observation~altitude*secteur, observations1)
ggplot(observations1, aes(x=altitude, y=obs_secondaire_nb_observation, fill=secteur))+
  geom_point(aes(colour=secteur))+
  geom_smooth(method='lm', formula= y~x)
#

sd_all <- aggregate(obs_secondaire_nb_observation~numero_maille, observations1, sd)
colnames(sd_all)[2] <- "sd"
mean_all <- aggregate(obs_secondaire_nb_observation~numero_maille, observations1, mean)
colnames(mean_all)[2] <- "mean"
all <- cbind(sd_all, mean=mean_all[,"mean"])
all$CV <- all$sd/all$mean
ggplot(all, aes(y=sd, x=mean))+
  geom_point()+
  ylab("Ecart-type")+
  xlab("Nombre de marmottes moyen par point de comptage")
min_all <- aggregate(obs_secondaire_nb_observation~numero_maille, observations1, min)
colnames(min_all)[2] <- "min"
max_all <- aggregate(obs_secondaire_nb_observation~numero_maille, observations1, max)
colnames(max_all)[2] <- "max"
min_max <- cbind(max_all, min=min_all[,"min"], mean=mean_all[,"mean"])
min_max$diff <- min_max$max-min_max$min
plot(diff~mean, min_max)
ggplot(min_max)+
  geom_point(aes(x=mean, y=diff))+
  ylab("Ecart entre le nb maximum et minimum de marmottes par point de comptage")+
  xlab("Nb de marmottes moyen par point de comptage")

ggplot(min_max)+
  geom_point(aes(x=mean, y=min))+
  geom_point(aes(x=mean, y=max))+
  ylab("Nb de marmottes min & max par point de comptage")+
  xlab("Nb de marmottes moyen par point de comptage")

# Dérangement ~ altitude
ggplot(observations1)+
  geom_boxplot(aes(x=derangement_new, y=altitude, fill=secteur))+
  labs(fill="Sites d'étude")+
  xlab("")+
  ylab("Altitude")

################################################################################
# Explorations Unmarked ---------------------------------------------------
################################################################################

observations1$prim <- observations1$secteur
observations1[is.element(observations1$secteur, "Bessans"),]$prim <- sample(c("a","b"), dim(observations1[is.element(observations1$secteur, "Bessans"),])[1], replace=T)
observations1[is.element(observations1$secteur, "Les Menuires"),]$prim <- sample(c("c","d"), dim(observations1[is.element(observations1$secteur, "Les Menuires"),])[1], replace=T)
observations1$sec <- observations1$secteur
observations1[is.element(observations1$secteur, "Bessans"),]$sec <- ifelse(observations1[is.element(observations1$secteur, "Bessans"),]$prim == "a", "b", "a")
observations1[is.element(observations1$secteur, "Les Menuires"),]$sec <- ifelse(observations1[is.element(observations1$secteur, "Les Menuires"),]$prim == "c", "d", "c")
# observations1$prim_bis <- observations1$prim
# observations1[is.element(observations1$secteur, "Bessans"),]$prim_bis <- ifelse(observations1[is.element(observations1$secteur, "Bessans"),]$prim=="a", 1, 2)
# observations1[is.element(observations1$secteur, "Les Menuires"),]$prim_bis <- ifelse(observations1[is.element(observations1$secteur, "Les Menuires"),]$prim=="c", 3, 4)
# observations1$sec_bis <- observations1$sec
# observations1[is.element(observations1$secteur, "Bessans"),]$sec_bis <- ifelse(observations1[is.element(observations1$secteur, "Bessans"),]$sec=="a", 1, 2)
# observations1[is.element(observations1$secteur, "Les Menuires"),]$sec_bis <- ifelse(observations1[is.element(observations1$secteur, "Les Menuires"),]$sec=="c", 3, 4)
# prim_bis_obs1 <- as.numeric(observations1$prim_bis)
# sec_bis_obs1 <- as.numeric(observations1$sec_bis)

# observations1$prim <- observations1$secteur
# observations1[is.element(observations1$secteur, "Bessans"),]$prim <- rep("a", dim(observations1[is.element(observations1$secteur, "Bessans"),])[1])
# observations1[is.element(observations1$secteur, "Les Menuires"),]$prim <- rep("c", dim(observations1[is.element(observations1$secteur, "Les Menuires"),])[1])
# observations1$sec <- observations1$secteur
# observations1[is.element(observations1$secteur, "Bessans"),]$sec <- rep("b", dim(observations1[is.element(observations1$secteur, "Bessans"),])[1])
# observations1[is.element(observations1$secteur, "Les Menuires"),]$sec <- rep("d", dim(observations1[is.element(observations1$secteur, "Les Menuires"),])[1])

# Merging some habitats
head(observations1)
table(observations1$typologie_habitat)
observations1$typologie_habitat_new <- observations1$typologie_habitat
observations1[is.element(observations1$typologie_habitat, "prairie_fauche"),]$typologie_habitat_new <- "alpage_herbace"
observations1[is.element(observations1$typologie_habitat, "foret_claire"),]$typologie_habitat_new <- "pre_bois"
table(observations1$typologie_habitat_new)

# Disturbance 0/1
table(observations1$derangement)
observations1$derangement_new <- observations1$derangement
observations1[is.element(observations1$derangement, c("fauche_TA","paturage","randonneur","travaux_divers")),]$derangement_new <- "derangement"
table(observations1$derangement_new)


## Unmarked Double-Observer sur toutes les obs
y_tot <- matrix(NA, dim(observations1), 2, byrow=T)
y_tot[,1] <- observations1$obs_primaire_nb_observation
y_tot[,2] <- observations1$obs_secondaire_nb_observation-observations1$obs_primaire_nb_observation
site.covs_tot <- data.frame(observations1[,c("secteur","condition_meteo","derangement",
                                             "derangement_new",
                                             "typologie_habitat","typologie_habitat_new",
                                             "orientation","x","y","julian_day",
                                             "periode","altitude","numero_maille")])
observer_m    <- as.matrix(
  observations1[,c("prim","sec")]
)


# y_tot <- matrix(NA, dim(observations1), 3, byrow=T)
# y_tot[,1] <- 0
# y_tot[,2] <- observations1$obs_secondaire_nb_observation-observations1$obs_primaire_nb_observation
# y_tot[,3] <- observations1$obs_secondaire_nb_observation
# site.covs_tot <- data.frame(observations1[,c("secteur","condition_meteo","derangement",
#                                              "typologie_habitat","orientation","x","y","julian_day",
#                                              "periode","elevation","numero_maille")])
# 

# Dependant double-observer method : The dependent double-observer point count
# method uses two observers communicating with each other during the count
# (dependent on each other). The primary observer identifies and verbally calls
# out (and points to) all birds seen or heard within defined time and distance
# intervals at the point count location. The secondary observer records all the
# observations dictated by the primary observer; the secondary observer also
# records birds missed by the primary observer.
umf_tot <- unmarkedFrameMPois(y=y_tot, 
                              siteCovs=site.covs_tot,
                              obsCovs = list(observer_m=observer_m),
                              type="depDouble")

# umf_tot <- unmarkedFrameMPois(y=y_tot, 
#                               siteCovs=site.covs_tot,
#                               obsCovs = list(observer_m=observer_m),
#                               type="double")

# Variables qui peuvent influencer la probabilité de détection : observateur, période (matin/soir), météo, habitat, secteur, durée de sortie de la marmotte (à intégrer avec les données de paramétrage)
# Variables qui peuvent influencer l'abondance : secteur, habitat, dérangement

# ~detection ~abondance
fm1_tot <- multinomPois(~1 ~1, umf_tot)
fm2_tot <- multinomPois(~secteur ~1, umf_tot)
fm3_tot <- multinomPois(~secteur+periode ~1, umf_tot)
fm4_tot <- multinomPois(~secteur+condition_meteo ~1, umf_tot)
fm5_tot <- multinomPois(~secteur+derangement_new ~1, umf_tot)
fm13_tot <- multinomPois(~observer_m ~1, umf_tot)

ms_tot_det <- fitList(fm1_tot=fm1_tot, fm2_tot=fm2_tot, fm3_tot=fm3_tot,
                  fm4_tot=fm4_tot, fm5_tot=fm5_tot, 
                  fm13_tot=fm13_tot)

modSel(ms_tot_det)

fm6_tot <- multinomPois(~secteur+derangement_new ~1, umf_tot)
fm7_tot <- multinomPois(~secteur+derangement_new ~secteur, umf_tot)
fm8_tot <- multinomPois(~secteur+derangement_new ~secteur+altitude, umf_tot)
fm9_tot <- multinomPois(~secteur+derangement_new ~secteur+typologie_habitat_new , umf_tot)
fm10_tot <- multinomPois(~observer_m-1 ~secteur+typologie_habitat_new, umf_tot)
fm11_tot <- multinomPois(~observer_m-1 ~secteur-1, umf_tot)

ms_tot_abo <- fitList(fm6_tot=fm6_tot, fm7_tot=fm7_tot, fm8_tot=fm8_tot, fm9_tot=fm9_tot,
                  fm10_tot=fm10_tot, fm11_tot=fm11_tot)
modSel(ms_tot_abo)
# Production de NaN dans fm8

## le problème survient quand la même variable "typologie_habitat" se retrouve à
## la fois sur l'abondance et la détection. Il doit y avoir un problème de
## colinéarité, on ne peut pas séparer l'effet de cette variable sur l'abondance
## et la recapture.

# Prédictions de l'abondance et de la probabilité de détection
newdata_tot <- data.frame(expand_grid(secteur=c("Bessans","Les Menuires"),
                                  typologie_habitat_new=unique(observations1$typologie_habitat_new),
                                  derangement_new=unique(observations1$derangement_new)))

pred_abundance_tot <- predict(fm9_tot, type="state", newdata=newdata_tot)
pred_abundance_tot1 <- cbind(newdata_tot, pred_abundance_tot)

ggplot(pred_abundance_tot1, aes(y=Predicted, x=typologie_habitat_new))+
  geom_bar(stat="identity")+
  facet_grid(.~secteur)+
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1, size=12))+
  ylab("Abondance prédite")+
  xlab("Habitat")

##
newdata_tot1 <- data.frame(expand_grid(secteur=c("Bessans","Les Menuires"),
                                      derangement_new=unique(observations1$derangement_new)))

pred_detection_tot <- predict(fm9_tot, type="det", newdata=newdata_tot1) # Probabilité de détection
pred_detection_tot1 <- cbind(newdata_tot1, pred_detection_tot) # Meilleure détection aux Menuires
pred_detection_tot2 <- aggregate(Predicted~secteur, pred_detection_tot1, mean)

ggplot(pred_detection_tot1, aes(y=Predicted, x=derangement_new))+
  facet_grid(.~secteur)+
  geom_bar(stat="identity", width=0.5)+
  ylab("Probabilité de détection")+
  xlab("Secteur")

## plot(ranef(fm9_tot, K=30), layout=c(10,7), xlim=c(-1, 10))
## plot(ranef(fm9_tot, K=10)@post[, 1,1], ranef(fm12_tot, K=30)@post[, 1,1])
## abline(0, 1)

## 
exp(coef(fm11_tot)[1:2])
# Nb/quart de maille

# Probabilité de détection de chaque observateur 
plogis(coef(fm11_tot)[3:6])

## Résidus
residuals <- as.data.frame(cbind(resid=residuals(fm9_tot)*residuals(fm9_tot), 
                                 x=site.covs_tot$x, 
                                 y=site.covs_tot$y))

residuals_menuires <- residuals[which(residuals$x<6.8),]
residuals_bessans <- residuals[which(residuals$x>6.8),]

ggmap::register_google(key="AIzaSyA53J3oEn4CEPw1xB7Grb2Ei_-AYYdcXes")
mapmarmot_menuires_unmkd <- get_map(location = c(lon = mean(residuals_menuires$x), lat = mean(residuals_menuires$y)), zoom = 13,
                                     maptype = "terrain", scale = 2)
ggmap(mapmarmot_menuires_unmkd) +
  geom_point(data = residuals_menuires, aes(x = x, y = y, size=V1))+
  ggtitle("Menuires")

mapmarmot_bessans_unmkd <- get_map(location = c(lon = mean(residuals_bessans$x), lat = mean(residuals_bessans$y)), zoom = 13,
                             maptype = "terrain", scale = 2)
ggmap(mapmarmot_bessans_unmkd) +
  geom_point(data = residuals_bessans, aes(x = x, y = y, size=V1))+
  ggtitle("Bessans")



## Questions
# - Ajouter effet aléatoire du secteur et de la maille

## Pour le point d'ajouter les effets aléatoires oui, il va bien falloir s'en
## occuper. Par contre, et c'est là que l'on va avoir besoin de JAGS, c'est que
## ces effets aléatoires devraient en être spatialement corrélés. Il y a
## certainement une corrélation spatiale qu'il conviendrait de prendre en
## compte. Voir par exemple:

## https://www.ecography.org/sites/ecography.org/files/appendix/e7853.pdf
## https://www.r-bloggers.com/2014/02/spatial-autocorrelation-of-errors-in-jags/

# - Zero-inflated double observer models ?

## Difficile de répondre pour le moment, il faudrait en passer par un test
## GOF. Je ne suis pas convaincu de la nécessité de complexifier trop le modèle
## (on se retrouve à modéliser le lambda et la proportion de zéro), là où il me
## semble que la strucuration spatiale est davantage critique. Donc pour le
## moment on reste sur du classique, et je préfère que l'on passe plus de temps
## à construire le modèle spatialement explicite.


