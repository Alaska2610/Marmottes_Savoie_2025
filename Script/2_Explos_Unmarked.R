if(CB == TRUE){setwd("../Dataset")}else
{setwd("G:/FDC_Savoie/Marmottes_Savoie_2025/")}

################################################################################
# Visualisation des données -----------------------------------------------
################################################################################
  ## Observations
# Création des tables pour les graphiques ci-après
tab_secteur <- table(observations1$secteur)
tab_maille_secteur <- aggregate(maille~secteur, observations1, function(x){length(unique(x))})
# 35 mailles sur chaque secteur

tab_meteo_secteur <- as.data.frame(table(observations1$meteo, observations1$secteur))
tab_habitat_secteur <- as.data.frame(table(observations1$habitat, observations1$secteur))
tab_habitat_derangement <- as.data.frame(table(observations1$derangement, observations1$secteur))

moy_obs_prim <- aggregate(nb_obs_prim~secteur, observations1, ci)
moy_obs_sec <- aggregate(nb_obs_sec~secteur, observations1, ci)

nb_obs_prim_sum <- aggregate(nb_obs_prim~secteur, observations1, sum)
# 42 à Beaufort, 119 à Lanslevillard
colnames(nb_obs_prim_sum)[2] <- "nb_obs"
nb_obs_prim_sum$observateur <- "primaire"
nb_obs_sec_sum <- aggregate(nb_obs_sec~secteur, observations1, sum)
# 47 à Beaufort, 123 à Lanslevillard
colnames(nb_obs_sec_sum)[2] <- "nb_obs"
nb_obs_sec_sum$observateur <- "secondaire"
nb_obs_prim_sec <- rbind(nb_obs_prim_sum, nb_obs_sec_sum)

nb_obs_max_prim <- aggregate(nb_obs_prim~secteur, observations1, max) # Max 6 ou 4 marmottes
nb_obs_min_prim <- aggregate(nb_obs_prim~secteur, observations1, min) # Min 0 marmottes

nb_periode_secteur <- aggregate(nb_obs_sec~secteur*session, observations1, sum)

nb_obs_mailles <- aggregate(nb_obs_sec~secteur*session*maille, observations1, sum)
nb_obs_mailles1 <- aggregate(nb_obs_sec~secteur*maille, nb_obs_mailles, mean)
nb_obs_mailles1$densite <- nb_obs_mailles1$nb_obs_sec/(pi*150*150/10000)
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

ggplot(observations1, aes(nb_obs_sec))+
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
lm_alt <- lm(nb_obs_sec~altitude*secteur, observations1)
ggplot(observations1, aes(x=altitude, y=nb_obs_sec, fill=secteur))+
  geom_point(aes(colour=secteur))
#

# Nb d'animaux comptés selon l'habitat
ggplot(observations1, aes(x=habitat, y=nb_obs_sec, fill=secteur))+
  geom_point(aes(colour=secteur))
#

# Dérangement ~ altitude
ggplot(observations1)+
  geom_boxplot(aes(x=derangement, y=altitude, fill=secteur))+
  labs(fill="Sites d'étude")+
  xlab("")+
  ylab("Altitude")

################################################################################
# Explorations Unmarked ---------------------------------------------------
################################################################################

# Disturbance 0/1
table(observations1$derangement)
observations1$derangement_new <- observations1$derangement
observations1[is.element(observations1$derangement, c("Agricole","Autre","Touristique")),]$derangement_new <- "Derangement"
table(observations1$derangement_new)


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

# Dependant double-observer method
umf_tot <- unmarkedFrameMPois(y=y_tot, 
                              siteCovs=site.covs_tot,
                              obsCovs = list(observer_m=observer_m),
                              type="depDouble")

# ~detection ~abondance
fm1_tot <- multinomPois(~1 ~1, umf_tot)
fm2_tot <- multinomPois(~secteur ~1, umf_tot)
fm3_tot <- multinomPois(~secteur+session ~1, umf_tot)
fm4_tot <- multinomPois(~secteur+meteo ~1, umf_tot)
fm5_tot <- multinomPois(~secteur+derangement_new ~1, umf_tot)
fm13_tot <- multinomPois(~observer_m ~1, umf_tot)

ms_tot_det <- fitList(fm1_tot=fm1_tot, fm2_tot=fm2_tot, fm3_tot=fm3_tot,
                  fm4_tot=fm4_tot, fm5_tot=fm5_tot, 
                  fm13_tot=fm13_tot)

modSel(ms_tot_det)

fm6_tot <- multinomPois(~secteur+derangement_new ~1, umf_tot)
fm7_tot <- multinomPois(~secteur+derangement_new ~secteur, umf_tot)
fm8_tot <- multinomPois(~secteur+derangement_new ~secteur+altitude, umf_tot)
fm9_tot <- multinomPois(~secteur+derangement_new ~secteur+habitat , umf_tot)
fm10_tot <- multinomPois(~observer_m-1 ~secteur+habitat, umf_tot)
fm11_tot <- multinomPois(~observer_m-1 ~secteur-1, umf_tot)

ms_tot_abo <- fitList(fm6_tot=fm6_tot, fm7_tot=fm7_tot, fm8_tot=fm8_tot, fm9_tot=fm9_tot,
                  fm10_tot=fm10_tot, fm11_tot=fm11_tot)
modSel(ms_tot_abo)
# Production de NaN dans fm8

# Prédictions de l'abondance et de la probabilité de détection
newdata_tot <- data.frame(expand_grid(secteur=c("Beaufort","Lanslevillard"),
                                  habitat=unique(observations1$habitat),
                                  derangement_new=unique(observations1$derangement_new)))

pred_abundance_tot <- predict(fm9_tot, type="state", newdata=newdata_tot)
pred_abundance_tot1 <- cbind(newdata_tot, pred_abundance_tot)

ggplot(pred_abundance_tot1, aes(y=Predicted, x=habitat))+
  geom_bar(stat="identity")+
  facet_grid(.~secteur)+
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1, size=12))+
  ylab("Abondance prédite")+
  xlab("Habitat")

##
newdata_tot1 <- data.frame(expand_grid(secteur=c("Beaufort","Lanslevillard"),
                                      derangement_new=unique(observations1$derangement_new)))

pred_detection_tot <- predict(fm9_tot, type="det", newdata=newdata_tot1) # Probabilité de détection
pred_detection_tot1 <- cbind(newdata_tot1, pred_detection_tot) # Meilleure détection aux Menuires
pred_detection_tot2 <- aggregate(Predicted~secteur, pred_detection_tot1, mean)

ggplot(pred_detection_tot1, aes(y=Predicted, x=derangement_new))+
  facet_grid(.~secteur)+
  geom_bar(stat="identity", width=0.5)+
  ylab("Probabilité de détection")+
  xlab("Secteur")

## 
exp(coef(fm11_tot)[1:2])
# Nb/quart de maille

# Probabilité de détection de chaque observateur 
plogis(coef(fm11_tot)[3:6])


