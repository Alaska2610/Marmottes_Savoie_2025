library(lubridate)
library(elevatr)
library(ggplot2)
library(gmodels)
library(unmarked)
library(dplyr)
library(tidyr)
library(terra)
library(sf)
library(rio)

if(CB == TRUE){setwd("../Dataset")}else
{setwd("G:/FDC_Savoie/Marmottes_Savoie_2025/Dataset")}

# Import des bases de données
observations <- import("Dataset/marmottes2025_renamedMB.xlsx", which="2025_Marmotte_Base_de_Donnees")
parametrage <- import("Dataset/marmottes2025_renamedMB.xlsx", which="Parametrage_Both_Sites")

# Ajout de colonnes (jour julien, altitude)
observations$date <- as.Date(observations$date_time, format="%d/%m/%Y")
observations$julian_day <- yday(observations$date)
pt <- st_as_sf(observations, coords = c("x", "y"), crs = 4326)
mnt <- get_elev_raster(pt, z=13)
observations$altitude <- extract(mnt,pt)
observations1 <- observations
observations1$diff_obs <- observations1$obs_secondaire_nb_observation-observations1$obs_primaire_nb_observation

# 1 NA dans la colonne "dérangement" -> on réattribue "aucun"
observations1[is.na(observations1$derangement),]$derangement <- "aucun"

# On "log" la colonne "durée" dans paramétrage
parametrage$log_duree <- log(parametrage$duree)

# On attribue aléatoirement les observateurs primaires et secondaires
prim_bis_obs1 <- sample(c(1:2), dim(observations1)[1], replace=T)
sec_bis_obs1 <- ifelse(prim_bis_obs1==1, 2, 1)

