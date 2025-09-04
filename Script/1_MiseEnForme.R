library(lubridate)
library(elevatr)
library(ggplot2)
library(gmodels)
library(unmarked)
library(dplyr)
library(tidyr)
library(terra)
library(sf)

if(CB == TRUE){setwd("../BDD")}else
{setwd("G:/FDC_Savoie/Marmottes_Savoie/BDD")}

# Import des bases de données
observations <- read.csv("Observations_renamedMB.csv", h=T, sep=";")
parametrage <- read.csv("Parametrage_marmottes.csv", h=T, sep=";", dec=",")

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

