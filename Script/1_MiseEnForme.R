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

if(CB == TRUE){setwd("../")}else
{setwd("G:/FDC_Savoie/Marmottes_Savoie_2025/")}

# Import des bases de données
observations <- import("Dataset/marmottes2025_renamedMB.xlsx", which="2025_Marmotte_Base_de_Donnees")
parametrage <- import("Dataset/marmottes2025_renamedMB.xlsx", which="Parametrage_Both_Sites")

# Ajout de colonnes (jour julien, altitude)
observations$julian_day <- yday(observations$date)
pt <- st_as_sf(observations, coords = c("longitude", "latitude"), crs = 4326)
mnt <- get_elev_raster(pt, z=13)
observations$altitude <- extract(mnt,pt)
observations1 <- observations
observations1$diff_obs <- observations1$nb_obs_sec-observations1$nb_obs_prim

save(observations1, file="Output/observations1.Rdata")
