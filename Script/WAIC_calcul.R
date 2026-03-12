# https://gist.github.com/oliviergimenez/68ad17910a62635ff6a062f8ec34292f

setwd("G:/FDC_Savoie/Marmottes_Savoie/BDD/Rdata")
load("samples.mnull1.Rdata")
load("samples.mnull2.Rdata")
load("samples.msecteurP.Rdata")
load("samples.mderangementP.Rdata")
load("samples.mmeteoP.Rdata")
load("samples.mperiodeP.Rdata")
load("samples.morientationP.Rdata")
load("samples.mderangement.meteoP.Rdata")
load("samples.maltitudeN.Rdata")
load("samples.mhabitatN.Rdata")
load("samples.mderangementP.altitudeN.Rdata")
load("samples.mderangementP.habitatN.Rdata")
load("samples.mderangementP.secteur.altitudeN.Rdata")

samples.mnull1$p_waic <- samples.mnull1$WAIC
samples.mnull1$waic <- samples.mnull1$deviance + samples.mnull1$p_waic
tmp <- sapply(samples.mnull1, sum)
waic.mnull1 <- round(c(waic = tmp[["waic"]], p_waic = tmp[["p_waic"]]),1)
waic.mnull1
# waic : 693.5 ; p_waic : 54.7

samples.mnull2$p_waic <- samples.mnull2$WAIC
samples.mnull2$waic <- samples.mnull2$deviance + samples.mnull2$p_waic
tmp <- sapply(samples.mnull2, sum)
waic.mnull2 <- round(c(waic = tmp[["waic"]], p_waic = tmp[["p_waic"]]),1)
waic.mnull2
# waic : 1131.9 ; p_waic : 5.8

samples.msecteurP$p_waic <- samples.msecteurP$WAIC
samples.msecteurP$waic <- samples.msecteurP$deviance + samples.msecteurP$p_waic
tmp <- sapply(samples.msecteurP, sum)
waic.msecteurP <- round(c(waic = tmp[["waic"]], p_waic = tmp[["p_waic"]]),1)
waic.msecteurP
# waic : 1099.4 ; p_waic : 23.8

samples.mderangementP$p_waic <- samples.mderangementP$WAIC
samples.mderangementP$waic <- samples.mderangementP$deviance + samples.mderangementP$p_waic
tmp <- sapply(samples.mderangementP, sum)
waic.mderangementP <- round(c(waic = tmp[["waic"]], p_waic = tmp[["p_waic"]]),1)
waic.mderangementP
# waic : 1080.2 ; p_waic : 7.3

samples.mmeteoP$p_waic <- samples.mmeteoP$WAIC
samples.mmeteoP$waic <- samples.mmeteoP$deviance + samples.mmeteoP$p_waic
tmp <- sapply(samples.mmeteoP, sum)
waic.mmeteoP <- round(c(waic = tmp[["waic"]], p_waic = tmp[["p_waic"]]),1)
waic.mmeteoP
# waic : 1116 ; p_waic : 39

samples.morientationP$p_waic <- samples.morientationP$WAIC
samples.morientationP$waic <- samples.morientationP$deviance + samples.morientationP$p_waic
tmp <- sapply(samples.morientationP, sum)
waic.morientationP <- round(c(waic = tmp[["waic"]], p_waic = tmp[["p_waic"]]),1)
waic.morientationP
# waic : 1087.1 ; p_waic : 13.3

samples.maltitudeN$p_waic <- samples.maltitudeN$WAIC
samples.maltitudeN$waic <- samples.maltitudeN$deviance + samples.maltitudeN$p_waic
tmp <- sapply(samples.maltitudeN, sum)
waic.maltitudeN <- round(c(waic = tmp[["waic"]], p_waic = tmp[["p_waic"]]),1)
waic.maltitudeN
# waic : 992.6 ; p_waic : 8.8

samples.mhabitatN$p_waic <- samples.mhabitatN$WAIC
samples.mhabitatN$waic <- samples.mhabitatN$deviance + samples.mhabitatN$p_waic
tmp <- sapply(samples.mhabitatN, sum)
waic.mhabitatN <- round(c(waic = tmp[["waic"]], p_waic = tmp[["p_waic"]]),1)
waic.mhabitatN
# waic : 954.7 ; p_waic : 15.1

samples.mderangement.meteoP$p_waic <- samples.mderangement.meteoP$WAIC
samples.mderangement.meteoP$waic <- samples.mderangement.meteoP$deviance + samples.mderangement.meteoP$p_waic
tmp <- sapply(samples.mderangement.meteoP, sum)
waic.mderangement.meteoP <- round(c(waic = tmp[["waic"]], p_waic = tmp[["p_waic"]]),1)
waic.mderangement.meteoP
# waic : 1101.1 ; p_waic: 41.7

samples.mderangementP.altitudeN$p_waic <- samples.mderangementP.altitudeN$WAIC
samples.mderangementP.altitudeN$waic <- samples.mderangementP.altitudeN$deviance + samples.mderangementP.altitudeN$p_waic
tmp <- sapply(samples.mderangementP.altitudeN, sum)
waic.mderangementP.altitudeN <- round(c(waic = tmp[["waic"]], p_waic = tmp[["p_waic"]]),1)
waic.mderangementP.altitudeN
# waic : 992.2 ; p_waic : 10.4

samples.mderangementP.habitatN$p_waic <- samples.mderangementP.habitatN$WAIC
samples.mderangementP.habitatN$waic <- samples.mderangementP.habitatN$deviance + samples.mderangementP.habitatN$p_waic
tmp <- sapply(samples.mderangementP.habitatN, sum)
waic.mderangementP.habitatN <- round(c(waic = tmp[["waic"]], p_waic = tmp[["p_waic"]]),1)
waic.mderangementP.habitatN
# waic : 953.5 ; p_waic : 16.8

samples.mderangementP.secteur.altitudeN$p_waic <- samples.mderangementP.secteur.altitudeN$WAIC
samples.mderangementP.secteur.altitudeN$waic <- samples.mderangementP.secteur.altitudeN$deviance + samples.mderangementP.secteur.altitudeN$p_waic
tmp <- sapply(samples.mderangementP.secteur.altitudeN, sum)
waic.mderangementP.secteur.altitudeN <- round(c(waic = tmp[["waic"]], p_waic = tmp[["p_waic"]]),1)
waic.mderangementP.secteur.altitudeN
# waic : 991.5 ; p_waic : 10

samples.mderangementP.secteur.altitudeN.ssCT$p_waic <- samples.mderangementP.secteur.altitudeN.ssCT$WAIC
samples.mderangementP.secteur.altitudeN.ssCT$waic <- samples.mderangementP.secteur.altitudeN.ssCT$deviance + samples.mderangementP.secteur.altitudeN.ssCT$p_waic
tmp <- sapply(samples.mderangementP.secteur.altitudeN.ssCT, sum)
waic.mderangementP.secteur.altitudeN.ssCT <- round(c(waic = tmp[["waic"]], p_waic = tmp[["p_waic"]]),1)
waic.mderangementP.secteur.altitudeN.ssCT
# waic : 985 ; p_waic : 10.6

