# https://gist.github.com/oliviergimenez/68ad17910a62635ff6a062f8ec34292f

setwd("/Volumes/Boulot/FDC_Savoie/2025")
load("samples.m1.Rdata")
load("samples.m2.Rdata")
load("samples.m3.Rdata")

samples.m1$p_waic <- samples.m1$WAIC
samples.m1$waic <- samples.m1$deviance + samples.m1$p_waic
tmp <- sapply(samples.m1, sum)
waic.m1 <- round(c(waic = tmp[["waic"]], p_waic = tmp[["p_waic"]]),1)
waic.m1
# waic : 823.1

samples.m2$p_waic <- samples.m2$WAIC
samples.m2$waic <- samples.m2$deviance + samples.m2$p_waic
tmp <- sapply(samples.m2, sum)
waic.m2 <- round(c(waic = tmp[["waic"]], p_waic = tmp[["p_waic"]]),1)
waic.m2
# waic : 821.2

samples.m3$p_waic <- samples.m3$WAIC
samples.m3$waic <- samples.m3$deviance + samples.m3$p_waic
tmp <- sapply(samples.m3, sum)
waic.m3 <- round(c(waic = tmp[["waic"]], p_waic = tmp[["p_waic"]]),1)
waic.m3
# waic : 772.9
