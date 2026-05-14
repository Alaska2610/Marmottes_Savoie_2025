library(fitdistrplus)
library(survival)
library(survminer)
library(rjags)
library(bayesplot)

CB <- TRUE
if(CB == TRUE){setwd("./Scripts")}else
{setwd("G:/FDC_Savoie/Marmottes_Savoie/BDD")}

load("Output/parametrage.Rdata")
head(parametrage)
parametrage$status <- 1
parametrage$duree_minutes <- (hour(parametrage$duree)*3600 +
  minute(parametrage$duree)*60 +
  second(parametrage$duree))/60
mean(parametrage$duree_minutes)
median(parametrage$duree_minutes)
#param$pourc <- param$duree/30
min(parametrage$duree_minutes)
max(parametrage$duree_minutes)
colnames(parametrage)[5] <- "secteur"

par(mfrow=c(1,2))
hist(parametrage[parametrage$secteur=="Beaufort",]$duree_minutes, breaks=5, 
     main="Beaufort", xlab="Durée d'observation", ylab="Fréquence",
     xlim=c(0,70), ylim=c(0,25))
hist(parametrage[parametrage$secteur=="Lanslevillard",]$duree_minutes, breaks=10, 
     main="Lanslevillard", xlab="Durée d'observation", ylab="Fréquence",
     xlim=c(0,70), ylim=c(0,25))

t.test(duree_minutes~secteur, parametrage)

aggregate(duree_minutes~secteur, parametrage, ci)
aggregate(duree_minutes~secteur, parametrage, median)

# On approxime le pourcentage moyen de temps passé dehors
# par l'intégrale de la courbe de survie entre 0 et 10 min 
# (aire sous la courbe)

## JAGS 
time <- seq(0.001, 10, by=0.001)
#time <- seq(0.001, 30, by=0.001) #pour faire la courbe d'évolution
h <- 0.001
parametrage_marmottes.string <-"
          model {

                # Extraction des parametrageètres de forme et d'échelle
                global.shape ~ dunif(0, 4)
                global.scale ~ dunif(0, 100)
                for(i in 1:n){
                  parametrage_obs[i] ~ dweib(global.shape, lambda[i])
                  lambda[i] <- pow((1/global.scale), global.shape)
                }
                
                # Modèle de survie
                for(j in 1:length(time)){
                  survie[j] <- exp(-pow((time[j]/global.scale),global.shape))
                }
                
                final <- h/2*(1+survie[10]+2*sum(survie[]))/10
                
                }
"

stops_marmottes_parametrage <- length(parametrage[,1])

parametrage_marmottes.string_model <- textConnection(parametrage_marmottes.string)
jags_marmottes_parametrage <- jags.model(
  parametrage_marmottes.string_model,
  data = list(
    parametrage_obs = parametrage$duree_minutes,
    n = stops_marmottes_parametrage, 
    time = time, 
    h = h
  ),
  inits = list(
    global.shape = 0.05,
    global.scale =  0.05
  ),
  n.chains = 4,
  n.adapt = 1000
)

## Once the model compile, run for a large number of iterations
update(jags_marmottes_parametrage, 5000)

## Save MCMC posteriors for estimated parametrageeters
out_marmottes_parametrage <- coda.samples(
  model = jags_marmottes_parametrage,
  variable.names =
    c(
      #"global.shape",
      #"global.scale",
      #"survie",
      "final"
    ),
  n.iter = 50000
)
save(out_marmottes_parametrage, file="Output/out_marmottes_parametrage.Rdata")

par(mfrow =  c(2, 1))
## Plot histograms
plot(out_marmottes_parametrage, trace = FALSE)
## Plot trace
plot(out_marmottes_parametrage, density = FALSE)
summary(out_marmottes_parametrage)

mean_survie <- summary(out_marmottes_parametrage)$statistics
CI_survie <- summary(out_marmottes_parametrage)$quantiles

save(mean_survie, file="Output/mean_survie.Rdata")
save(CI_survie, file="Output/CI_survie.Rdata")

## Get more parametrageeters
out_marmottes_parametrage1 <- jags.samples(
  model = jags_marmottes_parametrage,
  variable.names =
    c(
      "global.shape",
      "global.scale",
      "survie", 
      "final"
    ),
  n.iter = 10000
)
save(out_marmottes_parametrage1, file="Output/out_marmottes_parametrage1.Rdata")

out_marmottes_parametrage1$global.scale
out_marmottes_parametrage1$global.shape
quantile(out_marmottes_parametrage1$final, c(0.025,0.975))
mean(out_marmottes_parametrage1$final)

