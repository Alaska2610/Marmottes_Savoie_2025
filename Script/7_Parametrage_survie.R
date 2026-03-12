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
     xlim=c(0,35), ylim=c(0,22))
hist(parametrage[parametrage$secteur=="Lanslevillard",]$duree_minutes, breaks=10, 
     main="Lanslevillard", xlab="Durée d'observation", ylab="Fréquence",
     xlim=c(0,35), ylim=c(0,22))

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
#cbind(apply(out_marmottes_parametrage1$survie, 1, mean), time)
quantile(out_marmottes_parametrage1$final, c(0.025,0.975))
mean(out_marmottes_parametrage1$final)

survie_mean_ci <- cbind(mean_survie[4:30003,], CI_survie[4:30003,], time)
colnames(survie_mean_ci)[5] <- "IC2.5"
colnames(survie_mean_ci)[9] <- "IC97.5"
ess <- cbind(apply(out_marmottes_parametrage1$survie, 1, mean), 
             time)
colnames(ess) <- c("survie","time")
ggplot()+
  geom_histogram(data=parametrage, aes(x=duree_minutes))+
  ylab("Fréquence") + xlab("Durée d'observation des marmottes")

ggplot(data=survie_mean_ci, aes(x=time))+
  geom_line(aes(y=Mean))+
  geom_ribbon(aes(ymin=IC2.5, ymax=IC97.5), alpha=0.2)+
  ylab("Probabilité de disponibilité") + xlab("Durée d'observation des marmottes")

ess1 <- as.data.frame(ess)
ess1$mult <- ess1$survie*ess1$time
0.001/2*(1+2*sum(ess1[1:10000,"survie"])+0.135)


## Vérifications
fit.all <- fitdist(parametrage$duree_minutes, "weibull") 
fit.beaufort <- fitdist(parametrage[parametrage$secteur=="Beaufort",]$duree_minutes, "weibull")
fit.lanslevillard <- fitdist(parametrage[parametrage$secteur=="Lanslevillard",]$duree_minutes, "weibull")

summary(fit.all) # shape=0.87 scale=4.3
summary(fit.beaufort) # shape=0.84 scale=3.77
summary(fit.lanslevillard) # shape=0.9 scale=4.85

plot(fit.all)
plot(fit.beaufort)
plot(fit.lanslevillard)

pweibull(2.27, shape=0.868, scale=4.292) # fonction de distribution
exp(-(2.27/4.292)^0.868) # survie. 57% de survie après 2.27 min.

1-exp(-(2.27/4.292)^0.868) # p de ne pas avoir de marmottes visible après 2.27 min
1-0.57^(10/2.27) 
1-0.36^(10/4.5) 

length(which(parametrage$duree_minutes<=10))/dim(parametrage)[1]

# Probability of detecting a marmot within 10 minutes
1-exp(10/60*log(0.43))

fit <- survfit(Surv(duree_minutes, status) ~ 1, data = parametrage)
print(fit) # 50% de survie à 2.27 min
summary(fit) # survie=0.1 à 10 minutes, 10% de chances qu'une 
# marmotte soit visible plus de 10 minutes

ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw())
#fun="event"

y <- matrix(nrow=length(seq(0,1,by=0.01)), ncol=2)
y[,2] <- seq(0,1,by=0.01)
for(i in 1:length(seq(0,1,by=0.01))){
  y[i,1] <- exp(10/60*log(seq(0,1,by=0.01)[i]))
}
plot(y[,1]~y[,2])



# Aire sous la courbe à 10 minutes
scale <- 9.48616   # lambda
shape <- 0.903223  # k
t_max <- 10         # minutes

# fonction de survie Weibull
S <- function(t) {
  exp(-(t/scale)^shape)
}

# intégration numérique
auc <- integrate(S, lower = 0, upper = t_max)$value
auc/t_max


