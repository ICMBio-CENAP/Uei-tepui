## Script to analyse data from Uei-tepui, Monte Roraima National Park
# Elildo Carvalho Jr @ ICMBio/CENAP, 2020-05-10


##----- 1 - Load libraries-----
# obs: check which ones are really needed, same for source files below
library(TeachingDemos)
library(lubridate)
library(ggplot2)
library(dplyr)
library(reshape2)
library(here)
library(tidyverse)
library(abind)


##----- 2 - Source files-----
source(here("src", "functions-to-prepare-data.R"))


##-----3 - Read and fix data-----
wiiTiwi$Camera.Trap.Name <- as.factor(substr(wiiTiwi$deployment_id, start = 1, stop = 7))

deployments <- read.csv("/media/elildojr/Dados/r/monte_roraima/deployments.csv")
#View(deployments)

df2 <- merge(wiiTiwi, deployments, "deployment_id", sort = TRUE)
df2$Start.Date <- as.Date(df2$start_date)
df2$End.Date <- as.Date(df2$end_date)

# remove unidentified species, fix name
df2 <- subset(df2, bin != "cf. Nasuella" & bin != "Mazama sp" & bin != "Cebus olivaceus")
df2$bin[df2$bin == "Didelphis "] <- "Didelphis sp"
df2$bin[df2$bin == "Didelphis imperfecta"] <- "Didelphis sp"
df2$bin <- factor(df2$bin)

unique(df2$Camera.Trap.Name)
df2$Camera.Trap.Name <- factor(df2$Camera.Trap.Name)

##-----4 - Prepare data for Dorazio et al. multimodel code -----

# First, read/prepare occurrence data:

# The detection/non-detection data must be in a three dimensional 
# array X where the first dimension, j, is the point; the second 
# dimension, k, is the rep; and the last dimension, i, is the species

x <- f.matrix.creator.elildo(df2)
#x <- f.matrix.creator.elildo(jamari2017)
nrow(x[[1]]); ncol(x[[2]]) # check dimensions of individual matrices, here 120 rows by 24 cols

# convert X into array
X <- array(unlist(x), dim = c(nrow(x[[1]]), ncol(x[[1]]), length(x)))
dimnames(X) <- list(c=rownames(x[[1]]), seq(1:ncol(x[[2]])), names(x)) # use a seq instead of dates for the reps
#names(X3) <- c(Point, Rep, Species)
str(X)
X[,,1] # to check individual matrices

#---------
#Create all zero encounter histories to add to the detection array X 
#as part of the data augmentation to account for additional 
#species (beyond the n observed species). 

#--------
# Do this with a part of the f.matrix.creator function:

#list of sanpling units
data <- df2
data$Sampling.Unit.Name <- data$Camera.Trap.Name
cams <- unique(data$Sampling.Unit.Name)
cams <- sort(cams)
rows <- length(cams)
species <- unique(data$bin)
#start and end dates of sampling periods
#data <- dplyr::filter(data, Sampling.Period == year)
min <- min(data$Start.Date)
max <- max(data$End.Date)
cols <- max - min + 1
#cols <- ncol(x[[2]]) # a failed test by Elildo
#sampling period
date.header <- seq(from=min,to=max, by="days")
#date.header <- seq(from=1, to=ncol(x[[2]]))  # a failed test by Elildo
mat<-matrix(NA,rows,cols,dimnames=list(cams,as.character(date.header)))
#for all cameras, determine the open and close date and mark in the matrix
start.dates<-tapply(as.character(data$Start.Date),data$Sampling.Unit.Name,unique)
#start.dates <- start.dates[start.dates != ""]
nms<-names(start.dates)
start.dates<-ymd(start.dates)
names(start.dates)<-nms
end.dates<-tapply(as.character(data$End.Date),data$Sampling.Unit.Name,unique)
#end.dates <- end.dates[end.dates != ""]
end.dates<-ymd(end.dates)
names(end.dates)<-nms
#outline the sampling periods for each camera j
for(j in 1:length(start.dates)){
  #for each camera beginning and end of sampling
  low<-which(date.header==start.dates[j])
  hi<-which(date.header==end.dates[j])
  if(length(low)+length(hi)>0){
    indx<-seq(from=low,to=hi)
    mat[j,indx]<-0
  }}
mat.nas<-is.na(mat)
sum.nas<-apply(mat.nas,2,sum)
indx.nas<-which(sum.nas==rows)
if(length(indx.nas)>0){
  mat<-mat[,-indx.nas]
}
  # reduce the size of the matrix
  mat <- f.shrink.matrix.elildo(mat)
  #res<-c(res,list(mat))
  #return the matrix to its original form
  #mat<-mat.template

X.zero <- mat
str(X.zero)

#--------

#nzeroes is the number of all zero encounter histories to be added
nzeroes = 2*length(species) #100-length(species)
#Xaug is the augmented version of X.  The first n species were actually observed
#and the n+1 through nzeroes species are all zero encounter histories  

# create Xaug with the abind package and a loop
#library(abind)
Xaug <- X
for(i in 1:nzeroes){
  Xaug <- abind(Xaug, X.zero, along = 3)
}
dimnames(Xaug) <- list(c=rownames(x[[1]]), seq(1:ncol(x[[2]])), c(names(x), seq(1:nzeroes))) # use a seq instead of dates for the reps
str(Xaug)


#--------


#Find the number of unique species
uspecies = as.character(names(x))
#n is the number of observed species
n=length(uspecies)


#Find the number of unique sampling locations
upoints = rownames(x[[1]])
#J is the number of sampled points
J=length(upoints)

#K is a vector of length J indicating the number of reps at each point j  
KK <- X.zero
a=which(KK==0); KK[a] <- 1
K=apply(KK,1,sum, na.rm=TRUE)
K=as.vector(K)

## Covariates

# Put covars in the same order as Xaug
elevation <- elevation %>% arrange(Camera.Trap.Name)

# filter by selecting only Cams of interest
elevation <- subset(elevation, deployment_id %in% unique(df2$deployment_id)); elevation$deployment_id <- factor(elevation$deployment_id)

# elevation covar
elevation <- elevation$elevation
m.elevation <- mean(elevation, na.rm=TRUE)
sd.elevation <- sqrt(var(elevation[1:length(elevation)], na.rm=TRUE))
elevation <- (elevation-m.elevation) /  sd.elevation

# dates covar (just the sequence of occasions as a detection covariate)
dates <- c(1:11)
dates <-  matrix(rep(dates, J), nrow=J, ncol=11, byrow=T)
dates <- as.matrix(dates)
m.dates <- mean(dates, na.rm=TRUE)
sd.dates <- sqrt(var(dates[1:length(dates)], na.rm=TRUE))
dates <- (dates-m.dates) /  sd.dates
dates <- as.matrix(dates)

# bundle data
sp.data = list(n=n, nzeroes=nzeroes, J=J, K=K, X=Xaug, elevation=elevation, dates=dates)

#----------------------------------------------------------------------------------------

# heatmap plot

df2$bin <- factor(df2$bin)
site.spp.Matrix <- with(df2, table(bin, Camera.Trap.Name))

df_site.spp.Matrix <- melt(site.spp.Matrix);
colnames(df_site.spp.Matrix)[1] <- "species";
colnames(df_site.spp.Matrix)[2] <- "site";
colnames(df_site.spp.Matrix)[3] <- "detections";
head(df_site.spp.Matrix);

detections_heatmap_plot <-
  ggplot(df_site.spp.Matrix, aes(site, species)) +
  geom_tile(aes(fill = detections), colour = "white") +
  scale_fill_gradient(low = "white", high = "black") +
  labs(x = "Sítios (1 to 17)", y = "Número da espécie") +
  scale_x_discrete(expand = c(0, -5), breaks=(5 * (1:4))) +
  #scale_y_discrete(expand = c(0, -5), breaks=(5 * (1:5))) +
  scale_y_discrete(labels= seq(1,length(uspecies))) +
  ggtitle("Detecções das espécies nos sítios ao longo das visitas")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5));

jpeg("/media/elildojr/Dados/r/monte_roraima/heatmap.jpg", width = 800, height = 600) # Open jpeg file
plot(detections_heatmap_plot)
dev.off()
#--------------------------

# check effort etc for Table 1 of paper
df2$effort <- df2$End.Date-df2$Start.Date
table(df2$effort)

checkCAMS <- tibble(Camera.Trap.Name=df2$Camera.Trap.Name, effort=df2$effort, Start.date=df2$Start.Date, End.date=df2$End.Date)
checkCAMS <- arrange(checkCAMS, effort, End.date, Start.date)
checkCAMS <- distinct(checkCAMS, Camera.Trap.Name, effort, Start.date, End.date)
print(checkCAMS, n=Inf)

mean(checkCAMS$effort)
range(checkCAMS$effort)
min(checkCAMS$Start.date)
max(checkCAMS$End.date)


# 


#----- 1 - Load libraries -----
library(R2jags)

#----- 2 - Write the model -----
# Specify the model in JAGS language

sink("/media/elildojr/Dados/r/monte_roraima/mod1.txt")
cat("
   model{

## Define prior distributions for community-level model parameters

omega ~ dunif(0,1)                                # inclusion probability

# coefficients
mu.alpha.psi ~ dnorm(0, 0.001)
tau.alpha.psi ~ dgamma(0.1,0.1)
mua1 ~ dnorm(0, 0.001)
tau.a1 ~ dgamma(0.1,0.1)

mu.alpha.p ~ dnorm(0, 0.001)
tau.alpha.p ~ dgamma(0.1,0.1)
mub1 ~ dnorm(0, 0.001)
tau.b1 ~ dgamma(0.1,0.1) 


for (i in 1:(n+nzeroes)) {

# Create priors for species i from the community level prior distributions

    w[i] ~ dbern(omega)                             # inclusion indicators
    alpha.psi[i] ~ dnorm(mu.alpha.psi, tau.alpha.psi)
    a1[i] ~ dnorm(mua1, tau.a1)                     # elevation on psi
    alpha.p[i] ~ dnorm(mu.alpha.p, tau.alpha.p)
    b1[i] ~ dnorm(mub1, tau.b1)                     # occasion on p

# Create a loop to estimate the Z matrix (true occurrence for species i at point j)
   for (j in 1:J) {
       logit(psi[j,i]) <- alpha.psi[i] + a1[i]*elevation[j]

      mu.psi[j,i] <- psi[j,i]*w[i]
      Z[j,i] ~ dbern(mu.psi[j,i])

# Create a loop to estimate detection for species i at point j during sampling period k.      
   for (k in 1:K[j]) {  
      logit(p[j,k,i]) <-  alpha.p[i] + b1[i]*dates[j,k]
      
       mu.p[j,k,i] <- p[j,k,i]*Z[j,i]              # can only be detected if Z=1
       X[j,k,i] ~ dbern(mu.p[j,k,i])

    }#k
  }#j
}#i

# Derived quantities:

# Sum all species observed (n) and unobserved species (n0) to find the total estimated richness
n0 <- sum(w[(n+1):(n+nzeroes)])
N <- n + n0
#N <- sum(w[])

# Create a loop to determine point level richness for the whole community and for subsets of interest
for(j in 1:J){
Nsite[j]<- inprod(Z[j,1:(n+nzeroes)],w[1:(n+nzeroes)])
  }
}", fill=TRUE); sink()


sp.params = c('alpha.psi', 'a1', 'alpha.p', 'b1', 'omega','N', 'Nsite')

# Specify the initial values
sp.inits = function() {
  omegaGuess = runif(1, 0, (n/(n+nzeroes)))/2 # divided by 2 to make it even smaller
  psi.meanGuess = runif(1, 0.1, .5)
  list(omega=omegaGuess, w=c(rep(1, n), rbinom(nzeroes, size=1, prob=omegaGuess)),
       Z = matrix(c(rep(1, n*J),rep(0, nzeroes*J)), nrow=J, ncol=(n+nzeroes)),
       alpha.psi=rnorm(n+nzeroes), a1=rnorm(n+nzeroes),
       alpha.p=rnorm(n+nzeroes), b1=rnorm(n+nzeroes) )
}

# run model in jags
fit <- jags(data=sp.data, inits=sp.inits, parameters.to.save=sp.params, n.chains=3, n.iter=10000, n.burnin=5000, n.thin=50, model.file="/media/elildojr/Dados/r/monte_roraima/mod1.txt")
saveRDS(fit, "/media/elildojr/Dados/r/monte_roraima/fit.RDS")

print(fit)
traceplot(fit) # not working with the rds

#-------------------------------------------

# Species richness: metacommunity
summary(fit$BUGSoutput$sims.list$N)

# Species richness: site-level
fit$BUGSoutput$sims.list$Nsite # this is the posterior for site level richness
apply(fit$BUGSoutput$sims.list$Nsite, 2, mean) # this is the mean richness per site


#---------------------------------------------------------------
# See estimates of total richness (N) and estimates of richness at each of the 
#J sampling locations (Nsite)
N <- fit$BUGSoutput$sims.list$N # jags
mean(N); summary(N); quantile(N, probs=c(0.025,0.975))
plot(table(N), xlab="No. mammal species", ylab="Frequency", las=1)

# save as jpeg
jpeg("/media/elildojr/Dados/r/monte_roraima/SpeciesRichness.jpg", width = 800, height = 600) # Open jpeg file
plot(table(N), xlab="No. espécis de mamíferos", ylab="Frequência", las=1)
dev.off()


# richness per site, mean and SD
Nsite = fit$BUGSoutput$sims.list$Nsite
site.richness.matrix = cbind(apply(Nsite,2,mean), apply(Nsite,2,sd))
rownames(site.richness.matrix) <- cams; colnames(site.richness.matrix) <- c("mean", "sd")
site.richness.matrix
mean(site.richness.matrix[,1]); summary(site.richness.matrix[,1]); quantile(site.richness.matrix[,1], probs=c(0.025,0.975))

# glm riqueza vs elevação
glmRiquezaElev <- glm(apply(Nsite,2,mean)~elevation)
summary(glmRiquezaElev)

# same plot with jitter, linear trend and confidence interval
datap1 <- data.frame(cbind(elevation, apply(Nsite,2,mean))); names(datap1) <- c('elevação', 'riqueza.sítio')
p1 <- ggplot(datap1, aes(x=elevation, y=riqueza.sítio)) + 
  geom_point(shape = 19, alpha = 0.5, size = 4) +
  geom_smooth(method=lm , color="black", fill="grey", se=TRUE) +
  theme(axis.title.x = element_text(size=20, margin=margin(t=10, r=0, b=0, l=0)),
        axis.title.y = element_text(size=20, margin=margin(t=0, r=10, b=0, l=0))) +
  theme(panel.grid.minor = element_blank(), panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
  theme(axis.text.x = element_text(size=18), axis.text.y = element_text(size=18)) +
  theme(legend.title=element_text(size=18), legend.text=element_text(size=18)) +
  xlab("Elevação (padronizado)") + ylab("Riqueza de espécies por câmera") +
  ylim(0, 20)
p1

# save as jpeg
jpeg("/media/elildojr/Dados/r/monte_roraima/RichnessVselevacao.jpg", width = 800, height = 600) # Open jpeg file
p1
dev.off()

#---------------------------------------------------------------
# Check effect da elevação sobre especies individuais

## a1 (elevation)
a1bySpecies <- fit$BUGSoutput$sims.list$a1[,1:15]; colnames(a1bySpecies) <- species
a1Mean <- apply(a1bySpecies, 2, mean)
a1Low <- apply(a1bySpecies, 2, quantile, probs=0.05)
a1Upper <- apply(a1bySpecies, 2, quantile, probs=0.95)
a1DF <- rbind(a1Mean, a1Low, a1Upper); colnames(a1bySpecies) <- species #rownames(a1DF) <- cams; colnames(a1DF) <- c("mean", "lowerCI", "upperCI")

# transpose a1DF to check if a1 was significant (confint does not overlap 0) for any sppecies
a1bySpecies <- as.data.frame(t(a1DF[,1:ncol(a1DF)])) # as.data.frame(t(a1bySpecies[,3:ncol(a1bySpecies)])) in the previous version
colnames(a1bySpecies) <- c("mean", "lowerCI", "upperCI")
a1bySpecies <- a1bySpecies[order(a1bySpecies$mean),] # reorder based on mean response


# plot mean and 95% CI for a1 to check visually
par(mar=c(9,3,1.5,1.5))
plot(0,0, xlim=range(c(1:nrow(a1bySpecies))), ylim=range(c(a1bySpecies$lowerCI-1,a1bySpecies$upperCI+1)), type = "n", xaxt="n",
     pch=19, xlab="", ylab="", las=2)
abline(h = 0,col = "grey")
arrows(x0 = 1:nrow(a1bySpecies), x1 = 1:nrow(a1bySpecies), y0=a1bySpecies$lowerCI, y1=a1bySpecies$upperCI, length=0, angle=90, code=3)
points(1:nrow(a1bySpecies), a1bySpecies$mean, pch = 20, cex = 1)
axis(side=1,at=1:nrow(a1bySpecies),labels=row.names(a1bySpecies), las=2, cex.axis=0.8)

# same plot rotated (write as a function to save jpeg in next step)
a1Plot <- function(a1bySpecies){
  par(mar=c(3,11,1.5,1.5))
  plot(0,0, ylim=range(c(1:nrow(a1bySpecies))), xlim=range(c(a1bySpecies$lowerCI-1,a1bySpecies$upperCI+1)), type = "n", yaxt="n",
       pch=19, xlab="", ylab="", las=1)
  abline(v = 0,col = "grey")
  arrows(y0 = 1:nrow(a1bySpecies), y1 = 1:nrow(a1bySpecies), x0=a1bySpecies$lowerCI, x1=a1bySpecies$upperCI, length=0, angle=90, code=3)
  points(a1bySpecies$mean, 1:nrow(a1bySpecies), pch = 20, cex = 1)
  axis(side=2,at=1:nrow(a1bySpecies),labels=row.names(a1bySpecies), las=2, cex.axis=1)
}
a1Plot(a1bySpecies)

jpeg("/media/elildojr/Dados/r/monte_roraima/a1.jpg", width = 450, height = 600) # Open jpeg file
a1Plot(a1bySpecies)
dev.off()

#---------------------------------------------------------------
# Check effect of dates on species detection

## beta1.p (dates)
beta1.pbySpecies <- fit$BUGSoutput$sims.list$b1[,1:15]; colnames(beta1.pbySpecies) <- species
beta1.pMean <- apply(beta1.pbySpecies, 2, mean)
beta1.pLow <- apply(beta1.pbySpecies, 2, quantile, probs=0.05)
beta1.pUpper <- apply(beta1.pbySpecies, 2, quantile, probs=0.95)
beta1.pDF <- rbind(beta1.pMean, beta1.pLow, beta1.pUpper); colnames(beta1.pbySpecies) <- species #rownames(beta1.pDF) <- cams; colnames(beta1.pDF) <- c("mean", "lowerCI", "upperCI")

# transpose beta1.pDF to check if beta1.p was significant (confint does not overlap 0) for any sppecies
beta1.pbySpecies <- as.data.frame(t(beta1.pDF[,1:ncol(beta1.pDF)])) # as.data.frame(t(beta1.pbySpecies[,3:ncol(beta1.pbySpecies)])) in the previous version
colnames(beta1.pbySpecies) <- c("mean", "lowerCI", "upperCI")
beta1.pbySpecies <- beta1.pbySpecies[order(beta1.pbySpecies$mean),] # reorder based on mean response

# plot mean and 95% CI for beta1.p to check visually
par(mar=c(9,3,1.5,1.5))
plot(0,0, xlim=range(c(1:nrow(beta1.pbySpecies))), ylim=range(c(beta1.pbySpecies$lowerCI-1,beta1.pbySpecies$upperCI+1)), type = "n", xaxt="n",
     pch=19, xlab="", ylab="", las=2)
abline(h = 0,col = "grey")
arrows(x0 = 1:nrow(beta1.pbySpecies), x1 = 1:nrow(beta1.pbySpecies), y0=beta1.pbySpecies$lowerCI, y1=beta1.pbySpecies$upperCI, length=0, angle=90, code=3)
points(1:nrow(beta1.pbySpecies), beta1.pbySpecies$mean, pch = 20, cex = 1)
axis(side=1,at=1:nrow(beta1.pbySpecies),labels=row.names(beta1.pbySpecies), las=2, cex.axis=0.8)

# same plot rotated (write as a function to save jpeg in next step)
beta1.pPlot <- function(beta1.pbySpecies){
  par(mar=c(3,11,1.5,1.5))
  plot(0,0, ylim=range(c(1:nrow(beta1.pbySpecies))), xlim=range(c(beta1.pbySpecies$lowerCI-1,beta1.pbySpecies$upperCI+1)), type = "n", yaxt="n",
       pch=19, xlab="", ylab="", las=1)
  abline(v = 0,col = "grey")
  arrows(y0 = 1:nrow(beta1.pbySpecies), y1 = 1:nrow(beta1.pbySpecies), x0=beta1.pbySpecies$lowerCI, x1=beta1.pbySpecies$upperCI, length=0, angle=90, code=3)
  points(beta1.pbySpecies$mean, 1:nrow(beta1.pbySpecies), pch = 20, cex = 1)
  axis(side=2,at=1:nrow(beta1.pbySpecies),labels=row.names(beta1.pbySpecies), las=2, cex.axis=1)
}
dev.off()
beta1.pPlot(beta1.pbySpecies)

jpeg("/media/elildojr/Dados/r/monte_roraima/b1.p.jpg", width = 450, height = 600) # Open jpeg file
beta1.pPlot(beta1.pbySpecies)
dev.off()


#--------------------
# Examine how mean species-specific occupancy changes by recovery at logged sites

x <- seq(min(elevation), max(elevation), by=0.1) # use max(recovery) instead of 0.06 if all sites are included
y <- (x*sd(elevation)) + mean(elevation)

#alpha.site.PLOT <- fit$BUGSoutput$sims.list$alpha.site
alpha.psi.PLOT <- fit$BUGSoutput$sims.list$alpha.psi
a1.PLOT <- fit$BUGSoutput$sims.list$a1

dev.off()
sp.responses.elevation <- function() {
  
  par(mar = c(5,5,1,1))
  plot(y,y, type="l", ylim=c(0,1), xlim=c(min(y),max(y)),
       col="white", ylab="Probabilidade de ocupaçãoo", xlab="Elevação (padronizado)", cex=1.5, cex.lab=1.5)
  
  for (i in 1:n) {
    occPlot <- mean(alpha.psi.PLOT[,i]) + mean(a1.PLOT[,i])*x 
    lines(y,plogis(occPlot), type="l", ylim=c(0,1), xlim=c(0,40), main=i)
  }
  # to highlight an individual species:
  #species.1 <- mean(alpha.block.PLOT[,1,]) + mean(eps.PLOT[,1,]) + mean(alpha.psi.PLOT[,1]) +
  #                      mean(beta1.psi.PLOT[,1])*mean(treatment) + mean(beta2.psi.PLOT[,1])*x
  #lines(y,plogis(species.1), type="l", ylim=c(0,1), xlim=c(min(y),max(y)), lwd=4) # lwd for line thickness
}
sp.responses.elevation()

jpeg("/media/elildojr/Dados/r/monte_roraima/sp_responses_recovery.jpg", width = 750, height = 600) # Open jpeg file
sp.responses.recovery()
dev.off()


# Check individual species responses
species # check names 
species.check <- function(i) {
  plot(y,y, type="l", ylim=c(0,1), xlim=c(min(y),max(y)), main=paste(species[i], "response to recovery time"), cex.main=1, 
       col="white", ylab="Occupancy probability", xlab="Recovery time (standardized)")
  occPlot <- mean(alpha.block.PLOT[,i,]) + mean(eps.PLOT[,i,]) + mean(alpha.psi.PLOT[,i]) +
    mean(beta1.psi.PLOT[,i])*x
  lines(y,plogis(occPlot), type="l", ylim=c(0,1), xlim=c(0,40), main=i)
}
species.check(16) # mazama
species.check(2) # cuniculus
#etc

#----------------------
## plot relationships for a given species showing the uncertainty of the estimate

# predict effect of recovery with uncertainty
# write as a function

uncertainty.plot <- function(i){
  mcmc.sample <- fit$BUGSoutput$n.sims
  original.elevation <- seq(min(elevation), max(elevation), by=0.1)
  elevation.pred <- (original.elevation-mean(elevation))/sd(elevation)
  psi.pred.elevation <- plogis( mean(alpha.psi.PLOT[,i]) + mean(a1.PLOT[,i])*elevation.pred )
  
  array.psi.pred.elevation <- array(NA, dim=c(length(elevation.pred), mcmc.sample))
  
  for(j in 1:mcmc.sample) {
    array.psi.pred.elevation[,j] <- plogis(mean(fit$BUGSoutput$sims.list$alpha.psi[j,i]) + fit$BUGSoutput$sims.list$a1[j,i]*elevation.pred )
  }
  
  # plot for a subsample of mcmc draws
  sub.set <- sort(sample(1:mcmc.sample, size=200))
  plot(original.elevation, psi.pred.elevation, ylab="Occupancy probability", xlab="elevation time (standardized)", ylim=c(0,1),
       type="l", lwd=3, frame.plot=FALSE, main=species[i], cex.main=0.8)
  for(k in 1:length(sub.set)){
    lines(original.elevation, array.psi.pred.elevation[,k], type="l", lwd=1, col="azure3")
  }
  lines(original.elevation, psi.pred.elevation, type="l", lwd=2, col="blue")
}

uncertainty.plot(11) # Mazama spp
  uncertainty.plot(2) # Cuniculus
# etc.


# save multiplot as jpeg
jpeg("/media/elildojr/Dados/r/monte_roraima/multiplotOccVsElevPlusUncertainty.jpg", width = 800, height = 1000) # Open jpeg file
par(mfrow=c(3,5))
for(i in 1:15) {
  uncertainty.plot(i)
}
dev.off()

