## Camera trap Analytics Course - ICCB, Cartagena, JUly 22, 2017
## Jorge Ahumada @ Conservation International

# pegando diretório do script executado (Necessário usar ctrl/command-shift-s)
setwd(dirname(parent.frame(2)$ofile))

## ----Install these packages if you have not already---------------------------
install.packages(c("TeachingDemos","lubridate","unmarked","ggplot2","dplyr","chron","vegan", "activity", "ggmap"), repos="https://cloud.r-project.org")

## ----Load libraries------------------------
library(TeachingDemos)
library(lubridate)
library(unmarked)
library(ggplot2)
library(dplyr)
library(chron)
library(vegan)
library(activity)
library(ggmap)



## ----Source this file--------
source("camera trap analysis code-WILDID-09-20-17.R")

## ----Load data from Yanachaga,Peru-------
dataRBG <- f.readin.fix.data("rbg_tudo_set2017.csv")

# How many rows and columns?
dim(dataRBG)

#Look at the first 6 rows of data
head(dataRBG)
View(dataRBG)

## ----Subset the data-----------------
# Which photo types do we have?
table(dataRBG$Photo.Type)

# Only work with the animal photos
dataRBG <- filter(dataRBG, Photo.Type == "Animal")

## ----Number of records per camera and sampling regime--------------------
# Number of images
dim(dataRBG)

#Number of images per camera trap
imgsPerCT <- dataRBG %>% group_by(Camera.Trap.Name, Latitude, Longitude) %>% summarize(n = n()) %>% arrange(desc(n))
imgsPerCT

## ----make some plots-----------
#plot
p <- ggplot(imgsPerCT, aes(x = reorder(Camera.Trap.Name, -n), y = n)) + geom_col() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + xlab("Camera Trap Name") + ylab("Number of images")
p

## ----plot deployment in time---------------------------------------------
#Look at how camera traps were deployed
# First create a table
tabCTN_PD <- with(dataRBG, table(Camera.Trap.Name, Photo.Date))
head(tabCTN_PD)
#Get it ready for ggplot
tabCTN_PD <- melt(tabCTN_PD)

#Plot it
p <- ggplot(tabCTN_PD, aes(y = Camera.Trap.Name, x = Photo.Date)) + geom_raster(aes(fill=value), alpha = 0.8) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p


# Separate data into events for different types of analysis ---------------

# Extract the data from species Cuniculus paca
Daspry <- filter(dataRBG, bin == "Dasyprocta prymnolopha")

# Group by events that are two minutes apart
Daspry <- f.separate.events(Daspry, 2)

# Order the data sequentially
Daspry <- f.order.data(Daspry)

# Calculate for each event the time
# First get the median time for each event
acTimeDaspry <- Daspry %>% group_by(grp) %>% summarize(median = median(td.photo))

# Then extract just the hh:mm:ss information
acTimeDaspry <- f.extracthms(acTimeDaspry$median)

# Extract just the hour information - we will use later
acTimeDaspry_hour <- data.frame(seq = 1: length(acTimeDaspry), hour = hour(acTimeDaspry))

# Convert to radian time
acTimeDaspry <- f.timeformatfunc(acTimeDaspry)

actMod_Daspry <- fitact(acTimeDaspry, sample = "model", reps=100)




# Plot the result ---------------------------------------------------------

plot(actMod_Daspry, main = "Dasyprocta prymnolopha")
# Do a circular plot of activity data


# Plot results in a circular plot
ggplot(acTimeDaspry_hour, aes(x = hour)) + geom_histogram(breaks = seq(0, 24), colour = "grey") + coord_polar(start = 0) + theme_minimal() + ggtitle("Dasyprocta prymnolopha") + scale_x_continuous("", limits = c(0, 24), breaks = seq(0, 24), labels = seq(0, 24))

# Write a couple of functions to do this

calcActivity <- function(dataset, speciesName, threshold){
        # Extract the data from species
        spdata <- filter(dataset, bin == speciesName)
        
        # Group by events that are two minutes apart
        spdata <- f.separate.events(spdata, threshold)
        
        # Order the data sequentially
        spdata <- f.order.data(spdata)
        
        # Calculate for each event the time
        # First get the median time for each event
        acTime <- spdata %>% group_by(grp) %>% summarize(median = median(td.photo))
        # Then extract just the hh:mm:ss information
        acTime <- f.extracthms(acTime$median)
        # Extract the hour information for the circular plot
        acTime_hour <- data.frame(seq = 1:length(acTime), hour = hour(acTime))
        # Convert to radian time
        acTime <- f.timeformatfunc(acTime)
        # fit an activity kernel
        actModel <- fitact(acTime, sample = "model", reps=100)
        # Plot it in a couple ways
        
        plot(actModel, main = speciesName)
        
        list(actModel, acTime_hour)
}

plot.circular.activity <- function(actdata, speciesName) {
        ggplot(actdata, aes(x = hour)) + geom_histogram(breaks = seq(0, 24), colour = "grey") + coord_polar(start = 0)  + theme_minimal() + ggtitle(speciesName) + scale_x_continuous("", limits = c(0, 24), breaks = seq(0, 24), labels = seq(0, 24))
}

# Run it for other species
actMod_Mazame <- calcActivity(dataRBG, "Mazama americana", 2)

# Do a ciruclar plot of activity
plot.circular.activity(actMod_Mazame[[2]], "Mazama americana")

actMod_Mittub <- calcActivity(dataRBG, "Mitu tuberosum", 2)

plot.circular.activity(actMod_Mittub[[2]], "Mitu tuberosum")

# Compare the activity of two species
compareAct(list(actMod_Mazame[[1]],actMod_Daspry))


# Generate spatial distributions ----------------------------------

# Start with provide the lon/lat range of the data
lon <- range(dataRBG$Longitude)
lat <- range(dataRBG$Latitude)

# Extract the unique lat/lons and put them on a data frame
locationsRBG <- unique(cbind(as.character(dataRBG$Camera.Trap.Name), dataRBG$Latitude,dataRBG$Longitude))

locationsRBG <- data.frame(Camera.Trap.Name = locationsRBG[,1], Latitude = as.numeric(locationsRBG[,2]), Longitude = as.numeric(locationsRBG[,3]))

locationsRBG <- dplyr::arrange(locationsRBG, Camera.Trap.Name)

#locationsRBG <- locationsRBG[1:19,]

# If you have internet: Download the map from google
map <- get_map(location = c(c(lon[1],lat[1]),c(lon[2],lat[2])), zoom = 10, source = "google", maptype = "terrain")

# If you don't have internet
map <- readRDS("Gurupi_map.rds")

# Plot the locations of the camera traps
ggmap(map, extent = "normal", maprange = T) + geom_point(data=locationsRBG, aes(x = Longitude, y = Latitude), colour="red", size = 0.1)

# Plot the number of images per camera trap point
ggmap(map, extent = "normal", maprange = T) + geom_point(data = imgsPerCT, aes(x = Longitude, y = Latitude, color = n), size = 0.5)

# Plot as a surface
ggmap(map, extent = "device", legend = "topleft")  + stat_density2d(aes(x = Longitude, y = Latitude, fill = ..level..), data = dataRBG, geom = "polygon", size = 2, bins = 10, alpha = 0.5)

# Do it for a Dasyprocta
# First extract the number of photographic events by camera trap
sevPerCT_Daspry <- Daspry %>% group_by(Camera.Trap.Name, Latitude, Longitude) %>% summarize(n = length(unique(grp)))

# Then put them in a map
ggmap(map) + geom_point(data = sevPerCT_Daspry, aes(x = Longitude, y = Latitude, color = n), size = 0.5) + ggtitle("Dasyprocta prymnolopha")

# Can also be expressed as a relative abundance index - but CAREFUL - uncorrected for detection

ggmap(map) + geom_point(data = sevPerCT_Daspry, aes(x = Longitude, y = Latitude, color = n/sum(n)), size = 0.5) + ggtitle("Dasyprocta prymnolopha")

##---- Occupancy analysis ------------------------------------------------------
# First load some covariate data for Yanachaga

# Covariates --------------------------------------------------------------


covsRBG <- read.csv("team_gurupi_covars.csv", h = T)
head(covsRBG)

# We just need the cameras operating in 2016. Which ones are they?
suRBG <- unique(dataRBG$Camera.Trap.Name)

# filter our covariate file 
covsRBG <- filter(covsRBG, Camera.Trap.Name %in% suRBG)

# Normalize the two covariates - altitude and slope
covsRBG <- mutate(covsRBG, norm.altitude = (altitude.srtm - mean(altitude.srtm))/sd(altitude.srtm), norm.slope = (slope-mean(slope))/sd(slope))
head(covsRBG)

#Sort by camera trap unit name for analysis

covsRBG <- arrange(covsRBG, Camera.Trap.Name)
head(covsRBG)

# Create a matrix of presence/absence of each species
# rows are camera trap points (sampling units) and columns are dates
# Use function f.matrix.creator2 to do this

paMatsRBG <- f.matrix.creator2(dataRBG)

# This creates a list were each element of the list is a presence/absence matrix for a species
summary(paMatsRBG)
names(paMatsRBG)
# Look at the matrix for the third species D. prymnolopha

View(paMatsRBG[[3]])

# Convert these matrices into a special format to analyze them in unmarked

umRBG_Daspry <- unmarkedFrameOccu(y = paMatsRBG[[3]], siteCovs = covsRBG)
summary(umRBG_Daspry)
# fit a single season occupancy model with no covariates
occMod0_Daspry <- occu(~1 ~1, umRBG_Daspry)

# Look at the model results
occMod0_Daspry

# Transform the estimates from log to linear
backTransform(occMod0_Daspry, "state")
backTransform(occMod0_Daspry, "det")

occMod1_Daspry <- occu(~norm.altitude ~1, umRBG_Daspry)
occMod1_Daspry
backTransform(occMod1_Daspry, "state")

# We cannot use backTransform for p because it is dependent on elevation
# First create a dataframe with the average value of elevation (0)
newdata <- data.frame(norm.altitude = 0)
# Use predict to get the values
predict(occMod1_Daspry, type="det", newdata=newdata)

# How is detection varying with elevation
newdata <- data.frame(norm.altitude = seq(-1.5, 2.7, 0.05))
mod_pred <- predict(occMod1_Daspry, type="det", newdata=newdata)
mod_pred <- data.frame(mod_pred, elevation = newdata)

# Transform elevation back to original scale
mod_pred <- mutate(mod_pred, elev = sd(covsRBG$altitude.srtm)*norm.altitude + mean(covsRBG$altitude.srtm))
mod_pred
ggplot(mod_pred, aes(x = elev, y = Predicted)) + geom_line() + geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.4, fill = "blue") + xlab("Elevation") + ylab("Detection probability")

# Do another model where detection is a function of elevation and occupancy is a function of slope
occMod2_Daspry <- occu(~norm.altitude ~norm.slope, umRBG_Daspry)
occMod2_Daspry

# Compare all these models
modelList <- fitList(occMod0_Daspry, occMod1_Daspry, occMod2_Daspry)
modSel(modelList)
# Model 2 is better

# Lets do Dasyprocta fuliginosa
# umRBG_Dasful <- unmarkedFrameOccu(y = paMatsRBG[[15]], siteCovs = covsRBG)
# 
# 
# 
# occMod0_Dasful <- occu(~1 ~1, umRBG_Dasful)
# backTransform(occMod0_Dasful, "state")
# backTransform(occMod0_Dasful, "det")
# 
# # Lets fit some other models
# occMod1_Dasful <- occu(~1 ~norm.altitude, umRBG_Dasful)
# occMod1_Dasful
# confint(occMod1_Dasful,type = "state")
# 
# occMod2_Dasful <- occu(~ 1 ~ norm.slope, umRBG_Dasful)
# occMod2_Dasful
# confint(occMod2_Dasful,type = "state")
# 
# occMod3_Dasful <- occu(~ norm.altitude ~ norm.slope, umRBG_Dasful)
# occMod3_Dasful
# confint(occMod3_Dasful,type = "state")
# confint(occMod3_Dasful,type = "det")
# 
# occMod4_Dasful <- occu(~ norm.altitude + norm.slope ~ 1, umRBG_Dasful)
# occMod4_Dasful
# 
# # Compare all these models
# modelList <- fitList(occMod0_Dasful, occMod1_Dasful, occMod2_Dasful, occMod3_Dasful, occMod4_Dasful)
# modSel(modelList)

# Do it for Tapirus terrestris

umRBG_Tapter <- unmarkedFrameOccu(y = paMatsRBG[[14]], siteCovs = covsRBG)

# Fit the null model
occMod0_Tapter <- occu(~1 ~1, umRBG_Tapter)
backTransform(occMod0_Tapter, "state")
backTransform(occMod0_Tapter, "det")

# Fit a model were occupancy is changing with elevation
occMod1_Tapter <- occu(~1 ~norm.altitude, umRBG_Tapter)
occMod1_Tapter

# Fit a model were occupancy is changing with elevation and detection is changing with distance to edge
occMod2_Tapter <- occu(~norm.slope ~norm.altitude, umRBG_Tapter)
occMod2_Tapter

# Lets fit one more were we add distance to edge as a covariate with occupancy
occMod3_Tapter <- occu(~norm.slope ~norm.altitude + norm.slope, umRBG_Tapter)
occMod3_Tapter

# Compare these models
modelList <- fitList('psi()p()' = occMod0_Tapter, 'psi(slope)p()' = occMod1_Tapter, 'psi(slope)p(elev)' = occMod2_Tapter, 'psi(elev+slope)p(slope)' = occMod3_Tapter)

modSel(modelList)
# psi(elev)p(edge) & psi(elev+edge)p(edge) are better than the rest
# but the first one has less parameters

# Let's examine psi(slope)p(elev) in more detail
# First with psi 
newdata <- data.frame(norm.altitude = seq(-1.5, 2.7, 0.05))
mod_pred <- predict(occMod2_Tapter, type="state", newdata=newdata)
mod_pred <- data.frame(mod_pred, norm.altitude = newdata)
# Transform elevation back to original scale
mod_pred <- mutate(mod_pred, altitude.srtm = sd(covsRBG$altitude.srtm)*norm.altitude + mean(covsRBG$altitude.srtm))

# Plot
ggplot(mod_pred, aes(x = altitude.srtm, y = Predicted)) + geom_line() + geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.4, fill = "blue") + scale_y_continuous(limits = c(0,1)) + xlab("Elevation") + ylab("Occupancy") + ggtitle("Tapirus terrestris - RBG")
# Occupancy decreases with elevation - as expected by a lowland species

# look at detection probability
low_edge <- min(covsRBG$norm.slope)

hi_edge <- max(covsRBG$norm.slope)

newdata <- data.frame(norm.slope = seq(low_edge, hi_edge, 0.05))

mod_pred <- predict(occMod2_Tapter, type="det", newdata=newdata)

mod_pred <- data.frame(mod_pred, norm.slope = newdata)

# Transform elevation back to original scale
mod_pred <- mutate(mod_pred, slope = sd(covsRBG$slope)*norm.slope + mean(covsRBG$slope))

# Plot
ggplot(mod_pred, aes(x = slope, y = Predicted)) + geom_line() + geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.4, fill = "blue") + scale_y_continuous(limits = c(0,1)) + xlab("slope (%)") + ylab("Detection probability") + ggtitle("Tapirus terrestris - Yanachaga")
# Relatively minor effect compared to occupancy

# Put all this information in a spatial context
# First fill in the model values at each location for psi

newdata <- data.frame(norm.altitude = covsRBG$norm.altitude)

mod_pred <- predict(occMod2_Tapter, type="state", newdata=newdata)

mod_pred <- data.frame(mod_pred, norm.altitude = newdata)

# Transform elevation back to original scale
mod_pred <- mutate(mod_pred, 
                   elev = sd(covsRBG$altitude.srtm)*norm.altitude + mean(covsRBG$altitude.srtm),
                   Camera.Trap.Name = covsRBG$Camera.Trap.Name,
                   Latitude = locationsRBG$Latitude,
                   Longitude = locationsRBG$Longitude)

names(mod_pred)[1] <- "Occupancy"
# Put on a map
ggmap(map, extent = "normal", maprange = T) + geom_point(data = mod_pred, aes(x = Longitude, y = Latitude, color = Occupancy), size=0.5)


##---- Species richness analysis---------------------------
unique(dataRBG$bin)

#Remove the blank species
# Let's do this in a copy of dataRBG
dataRBG.copy <- dataRBG

dataRBG.copy$bin <- droplevels(dataRBG.copy$bin, exclude = " ")
unique(dataRBG.copy$bin)


# How often they show up in the camera traps

imgsPerSp <- dataRBG.copy %>% group_by(bin) %>% summarize(n = n()) %>% arrange(desc(n))
imgsPerSp

#plot
ggplot(imgsPerSp, aes(x = reorder(bin, -n), y = n)) + geom_col() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + xlab("Species") + ylab("Number of images")


# Build a simple species accumulation curve

# Put the data in matrix form
spMatrix <- with(dataRBG.copy, table(Camera.Trap.Name,bin))
spMatrix [1:10, ]

# Using function speaccum from package vegan
sp1 <- specaccum(spMatrix)
sp1
sp2 <- specaccum(spMatrix, "random")
sp2


# Plot accumulated species richness + confidence limits
plot(sp1, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="green")
boxplot(sp2, col="yellow", add=TRUE, pch="+")

# Exercise: do this for the whole Yanachaga data set and compare the results

##---- See Dorazio's (2006) approach to calculate species richness
# that takes into account imperfect detection




