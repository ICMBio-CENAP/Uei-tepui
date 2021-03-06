# Put wiiTiwi data in the same format as Wild.ID exports

# carregar bibliotecas
library(TeachingDemos)
library(lubridate)
library(unmarked)
library(ggplot2)
library(dplyr)
library(reshape2)
library(chron)
library(vegan)
library(activity)
library(ggmap)
library(here)

# carregar esse arquivo
source(here("bin", "camera trap analysis code-WILDID-09-20-17.R"))
source(here("bin", "time-lag.R"))

# carregar dados
wiiTiwi <- read.csv(here("data", "images.csv"))
deployments <- read.csv(here("data", "deployments.csv"))
elevation <- read.csv(here("data", "elevation_uei_tepui.csv"))

# change names
names(wiiTiwi) <- c("Project.Name", "Camera.Trap.Name", "image_id", "location", "is_blank", "Person.Identifying.the.Photo",
                    "wi_taxon_id", "Class", "Order", "Family", "Genus", "Species", "commom_name", "uncertainty", "timestamp", "age",
                    "sex", "animal_recognizable", "individual_id", "Number.of.Animals", "individual_animal_notes",
                    "highlighted", "color", "licence")

# new columns
wiiTiwi$ID <- seq(1:nrow(wiiTiwi))
wiiTiwi$Camera.Trap.Name <-  substr(wiiTiwi$Camera.Trap.Name, start = 1, stop = 7)
wiiTiwi$Photo.Date <- substr(wiiTiwi$timestamp, start = 1, stop = 10)
wiiTiwi$Photo.time <- substr(wiiTiwi$timestamp, start = 12, stop = nchar(wiiTiwi$timestamp))
wiiTiwi$Camera.Manufacturer <- "Bushnell"
wiiTiwi$Camera.Model <- "Tropy Cam HD"
wiiTiwi$Camera.Serial.Number <- NA
wiiTiwi$Sequence.Info <- NA
wiiTiwi$Moon.Phase <- NA
wiiTiwi$Temperature <- NA
wiiTiwi$Flash <- NA
wiiTiwi$Organization.Name <- "ICMBio/CENAP"
wiiTiwi$Sampling.Event <- 2018
wiiTiwi$Latitude <- NA
wiiTiwi$Longitude <- NA
wiiTiwi$Raw.Name <- NA
wiiTiwi$Photo.Type <- NA
wiiTiwi$Person.picking.up.the.Camera <- "Dionisio Ingariko"
wiiTiwi$Person.setting.up.the.Camera <- "Elildo Carvalho Jr"


# add latitude and longitude and start and end dates
df1 <- distinct(deployments, placename, latitude, longitude, start_date, end_date)
names(df1) <- c("Camera.Trap.Name", "Longitude", "Latitude", "Camera.Start.Date", "Camera.End.Date")
df1$Camera.Trap.Name <- as.factor(df1$Camera.Trap.Name)

# add unify coordinates for all sites-years in the dataset
df2 <- wiiTiwi
newcols <- c("Latitude", "Longitude", "Camera.Start.Date", "Camera.End.Date")
df2[newcols] <- lapply(newcols, function(x) df1[[x]][match(df2$Camera.Trap.Name, df1$Camera.Trap.Name)])
#View(df2)
wiiTiwi <- df2

# remove some columns
names(wiiTiwi)
wiiTiwi$image_id <- NULL
wiiTiwi$location <- NULL
wiiTiwi$is_blank <- NULL
wiiTiwi$image_id <- NULL
wiiTiwi$wi_taxon_id <- NULL
wiiTiwi$uncertainty <- NULL
wiiTiwi$timestamp <- NULL
wiiTiwi$age <- NULL
wiiTiwi$sex <- NULL
wiiTiwi$animal_recognizable <- NULL
wiiTiwi$individual_id <- NULL
wiiTiwi$individual_animal_notes <- NULL
wiiTiwi$commom_name <- NULL
wiiTiwi$highlighted <- NULL
wiiTiwi$color <- NULL
wiiTiwi$licence <- NULL

# change higher taxonomic levels to uppercase
wiiTiwi$Class <- toupper(wiiTiwi$Class)
wiiTiwi$Class[wiiTiwi$Class == "NO CV RESULT"] <- NA
wiiTiwi$Class[wiiTiwi$Class == "UNKNOWN"] <- NA
wiiTiwi$Order <- toupper(wiiTiwi$Order)
wiiTiwi$Order[wiiTiwi$Order == "NO CV RESULT"] <- NA
wiiTiwi$Order[wiiTiwi$Order == "UNKNOWN"] <- NA
wiiTiwi$Family <- toupper(wiiTiwi$Family)
wiiTiwi$Family[wiiTiwi$Family == "NO CV RESULT"] <- NA
wiiTiwi$Family[wiiTiwi$Family == "UNKNOWN"] <- NA

# reorder columns
col_order <- c("ID",	"Project.Name",	"Camera.Trap.Name",	"Latitude",	"Longitude",	
               "Sampling.Event",	"Photo.Type",	"Photo.Date",	"Photo.time",	"Raw.Name",	
               "Class",	"Order",	"Family",	"Genus",	"Species",	"Number.of.Animals",	
               "Person.Identifying.the.Photo",	"Camera.Serial.Number",	"Camera.Start.Date",	"Camera.End.Date",	"Person.setting.up.the.Camera",	"Person.picking.up.the.Camera",	"Camera.Manufacturer",	"Camera.Model",	"Sequence.Info",	"Moon.Phase",	"Temperature",	"Flash",	"Organization.Name")
wiiTiwi <- wiiTiwi[,col_order]

# check dates
sort(unique((wiiTiwi$Camera.Start.Date)))
sort(unique((wiiTiwi$Camera.End.Date)))
sort(unique((wiiTiwi$Photo.Date)))
wiiTiwi <- subset(wiiTiwi, Photo.Date >= "2019-10-30") # remove wrong dates


# check timelag
wiiTiwi$Camera.Start.Date <- as.Date(wiiTiwi$Camera.Start.Date)
wiiTiwi$Camera.End.Date <- as.Date(wiiTiwi$Camera.End.Date)
wiiTiwi$Photo.Date <- as.Date(wiiTiwi$Photo.Date)
time.lag(wiiTiwi)

# fix start:
wiiTiwi[wiiTiwi$Camera.Trap.Name=="MTRR-18",]$Camera.Start.Date <- min(subset(wiiTiwi, Camera.Trap.Name=="MTRR-18")$Photo.Date)
wiiTiwi[wiiTiwi$Camera.Trap.Name=="MTRR-16",]$Camera.Start.Date <- min(subset(wiiTiwi, Camera.Trap.Name=="MTRR-16")$Photo.Date)
wiiTiwi[wiiTiwi$Camera.Trap.Name=="MTRR-15",]$Camera.Start.Date <- min(subset(wiiTiwi, Camera.Trap.Name=="MTRR-15")$Photo.Date)
# fix end:
wiiTiwi[wiiTiwi$Camera.Trap.Name=="MTRR-22",]$Camera.End.Date <- max(subset(wiiTiwi, Camera.Trap.Name=="MTRR-22")$Photo.Date)
wiiTiwi[wiiTiwi$Camera.Trap.Name=="MTRR-21",]$Camera.End.Date <- max(subset(wiiTiwi, Camera.Trap.Name=="MTRR-21")$Photo.Date)
wiiTiwi[wiiTiwi$Camera.Trap.Name=="MTRR-20",]$Camera.End.Date <- max(subset(wiiTiwi, Camera.Trap.Name=="MTRR-20")$Photo.Date)
wiiTiwi[wiiTiwi$Camera.Trap.Name=="MTRR-19",]$Camera.End.Date <- max(subset(wiiTiwi, Camera.Trap.Name=="MTRR-19")$Photo.Date)
wiiTiwi[wiiTiwi$Camera.Trap.Name=="MTRR-16",]$Camera.End.Date <- max(subset(wiiTiwi, Camera.Trap.Name=="MTRR-16")$Photo.Date)
wiiTiwi[wiiTiwi$Camera.Trap.Name=="MTRR-11",]$Camera.End.Date <- max(subset(wiiTiwi, Camera.Trap.Name=="MTRR-11")$Photo.Date)
wiiTiwi[wiiTiwi$Camera.Trap.Name=="MTRR-10",]$Camera.End.Date <- max(subset(wiiTiwi, Camera.Trap.Name=="MTRR-10")$Photo.Date)
wiiTiwi[wiiTiwi$Camera.Trap.Name=="MTRR-05",]$Camera.End.Date <- max(subset(wiiTiwi, Camera.Trap.Name=="MTRR-05")$Photo.Date)
time.lag(wiiTiwi) # check

# save as csv
write.csv(wiiTiwi, here("data", "Wild_ID_UeiTepui.csv"), row.names = FALSE)


