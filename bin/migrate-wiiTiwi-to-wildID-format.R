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

# carregar dados
wiiTiwi <- read.csv(here("data", "images.csv"))
deployments <- read.csv(here("data", "deployments.csv"))
elevation <- read.csv(here("data", "elevation_uei_tepui.csv"))

# nova coluna bin
wiiTiwi$bin <- paste(wiiTiwi$genus, wiiTiwi$species, sep=" ")
# nova coluna Photo.Date
wiiTiwi$Photo.Date <- as.Date(substr(wiiTiwi$timestamp, start = 1, stop = 10))
# nova coluna td.photo
wiiTiwi$td.photo <- as.POSIXct(wiiTiwi$timestamp)
#Remove Homo and Canis
wiiTiwi <- subset(wiiTiwi, bin != "Homo sapiens" & bin != "Canis " & bin != "Canis lupus familiaris" & bin != "Proechimys ")
# fix some species names
wiiTiwi$bin[wiiTiwi$bin == "Cebus "] <- "Cebus olivaceus"
wiiTiwi$bin[wiiTiwi$bin == "Nasua "] <- "Nasua nasua"
wiiTiwi$bin[wiiTiwi$bin == "Mazama "] <- "Mazama sp"
#wiiTiwi$bin[wiiTiwi$bin == "Nasuella "] <- "cf. Nasuella"
wiiTiwi$bin[wiiTiwi$bin == "Nasuella "] <- "Nasua nasua"
wiiTiwi$bin[wiiTiwi$bin == "Didelphis "] <- "Didelphis imperfecta"


wiiTiwi$ID <- seq(1:nrow(wiiTiwi))
wiiTiwi$Photo.Date <- substr(wiiTiwi$date.time, start = 1, stop = 10)
wiiTiwi$Photo.time <- substr(wiiTiwi$date.time, start = 12, stop = nchar(wiiTiwi$date.time))
wiiTiwi$Class <- NA
wiiTiwi$Order <- NA
wiiTiwi$Family <- NA
wiiTiwi$Genus <- gsub( " .*$", "", wiiTiwi$species)
wiiTiwi$Species <- sub("^\\S+\\s+", '', wiiTiwi$species)
wiiTiwi$Camera.Manufacturer <- "Bushnell"
wiiTiwi$Camera.Model <- "Tropy Cam HD"
wiiTiwi$Sequence.Info <- NA
wiiTiwi$Moon.Phase <- NA
wiiTiwi$Temperature <- NA
wiiTiwi$Flash <- NA
wiiTiwi$Organization.Name <- "ICMBio/CENAP"
wiiTiwi$Sampling.Event <- 2018

# rename columns
names(wiiTiwi) <- c("Project.Name", "Camera.Trap.Name", "Latitude", "Longitude", "Raw.Name", "Photo.Time", "Photo.Type", "bin", "Number.of.Animals", "animal.recognizable", "uncertainty", "Person.Identifying.the.Photo", "Camera.Serial.Number", "Camera.Start.Date", "Camera.End.Date", "Person.setting.up.the.Camera", "Person.picking.up.the.Camera", "notes", "to.be.modified", "ID", "Photo.Date", "Photo.time", "Class", "Order", "Family", "Genus", "Species", "Camera.Manufacturer"	,"Camera.Model", "Sequence.Info",	"Moon.Phase",	"Temperature",	"Flash",	"Organization.Name", "Sampling.Event")

# remove columns
wiiTiwi$Photo.Time <- NULL
wiiTiwi$bin <- NULL
wiiTiwi$animal.recognizable <- NULL
wiiTiwi$uncertainty <- NULL
wiiTiwi$notes <- NULL
wiiTiwi$to.be.modified <- NULL

# reorder columns
col_order <- c("ID",	"Project.Name",	"Camera.Trap.Name",	"Latitude",	"Longitude",	"Sampling.Event",	"Photo.Type",	"Photo.Date",	"Photo.time",	"Raw.Name",	"Class",	"Order",	"Family",	"Genus",	"Species",	"Number.of.Animals",	"Person.Identifying.the.Photo",	"Camera.Serial.Number",	"Camera.Start.Date",	"Camera.End.Date",	"Person.setting.up.the.Camera",	"Person.picking.up.the.Camera",	"Camera.Manufacturer",	"Camera.Model",	"Sequence.Info",	"Moon.Phase",	"Temperature",	"Flash",	"Organization.Name")
wiiTiwi <- wiiTiwi[,col_order]
#View(wiiTiwi)
