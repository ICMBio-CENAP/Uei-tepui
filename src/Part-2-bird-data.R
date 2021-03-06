# monte roraima

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
source(here("bin", "funcao-para-grafico-composto.R"))

# carregar dados
#wiiTiwi <- read.csv(here("data", "images.csv"))
#deployments <- read.csv(here("data", "deployments.csv"))
#elevation <- read.csv(here("data", "elevation_uei_tepui.csv"))
wiiTiwi <- f.readin.fix.data(here("data", "Wild_ID_UeiTepui.csv"))


# Only work with bird photos
wiiTiwi <- filter(wiiTiwi, Class == "AVES")
photos <- as.data.frame(table(wiiTiwi$bin))

# Group by events that are 30 minutes apart
wiiTiwi <- f.separate.events(wiiTiwi, 30)
wiiTiwi <- distinct(wiiTiwi, Camera.Trap.Name, bin, grp, .keep_all = TRUE)
records <- as.data.frame(table(wiiTiwi$bin))

# create table for appendix of the report:
tableAppendix <- cbind(photos[,1:2], records[,2]); names(tableAppendix) <- c("bin", "N.fotos", "N.registros")
tableAppendix <- merge(tableAppendix, wiiTiwi[,c(12:13,32)], by="bin") # incluir Classe Ordem e Familia
tableAppendix <- distinct(tableAppendix, .keep_all = TRUE)
tableAppendix <- tableAppendix[,c(4:5,1:3)]
names(tableAppendix) <- c("Ordem", "Família", "Espécie", "N.fotos", "N.registros")
tableAppendix
write.csv(tableAppendix, here("results", "tableAppendix.csv"), row.names = FALSE)

# Number of images per camera trap
imgsPerCT <- wiiTiwi %>% group_by(Camera.Trap.Name) %>% summarize(n = n()) %>% arrange(desc(n))
imgsPerCT

# plot n. photos versus elevation
elevation <- read.csv(here("data", "elevation_uei_tepui.csv")) # load elevation
elevation$Camera.Trap.Name <- substr(elevation$deployment_id, 1, 7)
df1 <- merge(imgsPerCT, elevation, "Camera.Trap.Name", sort = TRUE)
with(df1, plot(elevation, n))

# Effort per camera trap
#df2 <- merge(wiiTiwi, deployments, "Camera.Trap.Name", sort = TRUE)
#df2$effort = as.numeric( as.Date(df2$end_date) - as.Date(df2$start_date) )
df2 <- wiiTiwi
df2$effort = as.numeric( as.Date(df2$End.Date) - as.Date(df2$Start.Date) )
effortPerCT <- distinct(df2, Camera.Trap.Name, effort, .keep_all = FALSE)
effortPerCT <- merge(imgsPerCT, effortPerCT, "Camera.Trap.Name", sort = TRUE)
effortPerCT$n.corrected <- effortPerCT$n / effortPerCT$effort

df3 <- merge(df1, effortPerCT, "Camera.Trap.Name", sort = TRUE)

#------------------
mod1 <- glm(n.corrected~elevation, data=df3)
summary(mod1)

# same plot with jitter, linear trend and confidence interval
p1 <- ggplot(df3, aes(x=elevation, y=n.corrected)) + 
  geom_point(position=position_jitter(h=0, w=0.1), shape = 19, alpha = 0.5, size = 4) +
  geom_smooth(method=lm , color="black", fill="grey", se=TRUE) +
  theme(axis.title.x = element_text(size=20, margin=margin(t=10, r=0, b=0, l=0)),
        axis.title.y = element_text(size=20, margin=margin(t=0, r=10, b=0, l=0))) +
  theme(panel.grid.minor = element_blank(), panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
  theme(axis.text.x = element_text(size=18), axis.text.y = element_text(size=18)) +
  theme(legend.title=element_text(size=18), legend.text=element_text(size=18)) +
  xlab("Elevação (m)") + ylab("Número de registros / dia")
p1

# save as jpeg
jpeg(here("results", "registrosVselevacaoAVES.jpg"), width = 800, height = 600) # Open jpeg file
p1
dev.off()


#-----------------
#plot it photos per camera
p <- ggplot(df1, aes(x = reorder(Camera.Trap.Name, -n), y = n)) + geom_col() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + xlab("Camera Trap") + ylab("Numero de registros") +
  theme(panel.grid.minor = element_blank(), panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
  theme(axis.text.x = element_text(size=18), axis.text.y = element_text(size=18))
p

# save as jpeg
jpeg(here("results", "imgsPerCTAVES.jpg"), width = 800, height = 600) # Open jpeg file
p
dev.off()



##---- Species richness analysis---------------------------
unique(wiiTiwi$bin)

# How often they show up in the camera traps
imgsPerSp <- wiiTiwi %>% group_by(bin) %>% summarize(n = n()) %>% arrange(desc(n))
imgsPerSp
imgsPerSp_No_Leptotilla <- imgsPerSp[-1,]

write.csv(imgsPerSp, here("results", "imgsPerSpAVES.csv"))

# how many sites each species appeared
spPerSite <- data.matrix(table(wiiTiwi$bin, wiiTiwi$Camera.Trap.Name))
spPerSite[spPerSite > 0] <- 1
spPerSite <- cbind(rownames(spPerSite), rowSums(spPerSite))

# plot imgsPerSp
p2 <- ggplot(imgsPerSp, aes(x = reorder(bin, -n), y = n)) + geom_col() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + xlab("") + ylab("Número de registros") +
  theme(panel.grid.minor = element_blank(), panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
  theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14)) +
  theme(axis.title.x = element_text(size=20, margin=margin(t=10, r=0, b=0, l=0)),
        axis.title.y = element_text(size=20, margin=margin(t=0, r=10, b=0, l=0)))
p2

# save as jpeg
jpeg(here("results", "imgsPerSPAVES.jpg"), width = 800, height = 600) # Open jpeg file
p2
dev.off()

# plot imgsPerSp no Leptotilla
p3 <- ggplot(imgsPerSp_No_Leptotilla, aes(x = reorder(bin, -n), y = n)) + geom_col() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + xlab("") + ylab("Número de registros") +
  theme(panel.grid.minor = element_blank(), panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
  theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14)) +
  theme(axis.title.x = element_text(size=20, margin=margin(t=10, r=0, b=0, l=0)),
        axis.title.y = element_text(size=20, margin=margin(t=0, r=10, b=0, l=0)))
p3

# save as jpeg
jpeg(here("results", "imgsPerSp_No_Leptotilla.jpg"), width = 800, height = 600) # Open jpeg file
p3
dev.off()


# which spp were more recorded
sort(imgsPerSp$n)
sum(imgsPerSp$n)

# Build a simple species accumulation curve

# Put the data in matrix form
spMatrix <- with(wiiTiwi, table(Photo.Date,bin))
spMatrix [1:10, ]


# Using function speaccum from package vegan
sp1 <- specaccum(spMatrix, "random")
sp1

# Plot accumulated species richness + confidence limits
plot(sp1, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="grey", ylab="Riqueza de espécies", xlab="Esforço (dias)")
#boxplot(sp1, col="yellow", add=TRUE, pch="+")

# save as jpeg
jpeg(here("results", "spAccumAVES.jpg"), width = 800, height = 600) # Open jpeg file
par(mar=c(6,6,1,1))
plot(sp1, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="grey", ylab="Riqueza de espécies", xlab="Esforço (dias)", cex.axis=1.5, cex.lab=2, las=1)
dev.off()


## gradient figure

# matrix by site
spMatrixBySite <- as.data.frame.matrix(with(wiiTiwi, table(Camera.Trap.Name,bin)))
spMatrixBySite <- spMatrixBySite[order(rownames(spMatrixBySite)),] # order by camera trap name
names(spMatrixBySite)
spMatrixBySite[1:19, 1:2]

# We just need cameras that worked. Which ones are they?
suWiiTiwi <- unique(wiiTiwi$Camera.Trap.Name)

# filter our covariate (elevation) file
elevation2 <- filter(elevation, Camera.Trap.Name %in% suWiiTiwi)
rownames(elevation2) <- elevation2$Camera.Trap.Name
elevation2 <- elevation2[order(rownames(elevation2)),] # order by camera trap name

# use "generico" function
generico(spMatrixBySite, elevation2$elevation, 24, "Elevação","Espécie","Elevação")

# save as jpeg
jpeg(here("results", "grad_elevacaoAVES.jpg")) # Open jpeg file
generico(spMatrixBySite, elevation2$elevation, 24, "Elevação","Espécie","Elevação")
dev.off()
