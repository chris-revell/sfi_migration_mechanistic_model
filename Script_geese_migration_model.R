library(maps)
library(mapdata)
library(raster)
library(sp)
library(rgdal)

## Load environmental data

NDVI.2006.feb.1 <- raster("/Users/mariussomveille/Desktop/NDVI_data/2006/February/06feb15a.n17-VIg_data.tif")
NDVI.2006.feb.1 <- crop(NDVI.2006.feb.1, extent(NDVI.2006.feb.1, 0, 700, 2000, 4000))
NDVI.2006.feb.1 <- calc(NDVI.2006.feb.1, fun=function(x){ifelse(x < -5000, NA, x)})
NDVI.2006.feb.2 <- raster("/Users/mariussomveille/Desktop/NDVI_data/2006/February/06feb15b.n17-VIg_data.tif")
NDVI.2006.feb.2 <- crop(NDVI.2006.feb.2, extent(NDVI.2006.feb.2, 0, 700, 2000, 4000))
NDVI.2006.feb.2 <- calc(NDVI.2006.feb.2, fun=function(x){ifelse(x < -5000, NA, x)})
NDVI.2006.mar.1 <- raster("/Users/mariussomveille/Desktop/NDVI_data/2006/March/06mar15a.n17-VIg_data.tif")
NDVI.2006.mar.1 <- crop(NDVI.2006.mar.1, extent(NDVI.2006.mar.1, 0, 700, 2000, 4000))
NDVI.2006.mar.1 <- calc(NDVI.2006.mar.1, fun=function(x){ifelse(x < -5000, NA, x)})
NDVI.2006.mar.2 <- raster("/Users/mariussomveille/Desktop/NDVI_data/2006/March/06mar15b.n17-VIg_data.tif")
NDVI.2006.mar.2 <- crop(NDVI.2006.mar.2, extent(NDVI.2006.mar.2, 0, 700, 2000, 4000))
NDVI.2006.mar.2 <- calc(NDVI.2006.mar.2, fun=function(x){ifelse(x < -5000, NA, x)})
NDVI.2006.apr.1 <- raster("/Users/mariussomveille/Desktop/NDVI_data/2006/April/06apr15a.n17-VIg_data.tif")
NDVI.2006.apr.1 <- crop(NDVI.2006.apr.1, extent(NDVI.2006.apr.1, 0, 700, 2000, 4000))
NDVI.2006.apr.1 <- calc(NDVI.2006.apr.1, fun=function(x){ifelse(x < -5000, NA, x)})
NDVI.2006.apr.2 <- raster("/Users/mariussomveille/Desktop/NDVI_data/2006/April/06apr15b.n17-VIg_data.tif")
NDVI.2006.apr.2 <- crop(NDVI.2006.apr.2, extent(NDVI.2006.apr.2, 0, 700, 2000, 4000))
NDVI.2006.apr.2 <- calc(NDVI.2006.apr.2, fun=function(x){ifelse(x < -5000, NA, x)})
NDVI.2006.may.1 <- raster("/Users/mariussomveille/Desktop/NDVI_data/2006/May/06may15a.n17-VIg_data.tif")
NDVI.2006.may.1 <- crop(NDVI.2006.may.1, extent(NDVI.2006.may.1, 0, 700, 2000, 4000))
NDVI.2006.may.1 <- calc(NDVI.2006.may.1, fun=function(x){ifelse(x < -5000, NA, x)})
NDVI.2006.may.2 <- raster("/Users/mariussomveille/Desktop/NDVI_data/2006/May/06may15b.n17-VIg_data.tif")
NDVI.2006.may.2 <- crop(NDVI.2006.may.2, extent(NDVI.2006.may.2, 0, 700, 2000, 4000))
NDVI.2006.may.2 <- calc(NDVI.2006.may.2, fun=function(x){ifelse(x < -5000, NA, x)})
NDVI.2006.jun.1 <- raster("/Users/mariussomveille/Desktop/NDVI_data/2006/June/06jun15a.n17-VIg_data.tif")
NDVI.2006.jun.1 <- crop(NDVI.2006.jun.1, extent(NDVI.2006.jun.1, 0, 700, 2000, 4000))
NDVI.2006.jun.1 <- calc(NDVI.2006.jun.1, fun=function(x){ifelse(x < -5000, NA, x)})
NDVI.2006.jun.2 <- raster("/Users/mariussomveille/Desktop/NDVI_data/2006/June/06jun15b.n17-VIg_data.tif")
NDVI.2006.jun.2 <- crop(NDVI.2006.jun.2, extent(NDVI.2006.jun.2, 0, 700, 2000, 4000))
NDVI.2006.jun.2 <- calc(NDVI.2006.jun.2, fun=function(x){ifelse(x < -5000, NA, x)})
NDVI.2006.jul.1 <- raster("/Users/mariussomveille/Desktop/NDVI_data/2006/July/06jul15a.n17-VIg_data.tif")
NDVI.2006.jul.1 <- crop(NDVI.2006.jul.1, extent(NDVI.2006.jul.1, 0, 700, 2000, 4000))
NDVI.2006.jul.1 <- calc(NDVI.2006.jul.1, fun=function(x){ifelse(x < -5000, NA, x)})
NDVI.2006.jul.2 <- raster("/Users/mariussomveille/Desktop/NDVI_data/2006/July/06jul15b.n17-VIg_data.tif")
NDVI.2006.jul.2 <- crop(NDVI.2006.jul.2, extent(NDVI.2006.jul.2, 0, 700, 2000, 4000))
NDVI.2006.jul.2 <- calc(NDVI.2006.jul.2, fun=function(x){ifelse(x < -5000, NA, x)})
NDVI.2006.aug.1 <- raster("/Users/mariussomveille/Desktop/NDVI_data/2006/August/06aug15a.n17-VIg_data.tif")
NDVI.2006.aug.1 <- crop(NDVI.2006.aug.1, extent(NDVI.2006.aug.1, 0, 700, 2000, 4000))
NDVI.2006.aug.1 <- calc(NDVI.2006.aug.1, fun=function(x){ifelse(x < -5000, NA, x)})
NDVI.2006.aug.2 <- raster("/Users/mariussomveille/Desktop/NDVI_data/2006/August/06aug15b.n17-VIg_data.tif")
NDVI.2006.aug.2 <- crop(NDVI.2006.aug.2, extent(NDVI.2006.aug.2, 0, 700, 2000, 4000))
NDVI.2006.aug.2 <- calc(NDVI.2006.aug.2, fun=function(x){ifelse(x < -5000, NA, x)})
NDVI.2006.sep.1 <- raster("/Users/mariussomveille/Desktop/NDVI_data/2006/September/06sep15a.n17-VIg_data.tif")
NDVI.2006.sep.1 <- crop(NDVI.2006.sep.1, extent(NDVI.2006.sep.1, 0, 700, 2000, 4000))
NDVI.2006.sep.1 <- calc(NDVI.2006.sep.1, fun=function(x){ifelse(x < -5000, NA, x)})
NDVI.2006.sep.2 <- raster("/Users/mariussomveille/Desktop/NDVI_data/2006/September/06sep15b.n17-VIg_data.tif")
NDVI.2006.sep.2 <- crop(NDVI.2006.sep.2, extent(NDVI.2006.sep.2, 0, 700, 2000, 4000))
NDVI.2006.sep.2 <- calc(NDVI.2006.sep.2, fun=function(x){ifelse(x < -5000, NA, x)})
NDVI.2006.oct.1 <- raster("/Users/mariussomveille/Desktop/NDVI_data/2006/October/06oct15a.n17-VIg_data.tif")
NDVI.2006.oct.1 <- crop(NDVI.2006.oct.1, extent(NDVI.2006.oct.1, 0, 700, 2000, 4000))
NDVI.2006.oct.1 <- calc(NDVI.2006.oct.1, fun=function(x){ifelse(x < -5000, NA, x)})
NDVI.2006.oct.2 <- raster("/Users/mariussomveille/Desktop/NDVI_data/2006/October/06oct15b.n17-VIg_data.tif")
NDVI.2006.oct.2 <- crop(NDVI.2006.oct.2, extent(NDVI.2006.oct.2, 0, 700, 2000, 4000))
NDVI.2006.oct.2 <- calc(NDVI.2006.oct.2, fun=function(x){ifelse(x < -5000, NA, x)})
NDVI.2006.nov.1 <- raster("/Users/mariussomveille/Desktop/NDVI_data/2006/November/06nov15a.n17-VIg_data.tif")
NDVI.2006.nov.1 <- crop(NDVI.2006.nov.1, extent(NDVI.2006.nov.1, 0, 700, 2000, 4000))
NDVI.2006.nov.1 <- calc(NDVI.2006.nov.1, fun=function(x){ifelse(x < -5000, NA, x)})
NDVI.2006.nov.2 <- raster("/Users/mariussomveille/Desktop/NDVI_data/2006/November/06nov15b.n17-VIg_data.tif")
NDVI.2006.nov.2 <- crop(NDVI.2006.nov.2, extent(NDVI.2006.nov.2, 0, 700, 2000, 4000))
NDVI.2006.nov.2 <- calc(NDVI.2006.nov.2, fun=function(x){ifelse(x < -5000, NA, x)})
NDVI.2006.dec.1 <- raster("/Users/mariussomveille/Desktop/NDVI_data/2006/December/06dec15a.n17-VIg_data.tif")
NDVI.2006.dec.1 <- crop(NDVI.2006.dec.1, extent(NDVI.2006.dec.1, 0, 700, 2000, 4000))
NDVI.2006.dec.1 <- calc(NDVI.2006.dec.1, fun=function(x){ifelse(x < -5000, NA, x)})
NDVI.2006.dec.2 <- raster("/Users/mariussomveille/Desktop/NDVI_data/2006/December/06dec15b.n17-VIg_data.tif")
NDVI.2006.dec.2 <- crop(NDVI.2006.dec.2, extent(NDVI.2006.dec.2, 0, 700, 2000, 4000))
NDVI.2006.dec.2 <- calc(NDVI.2006.dec.2, fun=function(x){ifelse(x < -5000, NA, x)})




# load bird tracking data 

data <- read.csv("/Users/mariussomveille/Desktop/CSSS/Tracking_animals/tracking goose data")
temp <- lapply(data$timestamp, function(x) as.numeric(unlist(strsplit(strsplit(as.character(x), " ")[[1]][1], "-"))))
data[,15:17] <- matrix(unlist(temp), ncol = 3, byrow = TRUE)
colnames(data)[15:17] <- c("year", "month", "day")

bird <- unique(data$individual.local.identifier)



## Plot ndvi and migration

# NDVI value
par(mfrow=c(2,3), mar=c(1.5,1.5,0.5,0.5), mgp=c(0,0.5,0))
plot(NDVI.2006.mar.2)
for(i in 1:length(bird)){
points(data$location.long[which(data$individual.local.identifier == bird[i] & data$migration.stage.standard == "spring-migration" & data$month == 3 & data$day > 15)], data$location.lat[which(data$individual.local.identifier == bird[i] & data$migration.stage.standard == "spring-migration"& data$month == 3 & data$day > 15)], pch=20, type="l")
}
plot(NDVI.2006.apr.1)
for(i in 1:length(bird)){
points(data$location.long[which(data$individual.local.identifier == bird[i] & data$migration.stage.standard == "spring-migration" & data$month == 4 & data$day <= 15)], data$location.lat[which(data$individual.local.identifier == bird[i] & data$migration.stage.standard == "spring-migration"& data$month == 4 & data$day <= 15)], pch=20, type="l")
}
plot(NDVI.2006.apr.2)
for(i in 1:length(bird)){
points(data$location.long[which(data$individual.local.identifier == bird[i] & data$migration.stage.standard == "spring-migration" & data$month == 4 & data$day > 15)], data$location.lat[which(data$individual.local.identifier == bird[i] & data$migration.stage.standard == "spring-migration"& data$month == 4 & data$day > 15)], pch=20, type="l")
}
plot(NDVI.2006.may.1)
for(i in 1:length(bird)){
points(data$location.long[which(data$individual.local.identifier == bird[i] & data$migration.stage.standard == "spring-migration" & data$month == 5 & data$day <= 15)], data$location.lat[which(data$individual.local.identifier == bird[i] & data$migration.stage.standard == "spring-migration"& data$month == 5 & data$day <= 15)], pch=20, type="l")
}
plot(NDVI.2006.may.2)
for(i in 1:length(bird)){
points(data$location.long[which(data$individual.local.identifier == bird[i] & data$migration.stage.standard == "spring-migration" & data$month == 5 & data$day > 15)], data$location.lat[which(data$individual.local.identifier == bird[i] & data$migration.stage.standard == "spring-migration"& data$month == 5 & data$day > 15)], pch=20, type="l")
}
plot(NDVI.2006.jun.1)
for(i in 1:length(bird)){
points(data$location.long[which(data$individual.local.identifier == bird[i] & data$migration.stage.standard == "spring-migration" & data$month == 6 & data$day <= 15)], data$location.lat[which(data$individual.local.identifier == bird[i] & data$migration.stage.standard == "spring-migration"& data$month == 6 & data$day <= 15)], pch=20, type="l")
}






## Modelling the bird movement from winter to breeding (spring migration)

winter.start.location <- data[which(data$individual.local.identifier == bird[2] & data$migration.stage.standard == "spring-migration" & data$migration.stage == "winter"), 4:5] # I just took the winter location of an individual in the dataset for the moment

winter.start.time <- data[which(data$individual.local.identifier == bird[2] & data$migration.stage.standard == "spring-migration" & data$migration.stage == "winter"), 15:17]

breeding.arrival.location <- data[which(data$individual.local.identifier == bird[2] & data$migration.stage.standard == "spring-migration" & data$migration.stage == "breeding"), 4:5] 

breeding.arrival.time <- data[which(data$individual.local.identifier == bird[2] & data$migration.stage.standard == "spring-migration" & data$migration.stage == "breeding"), 15:17]

# get the NDVI value of the focal location and its 8 neighbours (clockwise starting from north)
NDVI.conditions.start.point <- vector()
NDVI.conditions.start.point[1] <- extract(NDVI.2006.feb.2, cbind(48.775, 69.04783))
NDVI.conditions.start.point[2] <- NDVI.2006.feb.2[rowColFromCell(NDVI.2006.feb.2, cellFromXY(NDVI.2006.feb.2, c(48.775, 69.04783)))[1]+1, rowColFromCell(NDVI.2006.feb.2, cellFromXY(NDVI.2006.feb.2, c(48.775, 69.04783)))[2]]
NDVI.conditions.start.point[3] <- NDVI.2006.feb.2[rowColFromCell(NDVI.2006.feb.2, cellFromXY(NDVI.2006.feb.2, c(48.775, 69.04783)))[1]+1, rowColFromCell(NDVI.2006.feb.2, cellFromXY(NDVI.2006.feb.2, c(48.775, 69.04783)))[2]+1]
NDVI.conditions.start.point[4] <- NDVI.2006.feb.2[rowColFromCell(NDVI.2006.feb.2, cellFromXY(NDVI.2006.feb.2, c(48.775, 69.04783)))[1], rowColFromCell(NDVI.2006.feb.2, cellFromXY(NDVI.2006.feb.2, c(48.775, 69.04783)))[2]+1]
NDVI.conditions.start.point[5] <- NDVI.2006.feb.2[rowColFromCell(NDVI.2006.feb.2, cellFromXY(NDVI.2006.feb.2, c(48.775, 69.04783)))[1]-1, rowColFromCell(NDVI.2006.feb.2, cellFromXY(NDVI.2006.feb.2, c(48.775, 69.04783)))[2]+1]
NDVI.conditions.start.point[6] <- NDVI.2006.feb.2[rowColFromCell(NDVI.2006.feb.2, cellFromXY(NDVI.2006.feb.2, c(48.775, 69.04783)))[1]-1, rowColFromCell(NDVI.2006.feb.2, cellFromXY(NDVI.2006.feb.2, c(48.775, 69.04783)))[2]]
NDVI.conditions.start.point[7] <- NDVI.2006.feb.2[rowColFromCell(NDVI.2006.feb.2, cellFromXY(NDVI.2006.feb.2, c(48.775, 69.04783)))[1]-1, rowColFromCell(NDVI.2006.feb.2, cellFromXY(NDVI.2006.feb.2, c(48.775, 69.04783)))[2]-1]
NDVI.conditions.start.point[8] <- NDVI.2006.feb.2[rowColFromCell(NDVI.2006.feb.2, cellFromXY(NDVI.2006.feb.2, c(48.775, 69.04783)))[1], rowColFromCell(NDVI.2006.feb.2, cellFromXY(NDVI.2006.feb.2, c(48.775, 69.04783)))[2]-1]
NDVI.conditions.start.point[9] <- NDVI.2006.feb.2[rowColFromCell(NDVI.2006.feb.2, cellFromXY(NDVI.2006.feb.2, c(48.775, 69.04783)))[1]+1, rowColFromCell(NDVI.2006.feb.2, cellFromXY(NDVI.2006.feb.2, c(48.775, 69.04783)))[2]-1]


