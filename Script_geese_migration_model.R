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


cellInfo <- matrix(ncol=4, nrow=1400700)
colnames(cellInfo) <- c("rowNumber", "colNumber", "longitude", "latitude")
for(i in 1:1400700){
	cellInfo[i,1] <- rowFromCell(NDVI.2006.feb.1, i)
	cellInfo[i,2] <- colFromCell(NDVI.2006.feb.1, i)
	cellInfo[i,3] <- xFromCell(NDVI.2006.feb.1, i)
	cellInfo[i,4] <- yFromCell(NDVI.2006.feb.1, i)
}



## Trajectory similarity

library(utils)

observed.path.1 <- cbind(data$location.long[which(data$individual.local.identifier == bird[2] & data$migration.stage.standard == "spring-migration" & data$month == 3 & data$day > 15)], data$location.lat[which(data$individual.local.identifier == bird[2] & data$migration.stage.standard == "spring-migration"& data$month == 3 & data$day > 15)])

observed.path.2 <- cbind(data$location.long[which(data$individual.local.identifier == bird[3] & data$migration.stage.standard == "spring-migration" & data$month == 3 & data$day > 15)], data$location.lat[which(data$individual.local.identifier == bird[3] & data$migration.stage.standard == "spring-migration"& data$month == 3 & data$day > 15)])


# Function for movement direction

mov.dir <- function(trajectory.point.1, trajectory.point.2){
	if((trajectory.point.2[1] - trajectory.point.1[1]) >= 0){
		dir <- atan((trajectory.point.2[2] - trajectory.point.1[2])/(trajectory.point.2[1] - trajectory.point.1[1]))
	}else if((trajectory.point.2[2] - trajectory.point.1[2]) <= 0 & (trajectory.point.2[1] - trajectory.point.1[1]) < 0){
		dir <- atan((trajectory.point.2[2] - trajectory.point.1[2])/(trajectory.point.2[1] - trajectory.point.1[1])) - pi
	}else{
		dir <- atan((trajectory.point.2[2] - trajectory.point.1[2])/(trajectory.point.2[1] - trajectory.point.1[1])) + pi
	}
	return(dir)
}

movement.dir <- function(trajectory){
	dir.seq <- vector()
	for(i in 1:(nrow(trajectory)-1)){
		dir.seq[i] <- mov.dir(trajectory[i,], trajectory[i+1,])
	}
	return(dir.seq)
}


# Function for distance ratio

total.distance <- function(trajectory){
	total.dist = 0
	for(i in 1:(nrow(trajectory)-1)){
		total.dist <- total.dist + sqrt((trajectory[i+1,1] - trajectory[i,1])^2 + (trajectory[i+1,2] - trajectory[i,2])^2)
	}
	return(total.dist)
}
dista <- function(trajectory.point.1, trajectory.point.2, totalDist){
	if(totalDist > 0){
		d = sqrt((trajectory.point.2[1] - trajectory.point.1[1])^2 + (trajectory.point.2[2] - trajectory.point.1[2])^2)
	}else{
		d = 0
	}
	return(d)
}
distance.ratio.fct <- function(trajectory, totalDist){
	dist.ratio <- vector()
	for(i in 1:(nrow(trajectory)-1)){
		dist.ratio[i] <- dista(trajectory[i,], trajectory[i+1,], totalDist) / totalDist
	}
	return(dist.ratio)
}


movement.direction <- movement.dir(observed.path.1)
total.movement.distance <- total.distance(observed.path.1)
distance.ratio <- distance.ratio.fct(observed.path.1, total.movement.distance)

movement.sequence.1 <- cbind(movement.direction, distance.ratio)


QuantizationMap <- function(moveDir, distRatio){
	if(moveDir < (-3*pi/4) & distRatio < 0.125){ letter = "A" }
	if(moveDir > (-3*pi/4) & moveDir < (-pi/2) & distRatio < 0.125){ letter = "B" }
	if(moveDir > (-pi/2) & moveDir < (-pi/4) & distRatio < 0.125){ letter = "C" }
	if(moveDir > (-pi/4) & moveDir < 0 & distRatio < 0.125){ letter = "D" }
	if(moveDir > 0 & moveDir < (pi/4) & distRatio < 0.125){ letter = "E" }
	if(moveDir > (pi/4) & moveDir < (pi/2) & distRatio < 0.125){ letter = "F" }
	if(moveDir > (pi/2) & moveDir < (3*pi/4) & distRatio < 0.125){ letter = "G" }
	if(moveDir > (3*pi/4) & distRatio < 0.125){ letter = "H" }
	if(moveDir < (-3*pi/4) & distRatio > 0.125 & distRatio < 0.25){ letter = "I" }
	if(moveDir > (-3*pi/4) & moveDir < (-pi/2) & distRatio > 0.125 & distRatio < 0.25){ letter = "J" }
	if(moveDir > (-pi/2) & moveDir < (-pi/4) & distRatio > 0.125 & distRatio < 0.25){ letter = "K" }
	if(moveDir > (-pi/4) & moveDir < 0 & distRatio > 0.125 & distRatio < 0.25){ letter = "L" }
	if(moveDir > 0 & moveDir < (pi/4) & distRatio > 0.125 & distRatio < 0.25){ letter = "M" }
	if(moveDir > (pi/4) & moveDir < (pi/2) & distRatio > 0.125 & distRatio < 0.25){ letter = "N" }
	if(moveDir > (pi/2) & moveDir < (3*pi/4) & distRatio > 0.125 & distRatio < 0.25){ letter = "O" }
	if(moveDir > (3*pi/4) & distRatio > 0.125 & distRatio < 0.25){ letter = "P" }
	if(moveDir < (-3*pi/4) & distRatio > 0.25 & distRatio < 0.375){ letter = "Q" }
	if(moveDir > (-3*pi/4) & moveDir < (-pi/2) & distRatio > 0.25 & distRatio < 0.375){ letter = "R" }
	if(moveDir > (-pi/2) & moveDir < (-pi/4) & distRatio > 0.25 & distRatio < 0.375){ letter = "S" }
	if(moveDir > (-pi/4) & moveDir < 0 & distRatio > 0.25 & distRatio < 0.375){ letter = "T" }
	if(moveDir > 0 & moveDir < (pi/4) & distRatio > 0.25 & distRatio < 0.375){ letter = "U" }
	if(moveDir > (pi/4) & moveDir < (pi/2) & distRatio > 0.25 & distRatio < 0.375){ letter = "V" }
	if(moveDir > (pi/2) & moveDir < (3*pi/4) & distRatio > 0.25 & distRatio < 0.375){ letter = "W" }
	if(moveDir > (3*pi/4) & distRatio > 0.25 & distRatio < 0.375){ letter = "X" }
	if(moveDir < (-3*pi/4) & distRatio > 0.375 & distRatio < 0.5){ letter = "Y" }
	if(moveDir > (-3*pi/4) & moveDir < (-pi/2) & distRatio > 0.375 & distRatio < 0.5){ letter = "Z" }
	if(moveDir > (-pi/2) & moveDir < (-pi/4) & distRatio > 0.375 & distRatio < 0.5){ letter = "a" }
	if(moveDir > (-pi/4) & moveDir < 0 & distRatio > 0.375 & distRatio < 0.5){ letter = "b" }
	if(moveDir > 0 & moveDir < (pi/4) & distRatio > 0.375 & distRatio < 0.5){ letter = "c" }
	if(moveDir > (pi/4) & moveDir < (pi/2) & distRatio > 0.375 & distRatio < 0.5){ letter = "d" }
	if(moveDir > (pi/2) & moveDir < (3*pi/4) & distRatio > 0.375 & distRatio < 0.5){ letter = "e" }
	if(moveDir > (3*pi/4) & distRatio > 0.375 & distRatio < 0.5){ letter = "f" }
	if(moveDir < (-3*pi/4) & distRatio > 0.5 & distRatio < 0.625){ letter = "g" }
	if(moveDir > (-3*pi/4) & moveDir < (-pi/2) & distRatio > 0.5 & distRatio < 0.625){ letter = "h" }
	if(moveDir > (-pi/2) & moveDir < (-pi/4) & distRatio > 0.5 & distRatio < 0.625){ letter = "i" }
	if(moveDir > (-pi/4) & moveDir < 0 & distRatio > 0.5 & distRatio < 0.625){ letter = "j" }
	if(moveDir > 0 & moveDir < (pi/4) & distRatio > 0.5 & distRatio < 0.625){ letter = "k" }
	if(moveDir > (pi/4) & moveDir < (pi/2) & distRatio > 0.5 & distRatio < 0.625){ letter = "l" }
	if(moveDir > (pi/2) & moveDir < (3*pi/4) & distRatio > 0.5 & distRatio < 0.625){ letter = "m" }
	if(moveDir > (3*pi/4) & distRatio > 0.5 & distRatio < 0.625){ letter = "n" }
	if(moveDir < (-3*pi/4) & distRatio > 0.625 & distRatio < 0.75){ letter = "o" }
	if(moveDir > (-3*pi/4) & moveDir < (-pi/2) & distRatio > 0.625 & distRatio < 0.75){ letter = "p" }
	if(moveDir > (-pi/2) & moveDir < (-pi/4) & distRatio > 0.625 & distRatio < 0.75){ letter = "q" }
	if(moveDir > (-pi/4) & moveDir < 0 & distRatio > 0.625 & distRatio < 0.75){ letter = "r" }
	if(moveDir > 0 & moveDir < (pi/4) & distRatio > 0.625 & distRatio < 0.75){ letter = "s" }
	if(moveDir > (pi/4) & moveDir < (pi/2) & distRatio > 0.625 & distRatio < 0.75){ letter = "t" }
	if(moveDir > (pi/2) & moveDir < (3*pi/4) & distRatio > 0.625 & distRatio < 0.75){ letter = "u" }
	if(moveDir > (3*pi/4) & distRatio > 0.625 & distRatio < 0.75){ letter = "v" }
	if(moveDir < (-3*pi/4) & distRatio > 0.75 & distRatio < 0.875){ letter = "w" }
	if(moveDir > (-3*pi/4) & moveDir < (-pi/2) & distRatio > 0.75 & distRatio < 0.875){ letter = "x" }
	if(moveDir > (-pi/2) & moveDir < (-pi/4) & distRatio > 0.75 & distRatio < 0.875){ letter = "y" }
	if(moveDir > (-pi/4) & moveDir < 0 & distRatio > 0.75 & distRatio < 0.875){ letter = "z" }
	if(moveDir > 0 & moveDir < (pi/4) & distRatio > 0.75 & distRatio < 0.875){ letter = "1" }
	if(moveDir > (pi/4) & moveDir < (pi/2) & distRatio > 0.75 & distRatio < 0.875){ letter = "2" }
	if(moveDir > (pi/2) & moveDir < (3*pi/4) & distRatio > 0.75 & distRatio < 0.875){ letter = "3" }
	if(moveDir > (3*pi/4) & distRatio > 0.75 & distRatio < 0.875){ letter = "4" }
	if(moveDir < (-3*pi/4) & distRatio > 0.875){ letter = "5" }
	if(moveDir > (-3*pi/4) & moveDir < (-pi/2) & distRatio > 0.875){ letter = "6" }
	if(moveDir > (-pi/2) & moveDir < (-pi/4) & distRatio > 0.875){ letter = "7" }
	if(moveDir > (-pi/4) & moveDir < 0 & distRatio > 0.875){ letter = "8" }
	if(moveDir > 0 & moveDir < (pi/4) & distRatio > 0.875){ letter = "9" }
	if(moveDir > (pi/4) & moveDir < (pi/2) & distRatio > 0.875){ letter = "0" }
	if(moveDir > (pi/2) & moveDir < (3*pi/4) & distRatio > 0.875){ letter = "+" }
	if(moveDir > (3*pi/4) & distRatio > 0.875){ letter = "-" }
	return(letter)
}

letter.seq <- vector()
for(i in 1:nrow(movement.sequence.1)){
	if(movement.sequence.1[i,1] != "NaN" & movement.sequence.1[i,2] != "NaN"){
		letter.seq[i] <- QuantizationMap(movement.sequence.1[i,1], movement.sequence.1[i,2])
	}
}


st1 <- 
st2 <-
NED <- adist(st1, st2) / max(nchar(st1), nchar(st2))


