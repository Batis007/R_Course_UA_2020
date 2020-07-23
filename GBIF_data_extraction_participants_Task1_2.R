####################################################################
# Eduardo Batista                                                  #
# Project MDR                                                      #
# PhD Program in Biology and Ecology of Global Change BEGC         #
# 04-11-2017                                                       #
####################################################################

#Install packages 
#install.packages(c("sp","raster","dismo","rgdal","maptools","rworldmap"))

#Load the installed packages
library(sp)
library(raster)
library(dismo)
library(rgdal)
library(maptools)
library(rworldmap)
#library(tidyverse)
#  define file paths
main_dir <- "C:/Users/Eduardo/Desktop/Curso_R_2020"
#   Working with occurrence data - GBIF----
help(gbif) # check the gbif function and download data for Diplodia 
data <- gbif()
colnames(data)
head(data)
#Organize the data set by location and remove duplicas + NA cells ----

pointsdatasp <- data.frame(lon = data$lon,lat = data$lat, species =data$species)
pointsdatasp <- pointsdatasp[complete.cases(pointsdatasp),] # remove NA cells
pointsdatasp <- unique(pointsdatasp)

#Extract the associated host by species and countries ----


#Plot WorldMapa with data points ----
data(wrld_simpl)
plot(wrld_simpl, main = "Records of Diplodia spp.", axes=TRUE, col="light yellow")
# restore the box around the map 
box()

# add the points
points(pointsdatasp$lon, pointsdatasp$lat, col='orange', pch=20, cex=0.75)
# plot points again to add a border, for better visibility 
points(pointsdatasp$lon, pointsdatasp$lat, col='red', cex=0.75)

#Homework make a new map with what you learn from Rafael Felix ----

################### TASK2 ----
#Extract from Bioclim data for each record
# make pointsdatasp spatial
points <- data.frame(lon = pointsdatasp$lon, lat = pointsdatasp$lat)
coordinates(points) <- ~lon + lat

# get the worldclim data----
worldclim <- raster::getData('worldclim', var='bio', res=2.5)
rp_wc <- extract(worldclim, points)
# Combine pointsdatasp with rp_wc in a new file rp_wc1

# Check the NA values for modeling----
rp_wc1 <- rp_wc1[complete.cases(rp_wc1),] # remove NA cells

#Analyse max, min and mean annual temperature for Diplodia spp.----

#Plot annual temperature and occurrence location for Diplodia spp.----

#Plot WorldMap with data points
# restore the box around the map
box()
# add the points
points(pointsdatasp$lon, pointsdatasp$lat, col='orange', pch=20, cex=0.75)
# plot points again to add a border, for better visibility 
points(pointsdatasp$lon, pointsdatasp$lat, col='red', cex=0.75)

#Save files----
#setwd(output)
dev.copy(jpeg,filename="worldmapplot.jpeg");
dev.off ()
write.csv(data, file = "data.csv")
write.csv(pointsdatasp, file = "pointsdatasp.csv")
write.csv(rp_wc1, file = "rp_wc1.csv")
