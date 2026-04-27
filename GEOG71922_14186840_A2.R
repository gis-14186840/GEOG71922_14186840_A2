# GEOG71922: Assessment2
# Author: 14186840

#load data
lcm_raster=rast("LCMUK_2000.tif")
beetle_env=read.csv("scot_beetle_env.csv")
beetle_comm=read.csv("scot_beetle_community.csv",row.names=1)

#create spatial points object 
beetle_sf=st_as_sf(beetle_env,coords=c("X","Y"),crs=27700)

#extract land-cover values to points
lcm_extract=terra::extract(LCM, vect(beetle_sf))
beetle_env$Landcover=as.factor(lcm_extract[,2])
