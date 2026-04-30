# GEOG71922: Assessment2
# Author: 14186840
# Topic: Beetles

#load the required libraries
library(terra)
library(sf)
library(vegan)
library(dplyr)

#load data
lcm_raster=rast("LCMUK_2000.tif")
beetle_env=read.csv("scot_beetle_env.csv")
beetle_comm=read.csv("scot_beetle_community.csv",row.names=1)

#create spatial points object 
beetle_sf=st_as_sf(beetle_env,coords=c("X","Y"),crs=27700)

#extract land-cover values to points
lcm_extract=terra::extract(lcm_raster, vect(beetle_sf))
beetle_env$Landcover=as.factor(lcm_extract[,2])


#calculate richness and shannon
beetle_env$Richness=rowSums(beetle_comm > 0)
beetle_env$Shannon=diversity(beetle_comm, index="shannon")

#custom function from week 8 to calculate LCBD
function_lcbd=function(spe1){
  ss_mat=spe1
  ss_mat[]=0
  for(i in 1:ncol(spe1)){
    sp.i=spe1[,i]
    col.i_mean=mean(sp.i)
    ss_mat[,i]=sapply(sp.i,function(x) (x-col.i_mean)^2)
  }
  return(rowSums(ss_mat)/sum(ss_mat))
}
beetle_env$LCBD=function_lcbd(beetle_comm)

#test result
#print(table(beetle_env$Landcover,useNA="ifany"))


#broadleaf and grassland cover extraction (250m, 500m, 1000m)
scales=c(250,500,1000)

#250m buffer extraction
buffer_250=st_buffer(beetle_sf, dist=250)
lcm_250=terra::extract(lcm_raster, vect(buffer_250))
names(lcm_250)[2]= "LC_Class"
land_250=lcm_250 %>% group_by(ID) %>% 
  summarise(Broadleaf_250m=mean(LC_Class==1, na.rm = TRUE),
            Grassland_250m=mean(LC_Class==4, na.rm=TRUE))

#500m buffer extraction
buffer_500 = st_buffer(beetle_sf, dist=500)
lcm_500 = terra::extract(lcm_raster, vect(buffer_500))
names(lcm_500)[2] = "LC_Class"
land_500= lcm_500 %>% group_by(ID) %>% 
  summarise(Broadleaf_500m = mean(LC_Class==1, na.rm = TRUE),
            Grassland_500m=mean(LC_Class==4, na.rm=TRUE))

#1000m buffer extraction
buffer_1000=st_buffer(beetle_sf, dist=1000)
lcm_1000=terra::extract(lcm_raster, vect(buffer_1000))
names(lcm_1000)[2]= "LC_Class"
land_1000= lcm_1000 %>% group_by(ID) %>% 
  summarise(Broadleaf_1000m= mean(LC_Class==1,na.rm=TRUE),
            Grassland_1000m=mean(LC_Class==4, na.rm=TRUE))

#merge extracted metrics back to environmental data
beetle_env$ID=1:nrow(beetle_env)
beetle_env=beetle_env %>%
  left_join(land_250, by="ID") %>%
  left_join(land_500, by="ID") %>%
  left_join(land_1000, by="ID")

#test results
print(head(beetle_env[, c("Sites", "Broadleaf_250m", "Broadleaf_500m", "Broadleaf_1000m",
                          "Grassland_250m", "Grassland_500m", "Grassland_1000m")]))
print(colMeans(beetle_env[, c("Broadleaf_250m", "Broadleaf_500m", "Broadleaf_1000m",
                              "Grassland_250m", "Grassland_500m", "Grassland_1000m")]))
