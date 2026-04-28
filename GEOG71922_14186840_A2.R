# GEOG71922: Assessment2
# Author: 14186840
# Topic: Beetles

#load the required libraries
library(terra)
library(sf)
library(vegan)

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
print(table(beetle_env$Landcover,useNA="ifany"))
