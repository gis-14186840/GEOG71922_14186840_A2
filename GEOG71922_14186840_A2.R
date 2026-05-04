# GEOG71922: Assessment2
# Author: 14186840
# Topic: Beetles

#load the required libraries
library(terra)
library(sf)
library(vegan)
library(dplyr)

library(spdep)
library(blockCV)

library(Hmsc)
library(corrplot)

#load data
lcm_raster=rast("LCMUK_2000.tif")

#remove exported row number column
beetle_env=read.csv("scot_beetle_env.csv",row.names=1)
beetle_comm=read.csv("scot_beetle_community.csv",row.names=1)

#keep only species columns
beetle_comm=beetle_comm %>% select(starts_with("sp"))

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
    ss_mat[,i]=sapply(sp.i,function(x) (x-col.i_mean)^2)}
  return(rowSums(ss_mat)/sum(ss_mat))}

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
#print(head(beetle_env[,c("Sites", "Broadleaf_250m", "Broadleaf_500m", "Broadleaf_1000m","Grassland_250m", "Grassland_500m", "Grassland_1000m")]))
#print(colMeans(beetle_env[,c("Broadleaf_250m", "Broadleaf_500m", "Broadleaf_1000m","Grassland_250m", "Grassland_500m", "Grassland_1000m")]))

#compare scales for landscape variables

#hellinger-transform the community data
comm_hel=decostand(beetle_comm, method="hellinger")

#create an dataframe
scale_compare=data.frame(Scale=paste0(scales, "m"), Variance_mean=NA, R2=NA, p_value=NA)

#loop through each scale
set.seed(123)

for(i in 1:length(scales)){
  #extract landscape variables for the current scale
  current_scale=paste0(scales[i], "m")
  var_names=c(paste0("Broadleaf_", current_scale), paste0("Grassland_", current_scale))
  scale_vars=beetle_env[, var_names]
  
  #calculate mean variance of landscape variables at this scale
  scale_compare$Variance_mean[i]=mean(apply(scale_vars, 2, var))
  
  #run PERMANOVA for this specific scale
  ad=adonis2(comm_hel~.,data=scale_vars, permutations=999, method="bray")
  
  #store R2 and p-value
  scale_compare$R2[i]=round(ad$R2[1], 3)
  scale_compare$p_value[i]=ad$`Pr(>F)`[1]}

print(scale_compare)
#result shows 250m has the highest R2 among the tested scales


#prepare explanatory matrices
local_env=beetle_env[, c("pH", "Moist", "Elevation", "Management")]
land_env=beetle_env[, c("Broadleaf_250m", "Grassland_250m")]
space_env=as.data.frame(scale(beetle_env[, c("X", "Y")]))
names(space_env)=c("X_scaled", "Y_scaled")

#combine all explanatory variables
all_env=cbind(local_env, land_env, space_env)

#NMDS with envfit
set.seed(123)
nmds_res=metaMDS(comm_hel, distance="bray", k=2, trymax=100, autotransform=FALSE, trace=FALSE)

#plot NMDS
plot(nmds_res,type="n",main=paste("NMDS of Beetle Community Structure (Stress =",round(nmds_res$stress, 3),")"))
points(nmds_res, display="sites", pch=21, bg="grey75", col="black", cex=1.2) 
ef=envfit(nmds_res, all_env, permutations=999, na.rm=TRUE)
plot(ef, p.max=0.05, col="blue", cex=0.8)

#marginal PERMANOVA
adonis_res=adonis2(comm_hel~pH+Moist+Elevation+Management+Broadleaf_250m+Grassland_250m+X_scaled+Y_scaled, 
                     data=all_env, permutations=999, method="bray", by="margin")

#print result
print(adonis_res)

#variation partitioning
vp=varpart(comm_hel, local_env, land_env, space_env)
plot(vp, Xnames=c("Local", "Landscape", "Space"), bg=c("cadetblue1", "lightpink", "lightgreen"))
title("Variation Partitioning of Beetle Community")


#poisson GLM for species richness
rich_model=glm(Richness~pH+Moist+Elevation+Management+Broadleaf_250m+Grassland_250m, 
               family=poisson(link="log"), data=beetle_env)
summary(rich_model)

#check overdispersion
dispersion=sum(residuals(rich_model,type="pearson")^2)/rich_model$df.residual
print(dispersion)

#Moran's I test on Pearson residuals
#using k=8 to handle identical overlapping coordinates in dataset
coords=cbind(beetle_env$X, beetle_env$Y)
nb=suppressWarnings(knn2nb(knearneigh(coords, k=8)))
lw=nb2listw(nb, style="W")
moran_test=moran.test(residuals(rich_model, type="pearson"), lw)
print(moran_test)

#spatial block cross-validation
set.seed(42)
spatial_blocks=cv_spatial(x=beetle_sf, k=5, hexagon=TRUE, selection="random",
                          plot=FALSE, progress=FALSE)

#run spatial CV for richness model
fold_ids=spatial_blocks$folds_ids
nfolds=length(unique(fold_ids))
cv_results=data.frame(fold=1:nfolds, RMSE=NA, MAE=NA)

for(i in 1:nfolds){
  train=beetle_env[fold_ids!=i, ]
  test=beetle_env[fold_ids==i, ]
  m_cv=glm(Richness~pH+Moist+Elevation+Management+Broadleaf_250m+Grassland_250m, 
           family=poisson(link="log"), data=train)
  pred=predict(m_cv, newdata=test, type="response")
  cv_results$RMSE[i]=sqrt(mean((test$Richness-pred)^2))
  cv_results$MAE[i]=mean(abs(test$Richness-pred))}

#print results
print(cv_results)
print(paste("Mean spatial CV RMSE:", round(mean(cv_results$RMSE), 3)))
print(paste("Mean spatial CV MAE:", round(mean(cv_results$MAE), 3)))

print(paste("Richness mean:", round(mean(beetle_env$Richness), 2)))
print(paste("Richness SD:", round(sd(beetle_env$Richness), 2)))
print(paste("Richness range:", paste(range(beetle_env$Richness), collapse="-")))

#prepare data for Hmsc (standardize predictors for MCMC convergence)
Y=as.matrix(beetle_comm)
storage.mode(Y)="numeric"
XData=as.data.frame(scale(cbind(local_env, land_env)))
XFormula=~pH+Moist+Elevation+Management+Broadleaf_250m+Grassland_250m

#check species prevalence (informs discussion of rare-species warnings)
prevalence=colSums(beetle_comm>0)
print(paste("Species occurring at <5 sites:", sum(prevalence<5)))
print(prevalence)

#define site as random effect to account for residual co-occurrence
studyDesign=data.frame(site=as.factor(beetle_env$Sites))
rL=HmscRandomLevel(units=studyDesign$site)

#construct JSDM with lognormal Poisson for abundance data
m=Hmsc(Y=Y, XData=XData, XFormula=XFormula, 
       studyDesign=studyDesign, ranLevels=list(site=rL), 
       distr="lognormal poisson")

#fit model with extended MCMC sampling for better convergence
set.seed(123)
m=sampleMcmc(m, samples=2000, thin=20, transient=1000, nChains=2, nParallel=2)

#save fitted model to disk to avoid re-running MCMC
saveRDS(m, "hmsc_model_refined.rds")

#MCMC convergence diagnostics
mpost=convertToCodaObject(m)
ess_beta=effectiveSize(mpost$Beta)
gd_beta=gelman.diag(mpost$Beta, multivariate=FALSE)$psrf[,1]
ess_omega=effectiveSize(mpost$Omega[[1]])
gd_omega=gelman.diag(mpost$Omega[[1]], multivariate=FALSE)$psrf[,1]

#print results
print(paste("Mean ESS (Beta):", round(mean(ess_beta), 1)))
print(paste("Mean Gelman PSRF (Beta):", round(mean(gd_beta), 3)))
print(paste("Mean ESS (Omega):", round(mean(ess_omega), 1)))
print(paste("Mean Gelman PSRF (Omega):", round(mean(gd_omega), 3)))



#explanatory power
preds_expl=computePredictedValues(m)
MF_expl=evaluateModelFit(hM=m, predY=preds_expl)
print(paste("Mean explanatory SR2:", round(mean(MF_expl$SR2, na.rm=TRUE), 3)))

#predictive power via spatial block CV
preds_pred=computePredictedValues(m, partition=spatial_blocks$folds_ids)
MF_pred=evaluateModelFit(hM=m, predY=preds_pred)
print(paste("Mean predictive SR2 (spatial CV):", round(mean(MF_pred$SR2, na.rm=TRUE), 3)))

#save predictions to avoid re-running spatial CV
saveRDS(list(MF_expl=MF_expl, MF_pred=MF_pred), "hmsc_modelfit.rds")

#variance partitioning across environmental groups
VP=computeVariancePartitioning(m, group=c(1,1,1,1,1,2,2), 
                               groupnames=c("Local (incl. intercept)","Landscape"))
plotVariancePartitioning(m, VP=VP)

#residual species co-occurrence matrix
OmegaCor=computeAssociations(m)
supportLevel=OmegaCor[[1]]$support
corMatrix=OmegaCor[[1]]$mean

#filter weak associations
toPlot=corMatrix
toPlot[supportLevel<0.97 & supportLevel>0.03]=0

#plot residual species associations
corrplot(toPlot, method="color", type="lower", tl.col="black", tl.cex=0.7, 
         col=colorRampPalette(c("blue","white","red"))(200), 
         title="Residual species associations", mar=c(0,0,2,0))