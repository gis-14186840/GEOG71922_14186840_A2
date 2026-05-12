# GEOG71922: Assessment2
# Author: 14186840
# Topic: Beetles
# note: place input files in the working directory

# 1.Setup and loading data
#load the required libraries
library(terra)  #raster data
library(sf)   #spatial vector data
library(vegan)  #diversity and ordination
library(dplyr)  #data manipulation
library(spdep)  #Moran's I
library(blockCV) #spatial cross-validation
library(Hmsc)  #joint species distribution modelling
library(corrplot) #association matrix plots
library(foreach) #parallel loops
library(doParallel)

#load data
lcm_raster=rast("LCMUK_2000.tif")

#remove exported row number column
beetle_env=read.csv("scot_beetle_env.csv",row.names=1)
beetle_comm=read.csv("scot_beetle_community.csv",row.names=1)

#keep only species columns
beetle_comm=beetle_comm %>% select(starts_with("sp"))

#create spatial points object 
beetle_sf=st_as_sf(beetle_env,coords=c("X","Y"),crs=27700)

#align points CRS to raster
beetle_sf=st_transform(beetle_sf, crs(lcm_raster))

# 2.Community diversity metrics
#calculate richness
beetle_env$Richness=rowSums(beetle_comm > 0)

#calculate LCBD (Legendre & De Caceres 2013)
#use squared deviations from species means
function_lcbd=function(spe){
  spe_mat=as.matrix(spe)
  ss=sweep(spe_mat, 2, colMeans(spe_mat), "-")^2
  rowSums(ss)/sum(ss)}

beetle_env$LCBD=function_lcbd(beetle_comm)

#plot richness distribution
hist(beetle_env$Richness, breaks=10, col="grey80",
     main="Distribution of beetle species richness",
     xlab="Species richness per site",
     xlim=c(0, max(beetle_env$Richness) + 1))
abline(v=mean(beetle_env$Richness), lwd=2, lty=2)

# 3.Landscape variable extraction
#extract broadleaf and grassland proportions within each buffer
extract_landscape=function(scale_m){
  buf=st_buffer(beetle_sf, dist=scale_m)
  ext=terra::extract(lcm_raster, vect(buf))
  names(ext)[2]="LC_Class"
  out=ext %>% group_by(ID) %>%
    summarise(Broadleaf=mean(LC_Class==1, na.rm=TRUE),
              Grassland=mean(LC_Class==4, na.rm=TRUE))
  names(out)[-1]=paste0(names(out)[-1], "_", scale_m, "m")
  out}

#add extracted metrics to environmental data
beetle_env$ID=1:nrow(beetle_env)
beetle_env=beetle_env %>%
  left_join(extract_landscape(250),  by="ID") %>%
  left_join(extract_landscape(500),  by="ID") %>%
  left_join(extract_landscape(1000), by="ID")

# 4.Scale comparison
#hellinger-transform community data
comm_hel=decostand(beetle_comm, method="hellinger")

#create dataframe for scale comparison
scales=c(250, 500, 1000)
scale_compare=data.frame(Scale=paste0(scales, "m"), Variance_mean=NA, R2=NA, p_value=NA)

#loop over the three buffer scales
set.seed(123)

for(i in 1:length(scales)){
  #get landscape variables for the current scale
  current_scale=paste0(scales[i], "m")
  var_names=c(paste0("Broadleaf_", current_scale), paste0("Grassland_", current_scale))
  scale_vars=beetle_env[, var_names]
  
  #calculate mean variance at this scale
  scale_compare$Variance_mean[i]=mean(apply(scale_vars, 2, var))
  
  #run PERMANOVA for this scale
  ad=adonis2(comm_hel~.,data=scale_vars, permutations=999, method="bray")
  
  scale_compare$R2[i]=round(ad$R2[1], 3)
  scale_compare$p_value[i]=ad$`Pr(>F)`[1]}

print(scale_compare)
#result shows 250m has the highest R2 among the tested scales

# 5.Ordination and PERMANOVA
#prepare explanatory matrices
local_env=beetle_env[, c("pH", "Moist", "Elevation", "Management")]
land_env=beetle_env[, c("Broadleaf_250m", "Grassland_250m")]
space_env=as.data.frame(scale(beetle_env[, c("X", "Y")]))
names(space_env)=c("X_scaled", "Y_scaled")

#combine predictors for envfit and PERMANOVA
all_env=cbind(local_env, land_env, space_env)

#NMDS with envfit
#use k=3 to keep stress lower
set.seed(123)
nmds_res=metaMDS(comm_hel, distance="bray", k=3, trymax=200, 
                 autotransform=FALSE, trace=FALSE)

#plot first two NMDS axes
plot(nmds_res, type="n", choices=c(1,2),
     main=paste("NMDS of Beetle Community (Stress =", round(nmds_res$stress, 3), ")"))
points(nmds_res, display="sites", pch=21, bg="grey75", col="black", cex=1.2)
ef=envfit(nmds_res, all_env, permutations=999, choices=c(1,2), na.rm=TRUE)
plot(ef, p.max=0.05, col="blue", cex=0.8)

#marginal PERMANOVA tests each predictor after accounting for the others
adonis_res=adonis2(comm_hel~pH+Moist+Elevation+Management+Broadleaf_250m+Grassland_250m+X_scaled+Y_scaled, 
                     data=all_env, permutations=999, method="bray", by="margin")

print(as.matrix(adonis_res)[seq_len(nrow(adonis_res)-2), c("R2","F","Pr(>F)")])

# 6.Variation partitioning
#partition community variation among local, landscape, and spatial groups
vp=varpart(comm_hel, local_env, land_env, space_env)
plot(vp, Xnames=c("Local", "Landscape", "Space"), bg=c("cadetblue1", "lightpink", "lightgreen"), 
     digits=2)
title("Variation Partitioning of Beetle Community")

# 7.GLMs and spatial autocorrelation
#create k-nearest-neighbour weights
#use k=8 to handle repeated coordinates
coords=cbind(beetle_env$X, beetle_env$Y)
nb=suppressWarnings(knn2nb(knearneigh(coords, k=8)))
lw=nb2listw(nb, style="W")

#logit-transform LCBD before Gaussian GLM
beetle_env$LCBD_logit=log(beetle_env$LCBD/(1-beetle_env$LCBD))
lcbd_model=glm(LCBD_logit~pH+Moist+Elevation+Management+Broadleaf_250m+Grassland_250m,
               family=gaussian(), data=beetle_env)
print(round(coef(summary(lcbd_model))[,c(1,4)],3))

#test residual spatial autocorrelation
moran_lcbd=moran.test(residuals(lcbd_model), lw)
print(paste("Moran I (LCBD): I =", round(moran_lcbd$estimate[1],3),
            "p =", round(moran_lcbd$p.value,3)))

#poisson GLM for species richness
rich_model=glm(Richness~pH+Moist+Elevation+Management+Broadleaf_250m+Grassland_250m, 
               family=poisson(link="log"), data=beetle_env)

#check Poisson dispersion
dispersion=sum(residuals(rich_model,type="pearson")^2)/rich_model$df.residual
print(paste("Dispersion ratio:", round(dispersion, 3)))

#refit with quasipoisson to account for non-Poisson dispersion
rich_model_qp=glm(Richness~pH+Moist+Elevation+Management+Broadleaf_250m+Grassland_250m,
                  family=quasipoisson(link="log"), data=beetle_env)
print(round(coef(summary(rich_model_qp))[,c(1,4)], 3))

#residual spatial autocorrelation
moran_rich=moran.test(residuals(rich_model, type="pearson"), lw)
print(paste("Moran I (Richness): I =", round(moran_rich$estimate[1],3),
            "p =", round(moran_rich$p.value,3)))

# 8.Spatial cross-validation
#create spatial CV folds for richness model evaluation
set.seed(42)
spatial_blocks=cv_spatial(x=beetle_sf, k=5, hexagon=TRUE, selection="random",
                          plot=FALSE, progress=FALSE)
fold_ids=spatial_blocks$folds_ids
nfolds=length(unique(fold_ids))

#visualize spatial folds
plot(lcm_raster, main="Spatial Cross-Validation Folds over Land Cover", legend=FALSE, axes=FALSE)
fold_colors=c("red", "blue", "green", "yellow", "purple")
plot(st_geometry(beetle_sf), 
     col=fold_colors[as.numeric(fold_ids)], 
     pch=19, cex=1.5, add=TRUE)
legend("topright", legend=paste("Fold", 1:5), pch=19, 
       col=fold_colors, cex=0.8, bg="white", xpd=TRUE)

#check number of unique coordinate clusters
coord_id=paste(beetle_env$X, beetle_env$Y, sep="_")
n_coord=length(unique(coord_id))
print(paste("Unique coordinate clusters:", n_coord))

#run spatial CV for richness model
cv_results=data.frame(fold=1:nfolds, RMSE=NA, MAE=NA)

for(i in 1:nfolds){
  train=beetle_env[fold_ids!=i, ]
  test=beetle_env[fold_ids==i, ]
  m_cv=glm(Richness~pH+Moist+Elevation+Management+Broadleaf_250m+Grassland_250m, 
           family=quasipoisson(link="log"), data=train)
  pred=predict(m_cv, newdata=test, type="response")
  cv_results$RMSE[i]=sqrt(mean((test$Richness-pred)^2))
  cv_results$MAE[i]=mean(abs(test$Richness-pred))}

print(cv_results)
print(paste("Mean spatial CV RMSE:", round(mean(cv_results$RMSE), 3)))
print(paste("Mean spatial CV MAE:", round(mean(cv_results$MAE), 3)))
print(paste("Richness mean:", round(mean(beetle_env$Richness), 2)))
print(paste("Richness SD:", round(sd(beetle_env$Richness), 2)))
print(paste("Richness range:", paste(range(beetle_env$Richness), collapse="-")))

#random-fold baseline for comparison (Roberts et al. 2017)
set.seed(42)
random_folds=sample(rep(1:nfolds, length.out=nrow(beetle_env)))
cv_random=sapply(1:nfolds, function(i){
  m_r=glm(Richness~pH+Moist+Elevation+Management+Broadleaf_250m+Grassland_250m,
          family=quasipoisson(link="log"), data=beetle_env[random_folds!=i,])
  pred=predict(m_r, newdata=beetle_env[random_folds==i,], type="response")
  sqrt(mean((beetle_env$Richness[random_folds==i]-pred)^2))})
cv_diff=(mean(cv_results$RMSE)/mean(cv_random)-1)*100
print(paste("Random CV RMSE:", round(mean(cv_random),3),
            "Spatial CV RMSE:", round(mean(cv_results$RMSE),3),
            "spatial CV is", abs(round(cv_diff,1)),ifelse(cv_diff>=0,"% higher","% lower")))

# 9.HMSC model setup
#prepare response and predictor data
Y=as.matrix(beetle_comm)
storage.mode(Y)="numeric"

#scale continuous variables for MCMC convergence
cont_vars=as.data.frame(scale(beetle_env[,c("pH","Moist","Elevation",
                                            "Broadleaf_250m","Grassland_250m")]))

#treat Management as a factor in the main HMSC model
XData=data.frame(cont_vars, Management=as.factor(beetle_env$Management))
XFormula=~pH+Moist+Elevation+Management+Broadleaf_250m+Grassland_250m

#check species prevalence
prevalence=colSums(beetle_comm>0)
print(paste("Species occurring at <5 sites:", sum(prevalence<5)))
print(prevalence)

#use site as a random effect
studyDesign=data.frame(site=as.factor(beetle_env$Sites))
rL=HmscRandomLevel(units=studyDesign$site)

#construct JSDM with lognormal Poisson for abundance data
m=Hmsc(Y=Y, XData=XData, XFormula=XFormula, 
       studyDesign=studyDesign, ranLevels=list(site=rL), 
       distr="lognormal poisson")

#fit model or load saved model
#MCMC settings follow Tikhonov et al. (2019)
model_file="hmsc_model_refined.rds"

if(file.exists(model_file)){
  m=readRDS(model_file)
  print("Loaded saved Hmsc model")
}else{
  set.seed(123)
   m=sampleMcmc(m, samples=1000, thin=100, transient=10000, nChains=4, nParallel=4)
   saveRDS(m, model_file)}

# 10.HMSC diagnostics and outputs
#MCMC convergence diagnostics
mpost=convertToCodaObject(m)
ess_beta=effectiveSize(mpost$Beta)
gd_beta=gelman.diag(mpost$Beta, multivariate=FALSE)$psrf[,1]
ess_omega=effectiveSize(mpost$Omega[[1]])
gd_omega=gelman.diag(mpost$Omega[[1]], multivariate=FALSE)$psrf[,1]

print(paste("Mean ESS (Beta):", round(mean(ess_beta), 1)))
print(paste("Mean Gelman PSRF (Beta):", round(mean(gd_beta), 3)))
print(paste("Mean ESS (Omega):", round(mean(ess_omega), 1)))
print(paste("Mean Gelman PSRF (Omega):", round(mean(gd_omega), 3)))

#calculate explanatory and predictive fit
fit_file="hmsc_modelfit.rds"

if(file.exists(fit_file)){
  mf=readRDS(fit_file)
  MF_expl=mf$MF_expl
  MF_pred=mf$MF_pred
  print("Loaded saved Hmsc model fit")
}else{

  #explanatory power
  preds_expl=computePredictedValues(m)
  MF_expl=suppressWarnings(evaluateModelFit(hM=m, predY=preds_expl))

  #predictive power via spatial block CV
  cl=makeCluster(5)
  registerDoParallel(cl)
  fold_predictions=foreach(i=1:nfolds, .packages=c("Hmsc")) %dopar% {
    computePredictedValues(m, partition=as.numeric(fold_ids==i), nParallel=1)}
  stopCluster(cl)
  
  #combine fold predictions
  preds_pred=array(NA, dim=dim(preds_expl))
  for(i in 1:nfolds){
    idx=fold_ids==i
    preds_pred[idx,,]=fold_predictions[[i]][idx,,]}

  MF_pred=suppressWarnings(evaluateModelFit(hM=m, predY=preds_pred))

  #save model fit results
  saveRDS(list(MF_expl=MF_expl, MF_pred=MF_pred),fit_file)}

print(paste("Mean explanatory SR2:", round(mean(MF_expl$SR2, na.rm=TRUE), 3)))
print(paste("Mean predictive SR2 (spatial CV):", round(mean(MF_pred$SR2, na.rm=TRUE), 3)))

#partition HMSC variance between local and landscape predictors
#suppress repeated SD=0 warnings
VP=suppressWarnings(computeVariancePartitioning(m, group=c(1,1,1,1,1,2,2),
                                                groupnames=c("Local (incl. intercept)","Landscape")))

#plot per-species variance partitioning bars
plotVariancePartitioning(m,VP=VP,args.legend=list(x="bottomright",inset=c(0, 0.05),cex=0.8))

#residual species association matrix
OmegaCor=suppressWarnings(computeAssociations(m))
supportLevel=OmegaCor[[1]]$support
corMatrix=OmegaCor[[1]]$mean

#count supported associations at different thresholds
diag_idx=lower.tri(supportLevel)
n_strong=sum(supportLevel[diag_idx]>=0.95 | supportLevel[diag_idx]<=0.05)
n_mod=sum(supportLevel[diag_idx]>=0.90 | supportLevel[diag_idx]<=0.10)
n_weak=sum(supportLevel[diag_idx]>=0.85 | supportLevel[diag_idx]<=0.15)
print(paste("Associations:strong(>=0.95):",n_strong,"moderate(>=0.90):",n_mod,"weak(>=0.85):",n_weak))

#mask unsupported associations
toPlot=corMatrix
toPlot[supportLevel<0.95 & supportLevel>0.05]=NA
diag(toPlot)=NA

#plot supported residual species associations
corrplot(toPlot, method="color", type="lower", diag=FALSE,
         na.label=" ", tl.col="black", tl.cex=0.7, 
         col=colorRampPalette(c("blue","white","red"))(200), 
         title="Residual species associations", mar=c(0,0,2,0))

# 11.Sensitivity test
#optional check: treat Management as a continuous ordinal predictor
#set TRUE to rerun the sensitivity test
run_sensitivity=FALSE

if(run_sensitivity){
  
  #use continuous ordinal Management instead of factor
  XData_cont=data.frame(cont_vars,
                        Management=as.numeric(scale(beetle_env$Management)))
  
  m_cont=Hmsc(Y=Y, XData=XData_cont,
              XFormula=~pH+Moist+Elevation+Management+Broadleaf_250m+Grassland_250m,
              studyDesign=studyDesign, ranLevels=list(site=rL),
              distr="lognormal poisson")
  
  model_cont_file="hmsc_model_management_continuous.rds"
  
  if(file.exists(model_cont_file)){
    m_cont=readRDS(model_cont_file)
    print("Loaded saved continuous Management Hmsc model")
  }else{
    set.seed(123)
    m_cont=sampleMcmc(m_cont, samples=1000, thin=100, transient=10000,
                      nChains=4, nParallel=4)
    saveRDS(m_cont, model_cont_file)}
  
  #convergence diagnostics
  mp_cont=convertToCodaObject(m_cont)
  gd_beta_cont=gelman.diag(mp_cont$Beta, multivariate=FALSE)$psrf[,1]
  print(paste("Sensitivity Beta PSRF:", round(mean(gd_beta_cont), 3)))
  
  #explanatory power
  preds_cont=computePredictedValues(m_cont)
  MF_cont=suppressWarnings(evaluateModelFit(hM=m_cont, predY=preds_cont))
  print(paste("Sensitivity explanatory SR2:", round(mean(MF_cont$SR2, na.rm=TRUE), 3)))
  
  #count supported residual associations under the continuous Management model
  OC_cont=suppressWarnings(computeAssociations(m_cont))
  support_cont=OC_cont[[1]]$support
  idx=lower.tri(support_cont)
  n95_cont=sum(support_cont[idx]>=0.95 | support_cont[idx]<=0.05)
  print(paste("Sensitivity associations at >=0.95 support:", n95_cont))}