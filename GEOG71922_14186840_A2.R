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
library(foreach)
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

#align points CRS to raster to suppress on-the-fly
beetle_sf=st_transform(beetle_sf, crs(lcm_raster))

#calculate richness
beetle_env$Richness=rowSums(beetle_comm > 0)

#local contributions to beta diversity (Legendre & De Caceres 2013)
#squared deviations from species means, normalized by total SS
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


#broadleaf and grassland cover extraction (250m, 500m, 1000m)
extract_landscape=function(scale_m){
  buf=st_buffer(beetle_sf, dist=scale_m)
  ext=terra::extract(lcm_raster, vect(buf))
  names(ext)[2]="LC_Class"
  out=ext %>% group_by(ID) %>%
    summarise(Broadleaf=mean(LC_Class==1, na.rm=TRUE),
              Grassland=mean(LC_Class==4, na.rm=TRUE))
  names(out)[-1]=paste0(names(out)[-1], "_", scale_m, "m")
  out}

#merge extracted metrics back to environmental data
beetle_env$ID=1:nrow(beetle_env)
beetle_env=beetle_env %>%
  left_join(extract_landscape(250),  by="ID") %>%
  left_join(extract_landscape(500),  by="ID") %>%
  left_join(extract_landscape(1000), by="ID")


#compare scales for landscape variables
#hellinger-transform the community data
comm_hel=decostand(beetle_comm, method="hellinger")

#create an dataframe
scales=c(250, 500, 1000)
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
#use k=3 to reduce stress towards the Clarke (1993) acceptable range
#only the first two axes are visualized
set.seed(123)
nmds_res=metaMDS(comm_hel, distance="bray", k=3, trymax=200, 
                 autotransform=FALSE, trace=FALSE)

#plot first two NMDS axes
plot(nmds_res, type="n", choices=c(1,2),
     main=paste("NMDS of Beetle Community (Stress =", round(nmds_res$stress, 3), ")"))
points(nmds_res, display="sites", pch=21, bg="grey75", col="black", cex=1.2)
ef=envfit(nmds_res, all_env, permutations=999, choices=c(1,2), na.rm=TRUE)
plot(ef, p.max=0.05, col="blue", cex=0.8)

#marginal PERMANOVA
adonis_res=adonis2(comm_hel~pH+Moist+Elevation+Management+Broadleaf_250m+Grassland_250m+X_scaled+Y_scaled, 
                     data=all_env, permutations=999, method="bray", by="margin")

#print result
print(as.matrix(adonis_res)[1:8, c("R2","F","Pr(>F)")])

#variation partitioning
vp=varpart(comm_hel, local_env, land_env, space_env)
plot(vp, Xnames=c("Local", "Landscape", "Space"), bg=c("cadetblue1", "lightpink", "lightgreen"), 
     digits=2)
title("Variation Partitioning of Beetle Community")

#Moran's I test on Pearson residuals
#using k=8 to handle identical overlapping coordinates in dataset
coords=cbind(beetle_env$X, beetle_env$Y)
nb=suppressWarnings(knn2nb(knearneigh(coords, k=8)))
lw=nb2listw(nb, style="W")

#LCBD as response
#logit-transform LCBD (bounded [0,1]) for Gaussian regression
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

#under-dispersion detected (ratio<<1); refit with quasipoisson for valid SEs
rich_model_qp=glm(Richness~pH+Moist+Elevation+Management+Broadleaf_250m+Grassland_250m,
                  family=quasipoisson(link="log"), data=beetle_env)
print(round(coef(summary(rich_model_qp))[,c(1,4)], 3))

#residual spatial autocorrelation
moran_rich=moran.test(residuals(rich_model, type="pearson"), lw)
print(paste("Moran I (Richness): I =", round(moran_rich$estimate[1],3),
            "p =", round(moran_rich$p.value,3)))

#spatial block cross-validation
set.seed(42)
spatial_blocks=cv_spatial(x=beetle_sf, k=5, hexagon=TRUE, selection="random",
                          plot=FALSE, progress=FALSE)
fold_ids=spatial_blocks$folds_ids
nfolds=length(unique(fold_ids))

#visualize spatial folds
plot(lcm_raster, main="Spatial Cross-Validation Folds over Land Cover", legend=FALSE, axes=FALSE)
fold_colors = c("red", "blue", "green", "yellow", "purple")
plot(st_geometry(beetle_sf), 
     col=fold_colors[as.numeric(fold_ids)], 
     pch=19, cex=1.5, add=TRUE)
legend("topright", legend=paste("Fold", 1:5), pch=19, 
       col=fold_colors, cex=0.8, bg="white", xpd=TRUE)

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

#print results
print(cv_results)
print(paste("Mean spatial CV RMSE:", round(mean(cv_results$RMSE), 3)))
print(paste("Mean spatial CV MAE:", round(mean(cv_results$MAE), 3)))
print(paste("Richness mean:", round(mean(beetle_env$Richness), 2)))
print(paste("Richness SD:", round(sd(beetle_env$Richness), 2)))
print(paste("Richness range:", paste(range(beetle_env$Richness), collapse="-")))

#non-spatial CV baseline for comparison
#Roberts et al. (2017): random CV inflates predictive performance under spatial autocorrelation
set.seed(42)
random_folds=sample(rep(1:nfolds, length.out=nrow(beetle_env)))
cv_random=sapply(1:nfolds, function(i){
  m_r=glm(Richness~pH+Moist+Elevation+Management+Broadleaf_250m+Grassland_250m,
          family=quasipoisson(link="log"), data=beetle_env[random_folds!=i,])
  pred=predict(m_r, newdata=beetle_env[random_folds==i,], type="response")
  sqrt(mean((beetle_env$Richness[random_folds==i]-pred)^2))})
print(paste("Random CV RMSE:", round(mean(cv_random),3),
            "| Spatial CV RMSE:", round(mean(cv_results$RMSE),3),
            "| spatial is", round((mean(cv_results$RMSE)/mean(cv_random)-1)*100,1), "% higher"))

#prepare data for Hmsc
#standardize continuous predictors only
Y=as.matrix(beetle_comm)
storage.mode(Y)="numeric"

#scale continuous variables for MCMC convergence
cont_vars=as.data.frame(scale(beetle_env[,c("pH","Moist","Elevation",
                                            "Broadleaf_250m","Grassland_250m")]))

#combine with Management as factor
XData=data.frame(cont_vars, Management=as.factor(beetle_env$Management))
XFormula=~pH+Moist+Elevation+Management+Broadleaf_250m+Grassland_250m

#check species prevalence
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

#fit model with extended MCMC sampling
#following Tikhonov et al. (2019): 4 chains, thin=100, transient=10000
#total about 110000 iterations per chain, balances Omega convergence and compute
set.seed(123)
m=sampleMcmc(m, samples=1000, thin=100, transient=10000, nChains=4, nParallel=4)

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
cl=makeCluster(5)
registerDoParallel(cl)
fold_predictions=foreach(i=1:nfolds, .packages=c("Hmsc")) %dopar% {
  computePredictedValues(m, partition=as.numeric(fold_ids==i), nParallel=1)}
stopCluster(cl)

#combine predictions from all folds into single array
preds_pred=array(NA, dim=dim(preds_expl))
for(i in 1:nfolds){
  idx=fold_ids==i
  preds_pred[idx,,]=fold_predictions[[i]][idx,,]}

MF_pred=evaluateModelFit(hM=m, predY=preds_pred)
print(paste("Mean predictive SR2 (spatial CV):", round(mean(MF_pred$SR2, na.rm=TRUE), 3)))

#save predictions to avoid re-running spatial CV
saveRDS(list(MF_expl=MF_expl, MF_pred=MF_pred), "hmsc_modelfit.rds")

#variance partitioning across environmental groups
VP=computeVariancePartitioning(m, group=c(1,1,1,1,1,2,2), 
                               groupnames=c("Local (incl. intercept)","Landscape"))

#plot result
plotVariancePartitioning(m,VP=VP,args.legend=list(x="bottomright",inset=c(0, 0.05),cex=0.8))

#residual species co-occurrence matrix
OmegaCor=computeAssociations(m)
supportLevel=OmegaCor[[1]]$support
corMatrix=OmegaCor[[1]]$mean

#diagnostic across thresholds (used to justify the chosen support cut-off)
diag_idx=lower.tri(supportLevel)
n_strong=sum(supportLevel[diag_idx]>=0.95 | supportLevel[diag_idx]<=0.05)
n_mod=sum(supportLevel[diag_idx]>=0.90 | supportLevel[diag_idx]<=0.10)
n_weak=sum(supportLevel[diag_idx]>=0.85 | supportLevel[diag_idx]<=0.15)
print(paste("Associations:strong(>=0.95):",n_strong,"moderate(>=0.90):",n_mod,"weak(>=0.85):",n_weak))

#filter weak associations
toPlot=corMatrix
toPlot[supportLevel<0.95 & supportLevel>0.05]=0

#plot residual species associations
corrplot(toPlot, method="color", type="lower", tl.col="black", tl.cex=0.7, 
         col=colorRampPalette(c("blue","white","red"))(200), 
         title="Residual species associations", mar=c(0,0,2,0))