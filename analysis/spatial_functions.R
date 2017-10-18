#------------------------
#define functions 

ClassifySites= function(names, lon, lat, DistCut){
  
  new.names=names
  stat.coords= cbind(lon,lat)
  
  for(r in 1:length(lon) ){
    dists= spDistsN1(stat.coords, stat.coords[r,], longlat = TRUE) #find distances from focal site in km
    new.names[which(dists<DistCut)]= new.names[r] #rename sites within cutoff radius
  } #end loop sites
  
  return(new.names) #returns grouped site names
}

count=function(x) length(x)

#===================================================
# SPATIAL ANALYSIS

#FUNCTION TO RUN SPATIAL ANALYIS, error model, 40km neighborhood
spat.mod= function(absM){

#define spatial neighborhood
coords= as.matrix(absM[,c("Long","Lat")])
## Build neighborhood
m.nb40 <- dnearneigh(coords,0,40,longlat=TRUE)
# establish weights for the given list of neighbors, "B" = binary (all equal), "W" = row standardized
m.sw40 <- nb2listw(m.nb40, glist=NULL, style="W", zero.policy=TRUE)

#--------
#GREY MODELS
#Include interaction of J and temps as temps change through season
m_SAR40= spautolm(grey~PupalTave+ParentTave+J+Year+PupalTave*J, data=absM, listw= m.sw40, family = "SAR", method="eigen", na.action='na.pass')
#m_SAR40= spautolm(grey~PupalTave+ParentTave+J+Year+PupalTave*J+ParentTave*J, data=absM, listw= m.sw40, family = "SAR", method="eigen", na.action='na.pass')

#MODEL SELECTION
d_SAR40=dredge(m_SAR40)
# d_SAR40=dredge(m_SAR40, beta = "sd")

#extract best 5 models and weights
best.mods= d_SAR40[1:5,]

ma_SAR40= model.avg(d_SAR40)

#extract coefficients
name.ord= c("(Intercept)","J","Year","ParentTave","PupalTave","J:PupalTave","lambda")
ct= summary(ma_SAR40)$coefmat.full
match1= match(name.ord, rownames(ct))
coef.mods=ct[match1,]

#add importance values
imports= rep(NA, nrow(coef.mods))
match1=match(names(ma_SAR40$importance), row.names(coef.mods))
imports[match1]= ma_SAR40$importance
coef.mods= cbind(coef.mods, imports)

#full model as all terms have some support
m_x.errorSAR40 <- errorsarlm(grey~PupalTave+ParentTave+J+Year+PupalTave*J, data=absM, m.sw40, method="eigen", quiet=TRUE, zero.policy=TRUE, na.action=na.fail, interval=c(-0.5,0.5),tol.solve=1e-25)
#m_x.errorSAR40 <- errorsarlm(grey~PupalTave+ParentTave+J+Year+PupalTave*J+ParentTave*J, data=absM, m.sw40, method="eigen", quiet=TRUE, zero.policy=TRUE, na.action=na.fail, interval=c(-0.5,0.5),tol.solve=1e-25)

mod.sum= summary(m_x.errorSAR40, Nagelkerke=TRUE, Hausman=TRUE)

#extract Z values and pr(>|z|) 
zs=mod.sum$Coef[,1:4] #z, p

#Log-likelihood ratio tests
# only use LR.sarlm when nested. otherwise compare Loglikelihood values
mod.null <- errorsarlm(absM$grey~1, data=absM, m.sw40, method="eigen", quiet=TRUE, zero.policy=TRUE, na.action=na.fail, interval=c(-0.5,0.5),tol.solve=1e-25)
lr1= LR.sarlm(m_x.errorSAR40, mod.null)

mod.stats= c(mod.sum$NK,lr1$statistic, lr1$p.value )
names(mod.stats)[1]= "NKPseudoR2"

#---------------
#THORACIC FUR MODELS

#Include interaction of J and temps as temps change through season
m_SAR40= spautolm(Thorax~PupalTave+ParentTave+J+Year+PupalTave*J, data=absM, listw= m.sw40, family = "SAR", method="eigen", na.action='na.pass')
#m_SAR40= spautolm(Thorax~PupalTave+ParentTave+J+Year+PupalTave*J+ParentTave*J, data=absM, listw= m.sw40, family = "SAR", method="eigen", na.action='na.pass')

#MODEL SELECTION
d_SAR40=dredge(m_SAR40)

#extract best 5 models and weights
thorax.best.mods= d_SAR40[1:5,]

ma_SAR40= model.avg(d_SAR40)

#extract coefficients
ct= summary(ma_SAR40)$coefmat.full
match1= match(name.ord, rownames(ct))
thorax.coef.mods=ct[match1,]

#add importance values
imports= rep(NA, nrow(thorax.coef.mods))
match1=match(names(ma_SAR40$importance), row.names(thorax.coef.mods))
imports[match1]= ma_SAR40$importance
thorax.coef.mods= cbind(thorax.coef.mods, imports)

#full model as all terms have some support
m_x.errorSAR40 <- errorsarlm(Thorax~PupalTave+ParentTave+J+Year+PupalTave*J, data=absM, m.sw40, method="eigen", quiet=TRUE, zero.policy=TRUE, na.action=na.fail, interval=c(-0.5,0.5),tol.solve=1e-25)
#m_x.errorSAR40 <- errorsarlm(Thorax~PupalTave+ParentTave+J+Year+PupalTave*J+ParentTave*J, data=absM, m.sw40, method="eigen", quiet=TRUE, zero.policy=TRUE, na.action=na.fail, interval=c(-0.5,0.5),tol.solve=1e-25)

mod.sum= summary(m_x.errorSAR40, Nagelkerke=TRUE, Hausman=TRUE)

#extract Z values and pr(>|z|) 
thorax.zs=mod.sum$Coef[,1:4] #z, p

#Log-likelihood ratio tests
# only use LR.sarlm when nested. otherwise compare Loglikelihood values
mod.null <- errorsarlm(absM$Thorax~1, data=absM, m.sw40, method="eigen", quiet=TRUE, zero.policy=TRUE, na.action=na.fail, interval=c(-0.5,0.5),tol.solve=1e-25)
lr1= LR.sarlm(m_x.errorSAR40, mod.null)

thorax.mod.stats= c(mod.sum$NK,lr1$statistic, lr1$p.value )
names(thorax.mod.stats)[1]= "NKPseudoR2"

return(list(best.mods, coef.mods, zs, mod.stats, thorax.best.mods, thorax.coef.mods, thorax.zs, thorax.mod.stats) )

} #end spatial analysis function
