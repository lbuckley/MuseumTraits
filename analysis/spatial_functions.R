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

#=================================
absM2= na.omit(absM)
lon= absM2$Long
lat= absM2$Lat
dat= absM2[,c("Corr.Val","JJTave","doy") ]


#FUNCTION TO RUN SPATIAL ANALYIS, error model, 40km neighborhood
spat.mod.var= function(dat, yvar="Corr.Val", xvars=c("JJTave","doy","JJTave:doy"), lon, lat){
  
  #define spatial neighborhood
  coords= as.matrix(cbind(lon,lat))
  ## Build neighborhood
  m.nb40 <- dnearneigh(coords,0,40,longlat=TRUE)
  # establish weights for the given list of neighbors, "B" = binary (all equal), "W" = row standardized
  m.sw40 <- nb2listw(m.nb40, glist=NULL, style="W", zero.policy=TRUE)
  
  #--------
  #MODELS
  m_SAR40= spautolm(as.formula(paste(yvar, "~", paste(xvars, collapse="+"))), data=dat, listw= m.sw40, family = "SAR", method="eigen", na.action='na.pass')
  
  #MODEL SELECTION
  d_SAR40=dredge(m_SAR40)
  
  #extract best 5 models and weights
  best.mods= d_SAR40[1:5,]
  
  ma_SAR40= model.avg(d_SAR40)
  
  #extract coefficients
  name.ord= c("(Intercept)",xvars,"lambda")
  ct= summary(ma_SAR40)$coefmat.full
  match1= match(name.ord, rownames(ct))
  coef.mods=ct[match1,]
  
  #add importance values
  imports= rep(NA, nrow(coef.mods))
  match1=match(names(ma_SAR40$importance), row.names(coef.mods))
  imports[match1]= ma_SAR40$importance
  coef.mods= cbind(coef.mods, imports)
  
  #full model as all terms have some support
  m_x.errorSAR40 <- errorsarlm(as.formula(paste(yvar, "~", paste(xvars, collapse="+"))), data=absM, m.sw40, method="eigen", quiet=TRUE, zero.policy=TRUE, na.action=na.fail, interval=c(-0.5,0.5),tol.solve=1e-25)
 
  mod.sum= summary(m_x.errorSAR40, Nagelkerke=TRUE, Hausman=TRUE)
  
  #extract Z values and pr(>|z|) 
  zs=mod.sum$Coef[,1:4] #z, p
  
  #Log-likelihood ratio tests
  # only use LR.sarlm when nested. otherwise compare Loglikelihood values
  mod.null <- errorsarlm(yvar~1, data=absM, m.sw40, method="eigen", quiet=TRUE, zero.policy=TRUE, na.action=na.fail, interval=c(-0.5,0.5),tol.solve=1e-25)
  lr1= LR.sarlm(m_x.errorSAR40, mod.null)
  
  mod.stats= c(mod.sum$NK,lr1$statistic, lr1$p.value )
  names(mod.stats)[1]= "NKPseudoR2"
  
  return(list(best.mods, coef.mods, zs, mod.stats) )
  
} #end spatial analysis function


#=============================
my.prplot= function (g, i, xlabs)  #beautify plot
{
  xl <- xlabs[i]
  yl <- paste("")
  x <- model.matrix(g)[, i + 1]
  plot(x, g$coeff[i + 1] * x + g$res, xlab = xl, ylab = yl, col=rgb(0.2,0.2, 0.2, 0.5),pch=16)        #@I played with different plot characters to get different plots. Not sold
  abline(0, g$coeff[i + 1])
  invisible()
}
#--------------

#BOOTSTRAP AND DREDGE LINEAR MODEL
#set up data collection
out.mods= array(data=NA, dim=c(5,11,Nruns))
out.coefs=array(data=NA, dim=c(7,5,Nruns))
out.zs=array(data=NA, dim=c(6,4,Nruns))
out.stats= matrix(NA, 3, Nruns)

thorax.out.mods= array(data=NA, dim=c(5,11,Nruns))
thorax.out.coefs=array(data=NA, dim=c(7,5,Nruns))
thorax.out.zs=array(data=NA, dim=c(6,4,Nruns))
thorax.out.stats= matrix(NA, 3, Nruns)

#bootstrap
for(r in 1:Nruns){
  
  #sub sample
  z <- sapply(unique(absM$YrSite), FUN= function(x){ 
    sample(which(absM$YrSite==x), min(Nsamp, length(which(absM$YrSite==x))), FALSE)
  })
  absM.boot<- absM[unlist(z),]
  
  #run model
  mod1= lm(Corr.Val~PupalTave+ParentTave+J+Year+PupalTave*J, data=absM.boot, na.action=na.fail)
  
  #MODEL SELECTION
  d_SAR40=dredge(mod1)
  
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
  
  
}#end bootstrap

##UPDATE OUTPUT
#AVERAGE
coefs= apply(out.coefs, MARGIN=c(1,2), FUN="mean") 
zs= apply(out.zs, MARGIN=c(1,2), FUN="mean") 
stats=  apply(out.stats, MARGIN=c(1), FUN="mean") 

thorax.coefs= apply(thorax.out.coefs, MARGIN=c(1,2), FUN="mean") 
thorax.zs= apply(thorax.out.zs, MARGIN=c(1,2), FUN="mean") 
thorax.stats=  apply(thorax.out.stats, MARGIN=c(1), FUN="mean") 


