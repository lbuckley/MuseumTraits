########## INCLUDES OTHER DATASET

### ANALYSIS OF COLIAS ABSORPTIVITY DATA

#load libraries
library(sp)
library(Rmisc)
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
#---------------------------

#read data
setwd("C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\Butterflies\\Absorptivity\\data\\")
abs=read.table("Region1.SummerTemp.txt", header = T)
#abs=read.csv("region1.raw.newloc.lb.csv") #PREVIOUS DATASET

absM<- subset(abs,Sex =="M")

#group nearby localities
absM$NewLocation<- ClassifySites(absM$Location, absM$Long, absM$Lat, DistCut=5)    #radius of 5km

abs=absM
#-------------------------

#EXAMINE SINGLE SITES WITH ABUNDANT DATA

#look at counts across years
abs1= aggregate(abs, list(abs$Year, abs$NewLocation), FUN="mean", na.rm=TRUE)
colnames(abs1)[1:2]= c("Year", "NewLocation")
abs1.count= aggregate(abs, list(abs$NewLocation, abs$Year), FUN="count")
abs1.count= abs1.count[order(abs1.count$Group.1),]

#Select sites with good coverage
abs.sub= abs[abs$NewLocation=="Loveland_Pass",]
abs1.sub= abs1[abs1$NewLocation=="Loveland_Pass",]
#"4_July_Picnic_Grounds_Hessie" good past representation in old data, different name or cut from new?

#plot aggregated
plot(abs1.sub$Year, abs1.sub$Grey)
abline(lm(abs1.sub$Grey~abs1.sub$Year))

plot(abs1.sub$JulyMax, abs1.sub$Grey)
abline(lm(abs1.sub$Grey~abs1.sub$JulyMax))

plot(abs1.sub$JuneAv, abs1.sub$Grey)
abline(lm(abs1.sub$Grey~abs1.sub$JuneAv))

#plot without aggregation
plot(abs.sub$Year, abs.sub$Grey)
abline(lm(abs.sub$Grey~abs.sub$Year))

#LOAD DIFFERENT DATA SET, ALL REGIONS?
abs<-read.table ("ForMuseumAna2.txt", header=T)
abs$Grey<-as.numeric(1-abs$CV)
absM<- subset(abs,Sex =="M")

#group nearby localities
#omit NaN
absM<- absM[which(!is.nan(absM$Lat)),]
absM$NewLocation<- ClassifySites(absM$Location, absM$Long, absM$Lat, DistCut=5)    #radius of 5km
abs=absM

abs1= aggregate(abs, list(abs$Year, abs$Location), FUN="mean", na.rm=TRUE)
colnames(abs1)[1:2]= c("Year", "Location")
abs1.count= aggregate(abs, list(abs$Location, abs$Year), FUN="count")
abs1.count= abs1.count[order(abs1.count$Group.1),]

#Select sites with good coverage
abs.sub= abs[abs$NewLocation=="Beartooth_Pass",]
abs1.sub= abs1[abs1$NewLocation=="Beartooth_Pass",]

#plot aggregated
plot(abs1.sub$Year, abs1.sub$Grey)
abline(lm(abs1.sub$Grey~abs1.sub$Year))

plot(abs1.sub$JulyMax, abs1.sub$Grey)
abline(lm(abs1.sub$Grey~abs1.sub$JulyMax))

plot(abs1.sub$JuneAv, abs1.sub$Grey)
abline(lm(abs1.sub$Grey~abs1.sub$JuneAv))

#----------------------
#AGGREGATE BY COUNTY

abs1= aggregate(abs, list(abs$Year, abs$County), FUN="mean", na.rm=TRUE)
colnames(abs1)[1:2]= c("Year", "County")
abs1.count= aggregate(abs, list(abs$County, abs$Year), FUN="count")

abs1.count= abs1.count[order(abs1.count$Group.1),]

abs.sub= abs[abs$County=="ClearCreek",]
abs1.sub= abs1[abs1$County=="ClearCreek",]
abs.sub= abs[abs$County=="Park",]
abs1.sub= abs1[abs1$County=="Park",]
abs.sub= abs[abs$County=="Boulder",]
abs1.sub= abs1[abs1$County=="Boulder",]

plot(abs1.sub$Year, abs1.sub$Grey)
abline(lm(abs1.sub$Grey~abs1.sub$Year))

plot(abs1.sub$JulyMeanT, abs1.sub$Grey)
abline(lm(abs1.sub$Grey~abs1.sub$JulyMeanT))

#--------------------------------



 ###Linear Models + Mantel Test
 # grey~ year
 mod1R <- lme(grey~Year+estElevation, random = ~1|NewLocation, method = "ML", data = absM)
 anova(mod1R)
 summary(mod1R)
 anova(mod1R,test="Chi")

 dist1m<-dist(cbind(abs1$Lat, abs1$Long))
 dist0m<- dist(abs1$grey)
 dist1Rm<- dist(residuals(mod1R))

 mantel(dist1m, dist0m)
 mantel(dist1m, dist1Rm)

# grey ~ July 
 mod1JulyMa <- lme(grey~JulyMax+estElevation, random = ~1|NewLocation, method = "ML", data = absM)
 summary(mod1JulyMa)
 anova(mod1JulyMa,test="Chi")
  dist1m<-dist(cbind(abs1$Lat, abs1$Long))
 dist0m<- dist(abs1$grey)
 dist1JulyMa<- dist(residuals(mod1JulyMa))
 
 mantel(dist1m, dist0m)
 mantel(dist1m, dist1JulyMa)
 
 #grey ~ June
 mod1JuneR <- lme(grey~JuneAv+estElevation, random = ~1|NewLocation, method = "ML", data = absM)
 summary(mod1JuneR)
 anova(mod1JuneR,test="Chi")
  dist1m<-dist(cbind(abs1$Lat, abs1$Long))
 dist0m<- dist(abs1$grey)
 dist1JuneRm<- dist(residuals(mod1JuneR))
 
 mantel(dist1m, dist0m)
 mantel(dist1m, dist1JuneRm)

# Thoracic Fur
fabsM<-subset(absM,absM$Collection=="UF"|absM$Collection=="Yale"|absM$Collection=="MW")

 # TF~ year
 mod1R <- lme(Thorax~Year+estElevation, random = ~1|NewLocation, method = "ML", data = fabsM)
 anova(mod1R)
 summary(mod1R)
 anova(mod1R,test="Chi")

 dist1m<-dist(cbind(fabsM$Lat, abs1$Long))
 dist0m<- dist(fabsM$Thorax)
 dist1Rm<- dist(residuals(mod1R))

 mantel(dist1m, dist0m)
 mantel(dist1m, dist1Rm)

# TF ~ July 
 mod1JulyMa <- lme(Thorax~JulyMax+estElevation, random = ~1|NewLocation, method = "ML", data = fabsM)
 summary(mod1JulyMa)
 anova(mod1JulyMa,test="Chi")
 
 dist1m<-dist(cbind(fabsM$Lat, abs1$Long))
 dist0m<- dist(fabsM$Thorax)
 dist1JulyMa<- dist(residuals(mod1JulyMa))
 
 mantel(dist1m, dist0m)
 mantel(dist1m, dist1JulyMa)
 
 #TF ~ June
 mod1JuneR <- lme(Thorax~JuneAv+estElevation, random = ~1|NewLocation, method = "ML", data = fabsM)
 summary(mod1JuneR)
 anova(mod1JuneR,test="Chi")
 
 dist1m<-dist(cbind(fabsM$Lat, abs1$Long))
 dist0m<- dist(fabsM$Thorax)
 dist1JuneRm<- dist(residuals(mod1JuneR))
 
 mantel(dist1m, dist0m)
 mantel(dist1m, dist1JuneRm)


###############
# SPATIAL ANALYSIS
### References: http://geodacenter.asu.edu/drupal_files/spdepintro.pdf
# http://onlinelibrary.wiley.com/doi/10.1111/j.2007.0906-7590.05171.x/full
# http://onlinelibrary.wiley.com/doi/10.1111/j.1466-8238.2007.00334.x/full
# Code follows Buckley et al. 2008 Ecology
# adapted for C.meadii data


mod1<- lm(absM$grey~absM$Year )

absM$Lat=as.numeric(absM$Lat)
absM$Long=as.numeric(absM$Long)
coords= as.matrix(absM[,c("Long","Lat")])
## Build neighborhood sizes 
m.nb0 <- dnearneigh(coords,0,0,longlat=TRUE)
m.nb20 <- dnearneigh(coords,0,20,longlat=TRUE)
m.nb40 <- dnearneigh(coords,0,40,longlat=TRUE)
m.nb80 <- dnearneigh(coords,0,80,longlat=TRUE)

summary(m.nb20)
summary(m.nb40)
summary(m.nb80)

# Adjusted neighborhood size for biologically relevant regions in our analysis  
par(mfrow=c(2,2)  )
plot(m.nb0, coords,points=TRUE)         # just to see where the points are... 
plot(m.nb20, coords,points=TRUE)
plot(m.nb40, coords,points=TRUE)
plot(m.nb80, coords,points=TRUE)

# establish weights for the given list of neighbors, "B" = binary (all equal), "W" = row standardized

m.sw20 <- nb2listw(m.nb20, glist=NULL, style="W", zero.policy=TRUE)
m.sw40 <- nb2listw(m.nb40, glist=NULL, style="W", zero.policy=TRUE)
m.sw80 <- nb2listw(m.nb80, glist=NULL, style="W", zero.policy=TRUE)


## Use error models because error is inherent to the measurements and not to the trait... 

#grey ~ year 

m_x.errorSAR40.7 <- errorsarlm(absM$grey~absM$JulyMax+absM$estElevation, data=absM, m.sw40, method="eigen", quiet=FALSE, zero.policy=TRUE, interval=c(-0.5,0.5),tol.solve=1e-25)
m_x.errorSAR40.6<- errorsarlm(absM$grey~absM$JuneAv+absM$estElevation, data=absM, m.sw40, method="eigen", quiet=FALSE, zero.policy=TRUE, interval=c(-0.5,0.5),tol.solve=1e-25)
m_x.errorSAR40.Y <- errorsarlm(absM$grey~absM$Year+absM$estElevation, data=absM, m.sw40, method="eigen", quiet=FALSE, zero.policy=TRUE, interval=c(-0.5,0.5),tol.solve=1e-25)

summary(m_x.errorSAR40.Y)
summary(m_x.errorSAR40.6)
summary(m_x.errorSAR40.7)

#double check Moran Values
moran.test(residuals(m_x.errorSAR40.7), m.sw40, zero.policy = TRUE)
moran.test(residuals(m_x.errorSAR40.6), m.sw40, zero.policy = TRUE)
moran.test(residuals(m_x.errorSAR40.Y), m.sw40, zero.policy = TRUE)

#
par(mfrow=c(2,2))
cor.errorSAR40resid <- correlog(absM$Long,absM$Lat,z=residuals(m_x.errorSAR40.Y),increment=100,resamp=1,latlon=TRUE)
plot(cor.errorSAR40resid$mean.of.class[1:150], cor.errorSAR40resid$correlation[1:150],xlab="Distance class",ylab="Moran's I",type="b",main="spatial model- year",ylim=c(-0.4,0.6))
cor.errorSAR40resid.7 <- correlog(absM$Long,absM$Lat,z=residuals(m_x.errorSAR40.7),increment=100,resamp=1,latlon=TRUE)
plot(cor.errorSAR40resid.7$mean.of.class[1:150], cor.errorSAR40resid.7$correlation[1:150],xlab="Distance class",ylab="Moran's I",type="b",main="spatial model-july",ylim=c(-0.4,0.6))
cor.errorSAR40resid.6 <- correlog(absM$Long,absM$Lat,z=residuals(m_x.errorSAR40.6),increment=100,resamp=1,latlon=TRUE)
plot(cor.errorSAR40resid.6$mean.of.class[1:150], cor.errorSAR40resid.6$correlation[1:150],xlab="Distance class",ylab="Moran's I",type="b",main="spatial model-june",ylim=c(-0.4,0.6))


### Thorax Fur
mod1<- lm(fabsM$Thorax~fabsM$Year )

fabsM$Lat=as.numeric(fabsM$Lat)
fabsM$Long=as.numeric(fabsM$Long)
coords= as.matrix(fabsM[,c("Long","Lat")])
## Build neighborhood sizes 
m.nb0 <- dnearneigh(coords,0,0,longlat=TRUE)
m.nb20 <- dnearneigh(coords,0,20,longlat=TRUE)
m.nb40 <- dnearneigh(coords,0,40,longlat=TRUE)
m.nb80 <- dnearneigh(coords,0,80,longlat=TRUE)

summary(m.nb20)
summary(m.nb40)
summary(m.nb80)

# Adjusted neighborhood size for biologically relevant regions in our analysis  
par(mfrow=c(2,2)  )
plot(m.nb0, coords,points=TRUE)         # just to see where the points are... 
plot(m.nb20, coords,points=TRUE)
plot(m.nb40, coords,points=TRUE)
plot(m.nb80, coords,points=TRUE)

# establish weights for the given list of neighbors, "B" = binary (all equal), "W" = row standardized

m.sw20 <- nb2listw(m.nb20, glist=NULL, style="W", zero.policy=TRUE)
m.sw40 <- nb2listw(m.nb40, glist=NULL, style="W", zero.policy=TRUE)
m.sw80 <- nb2listw(m.nb80, glist=NULL, style="W", zero.policy=TRUE)


## Use error models because error is inherent to the measurements and not to the trait... 

#Thorax ~ year 

m_x.errorSAR40.7f <- errorsarlm(fabsM$Thorax~fabsM$JulyMax+fabsM$estElevation, data=fabsM, m.sw40, method="eigen", quiet=FALSE, zero.policy=TRUE, interval=c(-0.5,0.5),tol.solve=1e-25)
m_x.errorSAR40.6f<- errorsarlm(fabsM$Thorax~fabsM$JuneAv+fabsM$estElevation, data=fabsM, m.sw40, method="eigen", quiet=FALSE, zero.policy=TRUE, interval=c(-0.5,0.5),tol.solve=1e-25)
m_x.errorSAR40.Yf <- errorsarlm(fabsM$Thorax~fabsM$Year+fabsM$estElevation, data=fabsM, m.sw40, method="eigen", quiet=FALSE, zero.policy=TRUE, interval=c(-0.5,0.5),tol.solve=1e-25)

summary(m_x.errorSAR40.Yf)
summary(m_x.errorSAR40.6f)
summary(m_x.errorSAR40.7f)

#double check Moran Values
moran.test(residuals(m_x.errorSAR40.7f), m.sw40, zero.policy = TRUE)
moran.test(residuals(m_x.errorSAR40.6f), m.sw40, zero.policy = TRUE)
moran.test(residuals(m_x.errorSAR40.Yf), m.sw40, zero.policy = TRUE)

#
par(mfrow=c(2,2))
cor.errorSAR40resid <- correlog(fabsM$Long,fabsM$Lat,z=residuals(m_x.errorSAR40.Yf),increment=100,resamp=1,latlon=TRUE)
plot(cor.errorSAR40resid$mean.of.class[1:150], cor.errorSAR40resid$correlation[1:150],xlab="Distance class",ylab="Moran's I",type="b",main="spatial model- year",ylim=c(-0.4,0.6))
cor.errorSAR40resid.7 <- correlog(fabsM$Long,fabsM$Lat,z=residuals(m_x.errorSAR40.7f),increment=100,resamp=1,latlon=TRUE)
plot(cor.errorSAR40resid.7$mean.of.class[1:150], cor.errorSAR40resid.7$correlation[1:150],xlab="Distance class",ylab="Moran's I",type="b",main="spatial model-july",ylim=c(-0.4,0.6))
cor.errorSAR40resid.6 <- correlog(fabsM$Long,fabsM$Lat,z=residuals(m_x.errorSAR40.6f),increment=100,resamp=1,latlon=TRUE)
plot(cor.errorSAR40resid.6$mean.of.class[1:150], cor.errorSAR40resid.6$correlation[1:150],xlab="Distance class",ylab="Moran's I",type="b",main="spatial model-june",ylim=c(-0.4,0.6))



## Lag Models 
absM$Lat=as.numeric(absM$Lat)
absM$Long=as.numeric(absM$Long)
coords= as.matrix(absM[,c("Long","Lat")])
## Build neighborhood sizes 
m.nb40 <- dnearneigh(coords,0,40,longlat=TRUE)
m.sw40 <- nb2listw(m.nb40, glist=NULL, style="W", zero.policy=TRUE)
m_x.lagSAR40.7 <- lagsarlm(absM$grey~absM$JulyMax+absM$estElevation, data=absM, m.sw40, method="eigen", quiet=FALSE, zero.policy=TRUE, interval=c(-0.5,0.5))
summary(m_x.lagSAR40.7)
m_x.lagSAR40.Y <- lagsarlm(absM$grey~absM$Year+absM$estElevation, data=absM, m.sw40, method="eigen", quiet=FALSE, tol.solve=1e-12, zero.policy=TRUE, interval=c(-0.5,0.5))
summary(m_x.lagSAR40.Y)
m_x.lagSAR40.6 <- lagsarlm(absM$grey~absM$JuneAv+absM$estElevation, data=absM, m.sw40, method="eigen", quiet=FALSE, tol.solve=1e-12, zero.policy=TRUE, interval=c(-0.5,0.5))
summary(m_x.lagSAR40.6)


## Durbin Models
#Perform spatial simultaneous autoregressive lag + error (mixed) model (i.e. Durbin model)
m_x.mixedSAR40.7 <- lagsarlm(absM$grey~absM$JulyMax+absM$estElevation, data=absM, m.sw40, type="mixed",method="eigen", quiet=FALSE, zero.policy=TRUE, interval=c(-0.5,0.5),tol.solve=1e-25)
m_x.mixedSAR40.Y <- lagsarlm(absM$grey~absM$Year+absM$estElevation, data=absM, m.sw40, type="mixed",method="eigen", quiet=FALSE, zero.policy=TRUE, interval=c(-0.5,0.5),tol.solve=1e-25)
m_x.mixedSAR40.6 <- lagsarlm(absM$grey~absM$JuneAv+absM$estElevation, data=absM, m.sw40, type="mixed",method="eigen", quiet=FALSE, zero.policy=TRUE, interval=c(-0.5,0.5),tol.solve=1e-25)

summary(m_x.mixedSAR40.7)                
summary(m_x.mixedSAR40.Y)                
summary(m_x.mixedSAR40.6)

AIC(m_x.errorSAR40.7,m_x.errorSAR40.Y,m_x.errorSAR40.6,m_x.lagSAR40.7,m_x.lagSAR40.Y, m_x.lagSAR40.6, m_x.mixedSAR40.7,m_x.mixedSAR40.Y, m_x.mixedSAR40.6) 


## Thorax Lag Models 
fabsM$Lat=as.numeric(fabsM$Lat)
fabsM$Long=as.numeric(fabsM$Long)
coords= as.matrix(fabsM[,c("Long","Lat")])
## Build neighborhood sizes 
m.nb40f <- dnearneigh(coords,0,40,longlat=TRUE)
m.sw40f <- nb2listw(m.nb40f, glist=NULL, style="W", zero.policy=TRUE)
m_x.lagSAR40.7f <- lagsarlm(fabsM$Thorax~fabsM$JulyMax+fabsM$estElevation, data=fabsM, m.sw40f, method="eigen", quiet=FALSE, zero.policy=TRUE, interval=c(-0.5,0.5))
summary(m_x.lagSAR40.7f)
m_x.lagSAR40.Yf <- lagsarlm(fabsM$Thorax~fabsM$Year+fabsM$estElevation, data=fabsM, m.sw40f, method="eigen", quiet=FALSE, tol.solve=1e-12, zero.policy=TRUE, interval=c(-0.5,0.5))
summary(m_x.lagSAR40.Yf)
m_x.lagSAR40.6f <- lagsarlm(fabsM$Thorax~fabsM$JuneAv+fabsM$estElevation, data=fabsM, m.sw40f, method="eigen", quiet=FALSE, tol.solve=1e-12, zero.policy=TRUE, interval=c(-0.5,0.5))
summary(m_x.lagSAR40.6f)


## Durbin Models
#Perform spatial simultaneous autoregressive lag + error (mixed) model (i.e. Durbin model)
m_x.mixedSAR40.7f <- lagsarlm(fabsM$Thorax~fabsM$JulyMax+fabsM$estElevation, data=fabsM, m.sw40f, type="mixed",method="eigen", quiet=FALSE, zero.policy=TRUE, interval=c(-0.5,0.5),tol.solve=1e-25)
m_x.mixedSAR40.Yf <- lagsarlm(fabsM$Thorax~fabsM$Year+fabsM$estElevation, data=fabsM, m.sw40f, type="mixed",method="eigen", quiet=FALSE, zero.policy=TRUE, interval=c(-0.5,0.5),tol.solve=1e-25)
m_x.mixedSAR40.6f <- lagsarlm(fabsM$Thorax~fabsM$JuneAv+fabsM$estElevation, data=fabsM, m.sw40f, type="mixed",method="eigen", quiet=FALSE, zero.policy=TRUE, interval=c(-0.5,0.5),tol.solve=1e-25)

summary(m_x.mixedSAR40.7f)                
summary(m_x.mixedSAR40.Yf)                
summary(m_x.mixedSAR40.6f)

AIC(m_x.errorSAR40.7f,m_x.errorSAR40.Yf,m_x.errorSAR40.6f,m_x.lagSAR40.7f,m_x.lagSAR40.Yf, m_x.lagSAR40.6f, m_x.mixedSAR40.7f,m_x.mixedSAR40.Yf, m_x.mixedSAR40.6f) 


###############################
###############################
# Single location with good coverage
# Single year/location with good coverage
#all data... 


