### ANALYSIS OF COLIAS MUSEUM SPECIMENS

mydir= "C:\\Users\\lbuckley\\Documents\\MuseumTraits\\"
#mydir= "C:\\Users\\Buckley\\Google Drive\\ColiasEvolution\\HistoricalAnalysis\\"

#install packages
#install.packages(c("ecodist","spdep","faraway","MuMIn","lmtest"))

#load libraries
library(sp)
library(Rmisc)
library(nlme)
library(ecodist) #for Mantel ##HJM we used spdep for Mantel previously. 
library(spdep)
library(faraway) #for partial residual plot         
library(MuMIn) #for model averaging
require(lmtest) #for likelihood ratio test
require(dplyr)

library(ggplot2)
library(maptools)
library(tidyr)
library(plyr)
#for mapping
library(ggmap)
library(maps)
library(mapdata)

#daymetR, http://khufkens.github.io/daymetr/
#if(!require(devtools)){install.package(devtools)}
#devtools::install_github("khufkens/daymetr")
library(daymetr)

#test
download_daymet(site = "Oak Ridge National Laboratories",
                lat = 36.0133,
                lon = -84.2625,
                start = 1980,
                end = 2010,
                internal = TRUE)

#---------------------------
#READ AND MANIPULATE DATA

#read data
setwd(paste(mydir, "data\\", sep=""))
abs=read.csv("FullDataColiasMuseums.csv")

abs$lon= abs$Long
abs$lat= abs$Lat

#remove samples without locations
abs= subset(abs, !is.na(abs$Long)&!is.na(abs$Lat) )

#group nearby localities
ClassifySites= function(names, lon, lat, DistCut){
  
  new.names=names
  stat.coords= cbind(lon,lat)
  
  for(r in 1:length(lon) ){
    dists= spDistsN1(stat.coords, stat.coords[r,], longlat = TRUE) #find distances from focal site in km
    new.names[which(dists<DistCut)]= new.names[r] #rename sites within cutoff radius
  } #end loop sites
  
  return(new.names) #returns grouped site names
}

abs$NewLocation<- ClassifySites(abs$Location, abs$Long, abs$Lat, DistCut=5)    #radius of 5km

#calculate Julian
mdy=  paste(abs$Month,abs$Day,abs$Year, sep="-")
tmp <- as.POSIXlt(mdy, format = "%m-%d-%Y")
abs$doy=tmp$yday

#abbreviate columns to control for NAs
#absM= absM[,c("ID","NewLocation","Collection", "grey","Thorax","JulyAv","JulyMax","JuneAv","Year","estElevation","J", "Long","Lat")]
#abs= absM[,colnames(abs)!="BodyLength"]
# take out body length due to NAs

#make elevation numeric
abs$estElevation = as.numeric(as.character(abs$estElevation))

#specify regions
abs$region= 1
abs$region[which(abs$lat>42)]<-2
abs$region[which(abs$lat>48)]<-3

#fix spelling of Summit, check!
abs[which(abs$County=="Sumit") ,"County"]<-"Summit"

#GET RID OF TWO SPECIMENS WITH J=171
abs= abs[abs$doy>171,]

#Sex
absM<- subset(abs,Sex =="M")

#-----------------------
#DATA COUNTS
count=function(x) length(na.omit(x))

abs1.count= aggregate(abs, list(abs$State, abs$Year), FUN="count")
abs1.count= aggregate(absM, list(absM$NewLocation, absM$Year), FUN="count")

abs1.count= abs1.count[order(abs1.count$Group.1),1:3]

#-------------------------
#Subsample
# 
# #Assign site ID
# sites= unique(absM$NewLocation)
# absM$siteID= match(absM$NewLocation, sites)
# absM$YrSite= paste(absM$Year, absM$siteID, sep="")
# 
# Nruns= 50 #50 #number of bootstrapp runs                          
# Nsamp= 15 #max sample size of butterflies per site per year  
# 
# z <- sapply(unique(absM$YrSite), FUN= function(x){ 
#   sample(which(absM$YrSite==x), min(Nsamp, length(which(absM$YrSite==x))), FALSE)
# })
# absM.boot<- absM[unlist(z),]
# 
# absM<- absM.boot
#-------------------------

#CLIMATE DATA

#DaymetR, https://khufkens.github.io/daymetr/
#R prism
#rnoaa: R coop data  ghcnd

#Extract unique coordinates and years
locs= aggregate(absM, list(absM$NewLocation, absM$Year), FUN="count")

#CLIMAX
#ESTIMATE CLIMATE FOR ADULT, PUPAL, COMBINED 
setwd(paste(mydir, "data\\", sep=""))

clim= read.csv("ClimaxCOOP_5.4.15.csv", na.strings = "-9999")  

clim[,4:6]= clim[,4:6]/10   #convert to C
clim$Year= substr(clim$DATE,1,4)

#calculate Julian
tmp <- as.POSIXlt(as.character(clim$DATE), format = "%Y%m%d")
clim$J=tmp$yday
clim$JY= paste(clim$J, clim$Year, sep="")

#cut NAs 
#absM= na.omit(absM)

#match to temp
absM$AdultTmax=NA
absM$AdultTave=NA
absM$ParentTmax=NA
absM$ParentTave=NA
absM$PupalTmax=NA
absM$PupalTave=NA
absM$LifeTmax=NA
absM$LifeTave=NA

absM$J=absM$doy

for(i in 1:nrow(absM)){ #Lazily coding as a loop
  adult= paste( (absM$J[i]-5):(absM$J[i]+5), absM$Year[i], sep="")      # expanded to 11 day window
  # par.range= (absM$J[i]-10):(absM$J[i]+10)                             #@ Tried models with all permutations +-20, +- 10 , and using below percentile. 
  ## TRY FIXING PARENTAL RANGE
  par.range= 201:216  ## 25th and 75th percentile 
  ## OPTIONAL CODE TO CONSTRAIN PARENTAL RANGE  
  # x= length(which(par.range<171))
  #  if(x>0) par.range= (min(par.range)+x):(max(par.range)+x)
  #  x= length(which(par.range>225))
  # if(x>0) par.range= (min(par.range)-x):(max(par.range)-x)
  
  parent= paste( par.range, absM$Year[i]-1, sep="")  #YEAR BEFORE 
  pupal= paste( (absM$J[i]-26):(absM$J[i]-6), absM$Year[i], sep="")      
  
  adult.TMAX= clim[match(adult, clim$JY),"TMAX"]
  adult.TMIN= clim[match(adult, clim$JY),"TMIN"]
  parent.TMAX= clim[match(parent, clim$JY),"TMAX"]
  parent.TMIN= clim[match(parent, clim$JY),"TMIN"]
  pupal.TMAX= clim[match(pupal, clim$JY),"TMAX"]
  pupal.TMIN= clim[match(pupal, clim$JY),"TMIN"]
  
  absM$AdultTmax[i]=mean(adult.TMAX, na.rm=TRUE)
  absM$AdultTave[i]=mean(c(adult.TMAX, adult.TMIN), na.rm=TRUE)
  absM$ParentTmax[i]=mean(parent.TMAX, na.rm=TRUE)
  absM$ParentTave[i]=mean(c(parent.TMAX, adult.TMIN), na.rm=TRUE)
  absM$PupalTmax[i]=mean(pupal.TMAX, na.rm=TRUE)
  absM$PupalTave[i]=mean(c(pupal.TMAX, pupal.TMIN), na.rm=TRUE)
  absM$LifeTmax[i]=mean(c(pupal.TMAX,adult.TMAX), na.rm=TRUE)
  absM$LifeTave[i]=mean(c(pupal.TMAX,adult.TMAX,pupal.TMIN,adult.TMIN), na.rm=TRUE)
} #end loop rows

#cut NAs 
#absM= na.omit(absM)

#Add seasonal climate
#J= 152 to 212, June and July
clim.jj= clim[which(clim$J>151 & clim$J<213),]
clim.jj$TMEAN= (clim.jj$TMAX + clim.jj$TMIN)/2
clim.jj= aggregate(clim.jj, list(clim.jj$Year), FUN="mean", na.rm=TRUE)
clim.jj$Year= clim.jj$Group.1

match1= match(absM$Year, clim.jj$Year)
absM$JJTave= clim.jj$TMEAN[match1] 

#save
absM.all= absM

#-------
#For loveland pass
#Cabin creek 051186, http://climate.colostate.edu/data_access.html
#http://climatetrends.colostate.edu

#clim= read.csv("CabinCreek.csv", na.strings = "-9999")  #F and IN

#=======================
# Result 1. Maps and overview plots
#MAKE INITIAL MAPS

#set up map
bbox <- ggmap::make_bbox(lon, lat, absM, f = 0.1)
bbox[1]= bbox[1] -5
bbox[2]= bbox[2] -5
bbox[3]= bbox[3] +5
bbox[4]= bbox[4] +5

map_loc <- get_map(location = bbox, source = 'google', maptype = 'terrain')
map1=ggmap(map_loc, margins=FALSE)

#elevation
aper1.map<- map1 +geom_point(data=absM, aes(color=estElevation) ) + coord_cartesian()

#-------------------
#Elevation inset plot

#Lower elevation at higher lats
ggplot(data=absM, aes(x=lat, y = estElevation, color=doy ))+geom_point(alpha=0.8) +theme_bw()+
  geom_smooth(method="lm") #+ guides(color=FALSE)+labs(x = "",y="")

#==================================
#Result 2.	Temp change over time

#PLOT TEMP TREND
ggplot(data=absM.all, aes(x=Year, y = JJTave, color=Corr.Val ))+geom_point(alpha=0.8) +theme_bw()+
  geom_smooth(method="lm")

#! CHECK ONE VERY HOT YEAR

#==================================
#Result 3.	Shift in phenology as a function of temperature (Jun, July average; test others)
#Result 4.	Shift in phenotype as a function of phenology

#LOVELAND PASS
#"EisenhowerTunnel"/ Loveland pass ANALYSIS
abs.sub1= absM.all[absM.all$NewLocation =="EisenhowerTunnel", ]
mod1= lm(Corr.Val~doy*JJTave, data=abs.sub1)

#normalize
abs2<-scale(abs.sub1[,c("PupalTave","ParentTave","J","Year","JJTave")], center=TRUE, scale=FALSE)
abs.sub1[,c("PupalTave","ParentTave","J","Year","JJTave")]=abs2

absM<- na.omit(abs.sub1[,c("Corr.Val","ThoraxC","NewLocation","PupalTave","ParentTave","J","Year","JJTave")])

#----------
#phenology ##*
ggplot(data=abs.sub, aes(x=JJTave, y = J, color=Year ))+geom_point(alpha=0.8) +theme_bw()+
  geom_smooth(method="lm") 

ggplot(data=abs.sub, aes(x=J, y = Corr.Val, color=Year ))+geom_point(alpha=0.8) +theme_bw()+
  geom_smooth(method="lm") 

#Results in annual pattern
ggplot(data=abs.sub, aes(x=Year, y = Corr.Val, color=Year ))+geom_point(alpha=0.8) +theme_bw()+
  geom_smooth(method="lm") 

#------------
#REGIONS

#Phenology
ggplot(data=absM.all, aes(x=Year, y = doy, color=Year ))+geom_point(alpha=0.8) +theme_bw()+
  facet_wrap(~region)+geom_smooth(method="lm")
ggplot(data=absM.all, aes(x=JJTave, y = doy, color=Year ))+geom_point(alpha=0.8) +theme_bw()+
  facet_wrap(~region)+geom_smooth(method="lm")

#Abs
ggplot(data=absM.all, aes(x=doy, y = Corr.Val, color=Year ))+geom_point(alpha=0.8) +theme_bw()+
  facet_wrap(~region)+geom_smooth(method="lm")
#by year
ggplot(data=absM.all, aes(x=Year, y = Corr.Val, color=Year ))+geom_point(alpha=0.8) +theme_bw()+
  facet_wrap(~region)+geom_smooth(method="lm")

#-------------
#OTHER SUBSETS EXPLORED

#by county
abs1.count= aggregate(absM, list(absM$County), FUN="count")
counties= abs1.count[which(abs1.count$ID>75),"Group.1"]
absM2= absM[which(absM$County %in% counties) ,]

#ALberta 1980, n=144 
abs.sub= abs[abs$State=="Alberta" & abs$Year==1980,]

#Select sites with good coverage
abs.sub= absM[absM$NewLocation %in% c("Plateau_Mnt","EisenhowerTunnel","Clay_Butte_Beartooth_Plateau", "Libby_Flats" ), ]
#omit two 1920's specimens
abs.sub= abs.sub[which(abs.sub$Year>1930), ]
#body length to numeric
abs.sub$BodyLength= as.numeric(as.character(abs.sub$BodyLength))

#=====================
#STATISTICS

#BOOTSTRAP
#Assign site ID
sites= unique(absM$NewLocation)
absM$siteID= match(absM$NewLocation, sites)
absM$YrSite= paste(absM$Year, absM$siteID, sep="")

Nruns= 50 #50 #number of bootstrapp runs                          
Nsamp= 15 #max sample size of butterflies per site per year   

boot.lm<- function(x=absM$JJTave, y=absM$doy, sites= absM$YrSite, Nruns,Nsamp){
  out<- matrix(NA, nrow=Nruns, ncol=5)
  
  for(r in 1:Nruns){
    
    #sub sample
    z <- sapply(unique(sites), FUN= function(x){ 
      sample(which(sites==x), min(Nsamp, length(which(sites==x))), FALSE)
    })
    x.boot<- x[unlist(z)]
    y.boot<- y[unlist(z)]
    
    #run model
    mod1= lm(y.boot~x.boot)
    mod1$coefficients[2]
  
    out[r,]=c(summary(mod1)$coefficients[2,], summary(mod1)$r.squared )
  
  }# end loop
  
    #average
    return( colMeans(out) )
}

#-------
boot.lm(x=absM$JJTave, y=absM$doy, sites= absM$YrSite, Nruns,Nsamp)

#--------------------
#SPATIAL ANALYSIS  

#dat is matrix of y followed by predictor variables, assume interactions for now

#normalize
abs2<-scale(absM[,c("JJTave","doy","Year")], center=TRUE, scale=FALSE)
absM[,c("JJTave","doy","Year")]=abs2
#Remove NAs
dat= na.omit(absM[,c("Corr.Val","JJTave","doy","Year")])

boot.spatlm<- function(dat, sites= dat$YrSite, Nruns,Nsamp){
  out<- matrix(NA, nrow=Nruns, ncol=5)

#set up data collection
out.mods= array(data=NA, dim=c(5,11,Nruns))
out.coefs=array(data=NA, dim=c(7,5,Nruns))
out.zs=array(data=NA, dim=c(6,4,Nruns))
out.stats= matrix(NA, 3, Nruns)

#bootstrap
for(r in 1:Nruns){
  
  #sub sample
  z <- sapply(unique(sites), FUN= function(x){ 
    sample(which(sites==x), min(Nsamp, length(which(sites==x))), FALSE)
  })
  absM.boot<- absM[unlist(z),]
  
  #run spatial model
  out=spat.mod(absM.boot)
  
  #extract output
  bm=out[[1]]
  out.mods[,,r]=as.matrix(bm)
  out.coefs[,,r]=out[[2]]
  out.zs[,,r]=out[[3]]
  out.stats[,r]=out[[4]]
  
  bm=out[[5]]
  thorax.out.mods[,,r]=as.matrix(bm)
  thorax.out.coefs[,,r]=out[[6]]
  thorax.out.zs[,,r]=out[[7]]
  thorax.out.stats[,r]=out[[8]]
  
}#end bootstrap

#PROCESS OUTPUT
#ADD NAMES
dimnames(out.mods)[[2]]= colnames(out[[1]])
dimnames(out.coefs)[[1]]= rownames(out[[2]])
dimnames(out.coefs)[[2]]= colnames(out[[2]])
dimnames(out.zs)[[1]]= rownames(out[[3]])
dimnames(out.zs)[[2]]= colnames(out[[3]])
rownames(out.stats)= names(out[[4]])

dimnames(thorax.out.mods)[[2]]= colnames(out[[5]])
dimnames(thorax.out.coefs)[[1]]= rownames(out[[6]])
dimnames(thorax.out.coefs)[[2]]= colnames(out[[6]])
dimnames(thorax.out.zs)[[1]]= rownames(out[[7]])
dimnames(thorax.out.zs)[[2]]= colnames(out[[7]])
rownames(thorax.out.stats)= names(out[[8]])

#AVERAGE
coefs= apply(out.coefs, MARGIN=c(1,2), FUN="mean") 
zs= apply(out.zs, MARGIN=c(1,2), FUN="mean") 
stats=  apply(out.stats, MARGIN=c(1), FUN="mean") 

thorax.coefs= apply(thorax.out.coefs, MARGIN=c(1,2), FUN="mean") 
thorax.zs= apply(thorax.out.zs, MARGIN=c(1,2), FUN="mean") 
thorax.stats=  apply(thorax.out.stats, MARGIN=c(1), FUN="mean") 

#WRITE OUT
coefs # MODEL AVERAGING OUTPUT
zs #FULL MODEL OUTPUT
stats

thorax.coefs
thorax.zs
thorax.stats




#=================================

#PArtial residual plot
mod1= lm(absM$Corr.Val ~ absM$JJTave*absM$J) 
anova(mod1, type=2)
#Strong interaction influence phenology
prplot(mod1,1)
prplot(mod1,2)
prplot(mod1,3)
prplot(mod1,4)


#---------------------------------------------------

#=========================
#MODELS WHOLE DATA SET

#subset by region
absM.reg= subset(absM.all, absM.all$region==3)

absM1= na.omit(absM.reg)
absM1$estElevation = as.numeric(absM1$estElevation)

mod1R <- lme(Corr.Val~Year*doy, random = ~1|Location, method = "ML", data = absM1)
anova(mod1R)
summary(mod1R)
anova(mod1R,test="Chi")

#------------------------------

#===================================================
#PARTIAL REGRESSION PLOTS
#see http://webmail.dev411.com/p/r/r-sig-geo/1195dstqs8/partial-residual-plot-from-sar-model
#Partial residuals (residuals of regressing the response variable on the independent variables, but omitting the independent variable of interest) 
#Can depict low variance when predictors are correlated as is likely the case here, should look into.

z <- sapply(unique(absM$YrSite), FUN= function(x){ 
  sample(which(absM$YrSite==x), min(Nsamp, length(which(absM$YrSite==x))), FALSE)
})
#absM.boot<- absM[unlist(z),]
absM.boot= absM

#define spatial neighborhood
coords= as.matrix(absM.boot[,c("Long","Lat")])
## Build neighborhood
m.nb40 <- dnearneigh(coords,0,40,longlat=TRUE)
# establish weights for the given list of neighbors, "B" = binary (all equal), "W" = row standardized
m.sw40 <- nb2listw(m.nb40, glist=NULL, style="W", zero.policy=TRUE)

## PLOT ALL DATA
#GREY
m_x.errorSAR40 <- errorsarlm(grey~PupalTave+ParentTave+J+Year+PupalTave*J, data=absM.boot, m.sw40, method="eigen", quiet=TRUE, zero.policy=TRUE, na.action=na.fail, interval=c(-0.5,0.5),tol.solve=1e-25)
#m_x.errorSAR40 <- errorsarlm(grey~PupalTave+ParentTave+J+Year+PupalTave*J+ParentTave*J, data=absM.boot, m.sw40, method="eigen", quiet=TRUE, zero.policy=TRUE, na.action=na.fail, interval=c(-0.5,0.5),tol.solve=1e-25)
#! Take out interaction between flight season temperature and J

df0 <- as.data.frame(m_x.errorSAR40$tarX)
names(df0) <- names(m_x.errorSAR40$coefficients)
df <- cbind(y=m_x.errorSAR40$tary, df0[,-1])
str(df)
tarlm <- lm(absM.boot$grey~absM.boot$PupalTave+absM.boot$ParentTave+absM.boot$J+absM.boot$Year+absM.boot$PupalTave*absM.boot$J, data=df)
#tarlm <- lm(absM.boot$grey~absM.boot$PupalTave+absM.boot$ParentTave+absM.boot$J+absM.boot$Year+absM.boot$PupalTave*absM.boot$J+absM.boot$ParentTave*absM.boot$J, data=df)

#---------------------------------

xlabs= c("Pupal T (?C)","Flight season T (?C)", "Date of collection (J)", "Year", "Pupal T (?C): Date")

#@ xlabs= c("Pupal T (°C)","Parent T (°C)", "Date of collection (J)", "Year", "Pupal T(°C) : Date", "Parent T (°C): Date")  #@ Changed Labs to be capitalized. Doesn't really matter... 

#setwd(paste(mydir, "figs\\", sep=""))
#pdf("GreyPR_norm.pdf",height = 6, width = 10)

par(mfrow=c(2,3), cex=1.1, lwd=1, mar=c(3,2, 1, 1), mgp=c(1.3, 0.5, 0), oma=c(0,2,0,0), bty="l", cex.lab=1.2)
my.prplot(tarlm, 1, xlabs)
my.prplot(tarlm, 2, xlabs)
my.prplot(tarlm, 3, xlabs)
my.prplot(tarlm, 4, xlabs)
my.prplot(tarlm, 5, xlabs)
#my.prplot(tarlm, 6, xlabs)

mtext("Partial residual of wing melanism (grey level)", side=2, line = -0.5, cex=1.3, outer=TRUE)

#dev.off()
#----------------------------
#FUR
m_x.errorSAR40 <- errorsarlm(Thorax~PupalTave+ParentTave+J+Year+PupalTave*J, data=absM.boot, m.sw40, method="eigen", quiet=FALSE, zero.policy=TRUE, na.action=na.fail, interval=c(-0.5,0.5),tol.solve=1e-25)
#m_x.errorSAR40 <- errorsarlm(Thorax~PupalTave+ParentTave+J+Year+PupalTave*J+ParentTave*J, data=absM.boot, m.sw40, method="eigen", quiet=FALSE, zero.policy=TRUE, na.action=na.fail, interval=c(-0.5,0.5),tol.solve=1e-25)

df0 <- as.data.frame(m_x.errorSAR40$tarX)
names(df0) <- names(m_x.errorSAR40$coefficients)
df <- cbind(y=m_x.errorSAR40$tary, df0[,-1])
str(df)
tarlm <- lm(absM.boot$Thorax~absM.boot$PupalTave+absM.boot$ParentTave+absM.boot$J+absM.boot$Year+absM.boot$PupalTave*absM.boot$J, data=df)
#tarlm <- lm(absM.boot$Thorax~absM.boot$PupalTave+absM.boot$ParentTave+absM.boot$J+absM.boot$Year+absM.boot$PupalTave*absM.boot$J+absM.boot$ParentTave*absM.boot$J, data=df)

#pdf("FurPR_norm.pdf",height = 6, width = 10)

par(mfrow=c(2,3), cex=1.1, lwd=1, mar=c(3,2, 1, 1), mgp=c(1.3, 0.5, 0), oma=c(0,2,0,0), bty="l", cex.lab=1.2)
my.prplot(tarlm, 1, xlabs)
my.prplot(tarlm, 2, xlabs)
my.prplot(tarlm, 3, xlabs)
my.prplot(tarlm, 4, xlabs)
my.prplot(tarlm, 5, xlabs)
#my.prplot(tarlm, 6, xlabs)

mtext("Partial residual of setae length (mm)", side=2, line = -0.5, cex=1.3, outer=TRUE)

#dev.off()
#=========================================
# OTHER ANALYSES

plot(absM$Year, absM$J)
mod1= lm(absM$J ~ absM$Year)
abline(mod1)
## Collected later in later years
m_x.errorSAR40 <- errorsarlm(J~Year, data=absM.boot, m.sw40, method="eigen", quiet=TRUE, zero.policy=TRUE, na.action=na.fail, interval=c(-0.5,0.5),tol.solve=1e-25)

mod1= lm(absM$estElevation ~ absM$Year)
summary(mod1)

mod1= lm(absM$J ~ absM$Year * absM$estElevation)
anova(mod1, type=2)

#---------------------------------------
mod1 <- lm(absM$grey~absM$PupalTave+absM$ParentTave+absM$J+absM$Year+absM$PupalTave*absM$J+absM$ParentTave*absM$J, data=absM.boot)

#check dist all reasonably normal
hist(absM$grey)
hist(absM$PupalTave)
hist(absM$ParentTave)
hist(absM$J)
hist(absM$Year)




