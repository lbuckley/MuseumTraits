### ANALYSIS OF COLIAS MUSEUM SPECIMENS

mydir= "C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\MuseumTraits\\"
#mydir= "C:\\Users\\lbuckley\\Documents\\MuseumTraits\\"
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

#---------------------------
#READ AND MANIPULATE DATA

#read data
setwd(paste(mydir, "data\\", sep=""))
abs=read.csv("FullDataColiasMuseums.csv")

#make longitude numeric
abs$Long= as.numeric(as.character(abs$Long))

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
#CLIMATE DATA, use weather stations

#Region 1: Climax
#Region 2: MORAN 5WNW, WY USA, 	USC00486440, http://mco.cfc.umt.edu/ghcn/station/USC00486440.html
#Region 3: BLUEHILL LO, 1,950.70 m http://climate.weather.gc.ca/climate_data/

#Load
setwd(paste(mydir, "data\\", sep=""))
clim= read.csv("WeatherStationData_Regions.csv")

#calculate Julian
clim$DATE= paste(clim$YEAR,"/", clim$MONTH, "/", clim$DAY,sep="")
tmp <- as.POSIXlt(as.character(clim$DATE), format = "%Y/%m/%d")
clim$J=tmp$yday
clim$JYR= paste(clim$J, clim$YEAR,clim$region, sep="")
clim$MYR= paste(clim$MONTH, clim$YEAR,as.character(clim$region), sep="")

#---------------- 
#match to temp
absM$June=NA
absM$July=NA
absM$June15July15=NA
absM$doy162to202=NA

absM$J=absM$doy

for(i in 1:nrow(absM)){ #Lazily coding as a loop
  
  #June
  dev= paste( 6, absM$Year[i],absM$region[i], sep="")     
  cdat= clim[match(dev, clim$MYR),"TMEAN"]
  absM$June[i]=mean(cdat, na.rm=TRUE)
  
  #July
  dev= paste( 7, absM$Year[i],absM$region[i], sep="")     
  cdat= clim[match(dev, clim$MYR),"TMEAN"]
  absM$July[i]=mean(cdat, na.rm=TRUE)
  
  #June15 to July 15
  dev= paste( 152:181, absM$Year[i],absM$region[i], sep="")     
  cdat= clim[match(dev, clim$JYR),"TMEAN"]
  absM$June15July15[i]=mean(cdat, na.rm=TRUE)
  
  #doy160to202
  dev= paste( 162:202, absM$Year[i],absM$region[i], sep="")     
  cdat= clim[match(dev, clim$JYR),"TMEAN"]
  absM$doy162to202[i]=mean(cdat, na.rm=TRUE)
  
} #end loop rows

#JuneJuly
absM$JJTave= (absM$June + absM$July) / 2

#save
absM.all= absM

#set up bootstrap regressions
#BOOTSTRAP
#Assign site ID
sites= unique(absM.all$NewLocation)
absM.all$siteID= match(absM.all$NewLocation, sites)
absM$YrSite= paste(absM.all$Year, absM.all$siteID, sep="")

Nruns= 50 #50 #number of bootstrapp runs                          
Nsamp= 15 #max sample size of butterflies per site per year

#=======================
# Result 1. Maps and overview plots
#MAKE INITIAL MAPS

absM.all$Long= as.numeric(as.character(absM.all$Long))

#set up map
bbox <- ggmap::make_bbox(Long, Lat, absM.all, f = 0.1)
bbox[1]= bbox[1] -5
bbox[2]= bbox[2] -5
bbox[3]= bbox[3] +5
bbox[4]= bbox[4] +5

map_loc <- get_map(location = bbox, source = 'google', maptype = 'terrain')
map1=ggmap(map_loc, margins=FALSE) #

#elevation
absM$Long= as.numeric(as.character(absM$Long))

aper1.map<- map1 +geom_point(data=absM, aes(y=Lat, x=Long, color=estElevation) ) + coord_cartesian() + labs(x = "Longitude (°)",y="Latitude (°)", color="Elevation (m)") + theme(legend.position="bottom")+scale_color_gradientn(colours = topo.colors(5))

#-------------------
#Elevation inset plot

#Lower elevation at higher lats
elev.plot<-ggplot(data=absM, aes(x=lat, y = estElevation))+geom_point(alpha=0.8) +theme_classic()+
  labs(x = "Latitude (°)",y="Elevation (m)")
  # , color=doy #geom_smooth(method="lm")+

#-------------------
#Inset:	Temp change over time

#make region factor labels
clim$region.f= "1: Southern"
clim$region.f[which(clim$region==2)]<- "2: Northern"
clim$region.f[which(clim$region==3)]<- "3: Canadian"
clim$region= as.numeric(clim$region)

#calculate developmental temperatures each year
clim.dev= subset(clim, clim$J %in% 162:202)
#aggregate by year and region
clim.dev= aggregate(clim.dev, list(clim.dev$region.f, clim.dev$YEAR), FUN="mean", na.rm=TRUE)
clim.dev$region.f=clim.dev$Group.1

#time series for three regions
clim.plot= ggplot(data=clim.dev, aes(x=YEAR, y = TMEAN, color=region.f))+geom_line() +theme_classic()+geom_smooth(method="lm",se=FALSE)+xlim(1950,2013)+ylab("Developmental Temperature (°C)") +xlab("Year")+labs(color="Region") #+ theme(legend.position = c(0.2, 0.8))

#find years with samples
clim$region= as.factor(clim$region)
a.sub= subset(absM.all, absM.all$region==1)
yrs.reg= unique(a.sub$Year)
clim.dev1= subset(clim.dev, clim.dev$region==1 & clim.dev$YEAR %in% yrs.reg)

a.sub= subset(absM.all, absM.all$region==2)
yrs.reg= unique(a.sub$Year)
clim.dev2= subset(clim.dev, clim.dev$region==2 & clim.dev$YEAR %in% yrs.reg)

a.sub= subset(absM.all, absM.all$region==3)
yrs.reg= unique(a.sub$Year)
clim.dev3= subset(clim.dev, clim.dev$region==3 & clim.dev$YEAR %in% yrs.reg)

clim.devs= rbind(clim.dev1, clim.dev2, clim.dev3)

#add dots for years with samples
clim.plot= clim.plot+ geom_point(data=clim.devs, aes(x=YEAR, y = TMEAN, color=region.f))

#-------------------

#COMBINE
library(grid)

setwd(paste(mydir, "figures\\", sep=""))
pdf("ElevFig.pdf", height=8, width=8)

##open pdf
subvp.t<-viewport(width=.5,height=.38,x=.7,y=0.8)
subvp.e<-viewport(width=.4,height=.40,x=.29,y=0.34)
##Next, open the main graph which was stored in b by typing b at the prompt:
aper1.map
##Then, superimpose the graph stored in a on the viewport as:
print(elev.plot,vp=subvp.e)
print(clim.plot,vp=subvp.t)
dev.off()

#==================================
#Result 3.	Shift in phenology as a function of temperature (Jun, July average; test others)
#Result 4.	Shift in phenotype as a function of phenology

#LOVELAND PASS
#"EisenhowerTunnel"/ Loveland pass ANALYSIS
abs.sub1= absM.all[absM.all$NewLocation =="EisenhowerTunnel", ]

#-------------------------
#Statistics

abs.sub1$siteID= match(abs.sub1$NewLocation, sites)
abs.sub1$YrSite= paste(abs.sub1$Year, abs.sub1$siteID, sep="")

#Bootstrap
phen.mod= boot.lm(x=abs.sub1$doy162to202, y=abs.sub1$J, sites= abs.sub1$YrSite, Nruns,Nsamp)
plast.mod= boot.lm(x=abs.sub1$J, y = abs.sub1$Corr.Val, sites= abs.sub1$YrSite, Nruns,Nsamp)
year.mod= boot.lm(x=abs.sub1$Year, y = abs.sub1$Corr.Val, sites= abs.sub1$YrSite, Nruns,Nsamp)
  
#------------------------
#phenology 
fig2a=ggplot(data=abs.sub1, aes(x=doy162to202, y = J, color=Year ))+geom_point(alpha=0.8) +theme_classic()+
  xlab("Developmental Temperature (°C)") +ylab("Phenology (doy)") + theme(legend.position="bottom") +scale_color_gradientn(colours = topo.colors(5))+ theme(legend.position="bottom") 
#add trend
if(phen.mod["P"]<0.05) fig2a= fig2a + geom_abline( aes(slope=phen.mod["Estimate"],intercept=phen.mod["Intercept"]))

#plasticity
fig2b=ggplot(data=abs.sub1, aes(x=J, y = Corr.Val, color=Year ))+geom_point(alpha=0.8) +theme_classic() + xlab("Phenology (doy)") +ylab("Wing melanism (grey level)")+scale_color_gradientn(colours = topo.colors(5))+ theme(legend.position="none")
#add trend
if(plast.mod["P"]<0.05) fig2b= fig2b + geom_abline( aes(slope=plast.mod["Estimate"],intercept=plast.mod["Intercept"]))

#Results in annual pattern
fig2c=ggplot(data=abs.sub1, aes(x=Year, y = Corr.Val))+geom_point(alpha=0.8) +theme_classic() + theme(legend.position="bottom") + xlab("Year") +ylab("Wing melanism (grey level)")
#add trend
if(year.mod["P"]<0.05) fig2c= fig2c + geom_abline( aes(slope=year.mod["Estimate"],intercept=year.mod["Intercept"]))

#FIG 2
library(grid)

setwd(paste(mydir, "figures\\", sep=""))
pdf("Fig2_Loveland.pdf", height=8, width=4)

grid.newpage()
pushViewport(viewport(layout=grid.layout(3,1)))
vplayout<-function(x,y)
  viewport(layout.pos.row=x,layout.pos.col=y)
print(fig2a,vp=vplayout(1,1))
print(fig2b,vp=vplayout(2,1))
print(fig2c,vp=vplayout(3,1))

dev.off()

#------------
#REGIONS

#Statistics by region 
#region 1
areg= subset(absM.all, absM.all$region==1)
areg$siteID= match(areg$NewLocation, sites)
areg$YrSite= paste(areg$Year, areg$siteID, sep="")

#Bootstrap by region
phen.mod1= boot.lm(x=areg$doy162to202, y=areg$J, sites= areg$YrSite, Nruns,Nsamp)
plast.mod1= boot.lm(x=areg$J, y = areg$Corr.Val, sites= areg$YrSite, Nruns,Nsamp)
year.mod1= boot.lm(x=areg$Year, y = areg$Corr.Val, sites= areg$YrSite, Nruns,Nsamp)

#region 2
areg= subset(absM.all, absM.all$region==2)
areg$siteID= match(areg$NewLocation, sites)
areg$YrSite= paste(areg$Year, areg$siteID, sep="")

#Bootstrap by region
phen.mod2= boot.lm(x=areg$doy162to202, y=areg$J, sites= areg$YrSite, Nruns,Nsamp)
plast.mod2= boot.lm(x=areg$J, y = areg$Corr.Val, sites= areg$YrSite, Nruns,Nsamp)
year.mod2= boot.lm(x=areg$Year, y = areg$Corr.Val, sites= areg$YrSite, Nruns,Nsamp)

#region 3
areg= subset(absM.all, absM.all$region==3)
areg$siteID= match(areg$NewLocation, sites)
areg$YrSite= paste(areg$Year, areg$siteID, sep="")

#Bootstrap by region
phen.mod3= boot.lm(x=areg$doy162to202, y=areg$J, sites= areg$YrSite, Nruns,Nsamp)
plast.mod3= boot.lm(x=areg$J, y = areg$Corr.Val, sites= areg$YrSite, Nruns,Nsamp)
year.mod3= boot.lm(x=areg$Year, y = areg$Corr.Val, sites= areg$YrSite, Nruns,Nsamp)

#add slopes and intercepts
absM.all$phen.int= NA
if(phen.mod1["P"]<0.05) absM.all[which(absM.all$region.lab == "Region 1"),"phen.int"]= phen.mod1["Intercept"]
if(phen.mod2["P"]<0.05) absM.all[which(absM.all$region.lab == "Region 2"),"phen.int"]= phen.mod2["Intercept"]
if(phen.mod3["P"]<0.05) absM.all[which(absM.all$region.lab == "Region 3"),"phen.int"]= phen.mod3["Intercept"]

absM.all$phen.slope= NA
if(phen.mod1["P"]<0.05) absM.all[which(absM.all$region.lab == "Region 1"),"phen.slope"]= phen.mod1["Estimate"]
if(phen.mod2["P"]<0.05) absM.all[which(absM.all$region.lab == "Region 2"),"phen.slope"]= phen.mod2["Estimate"]
if(phen.mod3["P"]<0.05) absM.all[which(absM.all$region.lab == "Region 3"),"phen.slope"]= phen.mod3["Estimate"]

##ADD OTHER VARS

#-----------

#make region label
absM.all$region.lab<-"Region 1"
absM.all[which(absM.all$region==2), "region.lab"]<-"Region 2"
absM.all[which(absM.all$region==3), "region.lab"]<-"Region 3"

#Phenology 
fig3a<- ggplot(data=absM.all, aes(x=doy162to202, y = doy, color=Year ))+geom_point(alpha=0.8) +theme_classic()
#add trendlines
fig3a= fig3a+
  geom_abline(aes(slope=phen.slope,intercept=phen.int))+
 facet_wrap(~region.lab)+
  xlab("Developmental Temperature (°C)") +ylab("Phenology (doy)")+ theme(legend.position="bottom")+ scale_color_gradientn(colours = topo.colors(5))

#Abs
fig3b<- ggplot(data=absM.all, aes(x=doy, y = Corr.Val, color=Year ))+geom_point(alpha=0.8) +theme_classic() 
  #add trendlines
  fig3b= fig3b+
  geom_abline(aes(slope=plast.slope,intercept=plast.int))+
  facet_wrap(~region.lab)+ xlab("Phenology (doy)") +ylab("Wing melanism (grey level)")+scale_color_gradientn(colours = topo.colors(5))+ theme(legend.position="none")

#by year
fig3c<- ggplot(data=absM.all, aes(x=Year, y = Corr.Val))+geom_point(alpha=0.8) +theme_classic()+
  #add trendlines
  fig3c= fig3c+
  geom_abline(aes(slope=year.slope,intercept=year.int))+
  facet_wrap(~region.lab)+ xlab("Year") +ylab("Wing melanism (grey level)")

#---------
#Fig 3
setwd(paste(mydir, "figures\\", sep=""))
pdf("Fig3_Regions.pdf", height=8, width=8)

#grid.newpage()
pushViewport(viewport(layout=grid.layout(3,1)))
vplayout<-function(x,y)
  viewport(layout.pos.row=x,layout.pos.col=y)
print(fig3a,vp=vplayout(1,1))
print(fig3b,vp=vplayout(2,1))
print(fig3c,vp=vplayout(3,1))
dev.off()

#=====================
#STATISTICS
absM<- absM.all

   

#-------------------------------------
#Figure 2 Loveland Pass

abs.sub1$siteID= match(abs.sub1$NewLocation, sites)
abs.sub1$YrSite= paste(abs.sub1$Year, abs.sub1$siteID, sep="")

#Bootstrap
mod1= boot.lm(x=abs.sub1$doy162to202, y=abs.sub1$doy, sites= abs.sub1$YrSite, Nruns,Nsamp)

#---------------------------------------



#normalize
abs2<-scale(absM[,c("JJTave","doy","Year")], center=TRUE, scale=FALSE)
absM[,c("JJTave","doy","Year")]=abs2
#Remove NAs
dat= na.omit(absM[,c("Corr.Val","JJTave","doy","Year", "lon","lat","YrSite","region","NewLocation")])

#-------
boot.lm(x=absM$JJTave, y=absM$doy, sites= absM$YrSite, Nruns,Nsamp)

#--------------------
#SPATIAL ANALYSIS  

boot.spatlm(dat, sites= dat$YrSite, yvar="Corr.Val", xvars=c("JJTave","doy","doy:JJTave"), lon=dat$lon, lat=dat$lat, Nruns,Nsamp)

#All sites


mod1= boot.lm(x=dat$JJTave, y=dat$doy, sites= dat$YrSite, Nruns,Nsamp)




#=========================================

mod1= lm(Corr.Val~doy*JJTave, data=abs.sub1)

#normalize
abs2<-scale(abs.sub1[,c("J","Year","JJTave","JJTave.p")], center=TRUE, scale=FALSE)
abs.sub2= abs.sub1
abs.sub2[,c("J","Year","JJTave","JJTave.p")]=abs2

absM<- na.omit(abs.sub2[,c("Corr.Val","ThoraxC","NewLocation","J","Year","JJTave","JJTave.p")])

#----------

