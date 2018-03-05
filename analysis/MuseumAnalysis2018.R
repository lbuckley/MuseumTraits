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

setwd(paste(mydir, "analysis\\", sep=""))
source("spatial_functions_cleaned.R")

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
  
  for(r in 2:length(lon) ){
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

#After 1950?
absM<- subset(absM, Year >=1953 )

#Fix grey value
absM$Corr.Val= 1-absM$Corr.Val

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

#Region 2 weather station ~2000m too low, use ICAO environmental lapse rate of 6.5C/km to correct for purposed of plotting
#Makes too cold, move up 500m?
clim[which(clim$region==2),"TMEAN"]= clim[which(clim$region==2),"TMEAN"]-6.5/2

#calculate Julian
clim$DATE= paste(clim$YEAR,"/", clim$MONTH, "/", clim$DAY,sep="")
tmp <- as.POSIXlt(as.character(clim$DATE), format = "%Y/%m/%d")
clim$J=tmp$yday
clim$JYR= paste(clim$J, clim$YEAR,clim$region, sep="")
clim$MYR= paste(clim$MONTH, clim$YEAR,as.character(clim$region), sep="")

#----------------
#Plot temperature seasonality

clim1= subset(clim, region==1)
clim1= subset(clim1, YEAR>=1950)
clim1$time.per="pre 1975"
clim1$time.per[which(clim1$YEAR>=1975)]="post 1975"

max.plot= ggplot(clim1) + geom_smooth(aes(J, TMAX, group = YEAR, color = time.per), alpha = 0.2, se=FALSE) +xlim(160,240)+theme_classic()+ylab("Daily maximun temperature (°C)") +xlab("day of year")+labs(color="Time period")+ scale_colour_manual(values = c("gray","black"))+theme(legend.position=c(0.6,0.3))+ylim(8,20)

mean.plot= ggplot(clim1) + geom_smooth(aes(J, TMEAN, group = YEAR, color = time.per), alpha = 0.2, se=FALSE) +xlim(160,240)+theme_classic()+ylab("Daily mean temperature (°C)") +xlab("day of year")+labs(color="Time period")+ scale_colour_manual(values = c("gray","black"))+theme(legend.position=c(0.6,0.3))

min.plot= ggplot(clim1) + geom_smooth(aes(J, TMIN, group = YEAR, color = time.per), alpha = 0.2, se=FALSE) +xlim(160,240)+theme_classic()+ylab("Daily minimum temperature (°C)") +xlab("day of year")+labs(color="Time period")+ scale_colour_manual(values = c("gray","black"))+theme(legend.position=c(0.6,0.3))

#Fig S4
setwd(paste(mydir, "figures\\", sep=""))
pdf("FigS4_TempSeasonality.pdf", height=5, width=12)

plot_grid(max.plot, mean.plot, min.plot, align = "h", nrow = 1, rel_heights = c(1,1,1))

dev.off()

#---------------- 
#match to temp
absM$June=NA
absM$July=NA
absM$June15July15=NA
absM$doy162to202=NA
absM$Tpupal=NA

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
  
  #June 15 to July 15
  dev= paste( 152:181, absM$Year[i],absM$region[i], sep="")     
  cdat= clim[match(dev, clim$JYR),"TMEAN"]
  absM$June15July15[i]=mean(cdat, na.rm=TRUE)
  
  #doy162to202
  dev= paste( 162:202, absM$Year[i],absM$region[i], sep="")     
  cdat= clim[match(dev, clim$JYR),"TMEAN"]
  absM$doy162to202[i]=mean(cdat, na.rm=TRUE)
  
  #pupal
  pupal= paste( (absM$J[i]-26):(absM$J[i]-6), absM$Year[i],absM$region[i], sep="")      
  cdat= clim[match(pupal, clim$JYR),"TMEAN"]
  absM$Tpupal[i]=mean(cdat, na.rm=TRUE)
  
} #end loop rows

#JuneJuly
absM$JJTave= (absM$June + absM$July) / 2

#---------------------------
#DIVIDE BY TIME PERIOD
absM$time.per="pre 1980"
absM$time.per[which(absM$Year>1980) ]="post 1980"
absM$seas="early"
absM$seas[which(absM$doy>200) ]="late"
# absM.all$elev="low"
# absM.all$elev[which(absM.all$estElevation>3500) ]="high"
# absM.all$temp="cool"
# absM.all$temp[which(absM.all$doy162to202>10) ]="warm"
#-------------------------

#save
absM.all= absM

#set up bootstrap regressions
#BOOTSTRAP
#Assign site ID
sites= unique(absM.all$NewLocation)
absM.all$siteID= match(absM.all$NewLocation, sites)
absM.all$YrSite= paste(absM.all$Year, absM.all$siteID, sep="")

Nruns= 50 #50 #number of bootstrapp runs                          
Nsamp= 15 #max sample size of butterflies per site per year

#make region label
absM.all$region.lab<-"Region 1"
absM.all[which(absM.all$region==2), "region.lab"]<-"Region 2"
absM.all[which(absM.all$region==3), "region.lab"]<-"Region 3"

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
aper1.map<- map1 +geom_point(data=absM.all, aes(y=Lat, x=Long, color=estElevation) ) + coord_cartesian() + labs(x = "Longitude (°)",y="Latitude (°)", color="Elevation (m)") + theme(legend.position="bottom")+scale_color_gradientn(colours = topo.colors(5))+ theme(legend.key.width=unit(1,"cm"))

#-------------------
#Elevation inset plot

#Lower elevation at higher lats
elev.plot<-ggplot(data=absM.all, aes(x=lat, y = estElevation))+geom_point(alpha=0.8) +theme_classic()+
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

#cut post 1953
clim.dev= subset(clim.dev, YEAR>=1953 & YEAR<=2013)

#time series for three regions
clim.plot= ggplot(data=clim.dev, aes(x=YEAR, y = TMEAN, color=region.f))+geom_line() +theme_classic()+geom_smooth(method="lm",se=FALSE)+ylab("Developmental Temperature (°C)") +xlab("Year")+labs(color="Region") #+ theme(legend.position = c(0.2, 0.8))

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

#temp models
clim.sub= subset(clim.dev, region==3)
summary(lm(TMEAN~YEAR, data= clim.sub))

#-------------------

#COMBINE
library(grid)

setwd(paste(mydir, "figures\\", sep=""))
pdf("ElevFig.pdf", height=8, width=8)

##open pdf
subvp.t<-viewport(width=.5,height=.38,x=.74,y=0.8)
subvp.e<-viewport(width=.4,height=.40,x=.29,y=0.36)
##Next, open the main graph which was stored in b by typing b at the prompt:
aper1.map
##Then, superimpose the graph stored in a on the viewport as:
print(elev.plot,vp=subvp.e)
print(clim.plot,vp=subvp.t)
dev.off()

#============================================================

#LOVELAND PASS
#"EisenhowerTunnel"/ Loveland pass ANALYSIS
abs.sub1= absM.all[absM.all$NewLocation =="EisenhowerTunnel", ]

#-------------------------
#Statistics

abs.sub1$siteID= match(abs.sub1$NewLocation, sites)
abs.sub1$YrSite= paste(abs.sub1$Year, abs.sub1$siteID, sep="")

#Bootstrap
phen.mod= boot.lm(x=abs.sub1$Year, y=abs.sub1$J, sites= abs.sub1$YrSite, Nruns,Nsamp)
year.mod= boot.lm(x=abs.sub1$Year, y = abs.sub1$Corr.Val, sites= abs.sub1$YrSite, Nruns,Nsamp)


#------------------------
#phenology 
fig2a=ggplot(data=abs.sub1, aes(x=Year, y = doy, color=doy162to202 ))+geom_point(alpha=0.8) +theme_classic()+
  xlab("Year") +ylab("Phenology (doy)") + theme(legend.position="bottom") +scale_color_gradientn(colours = rev(heat.colors(5)))+ theme(legend.position="none")
#add trend
if(phen.mod["P"]<0.05) fig2a= fig2a + geom_abline( aes(slope=phen.mod["Estimate"],intercept=phen.mod["Intercept"]))

#Results in annual pattern
fig2c=ggplot(data=abs.sub1, aes(x=Year, y = Corr.Val, color=doy162to202))+geom_point(alpha=0.8) +theme_classic() + theme(legend.position="bottom")+ xlab("Year") +ylab("Wing melanism (grey level)")+scale_color_gradientn(colours = rev(heat.colors(5)))+ theme(legend.key.width=unit(1,"cm"))+labs(color="Developmental Temperature (°C)")
#add trend
if(year.mod["P"]<0.05) fig2c= fig2c + geom_abline( aes(slope=year.mod["Estimate"],intercept=year.mod["Intercept"]))

#FIG 2
library(grid)
library(cowplot)

setwd(paste(mydir, "figures\\", sep=""))
pdf("FigS1_Loveland.pdf", height=8, width=5)

plot_grid(fig2a,fig2c, align = "v", nrow = 2, rel_heights = c(1,1.4))
#plot_grid(fig2a, blank, fig2b,fig2bt, fig2c,blank, fig2d, fig2dt, align = "v", nrow = 4, rel_heights = c(1,1.4,1,1.4))

dev.off()

###################################################

#Figure 2: historical changes

#columns for slopes and intercepts
absM.all$tyear.int= NA
absM.all$tyear.slope= NA
absM.all$pyear.int= NA
absM.all$pyear.slope= NA
absM.all$myear.int= NA
absM.all$myear.slope= NA
absM.all$ryear.int= NA
absM.all$ryear.slope= NA

#STATS
for(reg in 1:3){
  areg= subset(absM.all, absM.all$region==reg)
  areg$siteID= match(areg$NewLocation, sites)
  areg$YrSite= paste(areg$Year, areg$siteID, sep="")
  areg= na.omit(areg[,c("ID","Corr.Val","doy","Tpupal","Year", "lon","lat","NewLocation","doy162to202","YrSite","time.per")])
  
  #match IDs
  match1= match(areg$ID, absM.all$ID)
  
  tyear.mod= boot.sar.lm(y=areg$doy162to202,x=areg$Year,lon=areg$lon,lat=areg$lat, sites= areg$YrSite, Nruns,Nsamp)
 
  pyear.mod= boot.sar.lm(y=areg$doy,x=areg$Year,lon=areg$lon,lat=areg$lat, sites= areg$YrSite, Nruns,Nsamp)
  
  myear.mod= boot.sar.lm(y=areg$Corr.Val,x=areg$Year,lon=areg$lon,lat=areg$lat, sites= areg$YrSite, Nruns,Nsamp)
  
  #save models
  if(reg==1){tyear.mod1=tyear.mod; pyear.mod1=pyear.mod; myear.mod1=myear.mod}
  if(reg==2){tyear.mod2=tyear.mod; pyear.mod2=pyear.mod; myear.mod2=myear.mod}
  if(reg==3){tyear.mod3=tyear.mod; pyear.mod3=pyear.mod; myear.mod3=myear.mod}
  
  #save slope and intercept
  inds= which(absM.all$region==reg)
  
  if(tyear.mod["P"]<0.05) absM.all[inds,"tyear.int"]= tyear.mod["Intercept"]
  if(tyear.mod["P"]<0.05) absM.all[inds,"tyear.slope"]= tyear.mod["Estimate"]
  
  if(pyear.mod["P"]<0.05) absM.all[inds,"pyear.int"]= pyear.mod["Intercept"]
  if(pyear.mod["P"]<0.05) absM.all[inds,"pyear.slope"]= pyear.mod["Estimate"]
  
  if(myear.mod["P"]<0.05) absM.all[inds,"myear.int"]= myear.mod["Intercept"]
  if(myear.mod["P"]<0.05) absM.all[inds,"myear.slope"]= myear.mod["Estimate"]
}

#---
tyear.mod3
pyear.mod3
myear.mod3

#--------------
#find unique
absM.all$RegYrDoy= paste(absM.all$region, absM.all$Year, absM.all$doy, sep="" )
dups= duplicated(absM.all$RegYrDoy)
abs.t<- absM.all[which(dups==FALSE),]

#Tdev 
fig2a<- ggplot(data=abs.t, aes(x=Year, y = doy162to202))+geom_point(alpha=0.8) +theme_classic()+
  ylab("Developmental Temperature (°C)") +xlab("year")+geom_smooth(method=lm, se=FALSE,color="black")
#add trendlines
fig2a= fig2a+
 # geom_abline(aes(slope=tyear.slope,intercept=tyear.int))+
  facet_wrap(~region.lab) 

#Phenology
fig2b<- ggplot(data=absM.all, aes(x=Year, y = doy, color=doy162to202 ))+geom_point(alpha=0.8) +theme_classic()+
  xlab("year") +ylab("Phenology (doy)")+ scale_color_gradientn(colours = rev(heat.colors(5)))+ theme(legend.position="none")+ theme(legend.key.width=unit(1,"cm"))
#add trendlines
fig2b= fig2b+
  geom_abline(aes(slope=pyear.slope,intercept=pyear.int))+
  facet_wrap(~region.lab)
#+ scale_color_gradientn(colours = topo.colors(5))

#melanism
#by year
fig2c<- ggplot(data=absM.all, aes(x=Year, y = Corr.Val, color=doy162to202))+geom_point(alpha=0.8) +theme_classic()+ xlab("year") +ylab("Wing melanism (gray level)")+ theme(legend.position="bottom")+scale_color_gradientn(colours = rev(heat.colors(5)))+ theme(legend.key.width=unit(1,"cm"))+labs(color="Developmental Temperature (°C)")

#add trendlines
fig2c= fig2c+
  geom_abline(aes(slope=myear.slope,intercept=myear.int))+
  facet_wrap(~region.lab)

#---------
setwd(paste(mydir, "figures\\", sep=""))
pdf("Fig2_hist.pdf", height=10, width=8)

plot_grid(fig2a, fig2b, fig2c, align = "v", nrow = 3, rel_heights = c(1,1,1.4))

dev.off()

#---------------------
#surface plots
library(akima)
library(cowplot)

#change Temp NaN to NA
absM.all$Tpupal[is.nan(absM.all$Tpupal)]<-NA
absM.all$doy162to202[is.nan(absM.all$doy162to202)]<-NA

#set limits
phen.lims= range(absM.all$doy)
Tpupal.lims= range(na.omit(absM.all$Tpupal))
Tdev.lims= range(na.omit(absM.all$doy162to202))
mel.lims= range(absM.all$Corr.Val)
year.lims= range(absM.all$Year)
  
for(reg in 1:3){

  areg= subset(absM.all, absM.all$region==reg)
  
#phen ~Tdev*year
  keep= which(!is.na(areg$doy162to202))
s=interp(x=areg$doy162to202[keep],y=areg$Year[keep],z=areg$doy[keep], duplicate="mean", nx=20, ny=30)

gdat <- interp2xyz(s, data.frame=TRUE)

plot.pty= ggplot(gdat) + 
  aes(x = x, y = y, z = z, fill = z) + 
  geom_tile() + 
  scale_fill_distiller(palette="Spectral", na.value="white", name="phenology (doy)", limits=phen.lims) +
  theme_bw(base_size=16)+xlab("developmental temperature")+ylab("year")+theme(legend.position="right")+xlim(Tdev.lims)+ylim(year.lims)

if(reg!=3) plot.pty= plot.pty +theme(legend.position="none")

#mel ~doy*year
s=interp(x=areg$doy,y=areg$Year,z=areg$Corr.Val, duplicate="mean", nx=20, ny=30)

gdat <- interp2xyz(s, data.frame=TRUE)

plot.mpy= ggplot(gdat) + 
  aes(x = x, y = y, z = z, fill = z) + 
  geom_tile() + 
  scale_fill_distiller(palette="Spectral", na.value="white", name="Wing melanism", limits=mel.lims) +
  theme_bw(base_size=16)+xlab("phenology (doy)")+ylab("year")+theme(legend.position="right")+xlim(phen.lims)+ylim(year.lims)

if(reg!=3) plot.mpy= plot.mpy +theme(legend.position="none")

#mel ~Tpup*year
#drop NAs for pupal temperature
keep= which(!is.na(areg$Tpupal))
s=interp(x=areg$Tpupal[keep],y=areg$Year[keep],z=areg$Corr.Val[keep], duplicate="mean", nx=20, ny=30)

gdat <- interp2xyz(s, data.frame=TRUE)

plot.mty= ggplot(gdat) + 
  aes(x = x, y = y, z = z, fill = z) + 
  geom_tile() + 
  scale_fill_distiller(palette="Spectral", na.value="white", name="Wing melanism", limits=mel.lims) +
  theme_bw(base_size=16)+xlab("pupal Temperature (°C)")+ylab("year")+theme(legend.position="right")+xlim(Tpupal.lims)+ylim(year.lims)

if(reg!=3) plot.mty= plot.mty +theme(legend.position="none")

#-----------------
#Assign region plots
if(reg==1){plot.pty1=plot.pty; plot.mpy1=plot.mpy; plot.mty1=plot.mty}
if(reg==2){plot.pty2=plot.pty; plot.mpy2=plot.mpy; plot.mty2=plot.mty}
if(reg==3){plot.pty3=plot.pty; plot.mpy3=plot.mpy; plot.mty3=plot.mty}

}
#end region loop

#------------
setwd(paste(mydir, "figures\\", sep=""))
pdf("Fig3_corr.pdf", height=10, width=10)

plot_grid(plot.pty1, plot.pty2, plot.pty3, plot.mpy1, plot.mpy2, plot.mpy3, plot.mty1, plot.mty2, plot.mty3, align = "v", nrow = 3, rel_widths = c(1,1,1.8))

dev.off()
--------------
#stats for correlations
   
for( reg in 1:3 ){
  areg.a= subset(absM.all, absM.all$region==reg)
  
  #phen ~Tdev*year
  areg= na.omit(areg.a[,c("doy", "doy162to202", "Year", "lon", "lat", "YrSite")] )
  areg[,1:3]<-scale(areg[,1:3], center=TRUE, scale=TRUE)
pty.mod= boot.sar.mult(y=areg$doy,x1=areg$doy162to202,x2=areg$Year,lon=areg$lon,lat=areg$lat, sites=areg$YrSite, Nruns,Nsamp)
    
  #mel ~doy*year
areg= na.omit(areg.a[,c("doy", "Corr.Val", "Year", "lon", "lat", "YrSite")] )
areg[,1:3]<-scale(areg[,1:3], center=TRUE, scale=TRUE)
mpy.mod= boot.sar.mult(y=areg$Corr.Val,x1=areg$doy,x2=areg$Year,lon=areg$lon,lat=areg$lat, sites=areg$YrSite, Nruns,Nsamp)  
  
  #mel ~Tpup*year
areg= na.omit(areg.a[,c("Corr.Val", "Tpupal", "Year", "lon", "lat", "YrSite")] )
areg[,1:3]<-scale(areg[,1:3], center=TRUE, scale=TRUE)
mty.mod= boot.sar.mult(y=areg$Corr.Val,x1=areg$Tpupal,x2=areg$Year,lon=areg$lon,lat=areg$lat, sites=areg$YrSite, Nruns,Nsamp)  

#mel ~Tpup*doy
areg= na.omit(areg.a[,c("Corr.Val", "Tpupal", "doy", "lon", "lat", "YrSite")] )
areg[,1:3]<-scale(areg[,1:3], center=TRUE, scale=TRUE)
mtp.mod= boot.sar.mult(y=areg$Corr.Val,x1=areg$Tpupal,x2=areg$doy,lon=areg$lon,lat=areg$lat, sites=areg$YrSite, Nruns,Nsamp)  

#Assign stats
if(reg==1){pty.mod1=pty.mod; mpy.mod1=mpy.mod; mty.mod1=mty.mod; mtp.mod1=mtp.mod}
if(reg==2){pty.mod2=pty.mod; mpy.mod2=mpy.mod; mty.mod2=mty.mod; mtp.mod2=mtp.mod}
if(reg==3){pty.mod3=pty.mod; mpy.mod3=mpy.mod; mty.mod3=mty.mod; mtp.mod3=mtp.mod}

}# loop regions

#----
corrs= rbind(pty.mod1, mpy.mod1,mty.mod1, mtp.mod1,pty.mod2, mpy.mod2,mty.mod2, mtp.mod2,pty.mod3, mpy.mod3,mty.mod3, mtp.mod3)
write.csv(corrs, "CorrelationStats.csv")

#======================
#Residual plots

absM.all$Corr.Resid= NA

for(reg in 1:3){
 if(reg==1) mtp.mod<- mtp.mod1
 if(reg==2) mtp.mod<- mtp.mod2
 if(reg==3) mtp.mod<- mtp.mod3
  
 inds= which(absM.all$region==reg)
 absM.all$Corr.Resid[inds]<- absM.all$Corr.Val[inds]-(mtp.mod["Estimate.x1"]*absM.all$Tpupal[inds] +mtp.mod["Estimate.x2"]*absM.all$doy[inds] +mtp.mod["Estimate.x1:x2"]*absM.all$Tpupal[inds]*absM.all$doy[inds] +mtp.mod["Intercept"])
 
 #regression stats
 areg= subset(absM.all, absM.all$region==reg)
 areg$siteID= match(areg$NewLocation, sites)
 areg$YrSite= paste(areg$Year, areg$siteID, sep="")
 areg= na.omit(areg[,c("ID","Corr.Resid","Year", "lon","lat","NewLocation","YrSite")])
 
 #match IDs
 match1= match(areg$ID, absM.all$ID)
 
 ryear.mod= boot.sar.lm(y=areg$Corr.Resid,x=areg$Year,lon=areg$lon,lat=areg$lat, sites= areg$YrSite, Nruns,Nsamp)
 
 #save models
 if(reg==1){ryear.mod1=ryear.mod}
 if(reg==2){ryear.mod2=ryear.mod}
 if(reg==3){ryear.mod3=ryear.mod}
 
 #save slope and intercept
 inds= which(absM.all$region==reg)
 
 if(ryear.mod["P"]<0.05) absM.all[inds,"ryear.int"]= ryear.mod["Intercept"]
 if(ryear.mod["P"]<0.05) absM.all[inds,"ryear.slope"]= ryear.mod["Estimate"]
 } #end region loop

#----------------
fig4<- ggplot(data=absM.all, aes(x=Year, y = Corr.Resid, color=doy162to202))+geom_point(alpha=0.8) +theme_classic()+ xlab("Year") +ylab("Residuals(Wing melanism) (gray level)")+ theme(legend.position="bottom")+scale_color_gradientn(colours = rev(heat.colors(5)))+labs(color="Developmental Temperature (°C)")
#add trendlines
fig4= fig4+
  geom_abline(aes(slope=ryear.slope,intercept=ryear.int))+
  facet_wrap(~region.lab)

#---------
setwd(paste(mydir, "figures\\", sep=""))
pdf("Fig4_resid.pdf", height=4, width=8)

fig4

dev.off()

#===========================================
#THORAX AND SIZE

#columns for slopes and intercepts
absM.all$syear.int= NA
absM.all$syear.slope= NA
absM.all$fyear.int= NA
absM.all$fyear.slope= NA

#STATS
for(reg in 1:3){
  areg= subset(absM.all, absM.all$region==reg)
  areg$siteID= match(areg$NewLocation, sites)
  areg$YrSite= paste(areg$Year, areg$siteID, sep="")
  areg= na.omit(areg[,c("ID","Thorax","FWL","Year", "lon","lat","NewLocation","doy162to202","YrSite","time.per")])
  
  #match IDs
  match1= match(areg$ID, absM.all$ID)
  
  syear.mod= boot.sar.lm(y=areg$FWL,x=areg$Year,lon=areg$lon,lat=areg$lat, sites= areg$YrSite, Nruns,Nsamp)
  fyear.mod= boot.sar.lm(y=areg$Thorax,x=areg$Year,lon=areg$lon,lat=areg$lat, sites= areg$YrSite, Nruns,Nsamp)
  
  #save models
  if(reg==1){syear.mod1=syear.mod; fyear.mod1=fyear.mod}
  if(reg==2){syear.mod2=syear.mod; fyear.mod2=fyear.mod}
  if(reg==3){syear.mod3=syear.mod; fyear.mod3=fyear.mod}
  
  #save slope and intercept
  inds= which(absM.all$region==reg)
  
  if(syear.mod["P"]<0.05) absM.all[inds,"syear.int"]= syear.mod["Intercept"]
  if(syear.mod["P"]<0.05) absM.all[inds,"syear.slope"]= syear.mod["Estimate"]
  
  if(fyear.mod["P"]<0.05) absM.all[inds,"fyear.int"]= fyear.mod["Intercept"]
  if(fyear.mod["P"]<0.05) absM.all[inds,"fyear.slope"]= fyear.mod["Estimate"]
}

#---
syear.mod3
fyear.mod3

#--------------
#thorax
#by year
figs2a<- ggplot(data=absM.all, aes(x=Year, y = FWL, color=doy162to202))+geom_point(alpha=0.8) +theme_classic()+ xlab("year") +ylab("Forewing length (mm)")+ theme(legend.position="none")+scale_color_gradientn(colours = rev(heat.colors(5)))+ theme(legend.key.width=unit(1,"cm"))+labs(color="Developmental Temperature (°C)")

#add trendlines
figs2a= figs2a+
  geom_abline(aes(slope=syear.slope,intercept=syear.int))+
  facet_wrap(~region.lab)

#-------------------
#setae
#by year
figs2b<- ggplot(data=absM.all, aes(x=Year, y = Thorax, color=doy162to202))+geom_point(alpha=0.8) +theme_classic()+ xlab("year") +ylab("Setae length (mm)")+ theme(legend.position="bottom")+scale_color_gradientn(colours = rev(heat.colors(5)))+ theme(legend.key.width=unit(1,"cm"))+labs(color="Developmental Temperature (°C)")

#add trendlines
figs2b= figs2b+
  geom_abline(aes(slope=fyear.slope,intercept=fyear.int))+
  facet_wrap(~region.lab)

#---------
setwd(paste(mydir, "figures\\", sep=""))
pdf("Figs2_sizefur.pdf", height=8, width=8)

plot_grid(figs2a, figs2b, align = "v", nrow = 2, rel_heights = c(1,1.4))

dev.off()

###############################################
#check data
areg= subset(absM.all, absM.all$region==1)

plot1=ggplot(data=areg, aes(x=Year, y = min_Standard))+geom_point()+geom_smooth(method="lm")
plot2=ggplot(data=areg, aes(x=Year, y = max_Sample))+geom_point()+geom_smooth(method="lm")
plot3=ggplot(data=areg, aes(x=Year, y = mean_Standard))+geom_point()+geom_smooth(method="lm")
plot_grid(plot1, plot2, plot3, align = "v", nrow = 1)


ggplot(data=areg, aes(x=Year, y = (areg$max_Sample- areg$min_Sample) ))+geom_point()+geom_smooth(method="lm")

areg$grey= 1-(areg$mean_Sample- areg$mean_Standard)/(255- areg$min_Standard  )

areg$grey= 1- (areg$mean_Sample- areg$min_Standard)/(255- areg$min_Standard)
areg$grey= 1- (areg$mean_Standard- areg$min_Standard)/(255- areg$min_Standard)

ggplot(data=areg, aes(x=Year, y = grey))+geom_point()+geom_smooth(method="lm")
ggplot(data=areg, aes(x=Year, y = Corr.Val))+geom_point()+geom_smooth(method="lm")

