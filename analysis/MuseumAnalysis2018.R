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

#daymetR, http://khufkens.github.io/daymetr/
#if(!require(devtools)){install.package(devtools)}
#devtools::install_github("khufkens/daymetr")
#library(daymetr)

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

#Region 1
setwd(paste(mydir, "data\\", sep=""))
climax=read.csv("ClimaxCOOP_F.csv")
#F to C
climax$TMAX= (climax$TMAX-32)*5/9
climax$TMIN= (climax$TMIN-32)*5/9
climax$TMEAN= (climax$TMAX + climax$TMIN)/2

#Region 2
setwd(paste(mydir, "data\\", sep=""))
moran=read.csv("Moran5WNWCOOP_F.csv")
#F to C
moran$TMAX= (moran$TMAX-32)*5/9
moran$TMIN= (moran$TMIN-32)*5/9
moran$TMEAN= (moran$TMAX + moran$TMIN)/2

#Regon 3
setwd(paste(mydir, "data\\BlueHillLO\\", sep=""))
bluehill=read.csv("bluehillLO_19622011.csv")
#many NAs after 2007
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

#counts
clim.counts= aggregate(clim.jj, list(clim.jj$Year), FUN="count")
clim.years= clim.counts[which(clim.counts$STATION>50) ,"Group.1"]
clim.years= as.numeric(as.character(clim.years))
#require 50 days of data #FIx 1980?
clim.jj= clim.jj[which(clim.jj$Year %in% clim.years),]

clim.jj= aggregate(clim.jj, list(clim.jj$Year), FUN="mean", na.rm=TRUE)
clim.jj$Year= clim.jj$Group.1
clim.jj$Year= as.numeric(as.character(clim.jj$Year))

match1= match(absM$Year, clim.jj$Year)
absM$JJTave= clim.jj$TMEAN[match1] 

#save
absM.all= absM

#-------
#ADD prism data
#http://ropensci.github.io/prism/
#https://rpubs.com/collnell/get_prism

#recent
setwd(paste(mydir, "data\\PRISMClimate\\", sep=""))

p.rec.june= read.csv("june_recent_tmean.csv", header=TRUE, row.names=1) 
p.rec.july= read.csv("july_recent_tmean.csv", header=TRUE, row.names=1)

#Historic data 
p.hist.june= read.csv("june_historic_tmean.csv", header=TRUE, row.names=1) 
p.hist.july= read.csv("july_historic_tmean.csv", header=TRUE, row.names=1)

#make names match
colnames(p.rec.june)[1]="years"
colnames(p.rec.july)[1]="years"
colnames(p.hist.june)[1]="years"
colnames(p.hist.july)[1]="years"

#-----------------
#Plot temperature trend

#combine
p.june= rbind(p.hist.june, p.rec.june)
p.july= rbind(p.hist.july, p.rec.july)

p.june.ave= rowMeans(p.june[,2:983], na.rm=TRUE)
p.july.ave= rowMeans(p.july[,2:983], na.rm=TRUE)

#ave of june and july
p.jj.ave= rowMeans(cbind(p.june.ave,p.july.ave))

p.june.ave= cbind(p.june[,1], p.june.ave, rep("june",73))
p.july.ave= cbind(p.june[,1], p.july.ave, rep("july",73))
p.jj.ave=  cbind(p.june[,1], p.jj.ave, rep("jj",73))

#plot
p.clim=p.jj.ave
#p.clim= rbind(p.june.ave, p.july.ave )
p.clim= as.data.frame(p.clim)
colnames(p.clim)=c("year","temp","month" )
p.clim$year= as.numeric(as.character(p.clim$year))
p.clim$temp= as.numeric(as.character(p.clim$temp))

#----------
#Add PRISM data
match.years= match(absM.all$Year, p.june[,1])

absM.all$Tjune.prism= NA
absM.all$Tjuly.prism= NA

#match data
for(k in 1:nrow(absM.all) ){
  absM.all$Tjune.prism[k]= p.june[match.years[k],k+1]
  absM.all$Tjuly.prism[k]= p.july[match.years[k],k+1]
}

#ave of june and july
absM.all$JJTave.p= rowMeans(absM.all[,c("Tjune.prism", "Tjuly.prism")])

#write out data including PRISM data
write.csv(absM.all, "MuseumData_wPRISM.csv")

#--------------------
#For loveland pass
#Cabin creek 051186, http://climate.colostate.edu/data_access.html
#http://climatetrends.colostate.edu
#clim= read.csv("CabinCreek.csv", na.strings = "-9999")  #F and IN

#Climax almost as close and more similar elevation

#---------
#Add Alberta Data
#https://sites.ualberta.ca/~ahamann/data.html

#read in
setwd(paste(mydir, "data\\", sep=""))
a.clim= read.csv("AlbertaClimate.csv" )

match1= match(a.clim$ID1, absM.all$ID)

alb.ind= which(absM.all$"State"=="Alberta")

#match data
for(k in 1:length(alb.ind) ){
  june.name=paste("june_", absM.all$Year[alb.ind[k]], sep="")
  july.name=paste("june_", absM.all$Year[alb.ind[k]], sep="")
  
  absM.all$Tjune.prism[alb.ind[k]]= a.clim[match1[k],june.name]
  absM.all$Tjuly.prism[alb.ind[k]]= a.clim[match1[k],july.name]
}

#ave of june and july
absM.all$JJTave.p= rowMeans(absM.all[,c("Tjune.prism", "Tjuly.prism")])

#=========================
#Use weather station data

#Region 1: Climax
#Region 2: MORAN 5WNW, WY USA, 	USC00486440, http://mco.cfc.umt.edu/ghcn/station/USC00486440.html
#Region 3: BLUEHILL LO, 1,950.70 m http://climate.weather.gc.ca/climate_data/

abs.co=subset(absM.all, absM.all$region==1)
summary(abs.co)


ggplot(data=abs.co, aes(x=JJTave, y = JJTave.p))+geom_point()+geom_smooth(method="lm")+theme_classic()


#Region 2
setwd(paste(mydir, "data\\", sep=""))
moran=read.csv("Moran5WNWCOOP_F.csv")

#Regon 3
setwd(paste(mydir, "data\\BlueHillLO\\", sep=""))
bluehill=read.csv("bluehillLO_19622011.csv")
#many NAs after 2007

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

#PLOT TEMP TREND PRISM DATA ACROSS ALL SITES IN ALL YEARS WITH SAMPLES
clim.plot= ggplot(data=p.clim, aes(x=year, y = temp))+geom_point() +theme_classic()+geom_smooth()+xlim(1950,2013)+ylab("June and July Temperature (°C)") +xlab("Year")

#-------------------

#COMBINE
library(grid)

setwd(paste(mydir, "figures\\", sep=""))
pdf("ElevFig.pdf", height=8, width=8)

##open pdf
subvp.t<-viewport(width=.38,height=.38,x=.7,y=0.8)
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
mod1= lm(Corr.Val~doy*JJTave, data=abs.sub1)

#normalize
abs2<-scale(abs.sub1[,c("PupalTave","ParentTave","J","Year","JJTave","JJTave.p")], center=TRUE, scale=FALSE)
abs.sub2= abs.sub1
abs.sub2[,c("PupalTave","ParentTave","J","Year","JJTave","JJTave.p")]=abs2

absM<- na.omit(abs.sub2[,c("Corr.Val","ThoraxC","NewLocation","PupalTave","ParentTave","J","Year","JJTave","JJTave.p")])

#----------
#phenology 
fig2a=ggplot(data=abs.sub1, aes(x=JJTave, y = J, color=Year ))+geom_point(alpha=0.8) +theme_classic()+
  xlab("June and July Temperature (°C)") +ylab("Phenology (doy)") +geom_smooth(method="lm")+ theme(legend.position="bottom") +scale_color_gradientn(colours = topo.colors(5))+ theme(legend.position="bottom") 

#plasticity
fig2b=ggplot(data=abs.sub1, aes(x=J, y = Corr.Val, color=Year ))+geom_point(alpha=0.8) +theme_classic()+  geom_smooth(method="lm") + xlab("Phenology (doy)") +ylab("Wing melanism (grey level)")+scale_color_gradientn(colours = topo.colors(5))+ theme(legend.position="none")

#Results in annual pattern
fig2c=ggplot(data=abs.sub1, aes(x=Year, y = Corr.Val))+geom_point(alpha=0.8) +theme_classic()+  geom_smooth(method="lm") + theme(legend.position="bottom") + xlab("Year") +ylab("Wing melanism (grey level)")

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

#make region label
absM.all$region.lab<-"Region 1"
absM.all[which(absM.all$region==2), "region.lab"]<-"Region 2"
absM.all[which(absM.all$region==3), "region.lab"]<-"Region 3"

#Phenology 
fig3a<- ggplot(data=absM.all, aes(x=JJTave, y = doy, color=Year ))+geom_point(alpha=0.8) +theme_classic()+ facet_wrap(~region.lab)+geom_smooth(method="lm")+
  xlab("June and July Temperature (°C)") +ylab("Phenology (doy)")+ theme(legend.position="bottom")+ scale_color_gradientn(colours = topo.colors(5))
##! PRISM data reverses trend. CHECK
#Tjune.prism

#Abs
fig3b<- ggplot(data=absM.all, aes(x=doy, y = Corr.Val, color=Year ))+geom_point(alpha=0.8) +theme_classic()+ facet_wrap(~region.lab)+geom_smooth(method="lm")+ xlab("Phenology (doy)") +ylab("Wing melanism (grey level)")+scale_color_gradientn(colours = topo.colors(5))+ theme(legend.position="none")

#by year
fig3c<- ggplot(data=absM.all, aes(x=Year, y = Corr.Val))+geom_point(alpha=0.8) +theme_classic()+facet_wrap(~region.lab)+geom_smooth(method="lm")+ xlab("Year") +ylab("Wing melanism (grey level)")

#---------
#Fig 3
setwd(paste(mydir, "figures\\", sep=""))
pdf("Fig3_Regions_climax.pdf", height=8, width=8)

#grid.newpage()
pushViewport(viewport(layout=grid.layout(3,1)))
vplayout<-function(x,y)
  viewport(layout.pos.row=x,layout.pos.col=y)
print(fig3a,vp=vplayout(1,1))
print(fig3b,vp=vplayout(2,1))
print(fig3c,vp=vplayout(3,1))
dev.off()

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
absM<- absM.all

#BOOTSTRAP
#Assign site ID
sites= unique(absM$NewLocation)
absM$siteID= match(absM$NewLocation, sites)
absM$YrSite= paste(absM$Year, absM$siteID, sep="")

Nruns= 50 #50 #number of bootstrapp runs                          
Nsamp= 15 #max sample size of butterflies per site per year   

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




