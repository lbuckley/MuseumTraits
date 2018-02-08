mydir= "C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\MuseumTraits\\"
#mydir= "C:\\Users\\lbuckley\\Documents\\MuseumTraits\\"
#mydir= "C:\\Users\\Buckley\\Google Drive\\ColiasEvolution\\HistoricalAnalysis\\"


#PROCESS Region 3 DATA

#Read and combine files
yrs= 1962:2011
setwd(paste(mydir, "data\\BlueHillLO\\", sep=""))

for(yr in yrs){
  file= paste("eng-daily-0101",yr,"-1231",yr,".csv", sep=""  )
  alb=read.csv(file, skip=25)
  
  if(yr==1962) bluehill<- alb[,c("Year","Month","Day","Mean.Temp..Â.C.")]
  if(yr>1962) bluehill<- rbind(bluehill, alb[,c("Year","Month","Day","Mean.Temp..Â.C.")])
}

#check NAs for june and july
bluehill.s= bluehill[which(bluehill$Month %in% c(5,6,7)),]

#write.out
write.csv(bluehill.s, "bluehillLO_19622011.csv")
#many NAs after 2007

#-----------------------------------
#COMBINE CLIMATE DATA, use weather stations

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
climax$region=1

#Region 2
setwd(paste(mydir, "data\\", sep=""))
moran=read.csv("Moran5WNWCOOP_F.csv")
#F to C
moran$TMAX= (moran$TMAX-32)*5/9
moran$TMIN= (moran$TMIN-32)*5/9
moran$TMEAN= (moran$TMAX + moran$TMIN)/2
moran$region=2

#Regon 3
setwd(paste(mydir, "data\\BlueHillLO\\", sep=""))
bluehill=read.csv("bluehillLO_19622007.csv")
bluehill$region=3
#many NAs after 2007

#combine
clim= rbind( climax[,c("YEAR","MONTH","DAY","TMEAN","region")], moran[,c("YEAR","MONTH","DAY","TMEAN","region")], bluehill[,c("YEAR","MONTH","DAY","TMEAN","region")])

#write out
write.csv(clim,"WeatherStationData_Regions.csv")

