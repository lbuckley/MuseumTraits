#Try prism data
#http://ropensci.github.io/prism/
#https://rpubs.com/collnell/get_prism

library(raster)

setwd("C:\\Users\\lbuckley\\Desktop\\Fall2017\\MuseumTraits\\PRISM_tmean_stable_4kmM2_198101_201705_bil\\")

#sort by year
absM= absM.all[sort(absM.all$Year),]
#years
years= unique(absM.all$Year)
years.rec= years[ which(years %in% 1981:2017) ]

p.dat= array(data=NA, dim=c(length(years.rec), 2, nrow(absM.all) ))

#extract data for each specimen (speed up by restricting to unique specimen)
pts<-SpatialPointsDataFrame(coords=absM.all[,c('lon','lat')], 
                            data=absM.all, proj4string = CRS("+proj=longlat +ellps=WGS84 +no_defs"))

#look at overall temporal temperature trend by examining all locations

for(k in 1:length(years.rec) ){
  
  #JUNE
  f.june= paste("PRISM_tmean_stable_4kmM2_",years.rec[k],"06_bil.bil",sep="")
  p<- raster(f.june)
  ##crop
  #ext= extent(t(matrix(c( -127,-98,30,61), nrow=2)))
  #p.crop=crop(p, ext)
  
  p.dat[k,1,]<-extract(p, pts, na.rm=FALSE)
  
  #JULY
  f.july= paste("PRISM_tmean_stable_4kmM2_",years.rec[k],"07_bil.bil",sep="")
  p<- raster(f.july)
  
  p.dat[k,2,]<-extract(p, pts, na.rm=FALSE)
  
}

p.rec.june= cbind(years.rec, p.dat[,1,])
p.rec.july= cbind(years.rec, p.dat[,2,])

#write out
setwd("C:\\Users\\lbuckley\\Desktop\\Fall2017\\MuseumTraits\\")
write.csv(cbind(years.rec,p.dat[,1,]), "june_recent_tmean.csv" )
write.csv(cbind(years.rec,p.dat[,2,]), "july_recent_tmean.csv" )

#------
#Historic data 

setwd("C:\\Users\\lbuckley\\Desktop\\Fall2017\\MuseumTraits\\PRISM_1895_1980\\")

years.hist= years[ which(years %in% 1895:1980) ]

p.dat= array(data=NA, dim=c(length(years.hist), 2, nrow(absM.all) ))

#look at overall temporal temperature trend by examining all locations

for(k in 1:length(years.hist) ){
  
  #JUNE
  f.june= paste("PRISM_tmean_stable_4kmM2_",years.hist[k],"06_bil.bil",sep="")
  p<- raster(f.june)
  ##crop
  #ext= extent(t(matrix(c( -127,-98,30,61), nrow=2)))
  #p.crop=crop(p, ext)
  
  p.dat[k,1,]<-extract(p, pts, na.rm=FALSE)
  
  #JULY
  f.july= paste("PRISM_tmean_stable_4kmM2_",years.hist[k],"07_bil.bil",sep="")
  p<- raster(f.july)
  
  p.dat[k,2,]<-extract(p, pts, na.rm=FALSE)
  
}

p.hist.june= cbind(years.hist, p.dat[,1,])
p.hist.july= cbind(years.hist, p.dat[,2,])

#write out
setwd("C:\\Users\\lbuckley\\Desktop\\Fall2017\\MuseumTraits\\")
write.csv(cbind(years.hist,p.dat[,1,]), "june_historic_tmean.csv" )
write.csv(cbind(years.hist,p.dat[,2,]), "july_historic_tmean.csv" )

#=====================================
#ADD ALBERTA DATA

#Add Alberta Data
#https://sites.ualberta.ca/~ahamann/data.html

alb.dat= subset(absM.all, absM.all$State=="Alberta")
alb.years= sort(unique(alb.dat$Year))

setwd(paste(mydir, "data\\AlbertaClimate\\", sep=""))

#combine files
a.clim.june= matrix(NA, 175, length(alb.years))
a.clim.july= matrix(NA, 175, length(alb.years))

for(yr in 1:length(alb.years))
{
file= paste("AlbertaSites_year_",alb.years[yr],"M.csv", sep="")
dat= read.csv(file, header=TRUE) 

a.clim.june[,yr]= dat$TAV06
a.clim.july[,yr]= dat$TAV07
}

colnames(a.clim.june)= paste("june_",alb.years, sep="")
colnames(a.clim.july)= paste("july_",alb.years, sep="")

#plot trend
a.ave.june=colMeans(a.clim.june)
a.ave.july=colMeans(a.clim.july)
a.ave= rbind(a.ave.june,a.ave.july )
a.ave= colMeans(a.ave)
plot(a.ave~alb.years)

#combine
a.clim= cbind( dat[,1:5], a.clim.june, a.clim.july)

#write out
setwd(paste(mydir, "data\\", sep=""))
write.csv(a.clim, "AlbertaClimate.csv" )

