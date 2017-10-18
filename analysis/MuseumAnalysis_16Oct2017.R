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

#---------------------------
#READ AND MANIPULATE DATA

#read data
setwd(paste(mydir, "data\\", sep=""))
abs=read.csv("FullDataColiasMuseums.csv")


absM<- subset(abs,Sex =="M")

#group nearby localities
absM$NewLocation<- ClassifySites(absM$Location, absM$Long, absM$Lat, DistCut=5)    #radius of 5km

#calculate Julian
mdy=  paste(absM$Month,absM$Day,absM$Year, sep="-")
tmp <- as.POSIXlt(mdy, format = "%m-%d-%Y")
absM$J=tmp$yday

#abbreviate columns to control for NAs
#absM= absM[,c("ID","NewLocation","Collection", "grey","Thorax","JulyAv","JulyMax","JuneAv","Year","estElevation","J", "Long","Lat")]
absM= absM[,colnames(absM)!="BodyLength"]
# take out body length due to NAs

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
absM= na.omit(absM)

#GET RID OF TWO SPECIMENS WITH J=171
absM= absM[absM$J>171,]

#match to temp
absM$AdultTmax=NA
absM$AdultTave=NA
absM$ParentTmax=NA
absM$ParentTave=NA
absM$PupalTmax=NA
absM$PupalTave=NA
absM$LifeTmax=NA
absM$LifeTave=NA


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
absM= na.omit(absM)

#normalize

abs2<-scale(absM[,c("PupalTave","ParentTave","J","Year")], center=TRUE, scale=FALSE)

#library(clusterSim)
#abs2= data.Normalization (absM[,c("PupalTave","ParentTave","J","Year")],type="n1",normalization="column")
absM[,c("PupalTave","ParentTave","J","Year")]=abs2

#@  absM<-subset(absM,absM$Collection=="UF"|abM3$Collection=="Yale"|absM$Collection=="MW")       #@ run subset for thorax data to confirm. I did this in a hacky way for the sake of time. 

#---------------------------------------------------
#BOOTSTRAP

#Assign site ID
sites= unique(absM$NewLocation)
absM$siteID= match(absM$NewLocation, sites)
absM$YrSite= paste(absM$Year, absM$siteID, sep="")

Nruns= 50 #50 #number of bootstrapp runs                           #@have not changed since we agreed. 
Nsamp= 15 #max sample size of butterflies per site per year    #@have not changed since we agreed. 

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

xlabs= c("Pupal T (?C)","Flight season T (?C)", "Date of collection (J)", "Year", "Pupal T (?C): Date")

#@ xlabs= c("Pupal T (°C)","Parent T (°C)", "Date of collection (J)", "Year", "Pupal T(°C) : Date", "Parent T (°C): Date")  #@ Changed Labs to be capitalized. Doesn't really matter... 

setwd(paste(mydir, "figs\\", sep=""))
pdf("GreyPR_norm.pdf",height = 6, width = 10)

par(mfrow=c(2,3), cex=1.1, lwd=1, mar=c(3,2, 1, 1), mgp=c(1.3, 0.5, 0), oma=c(0,2,0,0), bty="l", cex.lab=1.2)
my.prplot(tarlm, 1, xlabs)
my.prplot(tarlm, 2, xlabs)
my.prplot(tarlm, 3, xlabs)
my.prplot(tarlm, 4, xlabs)
my.prplot(tarlm, 5, xlabs)
#my.prplot(tarlm, 6, xlabs)

mtext("Partial residual of wing melanism (grey level)", side=2, line = -0.5, cex=1.3, outer=TRUE)

dev.off()
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

pdf("FurPR_norm.pdf",height = 6, width = 10)

par(mfrow=c(2,3), cex=1.1, lwd=1, mar=c(3,2, 1, 1), mgp=c(1.3, 0.5, 0), oma=c(0,2,0,0), bty="l", cex.lab=1.2)
my.prplot(tarlm, 1, xlabs)
my.prplot(tarlm, 2, xlabs)
my.prplot(tarlm, 3, xlabs)
my.prplot(tarlm, 4, xlabs)
my.prplot(tarlm, 5, xlabs)
#my.prplot(tarlm, 6, xlabs)

mtext("Partial residual of setae length (mm)", side=2, line = -0.5, cex=1.3, outer=TRUE)

dev.off()
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

#--------------------------
#LINEAR MODELS
mod1 <- lm(absM$grey~absM$PupalTave+absM$ParentTave+absM$J+absM$Year+absM$PupalTave*absM$J, data=absM.boot)
mod1.RE <- lme(grey~PupalTave+ParentTave+J+Year+PupalTave*J+ParentTave*J,random=~1|NewLocation, data=absM.boot, method="ML")      
 #@to include and serve the same purpose as before, should we not include the "new location" term as random and run it on the full data set? From my perspective, the point of the model comparison has changed since adding the bootstrap (which is why I thought we may want to leave it out). Is the question- do we get the same result in accounting for space in two different ways? OR is it- do we get the same qualitat  ive result when running the simplest linear model on a subset of the data?  
summary(mod1)
anova(mod1)
dist1m<-dist(cbind(absM.boot$Latitude, absM.boot$Long))
dist0m<- dist(absM.boot$grey)
dist1Rm<- dist(residuals(mod1.RE))
  #@ These models perform actually take care of the spatial autocorrelation as test with a mantel test AND approximately align with SAR significance 

mod2 <- lm(absM$Thorax~absM$PupalTave+absM$ParentTave+absM$J+absM$Year+absM$PupalTave*absM$J, data=absM.boot)

mod2.RE <- lme(Thorax~PupalTave+ParentTave+J+Year+PupalTave*J+ParentTave*J, random=~1|NewLocation, data=absM.boot, method="ML")        #@ see above logic
summary(mod2)
anova(mod2)

dist1m<-dist(cbind(absM.boot$Latitude, absM.boot$Long))
dist0m<- dist(absM.boot$Thorax)
dist1Rm<- dist(residuals(mod2))

mantel.rtest(dist1m, dist0m)
mantel.rtest(dist1m, dist1Rm)


#@ ran once with whole data set and once with the subset for TF. 

#=========================================
