setwd("~/PGY_crossvalidationchecking")

library(raster)
library(sp)
library(hydroGOF)
library(magrittr)
library(moments)

#####Loading cross validation matrix####
dat <- read.csv("G:\\My Drive\\PARAGUAY_COS\\Paraguay_Matriz_Validacion.csv")

#####Loading regression matrix#####
datrm <- read.csv("G:\\My Drive\\PARAGUAY_COS\\Paraguay_Matriz_Regresion.csv")
summary(datrm)
#####Loading final raster models and unit conversion to kg.m-2#####
RK <- raster("G:\\My Drive\\PARAGUAY_COS\\PRY_Mapa_COS_tnha_geo.tif")*0.1
RF <- raster("G:\\My Drive\\PARAGUAY_COS\\PRY1_COS_Rf_tnha_2a (1).tif")*0.1
RF <- resample(RF, RK, methods = 'ngb')
DIFF <- RK-RF
x11()
par(mfrow=c(1,3))
plot(RK,main="RK model")
plot(RF,main="RK model")
plot(DIFF,main="Difference")

#####Descriptive statistics for observed SOC stock-RK-RF#####
descr <- function(columna){
  MIN <- min(columna,na.rm=T)
  Q1 <- quantile(columna,0.25,na.rm=T)
  MED <- median(columna,na.rm=T)
  Q3 <- quantile(columna,0.75,na.rm=T)
  MAX <- max(columna,na.rm=T)
  PROM <- mean(columna,na.rm=T)
  DESVEST <- sd(columna,na.rm=T)
  CV <- DESVEST/PROM*100
  CURT <- kurtosis(columna,na.rm=T)
  SKEW <- skewness(columna,na.rm=T)
  return(data.frame(MIN,Q1,MED,Q3,MAX,PROM,DESVEST,CV,CURT,SKEW))%>% round(digits=2)
}
descr(datrm$DA_0.30)
descr(datrm$COS_0.30)
descr(datrm$COSkgm2) 
descr(values(RK))
descr(values(RF))


#####Checking external cross validation#####

#Stack RK, RF and raster of differences between the two models
mod <- stack(RK,RF,DIFF)
names(mod) <- c("RK","RF","DIFF")

# Converting dat to spatial points data frame, and extracting
dat_sp <- dat
names(dat)
coordinates(dat_sp) <- ~ LONGITUD + LATITUD

# Projecting spatial points data frame
dat_sp@proj4string <- CRS("+proj=longlat +datum=WGS84")
dat_sp <- spTransform (dat_sp, CRS=projection(mod))

#Extranting data from mod stack
dat <- cbind(dat, extract(mod, dat_sp))
summary(dat)
names(dat)

#COSkgm2
#RF
#RK
#DIFF

summary(dat$RK)
summary(dat$COSkgm2)
rmse(dat$RK,dat$COSkgm2)
mae(dat$RK,dat$COSkgm2)

summary(dat$RF)
rmse(dat$RF,dat$COSkgm2)
cor(dat$RF,dat$COSkgm2)
mae(dat$RF,dat$COSkgm2)

summary(dat$DIFF)
