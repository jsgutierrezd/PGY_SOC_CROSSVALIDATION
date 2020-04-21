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
summary(dat$RF)

ecv <- function(sim,pred){
  ME <- me(sim,pred)#mean error
  MAE <- mae(sim,pred)#mean absolute error
  MSE <- mse(sim,pred)#mean square error
  RMSE <- rmse(sim,pred)#root mean square error
  Q3 <- abs(sim - pred) %>% quantile(0.75, na.rm=TRUE)#Third quartile
  AVE <- 1 - sum((sim-pred)^2, na.rm=TRUE) / 
    sum( (pred - mean(pred, na.rm = TRUE))^2, 
         na.rm = TRUE)##Amount of variance explained
    R2 <- (cor(sim,pred)^2)#R2
  return(data.frame(ME,MAE,MSE,RMSE,Q3,R2,AVE))%>% round(digits=2)
}

ecv(dat$RK,dat$COSkgm2)
#   ME  MAE  MSE RMSE   Q3   R2 AVE
#-0.45 1.47 5.92 2.43 1.73 0.13 0.1


ecv(dat$RF,dat$COSkgm2)
#   ME  MAE  MSE  MSE   Q3   R2  AVE
#-0.45 1.46 5.90 2.43 1.68 0.14 0.11



