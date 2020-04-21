setwd("~/PGY_crossvalidationchecking")

#Descriptive statistics for observed SOC stock-RK-RF and external validation

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
RK <- raster("G:\\My Drive\\PARAGUAY_COS\\PRY_Mapa_COS_tnha_geo.tif")
RF <- raster("G:\\My Drive\\PARAGUAY_COS\\PRY1_COS_Rf_tnha_2a (1).tif")
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
# MIN   Q1  MED   Q3  MAX PROM DESVEST   CV CURT SKEW
#1.19 1.49 1.52 1.54 1.61 1.51    0.05 3.59 7.01 -1.8
descr(datrm$COS_0.30)
# MIN  Q1 MED    Q3   MAX  PROM DESVEST    CV  CURT SKEW
#0.88 7.8 9.4 12.65 61.55 11.01    5.63 51.12 15.23  2.6
descr(datrm$COSkgm2) 
# MIN  Q1  MED   Q3   MAX PROM DESVEST   CV  CURT SKEW
#0.42 3.6 4.29 5.65 21.98 4.89    2.17 44.5 11.04 2.11
descr(values(RK)*0.1)
# MIN   Q1  MED   Q3  MAX PROM DESVEST    CV CURT SKEW
#1.95 4.05 4.32 4.77 8.41 4.52    0.71 15.74 5.21  1.3
descr(values(RF)*0.1)
# MIN   Q1  MED   Q3  MAX PROM DESVEST    CV CURT SKEW
#1.4 4.14 4.36 4.68 8.32 4.52    0.63 14.01 5.69 1.42

#####Checking external validation#####

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

#Validation measurements
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



#####Updating plot and RF histogram#####
 
x11()
par(mar=c(5,5,5,0))
par(mgp=c(3,1,0))
hist(dat$COSkgm2*10, breaks = 50, main  = "Histogram of SOC stock measured",
     ylim=c(0,120),las=1,xlim=c(0,250),
     xlab= expression("SOC stock" ~ (t~ha^{-1}))) 
     
x11()
par(mar=c(5,5,5,0))
par(mgp=c(4,1,0))
hist(RK, breaks = 100, main  = "Histogram of SOC stock estimated values by RK",
     ylim=c(0,6000),las=1,
     xlab= expression("SOC stock" ~ (t~ha^{-1})))
     
x11()
par(mar=c(5,5,5,0))
par(mgp=c(4,1,0))
hist(RF, breaks = 100, main  = "Histogram of SOC stock estimated values by RF",
     ylim=c(0,35000),las=1,
     xlab= expression("SOC stock" ~ (t~ha^{-1})))


summary(dat$COSkgm2*10)






