setwd("~/PGY_crossvalidationchecking")

library(raster)
library(sp)
library(hydroGOF)
library(magrittr)
library(moments)

#Loading cross validation matrix
dat <- read.csv("G:\\My Drive\\PARAGUAY_COS\\Paraguay_Matriz_Validacion.csv")

#Loading regression matrix
datrm <- read.csv("G:\\My Drive\\PARAGUAY_COS\\Paraguay_Matriz_Regresion.csv")
summary(datrm)
#Summarizing observed SOC stock
names(datrm)

descr <- function(columna){
  MAX <- max(columna,na.rm=T)
  MIN <- min(columna,na.rm=T)
  PROM <- mean(columna,na.rm=T)
  DESVEST <- sd(columna,na.rm=T)
  CV <- DESVEST/PROM*100
  CURT <- kurtosis(columna,na.rm=T)
  SKEW <- skewness(columna,na.rm=T)
  return(data.frame(MAX, MIN, PROM,DESVEST,CV, CURT, SKEW))
}
descr(datrm$COSkgm2) %>% round(digits=2)

#             MAX   MIN  PROM DESVEST     CV   CURT  SKEW
# COSkgm2  21.976 0.423 4.887   2.175 44.499 11.037 2.112

#Loading final raster models and unit conversion to kg.m-2
RK <- raster("G:\\My Drive\\PARAGUAY_COS\\PRY_Mapa_COS_tnha_geo.tif")*0.1
RF <- raster("G:\\My Drive\\PARAGUAY_COS\\PRY1_COS_Rf_tnha_2a (1).tif")*0.1
RF <- resample(RF, RK, methods = 'ngb')
DIFF <- RK-RF

par(mfrow=c(1,3))
plot(RK,main="RK model")
plot(RF,main="RK model")
plot(DIFF,main="Difference")

#Cell stats for RK and RF
stats <- c("sum", "mean","min","max","sd","skew")

#RK model
df <- c()
for(i in 1:length(stats)){
  temp <- cellStats(RK,stats[i]) %>% as.numeric() %>% round(digits = 3)
  df <- c(df,temp)
}
df

#RF model
df1 <- c()
for(i in 1:length(stats)){
  temp <- cellStats(RF,stats[i]) %>% as.numeric() %>% round(digits = 3)
  df1 <- c(df1,temp)
}
df1


#DIFF
df2 <- c()
for(i in 1:length(stats)){
  temp <- cellStats(DIFF,stats[i]) %>% as.numeric() %>% round(digits = 3)
  df2 <- c(df2,temp)
}
df2

modStats <- rbind(df,df1,df2) 
row.names(modStats) <- c("RK","RF","DIFF")
colnames(modStats) <- stats
modStats

              #sum  mean    min   max    sd  skew
# RK   2280158.115 4.516  1.948 8.412 0.711 1.304
# RF   2286478.487 4.518  1.403 8.316 0.633 1.418
# DIFF     -88.038 0.000 -2.663 3.036 0.420 0.474

#Stack RK, RF and raster of differences between the two models
mod <- stack(RK,RF,DIFF)
names(mod) <- c("RK","RF","DIFF")

# Converting dat to spatial points data frame,  and extract
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
