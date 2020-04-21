setwd("~/PGY_crossvalidationchecking")
#

library(raster)
library(sp)
library(hydroGOF)
library(magrittr)


#Loading cross validation matrix
dat <- read.csv("G:\\My Drive\\PARAGUAY_COS\\Paraguay_Matriz_Validacion.csv")

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
  temp <- cellStats(RK,stats[i]) %>% as.numeric() %>% round(digits = 2)
  df <- c(df,temp)
}
df

#RF model
df1 <- c()
for(i in 1:length(stats)){
  temp <- cellStats(RF,stats[i]) %>% as.numeric() %>% round(digits = 2)
  df1 <- c(df1,temp)
}
df1


#DIFF
df2 <- c()
for(i in 1:length(stats)){
  temp <- cellStats(DIFF,stats[i]) %>% as.numeric() %>% round(digits = 2)
  df2 <- c(df2,temp)
}
df2

modStats <- rbind(df,df1,df2) 
row.names(modStats) <- c("RK","RF","DIFF")
colnames(modStats) <- stats
modStats

#           sum mean   min  max   sd skew
#RK   2280158.11 4.52  1.95 8.41 0.71 1.30
#RF   2286478.49 4.52  1.40 8.32 0.63 1.42
#DIFF     -88.04 0.00 -2.66 3.04 0.42 0.47

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
