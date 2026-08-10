###Predict from models onto spatial grids

set.seed <- 787

##Load packages from https://github.com/TobieSurette 

 library(gulf.data)
 library(gulf.stats)
library(gulf.graphics)
library(gulf.spatial)
#library(splines)
#library(MASS) ##for glm.nb function
library(sdmTMB)
library(sdmTMBextra)
library(dplyr)
library(ggplot2)
library(assertthat)
#library(rnaturalearthhires)
library(INLA)
library(sf)
library(reshape2)

##Load gulf package (DFO package to access and manipulate data)
library(gulf)
#library(DHARMa)

##this should get rid of scientific notation
options(scipen = 999)

## Set paths
main.fp <- getwd()

fig.path <- file.path(main.fp, "results/figures")
data.path <- file.path(main.fp, "data")
results.path <- file.path(main.fp, "results/tables")

path <- paste0(main.fp, "/programs/")
rdatapath <- paste0(main.fp, "/programs/RData/")

##needed for grid
source(paste0(path, "make.grid.GULF.R"))
source(paste0(path, "temperature_near.neigh.R"))
source(paste0(path, "get.bottom.temp.R"))
source(paste0(path, "get.substrate.grid.R"))

##set up a summary table
columns <- c("survey", "mod.survey", "mod.all" )
sum.table <- as.data.frame( matrix(data=NA, nrow =3, ncol=length(columns), byrow=TRUE))
colnames(sum.table )<- columns

##list surveys
sum.table$survey <- c("NS", "RV", "MI")

##define best model for each survey
sum.table$mod.survey [sum.table$survey =="NS"] <- "m15"
sum.table$mod.survey [sum.table$survey =="RV"] <- "m6"
sum.table$mod.survey [sum.table$survey =="MI"] <- "m7"

sum.table$mod.all <- "m13"

mod.label<- as.character(sum.table$mod.all[1])

##load general values from integrated dataset 

load( file = paste0(rdatapath, "sGSL_", "GeneralValues",   ".rdata")) 
label <- "commercial biomass"

##load models
##NS
load(file = paste0(rdatapath, "sGSL_",  sum.table$mod.survey [sum.table$survey =="NS"] , "_", label,  "_NS.rdata")) ##load model
modNS <- m

##RV
load(file = paste0(rdatapath, "sGSL_",  sum.table$mod.survey [sum.table$survey =="RV"], "_", label,  "_RV.rdata")) ##load model
modRV <- m

##MI
load(file = paste0(rdatapath, "sGSL_",  sum.table$mod.survey [sum.table$survey =="MI"], "_", label,  "_MI.rdata")) ##load model
modMI <- m

##ALL
load(file = paste0(rdatapath, "sGSL_",  sum.table$mod.all[1], "_", label,  "_ALL.rdata")) ##load model
modALL <- m

##labels for file names #3reset after loading in data files as they include other info. 
label <- "commercial biomass"

##load in base map
load(file = paste0(data.path, "/mapdata.rdata"))

##bring in data used in individual survey analyses and append to identify as that survey

load(paste0(rdatapath, "/sGSL_commercial biomass_NS.rdata"))

data_NS <- data
max.depth_NS <- max.depth
min.depth_NS <- min.depth


load(paste0(rdatapath, "/sGSL_commercial biomass_RV.rdata"))
data_RV <- data
max.depth_RV <- max.depth
min.depth_RV <- min.depth


load(paste0(rdatapath, "/sGSL_commercial biomass_MI.rdata"))
data_MI <- data
max.depth_MI <- max.depth
min.depth_MI <- min.depth


load(paste0(rdatapath, "/sGSL_commercial biomass_ALL.rdata"))
data_ALL <- data

xlimALL <- xlim
ylimALL <- ylim


## set number of simulations for SE

nsim <- 1000 

surveys <- c("NS", "RV", "MI")

##ns Survey

{
model.run <- "NS"
  
  width <- if (model.run == "NS") {
    15
  }else {
    12
  }
  
  height  <- if (model.run == "NS")  {
    11
  }else {
    11
  }
  

  res <- if (model.run == "MI") {
    0.5 
  }else{
    1
  } 

max.depth <- max.depth_NS

min.depth <-  min.depth_NS


pred.survey <- "NS_Fall" 
label <- "commercial biomass_NS"

grid <- make.grid.GULF (res = res, xlim = xlimNS , ylim = ylimNS,  mean.depth = mean.depth, sd.depth = sd.depth, mean.temp=mean.temp, sd.temp=sd.temp,
                       year=years, temp.month="september")


grid$survey <- pred.survey ##defined at top

##map out grid 

x <- grid [grid$ix==TRUE  & grid$ns.study.area,]

png(file = paste0(fig.path, "/sGSL_", "NS", "_grid_", res, "Kres_pts.png"), res = 500, units = "in", height = height, width = width)
plot(xlimNS,ylimNS, type = "n", xlab = "", ylab = "", 
     xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n")
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "white")
points(x$X, x$Y, pch=1,  cex=0.05, col="blue")
plot(map_data[1],xlim=xlimNS, ylim=ylimNS, col= "grey75", border="grey15" , main="", add=TRUE)
dev.off()


tmp <- NULL
tmp$x <- c(519, 519, 519+res, 519+res)
tmp$y <-c(5177, 5177+res, 5177+res, 5177) ##determines boundaries of one grid square in km in center of study area
area <- 1000000 * gulf.graphics::area( as.polygon(tmp$x, tmp$y)) #calculates the area of one grid square in m2 


grid <- as.data.frame(grid) 

results <- NULL
index.modNS <- NULL
index.modALL <- NULL

#define years
yr.figure <- if (model.run == "MI") {
  c(2001:2022) 
}else{
  c(2001: 2025) 
} 

for (i in min(yr.figure) :max(yr.figure)){
  print(i)
  
  grid.yr <- grid[grid$year==i,]
  
  ##prediction from NS model
  p <- predict(modNS, newdata = grid.yr[grid.yr$ix & grid.yr$ns.study.area, ], return_tmb_object = TRUE)
  ind.modNS <- get_index(p, area = area, bias_correct = TRUE)
  index.modNS <- rbind(ind.modNS, index.modNS)
  
  ##redo prediction without return_tmb_object
  p <- predict(modNS, newdata = grid.yr[grid.yr$ix & grid.yr$ns.study.area, ])
  grid.yr$pred_modNS [grid.yr$ix & grid.yr$ns.study.area] <- exp(p$est) * 1000000
  
  pp <- predict(modNS, newdata = grid.yr[grid.yr$ix & grid.yr$ns.study.area,], nsim = nsim)
  pp2 <- 1000000 * exp(pp) ##exponentiate to be working in real values * 1000 000 to get results in km2
  pp3 <- apply(pp2, 1, sd) ##get sd of nsim simulations
  grid.yr$pred_SE_modNS[grid.yr$ix & grid.yr$ns.study.area] <- pp3  ## assign results to points to then map
  
  grid.yr$pred_CV_modNS <-   grid.yr$pred_SE_modNS / grid.yr$pred_modNS
  
  ##prediction from integrated model
  
  p <- predict(modALL, newdata = grid.yr[grid.yr$ix & grid.yr$ns.study.area, ], return_tmb_object = TRUE)
  ind.modALL <- get_index(p, area = area, bias_correct = TRUE)
  index.modALL <- rbind(ind.modALL, index.modALL)
  
  ##redo prediction without return_tmb_object
  p <- predict(modALL, newdata = grid.yr[grid.yr$ix & grid.yr$ns.study.area, ], )
  grid.yr$pred_modALL [grid.yr$ix & grid.yr$ns.study.area] <- exp(p$est) * 1000000
  
  pp <- predict(modALL, newdata = grid.yr[grid.yr$ix & grid.yr$ns.study.area,], nsim = nsim)
  pp2 <- 1000000 * exp(pp) ##exponentiate to be working in real values * 1000 000 to get results in km2
  pp3 <- apply(pp2, 1, sd) ##get sd of nsim simulations
  grid.yr$pred_SE_modALL[grid.yr$ix & grid.yr$ns.study.area] <- pp3  ## assign results to points to then map
  grid.yr$pred_CV_modALL <-   grid.yr$pred_SE_modALL / grid.yr$pred_modALL

  grid.yr$SE_ratio <- grid.yr$pred_SE_modALL / grid.yr$pred_SE_modNS
  grid.yr$CV_ratio <- grid.yr$pred_CV_modALL / grid.yr$pred_CV_modNS
  
  results <- rbind(results, grid.yr)
  
}


vars <- c("est", "lwr", "upr")
log.vars <- c("log_est", "se")

index.modNS[vars] <- index.modNS[vars] / 1000000 ##divide by 1000000 to be in1000 s of tons
index.modALL[vars] <- index.modALL[vars] / 1000000

colnames(index.modNS) <- paste0(colnames(index.modNS), "_modNS")
colnames(index.modALL) <- paste0(colnames(index.modALL), "_modALL")

res.index.NS <- cbind(index.modNS, index.modALL)
res.index.NS$year <- res.index.NS$year_modNS
res.index.NS <- res.index.NS [,c("year",  "est_modNS"   ,   "lwr_modNS"    ,  "upr_modNS"    ,  "log_est_modNS" , "se_modNS" ,
                                 "est_modALL"  ,   "lwr_modALL"  ,   "upr_modALL"   ,  "log_est_modALL" ,"se_modALL" )]


results <- results[,!names(results) %in% c("geometry")]

write.csv(results, file=paste0(results.path,"/sGSL_",   min(years), "_", max(years),   "_", res, "Kres_", model.run, "_Predictions_", mod.label, "ALL.csv"), row.names = FALSE)
write.csv(res.index.NS, file=paste0(results.path,"/sGSL_",   min(years), "_", max(years),   "_", res, "Kres_", model.run, "_Index_", mod.label, "ALL.csv"), row.names = FALSE)

}

##RV Survey
{
  model.run <- "RV"
  
  width <- if (model.run == "NS") {
    15
  }else {
    12
  }
  
  height  <- if (model.run == "NS")  {
    11
  }else {
    11
  }
  
  res <- if (model.run == "MI") {
    0.5 
  }else{
    1
  } 
  
  max.depth <- max.depth_RV
  
  min.depth <-  min.depth_RV
  
  
  pred.survey <- "RV" 
  label <- "commercial biomass_RV"
  
  #grid <- make.grid.GULF (res = res, xlim = xlimRV , ylim = ylimRV,  mean.depth = mean.depth, sd.depth = sd.depth, mean.temp=mean.temp, sd.temp=sd.temp,
                     #     year=years, temp.month="september")
  
  #grid$survey <- pred.survey ##defined at top
  
  #save(grid, file = paste0(path, "sGSL_",label, "_grid_" ,  res, "Kres.rdata"))
  
  load(file = paste0(path, "sGSL_",label, "_grid_" ,  res, "Kres.rdata"))
  

    ##map out grid 
  
  x <- grid [grid$ix==TRUE  & grid$rv.study.area,]
  
  png(file = paste0(fig.path, "/sGSL_", "RV", "_grid_", res, "Kres_pts.png"), res = 500, units = "in", height = height, width = width)
  plot(xlimRV,ylimRV, type = "n", xlab = "", ylab = "", 
       xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n")
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "white")
  points(x$X, x$Y, pch=1,  cex=0.05, col="blue")
  plot(map_data[1],xlim=xlimRV, ylim=ylimRV, col= "grey75", border="grey15" , main="", add=TRUE)
  dev.off()
  
  
  tmp <- NULL
  tmp$x <- c(519, 519, 519+res, 519+res)
  tmp$y <-c(5177, 5177+res, 5177+res, 5177) ##determines boundaries of one grid square in km in center of study area. location more important when working in lats and longs
  area <- 1000000 * gulf.graphics::area( as.polygon(tmp$x, tmp$y)) #calculates the area of one grid square in m2 - with a res =1 area is 1 km2
  
  
  grid <- as.data.frame(grid) ##otherwise I get an error
  
  results <- NULL
  index.modRV <- NULL
  index.modALL <- NULL
  
  #define years
  yr.figure <- if (model.run == "MI") {
    c(2001:2022) 
  }else{
    c(2001: 2025) 
  } 
  
  for (i in min(yr.figure) :max(yr.figure)){
    print(i)
    
    grid.yr <- grid[grid$year==i,]
    
    ##prediction from RV model
    p <- predict(modRV, newdata = grid.yr[grid.yr$ix & grid.yr$rv.study.area, ], return_tmb_object = TRUE)
    ind.modRV <- get_index(p, area = area, bias_correct = TRUE)
    index.modRV <- rbind(ind.modRV, index.modRV)
    
    ##redo prediction without return_tmb_object
    p <- predict(modRV, newdata = grid.yr[grid.yr$ix & grid.yr$rv.study.area, ])
    grid.yr$pred_modRV [grid.yr$ix & grid.yr$rv.study.area] <- exp(p$est) * 1000000
    
    pp <- predict(modRV, newdata = grid.yr[grid.yr$ix & grid.yr$rv.study.area,], nsim = nsim)
    pp2 <- 1000000 * exp(pp) ##exponentiate to be working in real values * 1000 000 to get results in km2
    pp3 <- apply(pp2, 1, sd) ##get sd of nsim simulations
    grid.yr$pred_SE_modRV[grid.yr$ix & grid.yr$rv.study.area] <- pp3  ## assign results to points to then map
    
    grid.yr$pred_CV_modRV <-   grid.yr$pred_SE_modRV / grid.yr$pred_modRV
    
    ##prediction from integrated model
    
    p <- predict(modALL, newdata = grid.yr[grid.yr$ix & grid.yr$rv.study.area, ], return_tmb_object = TRUE)
    ind.modALL <- get_index(p, area = area, bias_correct = TRUE)
    index.modALL <- rbind(ind.modALL, index.modALL)
    
    ##redo prediction without return_tmb_object
    p <- predict(modALL, newdata = grid.yr[grid.yr$ix & grid.yr$rv.study.area, ], )
    grid.yr$pred_modALL [grid.yr$ix & grid.yr$rv.study.area] <- exp(p$est) * 1000000
    
    pp <- predict(modALL, newdata = grid.yr[grid.yr$ix & grid.yr$rv.study.area,], nsim = nsim)
    pp2 <- 1000000 * exp(pp) ##exponentiate to be working in real values * 1000 000 to get results in km2
    pp3 <- apply(pp2, 1, sd) ##get sd of nsim simulations
    grid.yr$pred_SE_modALL[grid.yr$ix & grid.yr$rv.study.area] <- pp3  ## assign results to points to then map
    grid.yr$pred_CV_modALL <-   grid.yr$pred_SE_modALL / grid.yr$pred_modALL
    
    grid.yr$SE_ratio <- grid.yr$pred_SE_modALL / grid.yr$pred_SE_modRV
    grid.yr$CV_ratio <- grid.yr$pred_CV_modALL / grid.yr$pred_CV_modRV
    
    results <- rbind(results, grid.yr)
    
  }
  
  
  vars <- c("est", "lwr", "upr")
  log.vars <- c("log_est", "se")
  
  index.modRV[vars] <- index.modRV[vars] / 1000000 ##divide by 1000000 to be in1000 s of tons
  index.modALL[vars] <- index.modALL[vars] / 1000000
  
  colnames(index.modRV) <- paste0(colnames(index.modRV), "_modRV")
  colnames(index.modALL) <- paste0(colnames(index.modALL), "_modALL")
  
  res.index.RV <- cbind(index.modRV, index.modALL)
  res.index.RV$year <- res.index.RV$year_modRV
  res.index.RV <- res.index.RV [,c("year",  "est_modRV"   ,   "lwr_modRV"    ,  "upr_modRV"    ,  "log_est_modRV" , "se_modRV" ,
                                   "est_modALL"  ,   "lwr_modALL"  ,   "upr_modALL"   ,  "log_est_modALL" ,"se_modALL" )]
  
  
  results <- results[,!names(results) %in% c("geometry")]
  
  write.csv(results, file=paste0(results.path,"/sGSL_",   min(years), "_", max(years),   "_", res, "Kres_", model.run, "_Predictions_", mod.label, "ALL.csv"), row.names = FALSE)
  write.csv(res.index.RV, file=paste0(results.path,"/sGSL_",   min(years), "_", max(years),   "_", res, "Kres_", model.run, "_Index_", mod.label, "ALL.csv"), row.names = FALSE)
  
}


##MI Survey

{
  model.run <- "MI"
  
  width <- if (model.run == "NS") {
    15
  }else {
    12
  }
  
  height  <- if (model.run == "NS")  {
    11
  }else {
    11
  }
  
  res <- if (model.run == "MI") {
    0.5 
  }else{
    1
  } 
  
  max.depth <- max.depth_MI
  
  min.depth <-  min.depth_MI
  
  
  pred.survey <- "MI" 
  label <- "commercial biomass_MI"
  
  grid <- make.grid.GULF (res = res, xlim = xlimMI , ylim = ylimMI,  mean.depth = mean.depth, sd.depth = sd.depth, mean.temp=mean.temp, sd.temp=sd.temp,
                          year=years, temp.month="september")
  
  grid$survey <- pred.survey ##defined at top
  
  grid$substrate <- factor(
    as.character(grid$substrate),
    levels = levels(factor(modALL$data$substrate))
  )
  
  ##map out grid 
  
  x <- grid [grid$ix==TRUE  & grid$mi.study.area,]
  
  png(file = paste0(fig.path, "/sGSL_", "MI", "_grid_", res, "Kres_pts.png"), res = 500, units = "in", height = height, width = width)
  plot(xlimMI,ylimMI, type = "n", xlab = "", ylab = "", 
       xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n")
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "white")
  points(x$X, x$Y, pch=1,  cex=0.05, col="blue")
  plot(map_data[1],xlim=xlimMI, ylim=ylimMI, col= "grey75", border="grey15" , main="", add=TRUE)
  dev.off()
  
  
  
  tmp <- NULL
  tmp$x <- c(519, 519, 519+res, 519+res)
  tmp$y <-c(5177, 5177+res, 5177+res, 5177) ##determines boundaries of one grid square in km in center of study area. location more important when working in lats and longs
  area <- 1000000 * gulf.graphics::area( as.polygon(tmp$x, tmp$y)) #calculates the area of one grid square in m2 - with a res =1 area is 1 km2
  
  
  grid <- as.data.frame(grid) ##otherwise I get an error
  
  results <- NULL
  index.modMI <- NULL
  index.modALL <- NULL
  
  #define years
  yr.figure <- if (model.run == "MI") {
    c(2001:2022) 
  }else{
    c(2001: 2025) 
  } 
  
  index.modALL <- NULL
  
  for (i in min(yr.figure) :max(yr.figure)){
    print(i)
    
    grid.yr <- grid[grid$year==i,]
    
    ##prediction from MI model
    p <- predict(modMI, newdata = grid.yr[grid.yr$ix & grid.yr$mi.study.area, ], return_tmb_object = TRUE)
    ind.modMI <- get_index(p, area = area, bias_correct = TRUE)
    index.modMI <- rbind(ind.modMI, index.modMI)
    
    ##redo prediction without return_tmb_object
    p <- predict(modMI, newdata = grid.yr[grid.yr$ix & grid.yr$mi.study.area, ])
    grid.yr$pred_modMI [grid.yr$ix & grid.yr$mi.study.area] <- exp(p$est) * 1000000
    
    pp <- predict(modMI, newdata = grid.yr[grid.yr$ix & grid.yr$mi.study.area,], nsim = nsim)
    pp2 <- 1000000 * exp(pp) ##exponentiate to be working in real values * 1000 000 to get results in km2
    pp3 <- apply(pp2, 1, sd) ##get sd of nsim simulations
    grid.yr$pred_SE_modMI[grid.yr$ix & grid.yr$mi.study.area] <- pp3  ## assign results to points to then map
    
    grid.yr$pred_CV_modMI <-   grid.yr$pred_SE_modMI / grid.yr$pred_modMI
    
    ##prediction from integrated model
    
    p <- predict(modALL, newdata = grid.yr[grid.yr$ix & grid.yr$mi.study.area, ], return_tmb_object = TRUE)
    ind.modALL <- get_index(p, area = area, bias_correct = TRUE)
    index.modALL <- rbind(ind.modALL, index.modALL)
    
    ##redo prediction without return_tmb_object
    p <- predict(modALL, newdata = grid.yr[grid.yr$ix & grid.yr$mi.study.area, ], )
    grid.yr$pred_modALL [grid.yr$ix & grid.yr$mi.study.area] <- exp(p$est) * 1000000
    
    pp <- predict(modALL, newdata = grid.yr[grid.yr$ix & grid.yr$mi.study.area,], nsim = nsim)
    pp2 <- 1000000 * exp(pp) ##exponentiate to be working in real values * 1000 000 to get results in km2
    pp3 <- apply(pp2, 1, sd) ##get sd of nsim simulations
    grid.yr$pred_SE_modALL[grid.yr$ix & grid.yr$mi.study.area] <- pp3  ## assign results to points to then map
    grid.yr$pred_CV_modALL <-   grid.yr$pred_SE_modALL / grid.yr$pred_modALL
    
    grid.yr$SE_ratio <- grid.yr$pred_SE_modALL / grid.yr$pred_SE_modMI
    grid.yr$CV_ratio <- grid.yr$pred_CV_modALL / grid.yr$pred_CV_modMI
    
    results <- rbind(results, grid.yr)
    
  }
  
  
  vars <- c("est", "lwr", "upr")
  log.vars <- c("log_est", "se")
  
  index.modMI[vars] <- index.modMI[vars] / 1000000 ##divide by 1000000 to be in1000 s of tons
  index.modALL[vars] <- index.modALL[vars] / 1000000
  
  colnames(index.modMI) <- paste0(colnames(index.modMI), "_modMI")
  colnames(index.modALL) <- paste0(colnames(index.modALL), "_modALL")
  
  res.index.MI <- cbind(index.modMI, index.modALL)
  res.index.MI$year <- res.index.MI$year_modMI
  res.index.MI <- res.index.MI [,c("year",  "est_modMI"   ,   "lwr_modMI"    ,  "upr_modMI"    ,  "log_est_modMI" , "se_modMI" ,
                                   "est_modALL"  ,   "lwr_modALL"  ,   "upr_modALL"   ,  "log_est_modALL" ,"se_modALL" )]
  
  
  results <- results[,!names(results) %in% c("geometry")]
  
  write.csv(results, file=paste0(results.path,"/sGSL_",   min(years), "_", max(years),   "_", res, "Kres_", model.run, "_Predictions_", mod.label, "ALL.csv"), row.names = FALSE)
  write.csv(res.index.MI, file=paste0(results.path,"/sGSL_",   min(years), "_", max(years),   "_", res, "Kres_", model.run, "_Index_", mod.label, "ALL.csv"), row.names = FALSE)
  
}

##ALL prediction, no ratios


{
  model.run <- "ALL"
  
  width <- if (model.run == "NS") {
    15
  }else {
    12
  }
  
  height  <- if (model.run == "NS")  {
    11
  }else {
    11
  }
  
  xlim <-  xlimALL
  
  ylim <-ylimALL
  
  
  max.depth <- max.depth
  
  min.depth <-  min.depth
  
  res <- if (model.run == "MI") {
    0.5 
  }else{
    1
  } 
  
  pred.survey <- "RV" 
  label <- "commercial biomass_ALL"
  
  grid <- make.grid.GULF (res = res, xlim = xlim , ylim = ylim,  mean.depth = mean.depth, sd.depth = sd.depth, mean.temp=mean.temp, sd.temp=sd.temp,
                          year=years, temp.month="september")
  
  grid$survey <- pred.survey ##defined at top
  
  save(grid, file = paste0(path, "sGSL_",label, "_grid_" ,  res, "Kres.rdata"))
  
  #load(file = paste0(path, "sGSL_",label, "_grid_" ,  res, "Kres.rdata"))
  
  
  grid <- as.data.frame(grid) ##otherwise I get an error
  
  
  
  
  
  ##map out grid 
  
  x <- grid [grid$ix==TRUE &  (grid$rv.study.area | grid$ns.study.area | grid$mi.study.area),]
  
  png(file = paste0(fig.path, "/sGSL_", "ALL", "_grid_", res, "Kres_pts.png"), res = 500, units = "in", height = height, width = width)
  plot(xlim,ylim, type = "n", xlab = "", ylab = "", 
       xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n")
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "white")
  points(x$X, x$Y, pch=1,  cex=0.05, col="blue")
  plot(map_data[1],xlim=xlim, ylim=ylim, col= "grey75", border="grey15" , main="", add=TRUE)
  dev.off()
  
  tmp <- NULL
  tmp$x <- c(519, 519, 519+res, 519+res)
  tmp$y <-c(5177, 5177+res, 5177+res, 5177) ##determines boundaries of one grid square in km in center of study area. location more important when working in lats and longs
  area <- 1000000 * gulf.graphics::area( as.polygon(tmp$x, tmp$y)) #calculates the area of one grid square in m2 - with a res =1 area is 1 km2
  
  
  results <- NULL
  index.modALL <- NULL
  
  
  #define years
  yr.figure <- if (model.run == "MI") {
    c(2001:2022) 
  }else{
    c(2001: 2025) 
  } 
  

  
  for (i in min(yr.figure) :max(yr.figure)){
    print(i)
    
    grid.yr <- grid[grid$year==i,]
    
    ##prediction from integrated model
    
    p <- predict(modALL, newdata = grid.yr[grid.yr$ix & (grid.yr$rv.study.area | grid.yr$ns.study.area | grid.yr$mi.study.area), ], return_tmb_object = TRUE)
    ind.modALL <- get_index(p, area = area, bias_correct = TRUE)
    index.modALL <- rbind(ind.modALL, index.modALL)
    
    ##redo prediction without return_tmb_object
    
    p <- predict(modALL, newdata = grid.yr[grid.yr$ix & (grid.yr$rv.study.area | grid.yr$ns.study.area | grid.yr$mi.study.area), ], )
    
    
    grid.yr$pred_modALL [grid.yr$ix & (grid.yr$rv.study.area | grid.yr$ns.study.area | grid.yr$mi.study.area)] <- exp(p$est) * 1000000
    
    pp <- predict(modALL, newdata = grid.yr[grid.yr$ix & (grid.yr$rv.study.area | grid.yr$ns.study.area | grid.yr$mi.study.area), ], nsim = nsim)
    pp2 <- 1000000 * exp(pp) ##exponentiate to be working in real values * 1000 000 to get results in km2
    pp3 <- apply(pp2, 1, sd) ##get sd of nsim simulations
    grid.yr$pred_SE_modALL[grid.yr$ix & (grid.yr$rv.study.area | grid.yr$ns.study.area | grid.yr$mi.study.area)] <- pp3  ## assign results to points to then map
    grid.yr$pred_CV_modALL <-   grid.yr$pred_SE_modALL / grid.yr$pred_modALL
    
    results <- rbind(results, grid.yr)
    
  }
  
  vars <- c("est", "lwr", "upr")
  log.vars <- c("log_est", "se")
  

  index.modALL[vars] <- index.modALL[vars] / 1000000 ##divide by 1000000 to be in1000 s of tons
  
 
  colnames(index.modALL) <- paste0(colnames(index.modALL), "_modALL")
  
  res.index.ALL <- index.modALL
  res.index.ALL$year <- res.index.ALL$year_modALL
  res.index.ALL <- res.index.ALL [,c("year",  
                                   "est_modALL"  ,   "lwr_modALL"  ,   "upr_modALL"   ,  "log_est_modALL" ,"se_modALL" )]
  
  results <- results[,!names(results) %in% c("geometry")]
  
  write.csv(results, file=paste0(results.path,"/sGSL_",   min(years), "_", max(years),   "_", res, "Kres_", model.run, "_Predictions_", mod.label, "ALL.csv"), row.names = FALSE)
  write.csv(res.index.ALL, file=paste0(results.path,"/sGSL_",   min(years), "_", max(years),   "_", res, "Kres_", model.run, "_Index_", mod.label, "ALL.csv"), row.names = FALSE)
  
  ### predict for sgsl LFAs


  p <- list()
  index <- list()

  p[[1]] <- predict(modALL, newdata = grid[which(grid$ix &  (grid$rv.study.area | grid$ns.study.area | grid$mi.study.area) & (grid$lfa %in% c( "23",  "24", "25", "26A", "26B"))), ], return_tmb_object = TRUE)
  index[[1]] <- get_index(p[[1]], area = area, bias_correct = TRUE)

 
  names(index) <- c("total" )

  tab <- NULL
  for (i in 1:length(index)){
    tmp <- index[[i]]
    tmp$area <- names(index)[i]
    tab <- rbind(tab, tmp)
  }
  vars <- c("est", "lwr", "upr")
  log.vars <- c("log_est", "se")
  tab <- tab[c("area", "year", vars, log.vars)]
  tab[vars] <- tab[vars] / 1000000

  # Export to file:

  write.csv(tab, file = paste0(results.path, "/sGSL_",  mod.label, "_", "CommercialBiomasssGSL_", " 1000 tons_", min(years), "_", max(years), "_",  res, "Kres.csv"), row.names = FALSE)


  
}

