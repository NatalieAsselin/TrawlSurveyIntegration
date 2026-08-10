##Script to run models for MIS

##load libraries

library(sdmTMB)
library(sdmTMBextra)

##these are used to make the mesh
library(dplyr)
library(ggplot2)
library(assertthat)
library(rnaturalearthhires)
library(INLA)

library(sf)




## set paths
main.fp <- getwd()


fig.path <- file.path(main.fp, "results/figures")
data.path <- file.path(main.fp, "data")
results.path <- file.path(main.fp, "results/tables")

path <- paste0(main.fp, "/programs/")
rdatapath <- paste0(main.fp, "/programs/RData/")


##label for file names
label <- "commercial biomass_MI"

source(paste0(path, "/make.mesh.GULF.R"))
source(paste0(path, "/BIO_functions.r"))

#load dataset

load(paste0(data.path, "/MI_Lobster data.rdata"))

##load general values from complete dataset, used to scale depth and temp 

load( file = paste0(rdatapath, "sGSL_", "GeneralValues",   ".rdata"))


# Define dataset:
data <- setMI

##number of folds for cross validation, 90% of data to predict 10%
k_folds <- as.numeric(length(unique(data$fold_id)))

##scale the depth and temp by mean and sd from integrated dataset

data$depth.s <- (data$depth - mean.depth)/sd.depth
data$depth.s2 <- data$depth.s ^ 2

data$bottom.temp.s <- (data$bottom.temp - mean.temp)/sd.temp
data$bottom.temp.s2 <- data$bottom.temp.s ^ 2

max.depthMI <- max(data$depth) ##get max.depth to limit grid to places where I have data
min.depthMI <- min(data$depth)

yr.model <- c(2001:2022)
 
a <- which(data$year %in% yr.model)

data <- data [a,]

##identify years in dataset

years <- c(min(data$year):max(data$year))

xlimMI = c(580, 630)
ylimMI = c(5230, 5275)

##change name of offset

data$log.effort <- data$log.swept.area



# Build mesh:

cutoff <- 1

mesh <- make.barrier.mesh (data=data, xlim = xlimMI , ylim= ylimMI)

 png(file = paste0(fig.path, "/sGSL_", label, "_", cutoff, "cutoff_abundance mesh.png"), res = 500, units = "in", height = 7.0, width = 10)
 plot(mesh, cex = 0.8, pch = 21, bg = "grey60")
 dev.off()

##Models to run
models <- paste("m", c(1:3,5:9,11,15), sep="")


##model selection table
##set up a blank table

columns <- c("Model", "Formula", "Time-varying", "Family", "AIC", "cAIC", "rho", "Matern range", "Spatial SD", "Hessian_positive", "Sum loglik", "pbias_train", "pbias_test", "MAE_train","MAE_test", "RMSE_train", "RMSE_test" )
mod.select <- as.data.frame( matrix(data=NA, nrow =0, ncol=length(columns), byrow=TRUE))
colnames(mod.select) <- columns


##model selection table function

mod.select.fn <- function (){
  
  c<- as.data.frame( matrix(data=NA, nrow =1, ncol=length(columns), byrow=TRUE))
  colnames(c) <- columns
  c$Model <- mod.label
  c$Formula <-m$formula [1]
  c$"Time-varying" <- ifelse(is.null (m$time_varying), NA, paste(m$time_varying[1], m$time_varying[2])  )
  c$"Family" <- ifelse(m$family[1]=="tweedie", paste0(m$family[1], "(link = ", m$family[2], ")"), m$family["clean_name"])
  
  
  ##spatial model 
  c$AIC <- AIC (m)
  c$rho <- m$sd_report[[1]]["rho"]
  c$`Matern range` <- m$sd_report[[1]]["range"]
  c$`Spatial SD` <- m$sd_report[[1]]["sigma_O"]
  c$"Hessian_positive" <- m$pos_def_hessian
  
  ##model validation 
  
  c$"Sum loglik" <- m_cv$sum_loglik
  
  m_cvTT = sdmTMBcv_tntpreds (m_cv)
  
  fitTT = dplyr::bind_rows(m_cvTT)
  fitTT$sqR = fitTT$commercial.weight - fitTT$pred
  
  
  c$pbias_test<-  with(fitTT[fitTT$tt=='test',],pbias(as.numeric(commercial.weight ),as.numeric(pred)))
  c$pbias_train<-  with(fitTT[fitTT$tt=='train',],pbias(as.numeric(commercial.weight ),as.numeric(pred))) 
  c$MAE_test<-  with(fitTT[fitTT$tt=='test',],mae(as.numeric(commercial.weight ),as.numeric(pred)))
  c$MAE_train<-  with(fitTT[fitTT$tt=='train',],mae(as.numeric(commercial.weight ),as.numeric(pred)))
  c$RMSE_test <- with(fitTT[fitTT$tt=='test',],rmse(as.numeric(commercial.weight ),as.numeric(pred)))
  c$RMSE_train <- with(fitTT[fitTT$tt=='train',],rmse(as.numeric(commercial.weight ),as.numeric(pred)))
  
  return(c)
  
}
 
##Models
 
if ("m1" %in% models) {
  
  mod.label <- "m1" 
  
  ##model
  
  m <- sdmTMB(
    data=data,
    formula = commercial.weight ~ 0 ,
    mesh = mesh,
    time_varying = ~ 1 ,
    time_varying_type = "rw0",
    spatial = "on",
    family =  tweedie(link = "log"),
    time = "year",
    spatiotemporal = "ar1",
    offset = data$log.effort
  )
  
  ##model cross-validation
  
  m_cv <- sdmTMB_cv(
    data=data,
    formula = commercial.weight ~ 0 ,
    mesh = mesh, 
    time_varying = ~ 1 ,
    time_varying_type = "rw0",
    spatial = "on",
    family =  tweedie(link = "log"),
    time = "year",
    spatiotemporal = "ar1",
    offset = "log.effort",
    fold_ids = 'fold_id',
    k_folds = k_folds
  )
  
  c <-mod.select.fn()
  mod.select <- rbind(mod.select, c)
  
  summary(m)
  
  # Save model results:
  save(m, file = paste0(rdatapath, "sGSL_",  mod.label, "_", label,  ".rdata"))
  
  m1 <- m
  m1_cv <- m_cv
  
  rm(m, m_cv)        
  mod.select2 <- as.data.frame(lapply(mod.select, as.character), stringsAsFactors=FALSE)        
  write.csv(mod.select2, file = paste0(results.path, "/sGSL_",  label, "_Model Selection table_m1_", mod.label, "_", min(years), "_", max(years), "_", cutoff, "cutoff.csv"), row.names = FALSE)    
  
}



if ("m2" %in% models) {
  
  mod.label <- "m2" 
  
  m <- sdmTMB(  
    data=data,
    formula = commercial.weight ~ 0 + depth.s , 
    mesh = mesh,
    time_varying = ~ 1 ,
    time_varying_type = "rw0",
    spatial = "on",
    family =  tweedie(link = "log"),
    time = "year",
    spatiotemporal = "ar1",
    offset = data$log.effort
  )
  
  m_cv <- sdmTMB_cv(  
    data=data,
    formula = commercial.weight ~ 0 + depth.s ,  
    mesh = mesh,
    time_varying = ~ 1 , 
    time_varying_type = "rw0",
    spatial = "on",
    family =  tweedie(link = "log"),
    time = "year",
    spatiotemporal = "ar1",
    offset = "log.effort",
    fold_ids = 'fold_id',
    k_folds = k_folds
  )
  
  c <-mod.select.fn()
  mod.select <- rbind(mod.select, c)
  
  summary(m)
  
  # Save model results:
  save(m, file = paste0(rdatapath, "sGSL_",  mod.label, "_", label,  ".rdata"))
  
  m2 <- m
  m2_cv <- m_cv
  
  rm(m, m_cv)        
  mod.select2 <- as.data.frame(lapply(mod.select, as.character), stringsAsFactors=FALSE)        
  write.csv(mod.select2, file = paste0(results.path, "/sGSL_",  label, "_Model Selection table_m1_", mod.label, "_", min(years), "_", max(years), "_", cutoff, "cutoff.csv"), row.names = FALSE)    
  
}


if ("m3" %in% models) {
  
  mod.label <- "m3" 
  
  m <- sdmTMB(  
    data=data,
    formula = commercial.weight ~ 0 + bottom.temp.s , 
    mesh = mesh,
    time_varying = ~ 1 ,
    time_varying_type = "rw0",
    spatial = "on",
    family =  tweedie(link = "log"),
    time = "year",
    spatiotemporal = "ar1",
    offset = data$log.effort
  )
  
  m_cv <- sdmTMB_cv(  
    data=data,
    formula = commercial.weight ~ 0 + bottom.temp.s ,  
    mesh = mesh,
    time_varying = ~ 1 , 
    time_varying_type = "rw0",
    spatial = "on",
    family =  tweedie(link = "log"),
    time = "year",
    spatiotemporal = "ar1",
    offset = "log.effort",
    fold_ids = 'fold_id',
    k_folds = k_folds
  )
  
  c <-mod.select.fn()
  mod.select <- rbind(mod.select, c)
  
  summary(m)
  
  # Save model results:
  save(m, file = paste0(rdatapath, "sGSL_",  mod.label, "_", label,  ".rdata"))
  
  m3 <- m
  m3_cv <- m_cv
  
  rm(m, m_cv)        
  mod.select2 <- as.data.frame(lapply(mod.select, as.character), stringsAsFactors=FALSE)        
  write.csv(mod.select2, file = paste0(results.path, "/sGSL_",  label, "_Model Selection table_m1_", mod.label, "_", min(years), "_", max(years), "_", cutoff, "cutoff.csv"), row.names = FALSE)    
  
}


if ("m5" %in% models) {
  
  mod.label <- "m5" 
  
  m <- sdmTMB(  
    data=data,
    formula = commercial.weight ~ 0  + depth.s + depth.s2 , 
    mesh = mesh,
    time_varying = ~ 1 ,
    time_varying_type = "rw0",
    spatial = "on",
    family =  tweedie(link = "log"),
    time = "year",
    spatiotemporal = "ar1",
    offset = data$log.effort
  )
  
  m_cv <- sdmTMB_cv(  
    data=data,
    formula = commercial.weight ~ 0  + depth.s + depth.s2 , 
    mesh = mesh,
    time_varying = ~ 1 ,
    time_varying_type = "rw0",
    spatial = "on",
    family =  tweedie(link = "log"),
    time = "year",
    spatiotemporal = "ar1",
    offset = "log.effort",
    fold_ids = 'fold_id',
    k_folds = k_folds
  )
  
  c <-mod.select.fn()
  mod.select <- rbind(mod.select, c)
  
  summary(m)
  
  # Save model results:
  save(m, file = paste0(rdatapath, "sGSL_",  mod.label, "_", label,  ".rdata"))
  
  m5 <- m
  m5_cv <- m_cv
  
  rm(m, m_cv)       
  mod.select2 <- as.data.frame(lapply(mod.select, as.character), stringsAsFactors=FALSE)       
  write.csv(mod.select2, file = paste0(results.path, "/sGSL_",  label, "_Model Selection table_m1_", mod.label, "_", min(years), "_", max(years), "_", cutoff, "cutoff.csv"), row.names = FALSE)    
  
}


if ("m6" %in% models) {
  
  mod.label <- "m6" 
  
  m <- sdmTMB(  
    data=data,
    formula = commercial.weight ~ 0  + bottom.temp.s + bottom.temp.s2 , 
    mesh = mesh,
    time_varying = ~ 1 ,
    time_varying_type = "rw0",
    spatial = "on",
    family =  tweedie(link = "log"),
    time = "year",
    spatiotemporal = "ar1",
    offset = data$log.effort
  )
  
  m_cv <- sdmTMB_cv(  
    data=data,
    formula = commercial.weight ~ 0  + bottom.temp.s + bottom.temp.s2 , 
    mesh = mesh,
    time_varying = ~ 1 ,
    time_varying_type = "rw0",
    spatial = "on",
    family =  tweedie(link = "log"),
    time = "year",
    spatiotemporal = "ar1",
    offset = "log.effort",
    fold_ids = 'fold_id',
    k_folds = k_folds
  )
  
  c <-mod.select.fn()
  mod.select <- rbind(mod.select, c)
  
  summary(m)
  
  # Save model results:
  save(m, file = paste0(rdatapath, "sGSL_",  mod.label, "_", label,  ".rdata"))
  
  m6 <- m
  m6_cv <- m_cv
  
  rm(m, m_cv)       
  mod.select2 <- as.data.frame(lapply(mod.select, as.character), stringsAsFactors=FALSE)       
  write.csv(mod.select2, file = paste0(results.path, "/sGSL_",  label, "_Model Selection table_m1_", mod.label, "_", min(years), "_", max(years), "_", cutoff, "cutoff.csv"), row.names = FALSE)    
  
}


if ("m7" %in% models) {
  
  mod.label <- "m7" 
  
  m <- sdmTMB(  
    data=data,
    formula = commercial.weight ~  0 , 
    mesh = mesh,
    time_varying = ~ 1 + depth.s,
    time_varying_type = "rw0",
    spatial = "on",
    family =  tweedie(link = "log"),
    time = "year",
    spatiotemporal = "ar1",
    offset = data$log.effort
  )
  
  m_cv <- sdmTMB_cv(  
    data=data,
    formula = commercial.weight ~  0 , 
    mesh = mesh,
    time_varying = ~ 1 + depth.s,
    time_varying_type = "rw0", 
    spatial = "on",
    family =  tweedie(link = "log"),
    time = "year",
    spatiotemporal = "ar1",
    offset = "log.effort",
    fold_ids = 'fold_id',
    k_folds = k_folds
  )
  
  c <-mod.select.fn()
  mod.select <- rbind(mod.select, c)
  
  summary(m)
  
  # Save model results:
  save(m, file = paste0(rdatapath, "sGSL_",  mod.label, "_", label,  ".rdata"))
  
  m7 <- m
  m7_cv <- m_cv
  
  rm(m, m_cv)        
  mod.select2 <- as.data.frame(lapply(mod.select, as.character), stringsAsFactors=FALSE)        
  write.csv(mod.select2, file = paste0(results.path, "/sGSL_",  label, "_Model Selection table_m1_", mod.label, "_", min(years), "_", max(years), "_", cutoff, "cutoff.csv"), row.names = FALSE)    
  
}

if ("m8" %in% models) {
  
  mod.label <- "m8" 
  
  m <- sdmTMB(  
    data=data,
    formula = commercial.weight ~ 0   , 
    mesh = mesh,
    time_varying = ~ 1 + depth.s + depth.s2,
    time_varying_type = "rw0",
    spatial = "on",
    family =  tweedie(link = "log"),
    time = "year",
    spatiotemporal = "ar1",
    offset = data$log.effort
  )
  
  m_cv <- sdmTMB_cv(  
    data=data,
    formula = commercial.weight ~ 0  , 
    mesh = mesh,
    time_varying = ~ 1 + depth.s + depth.s2, 
    time_varying_type = "rw0",
    spatial = "on",
    family =  tweedie(link = "log"),
    time = "year",
    spatiotemporal = "ar1",
    offset = "log.effort",
    fold_ids = 'fold_id',
    k_folds = k_folds
  )
  
  c <-mod.select.fn()
  mod.select <- rbind(mod.select, c)
  
  summary(m)
  
  # Save model results:
  save(m, file = paste0(rdatapath, "sGSL_",  mod.label, "_", label,  ".rdata"))
  
  m8 <- m
  m8_cv <- m_cv
  
  rm(m, m_cv)        
  mod.select2 <- as.data.frame(lapply(mod.select, as.character), stringsAsFactors=FALSE)        
  write.csv(mod.select2, file = paste0(results.path, "/sGSL_",  label, "_Model Selection table_m1_", mod.label, "_", min(years), "_", max(years), "_", cutoff, "cutoff.csv"), row.names = FALSE)    
  
}

if ("m9" %in% models) {
  
  mod.label <- "m9" 
  
  m <- sdmTMB(  
    data=data,
    formula = commercial.weight ~ 0 + bottom.temp.s  , 
    mesh = mesh,
    time_varying = ~ 1 + depth.s + depth.s2,
    time_varying_type = "rw0",
    spatial = "on",
    family =  tweedie(link = "log"),
    time = "year",
    spatiotemporal = "ar1",
    offset = data$log.effort
  )
  
  m_cv <- sdmTMB_cv(  
    data=data,
    formula = commercial.weight ~ 0 + bottom.temp.s  , 
    mesh = mesh,
    time_varying = ~ 1 + depth.s + depth.s2, 
    time_varying_type = "rw0",
    spatial = "on",
    family =  tweedie(link = "log"),
    time = "year",
    spatiotemporal = "ar1",
    offset = "log.effort",
    fold_ids = 'fold_id',
    k_folds = k_folds
  )
  
  
  c <-mod.select.fn()
  mod.select <- rbind(mod.select, c)
  
  summary(m)
  
  # Save model results:
  save(m, file = paste0(rdatapath, "sGSL_",  mod.label, "_", label,  ".rdata"))
  
  m9 <- m
  m9_cv <- m_cv
  
  rm(m, m_cv)        
  mod.select2 <- as.data.frame(lapply(mod.select, as.character), stringsAsFactors=FALSE)        
  write.csv(mod.select2, file = paste0(results.path, "/sGSL_",  label, "_Model Selection table_m1_", mod.label, "_", min(years), "_", max(years), "_", cutoff, "cutoff.csv"), row.names = FALSE)    
  
}





if ("m11" %in% models) {
  
  mod.label <- "m11" 
  
  m <- sdmTMB(  
    data=data,
    formula = commercial.weight ~ 0  , 
    mesh = mesh,
    time_varying = ~ 1 + depth.s + depth.s2 + bottom.temp.s,
    time_varying_type = "rw0",
    spatial = "on",
    family =  tweedie(link = "log"),
    time = "year",
    spatiotemporal = "ar1",
    offset = data$log.effort
  )
  
  m_cv <- sdmTMB_cv(  
    data=data,
    formula = commercial.weight ~ 0  , 
    mesh = mesh,
    time_varying = ~ 1 + depth.s + depth.s2 + bottom.temp.s,
    time_varying_type = "rw0",
    spatial = "on",
    family =  tweedie(link = "log"),
    time = "year",
    spatiotemporal = "ar1",
    offset = "log.effort",
    fold_ids = 'fold_id',
    k_folds = k_folds
  )
  
  c <-mod.select.fn()
  mod.select <- rbind(mod.select, c)
  
  summary(m)
  
  # Save model results:
  save(m, file = paste0(rdatapath, "sGSL_",  mod.label, "_", label,  ".rdata"))
  
  m11 <- m
  m11_cv <- m_cv
  
  rm(m, m_cv)        
  mod.select2 <- as.data.frame(lapply(mod.select, as.character), stringsAsFactors=FALSE)        
  write.csv(mod.select2, file = paste0(results.path, "/sGSL_",  label, "_Model Selection table_m1_", mod.label, "_", min(years), "_", max(years), "_", cutoff, "cutoff.csv"), row.names = FALSE)    
  
}






if ("m15" %in% models) {
  
  mod.label <- "m15" 
  
  m <- sdmTMB(  
    data=data,
    formula = commercial.weight ~ 0  , 
    mesh = mesh,
    time_varying = ~ 1 + depth.s + depth.s2 + bottom.temp.s + bottom.temp.s2,
    time_varying_type = "rw0",
    spatial = "on",
    family =  tweedie(link = "log"),
    time = "year",
    spatiotemporal = "ar1",
    offset = data$log.effort
  )
  
  m_cv <- sdmTMB_cv(  
    data=data,
    formula = commercial.weight ~ 0  , 
    mesh = mesh,
    time_varying = ~ 1 + depth.s + depth.s2 + bottom.temp.s + bottom.temp.s2,
    time_varying_type = "rw0",
    spatial = "on",
    family =  tweedie(link = "log"),
    time = "year",
    spatiotemporal = "ar1",
    offset = "log.effort",
    fold_ids = 'fold_id',
    k_folds = k_folds
  )
  
  c <-mod.select.fn()
  mod.select <- rbind(mod.select, c)
  
  summary(m)
  
  # Save model results:
  save(m, file = paste0(rdatapath, "sGSL_",  mod.label, "_", label,  ".rdata"))
  
  m15 <- m
  m15_cv <- m_cv
  
  rm(m, m_cv)        
  mod.select2 <- as.data.frame(lapply(mod.select, as.character), stringsAsFactors=FALSE)        
  write.csv(mod.select2, file = paste0(results.path, "/sGSL_",  label, "_Model Selection table_m1_", mod.label, "_", min(years), "_", max(years), "_", cutoff, "cutoff.csv"), row.names = FALSE)    
  
}


 # Save all the files needed for maps and diagnostics, so you don't need to run models again:
 save(mod.select, cutoff, label, data, max.depthMI,  min.depthMI,  years, xlimMI, ylimMI,   file = paste0(rdatapath, "sGSL_", label,  ".rdata"))
 
 