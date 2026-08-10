##Script to compare predictions from all survey model to individual survey models. 

##load libraries

library(sdmTMB)
library(mgcv)
library(dplyr)


## Set paths
main.fp <- getwd()

fig.path <- file.path(main.fp, "results/figures")
data.path <- file.path(main.fp, "data")
results.path <- file.path(main.fp, "results/tables")

path <- paste0(main.fp, "/programs/")
rdatapath <- paste0(main.fp, "/programs/RData/")

##set up a summary table
columns <- c("survey", "mod.survey", "mod.all" )
sum.table <- as.data.frame( matrix(data=NA, nrow =4, ncol=length(columns), byrow=TRUE))
colnames(sum.table )<- columns

##list surveys
sum.table$survey <- c("RV", "NS", "MI", "sGSLM")

##define best model for each survey
sum.table$mod.survey [sum.table$survey =="NS"] <- "m15"
sum.table$mod.survey [sum.table$survey =="RV"] <- "m6"
sum.table$mod.survey [sum.table$survey =="MI"] <- "m7"
sum.table$mod.survey [sum.table$survey == "sGSLM"] <- "m13"

sum.table$mod.all <- "m13"

mod.label<- as.character(sum.table$mod.all[1])

##labels for file names  
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

##bring in data used in individual survey analyses

load(paste0(rdatapath, "/sGSL_commercial biomass_NS.rdata"))
dataNS <- data

load(paste0(rdatapath, "/sGSL_commercial biomass_RV.rdata"))
dataRV <- data

load(paste0(rdatapath, "/sGSL_commercial biomass_MI.rdata"))
dataMI <- data

load(paste0(rdatapath, "/sGSL_commercial biomass_ALL.rdata"))
dataAll <- data

##load general values from integrated dataset 

load( file = paste0(rdatapath, "sGSL_", "GeneralValues",   ".rdata"))

##labels for file names #3reset after loading in data files as they include other info. 
label <- "commercial biomass"

##set up a summary table for results, to merge with overall tables
columns <- c( "survey", "pbias_survey", "pbias_all", "MAE_survey", "MAE_all", "RMSE_survey", "RMSE_all" )
d <- as.data.frame( matrix(data=NA, nrow =0, ncol=length(columns), byrow=TRUE))
colnames(d)<- columns

for (i in 1:length(sum.table$survey)) {
 s <- sum.table$survey[i]

 c <- NULL
 c$survey <- s
 
 ##load in survey specific data
 nd <- if (s=="NS") {
   dataNS 
 } else if(s=="RV")  {
   dataRV 
 } else if(s=="MI")  {
   dataMI 
 } else if(s=="sGSLM")  {
   dataAll 
 }
 
 ##load in survey specific model 
 
 mod <- if (s=="NS") {
   modNS 
 } else if(s=="RV")  {
   modRV 
 } else if(s=="MI")  {
   modMI 
 } else if(s=="sGSLM")  {
   modALL 
 }
 
 
##predict from survey specific model

nsim <- 5000
pp <- predict(mod, newdata = nd, se_fit = FALSE, re_form = NULL, offset=nd$log.effort, nsim = nsim)
pp2 <- exp(pp) 
nd$predicted_survey <- apply(pp2, 1, mean) ##get mean of nsim simulations

c$pbias_survey <- (sum(nd$predicted_survey - nd$commercial.weight)/ sum(nd$commercial.weight)) *100
c$MAE_survey <- mean(abs(nd$commercial.weight - nd$predicted_survey))
c$RMSE_survey <- sqrt(mean((nd$commercial.weight - nd$predicted_survey)^2))

if (s %in% c("RV", "NS", "sGSLM")) {

pp <- predict(modALL, newdata = nd, se_fit = FALSE, re_form = NULL, offset=nd$log.effort, nsim = nsim)

} else if (s %in% c("MI")) {

  
  ##need to do this as MI only has one level of substrate factors
  nd$substrate <- factor(
    as.character(nd$substrate),
    levels = levels(factor(modALL$data$substrate))
  )
    pp <- predict(modALL, newdata = nd, se_fit = FALSE, re_form = NULL, offset=nd$log.effort, nsim = nsim)
  
  
}


pp2 <- exp(pp) 
nd$predicted_ALL <- apply(pp2, 1, mean) ##get mean of nsim simulations

c$pbias_all <-  (sum(nd$predicted_ALL- nd$commercial.weight)/ sum(nd$commercial.weight)) * 100
c$MAE_all <- mean(abs(nd$commercial.weight - nd$predicted_ALL))
c$RMSE_all<- sqrt(mean((nd$commercial.weight - nd$predicted_ALL)^2))

d <- rbind(d, c)

if (s=="NS") {
  NS_nd <-   nd
} else if (s=="RV") {
  RV_nd <-   nd
} else if (s=="MI") {
  MI_nd <-   nd
} else if (s=="sGSLM") {
  ALL_nd <-   nd
}

}

sum.table <-merge(sum.table, d)


# Export to file:

write.csv(sum.table, file = paste0(results.path, "/sGSL_", label, "_MultiSurveyComparisonResults_",mod.label, ".csv"), row.names = FALSE)

save(sum.table, NS_nd, RV_nd, MI_nd, ALL_nd, file = paste0(rdatapath, "sGSL_", label, "_MultiSurveyComparisonResults_",mod.label,  ".rdata"))

##calculate pbias by year

surveys <- c("SRVS", "NSS", "MIS", "ID")

res_year <- NULL

for (i in 1:length(surveys)) {
  
  s <- surveys[i]
  
  data <- if (s=="NSS") {
    NS_nd
  } else if(s=="SRVS")  {
    RV_nd
  } else if(s=="MIS")  {
    MI_nd
  } else if(s=="ID")  {
    ALL_nd
  }


# Calculate bias by year
a <- data %>%
  group_by(year) %>%
  dplyr::summarise(
    pBIAS_survey = (sum(predicted_survey - commercial.weight)/sum(commercial.weight)) * 100
  )

b <- data %>%
  group_by(year) %>%
 dplyr:: summarise(
    pBIAS_All = (sum(predicted_ALL - commercial.weight)/sum(commercial.weight)) * 100
  )

c <- merge(a,b)
c$survey <- s
res_year <- rbind(res_year, c)
}

write.csv(res_year, file = paste0(results.path, "/sGSL_", label, "_PBIAS_By_Year_",mod.label, ".csv"), row.names = FALSE)

save(res_year, file = paste0(rdatapath, "sGSL_", label, "_PBIAS_By_Year_",mod.label,  ".rdata"))


