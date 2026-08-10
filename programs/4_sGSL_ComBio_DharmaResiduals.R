##script for Dharma residuals and associated plots

set.seed <- 787

##load libraries

library(gulf)
library(DHARMa)

## set paths
main.fp <- getwd()

fig.path <- file.path(main.fp, "results/figures")
data.path <- file.path(main.fp, "data")
results.path <- file.path(main.fp, "results/tables")

path <- paste0(main.fp, "/programs/")
rdatapath <- paste0(main.fp, "/programs/RData/")

##labels for file names
#label <- "commercial biomass_NS"
#label <- "commercial biomass_RV"
label <- "commercial biomass_MI"
#label <- "commercial biomass_ALL"

#model.lab <- "m15" ##NS
#model.lab <- "m6" ##RV
model.lab <- "m7" ##MI 
#model.lab <- "m13" ##ALL

##dimensions for multi-year maps

width <- if (label == "commercial biomass_NS") {
  15
}else if (label == "commercial biomass_RV") {
  12
}else {
  12
}

height  <- if (label == "commercial biomass_NS") {
  11
}else if (label == "commercial biomass_RV") {
  11
}else {
  11
}

##load model 

load( file = paste0(rdatapath, "sGSL_",  model.lab, "_", label,  ".rdata")) ##load model
load(file = paste0(rdatapath, "sGSL_", label,  ".rdata")) ##load extra data

##load general values from integrated dataset 

load( file = paste0(rdatapath, "sGSL_", "GeneralValues",   ".rdata"))

model <- m ##assign m to model


xlim_map <- if (label == "commercial biomass_NS") {
  xlimNS
}else if (label == "commercial biomass_RV") {
  xlimRV
}else if (label == "commercial biomass_MI") {
  xlimMI
}else if (label == "commercial biomass_ALL") {
  xlim
}


ylim_map <- if (label == "commercial biomass_NS") {
  ylimNS
}else if (label == "commercial biomass_RV") {
  ylimRV
}else if (label == "commercial biomass_MI") {
  ylimMI
}else if (label == "commercial biomass_ALL") {
  ylim
}


##bring in a base map
load(file = paste0(data.path, "/mapdata.rdata"))


##dharma residuals

sim.mod <- simulate(model, nsim=1000, type = "mle-mvn")

##check if data and simulated data have roughly equal proportion of 0s
sum(data$commercial.weight == 0) / length(data$commercial.weight)
sum(sim.mod == 0 )/ length(sim.mod)

r_mod <- dharma_residuals(sim.mod, model, return_DHARMa = TRUE)

data$dharma.res <- r_mod$scaledResiduals

png(file = paste0(fig.path, "/sGSL_", model.lab, "_", label,  "_dharma_residuals", cutoff, "cutoff.png"), res = 500, units = "in", height = 11, width = 8)

m <- kronecker(matrix(1:5, ncol = 1), matrix(1, ncol = 5, nrow = 5))

m <- rbind(0, cbind(0, m, 0), 0)
layout(m)
par(mar = c(2, 0, 5, 0))

testUniformity(r_mod) 
plotResiduals(r_mod) 
testDispersion(r_mod) 
testOutliers(r_mod, type="binomial") 
testZeroInflation(r_mod)

dev.off()


png(file = paste0(fig.path, "/sGSL_", model.lab, "_", label,  "_dharma_Predictor_residuals", cutoff, "cutoff.png"), res = 500, units = "in", height = 8, width =11)

m <- kronecker(matrix(1:4, ncol = 2), matrix(1, ncol = 5, nrow = 5))

m <- rbind(0, cbind(0, m, 0), 0)
layout(m)
par(mar = c(5, 2, 5, 2))

DHARMa::plotResiduals(r_mod, data$depth)
DHARMa::plotResiduals(r_mod, data$bottom.temp)
DHARMa::plotResiduals(r_mod, as.factor(data$substrate))
DHARMa::plotResiduals(r_mod, as.factor(data$survey))

dev.off()

##map out dharma residuals

clg()

png(file = paste0(fig.path, "/sGSL", "_", model.lab, "_", label, "_dharma_residuals_map_", min(years), "_", max(years), ".png"), 
    res = 500, units = "in", height = height, width = width)

##layout for 2001 to 2024
m <- kronecker(matrix(1:25, ncol = 5), matrix(1, ncol = 5, nrow = 5))
m <- rbind(0, cbind(0, m, 0), 0)
layout(m)
par(mar = c(0, 0, 0, 0))

cols <- colorRampPalette(c( "red", "white", "blue"))

breaks <- seq(from = 0, to = 1, by= 0.1)
leg.breaks <- seq(from = 0, to = 1, by= 0.1) ##same breaks in legend

yr.figure <- if (label == "commercial biomass_MI") {

 c(2001:2022) 
}else{
 c(2001:2025) 

}

data$col <- cols(length(breaks) -1 ) [as.numeric(cut(data$dharma.res, breaks=breaks))]

for (i in min(yr.figure) :max(yr.figure)){
  print(i)
  
  data.yr <- data[data$year==i,]
  
  # Plot map:
  plot(xlim_map, ylim_map, type = "n", xlab = "", ylab = "", 
       xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n")
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "white")
  
  points(data.yr$X, data.yr$Y, pch=20, col=data.yr$col, cex=1)
  
  plot(map_data[1],xlim=xlim_map, ylim=ylim_map, col= "grey75", border="grey15" , main="", add=TRUE)
  
  # Display survey year:  top right    
  text(0.85*(xlim_map[2]-xlim_map[1])+xlim_map[1],
       0.8 * (ylim_map[2]-ylim_map[1]) + ylim_map [1],
       i, cex = 2.5, pos = 3)
}

leg_title<- "Residuals"

box(col = "grey50")

dev.off()
