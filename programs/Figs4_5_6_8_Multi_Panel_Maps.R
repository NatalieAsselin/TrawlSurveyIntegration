##Script to make multi-panel figures of data and results
library(sf)

options(scipen = 999) ##this should take out scientific notation 


##set paths

main.fp <- getwd()

fig.path <- file.path(main.fp, "results/figures")
data.path <- file.path(main.fp, "data")
results.path <- file.path(main.fp, "results/tables")

path <- paste0(main.fp, "/programs/")
rdatapath <- paste0(main.fp, "/programs/RData/")

mod.label <- "m11"



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


##labels for file names #3reset after loading in data files as they include other info. 

##load general values from integrated dataset 
##includes knots for basis splines (k3 and k4)- defined for integrated dataset and imported here, SO THAT k3 AND K4 ARE CONSISTENT FOR ALL MODELS. 

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

##bring in data used in individual survey analyses and append to identify as that survey

load(paste0(rdatapath, "/sGSL_commercial biomass_NS.rdata"))

data_NS <- data
# max.depth_NS <- max.depth
# min.depth_NS <- min.depth


load(paste0(rdatapath, "/sGSL_commercial biomass_RV.rdata"))
data_RV <- data
# max.depth_RV <- max.depth
# min.depth_RV <- min.depth


load(paste0(rdatapath, "/sGSL_commercial biomass_MI.rdata"))
data_MI <- data
# max.depth_MI <- max.depth
# min.depth_MI <- min.depth


load(paste0(rdatapath, "/sGSL_commercial biomass_ALL.rdata"))
data_ALL <- data

##load  modelled dataset
glf <- read.csv(file = paste0(data.path, "/sGSL_SpatioTemporalModelDataset.csv"))

glf$survey <- ifelse (glf$survey %in% c("NS_Fall", "NS_Neph_Fall", "NS_Neph_Spring", "NS_Spring"), "NS", glf$survey)

glf$density <- (glf$commercial.weight / glf$swept.area) * 1000000


##load in individual survey results
SRVS <- read.csv (paste0(results.path, "/sGSL_2001_2025_1Kres_RV_Predictions_", mod.label, "ALL.csv"))
SRVS$survey <- "SRVS"
SRVS$CV_ratio_log <- log(SRVS$CV_ratio)

NS <- read.csv (paste0(results.path, "/sGSL_2001_2025_1Kres_NS_Predictions_", mod.label, "ALL.csv"))
NS$survey <- "NS"
NS$CV_ratio_log <- log(NS$CV_ratio)

MI <- read.csv (paste0(results.path, "/sGSL_2001_2025_0.5Kres_MI_Predictions_", mod.label, "ALL.csv"))
MI$survey <- "MI"
MI$CV_ratio_log <- log(MI$CV_ratio)

ALL <- read.csv (paste0(results.path, "/sGSL_2001_2025_1Kres_ALL_Predictions_", mod.label, "ALL.csv"))
ALL$survey <- "integrated"



survey <- c("NS", "RV", "MI")




load(file = paste0(data.path, "/mapdata.rdata"))

##multi-panel figures of the three surveys

for (i in 1:length(survey)) {
  
  survey[i]
  

  yr.fig  <- if (survey [i] == "NS") {
    c(2001, 2005, 2010, 2015, 2020, 2025)
  }else if (survey [i] == "RV") {
    c(2001, 2005, 2010, 2015, 2020, 2025)
  }else if (survey [i] == "MI") {
    c(2001, 2005, 2010, 2015, 2020, 2022)
  }else {
    c(2001, 2005, 2010, 2015, 2020, 2025)
  }

  ##looking at the last few years
  
  # yr.fig  <- if (survey [i] == "NS") {
  #   c(2018:2023)
  # }else if (survey [i] == "RV") {
  #   c(2018:2023)
  # }else if (survey [i] == "MI") {
  #   c(2017:2022)
  # }else {
  #   c(2018:2023)
  # }
  
  width <- 18.2 /2.54
  
  height <- if (survey [i] == "NS") {
   9.7/2.54
  }else if (survey [i] == "RV") {
   11.2/2.54
  }else if (survey [i] == "MI") {
    11.2/2.54
  }else {
    11.2/2.54
  }
  
#png(file = paste0(fig.path, "/sGSL_Multi_Panel_Results_", survey [i], ".png"), res = 500, units = "in", height = height, width = width)

#png(file = paste0(fig.path, "/sGSL_Multi_Panel_Results_", survey [i], "_", min(yr.fig), "_", max(yr.fig), ".png"), res = 500, units = "in", height = height, width = width)
grDevices::pdf(file = paste0(fig.path, "/sGSL_Multi_Panel_Results_", survey [i], "_", min(yr.fig), "_", max(yr.fig),  ".pdf"),  height = height, width = width, pointsize=7)


m <- kronecker(matrix(c(1:6, 8:13, 15:20, 22:27), ncol = 6, byrow=TRUE), matrix(1, ncol = 5, nrow = 5)) ##larger squares for the maps
n <-kronecker(matrix(c(7, 14, 21, 28), ncol = 1), matrix(1, ncol = 2, nrow = 5)) ##smaller squares for the maps

m <- rbind(0, cbind(0, m, n, 0), 0)
layout(m)
par(mar = c(0, 0, 0, 0))


##set up xlim and ylim

xlim.map <- if (survey[i]=="NS") {
  xlimNS 
}else if (survey[i]=="RV"){
  xlimRV
}else if (survey[i]=="MI"){
  xlimMI
}

ylim.map <- if (survey[i]=="NS") {
  ylimNS 
}else if (survey[i]=="RV"){
  ylimRV
}else if (survey[i]=="MI"){
  ylimMI
}

###annual points

x <- glf [glf$survey ==survey [i],]

for (j in 1:length(yr.fig)) {

  y <- x [x$year == yr.fig [j],  ]
  

plot(xlim.map, ylim.map, type = "n", xlab = "", ylab = "", 
     xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n")
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "white")
plot(map_data[1],xlim=xlim.map, ylim=ylim.map, col= "grey75", border="grey40" , main="", add=TRUE)

#points(y$X, y$Y, pch = 19, bg = "black", col = "black", cex = 0.5)

scale = 0.02    #play with scale for size of points here
index <- y$density == 0

points(y$X [index], y$Y[index], pch = 4, cex = 0.3, lwd = 0.5)  #zero points

col_circles <- grDevices::adjustcolor("black", alpha = 0.50)
points(y$X[!index], y$Y[!index], pch = 1, cex = scale*sqrt(y$density[!index]), col=col_circles)  #non-zero points adjusted to scale

v <- c(0, 500, 2500, 5000, 15000, 30000)  #vector of the number per hectare to identify in legend. keep same number of categories


# Display survey year:  top right    for NS and RV, bottom right for MI


if (survey[i]=="NS") {
  text(0.98*(xlim.map[2]-xlim.map[1])+xlim.map[1],
       0.9 * (ylim.map[2]-ylim.map[1]) + ylim.map [1],
       yr.fig [j], cex = 1.8, pos = 2)
}else if (survey[i]=="RV"){
  text(0.98*(xlim.map[2]-xlim.map[1])+xlim.map[1],
       0.9 * (ylim.map[2]-ylim.map[1]) + ylim.map [1],
       yr.fig [j], cex = 2, pos = 2)
}else if (survey[i]=="MI"){
  text(0.98*(xlim.map[2]-xlim.map[1])+xlim.map[1],
       0.2 * (ylim.map[2]-ylim.map[1]) + ylim.map [1],
       yr.fig [j], cex = 2, pos = 2)
}


##add label with number of sets

nsets <- paste0("n = ", nrow(y))

if (survey[i]=="NS") {
  text(0.98*(xlim.map[2]-xlim.map[1])+xlim.map[1],
       0.82 * (ylim.map[2]-ylim.map[1]) + ylim.map [1],
       nsets , cex = 1, pos = 2)
}else if (survey[i]=="RV"){
  text(0.98*(xlim.map[2]-xlim.map[1])+xlim.map[1],
       0.82 * (ylim.map[2]-ylim.map[1]) + ylim.map [1],
       nsets , cex =1, pos = 2)
}else if (survey[i]=="MI"){
      text(0.98*(xlim.map[2]-xlim.map[1])+xlim.map[1],
       0.13 * (ylim.map[2]-ylim.map[1]) + ylim.map [1],
       nsets , cex = 1, pos = 2)
}




##add probability of occurence 

p <- round((sum(y$commercial.weight > 0) )/ length(y$commercial.weight), 2)

p.occ = bquote(italic(P) == .(format(p, digits = 2)))##make label with italics for probability of occurrence
#p.occ<-  paste0("P = ", p)


if (survey[i]=="NS") {
  text(0.98*(xlim.map[2]-xlim.map[1])+xlim.map[1],
       0.75 * (ylim.map[2]-ylim.map[1]) + ylim.map [1],
       p.occ , cex = 1, pos = 2)
}else if (survey[i]=="RV"){
  text(0.98*(xlim.map[2]-xlim.map[1])+xlim.map[1],
       0.75 * (ylim.map[2]-ylim.map[1]) + ylim.map [1],
       p.occ , cex =1, pos = 2)
}else if (survey[i]=="MI"){
  text(0.98*(xlim.map[2]-xlim.map[1])+xlim.map[1],
       0.06 * (ylim.map[2]-ylim.map[1]) + ylim.map [1],
       p.occ , cex = 1, pos = 2)
}


box(col = "grey50")

}

##legend for points

leg_title <- expression(Kg/km^{2}) 
plot.new() 

vlabels <- format(v, big.mark = ",", scientific = FALSE)

legend("center",
       legend = vlabels,   #use vector v for legend
       pch = c(4, 1, 1, 1, 1, 1),  #symbols to use for the numbers (0, 1, 5, 25, 100)
       pt.cex = c(0.5, scale*sqrt(v[2:6])),  #size of the symbols with adjustment with scale for non-zero values
       bg = "white",     #background color of the legend
       bty="n",
       ncol=1,
       #text.width = c(0.3, 0.3),
       x.intersp = 1.7,
      # y.intersp= c(1, 2.1, 2.1, 2.1, 2.1, 2.1) ,
      y.intersp= c(0.9, 1.9, 1.9, 1.9, 1.9, 1.9) ,
     # y.intersp= c(0.7, 1.5, 1.5, 1.5, 1.5, 1.5) ,
      title = leg_title,
      #title = NULL,
       cex=0.8)



box(col = "grey50")

###predictions 

##color scheme for survey predictions
cols <- colorRampPalette(c("blue4", "blue", "mediumturquoise", "yellow", "orange", "red", "firebrick4"))
#cols <- colorRampPalette(c("lightsteelblue1", "skyblue1", "mediumturquoise", "yellow", "orange", "red", "firebrick4"))
#cols <- colorRampPalette(c("white", "lightyellow", "yellow", "orange", "darkorange", "red", "firebrick4"))

# cols <- colorRampPalette( c("aliceblue",
#                             "lightskyblue",
#                             "lemonchiffon",
#                             "gold",
#                             "darkorange",
#                             "red",
#                             "firebrick4"))

breaks <- c(0, 10,  50, 100, 500, 1000, 5000, 10000, 25000, 50000, 100000)


#breaks <- c(0, 10,  50, 100, 500, 1000, 2500, 5000, 10000, 25000, 50000)
x <- if (survey[i]=="NS") {
  NS 
}else if (survey[i]=="RV"){
  SRVS
}else if (survey[i]=="MI"){
  MI
}

###prediction from integrated survey model

for (j in 1:length(yr.fig)) {
  
  y <- x [x$year == yr.fig [j],  ]

  # Square out predictions: 
  zz <- if (survey[i]=="NS") {
    y$pred_modNS 
  }else if (survey[i]=="RV"){
    y$pred_modRV 
  }else if (survey[i]=="MI"){
    y$pred_modMI 
  }
  
  
  xgrid <- unique(y$X)
  ygrid <- unique(y$Y)
  dim(zz) <- c(length(xgrid), length(ygrid))
  

  
  # Plot map:

  plot(xlim.map, ylim.map, type = "n", xlab = "", ylab = "", 
       xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n")
  
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "white")
  
  image(xgrid, ygrid, zz, col = cols(length(breaks)-1), 
        
        xaxt = "n", yaxt = "n", breaks = breaks, 
        xlim = xlim.map, ylim = ylim.map, add = TRUE)
  
   plot(map_data[1],xlim=xlim.map, ylim=ylim.map, col= "grey75", border="grey15" , main="", add=TRUE)
   box(col = "grey50")
}

##legend for points

leg_title<- expression(Kg/km^{2}) 
plot.new() 

leg.breaks <- format(breaks, big.mark = ",", scientific = FALSE)
legend("center", legend = leg.breaks, pch = 22, pt.bg = cols(length(breaks)), 
       bty="n", pt.cex = 2, pt.lwd = 0.5, cex = 0.8, title =  leg_title)

box(col = "grey50")


###prediction from integrated survey model


for (j in 1:length(yr.fig)) {
  
  y <- x [x$year == yr.fig [j],  ]
  
  # Square out predictions: 
  zz <-  y$pred_modALL 
  xgrid <- unique(y$X)
  ygrid <- unique(y$Y)
  dim(zz) <- c(length(xgrid), length(ygrid))
  
  
  # Plot map:
  plot(xlim.map, ylim.map, type = "n", xlab = "", ylab = "", 
       xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n")
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "white")
  
  image(xgrid, ygrid, zz, col = cols(length(breaks)-1), 
        
        xaxt = "n", yaxt = "n", breaks = breaks, 
        xlim = xlim.map, ylim = ylim.map, add = TRUE)
  
  plot(map_data[1],xlim=xlim.map, ylim=ylim.map, col= "grey75", border="grey15" , main="", add=TRUE)
  box(col = "grey50")
}

##legend for points

leg_title<- expression(Kg/km^{2}) 
plot.new() 

leg.breaks <- format(breaks, big.mark = ",", scientific = FALSE)
legend("center", legend = leg.breaks, pch = 22, pt.bg = cols(length(breaks)), 
       bty="n", pt.cex = 2, pt.lwd = 0.5, cex = 0.8, title =  leg_title)

box(col = "grey50")
##CV ratio of integrated vs single survey model


for (j in 1:length(yr.fig)) {
  
  y <- x [x$year == yr.fig [j],  ]
  
  cols <-  colorRampPalette(c( "blue", "white", "red"))
# y$CV_ratio_map <- ifelse(y$CV_ratio>2, 2, y$CV_ratio)
  y$CV_ratio_log_map <- y$CV_ratio_log 
    
  y$CV_ratio_log_map <- ifelse(y$CV_ratio_log_map <= -1, -1 ,y$CV_ratio_log_map)
  y$CV_ratio_log_map <- ifelse(y$CV_ratio_log_map >= 1, 1, y$CV_ratio_log_map)

  
  breaks <-seq(from = -1, to = 1, by= 0.001)
  leg.breaks <- seq(from = -1, to = 1, by= 0.2)
 
  # leg.breaks.label <- c(seq(from = 0, to = 1.8, by= 0.2), expression(Value %>=% 2))

  
  # Square out predictions: 
  zz <-  y$CV_ratio_log_map
  xgrid <- unique(y$X)
  ygrid <- unique(y$Y)
  dim(zz) <- c(length(xgrid), length(ygrid))
  
# Plot map:
  plot(xlim.map, ylim.map, type = "n", xlab = "", ylab = "", 
       xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n")
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "white")
  
  image(xgrid, ygrid, zz, col = cols(length(breaks)-1), 
        
        xaxt = "n", yaxt = "n", breaks = breaks, 
        xlim = xlim.map, ylim = ylim.map, add = TRUE)
  
  plot(map_data[1],xlim=xlim.map, ylim=ylim.map, col= "grey75", border="grey15" , main="", add=TRUE)
  
  # ##block out areas of very low density
  # y$mask <- ifelse(y$pred_modALL > 10, NA, 1 )
  # 
  # zz <- y$mask
  # xgrid <- unique(y$X)
  # ygrid <- unique(y$Y)
  # dim(zz) <- c(length(xgrid), length(ygrid))
  # 
  # mask.col <- grDevices::adjustcolor("grey90", alpha = 0.50)
  # 
  # image(xgrid, ygrid, zz, col = mask.col,
  #     xaxt = "n", yaxt = "n",
  #       xlim = xlim.map, ylim = ylim.map, add = TRUE)
  
  
  
  
  plot(map_data[1],xlim=xlim.map, ylim=ylim.map, col= "grey75", border="grey15" , main="", add=TRUE)
  

  box(col = "grey50")
  

  
  }

##legend for points

leg_title<- "Log CV ratio"
plot.new() 
legend("center", legend =  c(expression(""<="-1"),seq(from = -0.8, to = 0.8, by= 0.2), expression("">="1") ), pch = 22, pt.bg = cols(length(leg.breaks)),
   bty="n",pt.cex = 2, pt.lwd = 0.5, cex = 0.8, title =  leg_title)


#legend("center", legend =  leg.breaks, pch = 22, pt.bg = cols(length(leg.breaks)),
 #      bty="n",pt.cex = 2, pt.lwd = 0.5, cex = 1, title =  leg_title)

box(col = "grey50")

dev.off()
}

#density check

x <- rev(sort(ALL$pred_modALL))[1:100]

y <- rev(sort(NS$pred_modALL ))[1:100]


##Multi-panel figure

##integrated survey model resuls. 


{
 yr.fig  <- c(2001, 2005, 2010, 2015, 2020, 2025)
  
  #yr.fig  <- c(2018:2023)
  
  width <- 18.2/2.54

  height <- 11.2/2.54
  
  
 # png(file = paste0(fig.path, "/sGSL_Multi_Panel_Results_", "ALL", "_", min(yr.fig), "_", max(yr.fig), ".png"), res = 500, units = "in", height = height, width = width)
  
  grDevices::pdf(file = paste0(fig.path, "/sGSL_Multi_Panel_Results_", "ALL", "_", min(yr.fig), "_", max(yr.fig), ".pdf"),  height = height, width = width, pointsize=7)
  
  
  m <- kronecker(matrix(c(1:6, 8:13, 15:20, 22:27), ncol = 6, byrow=TRUE), matrix(1, ncol = 5, nrow = 5)) ##larger squares for the maps
  n <-kronecker(matrix(c(7, 14, 21, 28), ncol = 1), matrix(1, ncol = 2, nrow = 5)) ##smaller squares for the legends
  
  m <- rbind(0, cbind(0, m, n, 0), 0)
  layout(m)
  par(mar = c(0, 0, 0, 0))
  
  
  ##set up xlim and ylim
  
  xlim.map <-  xlimRV
 
  
  ylim.map <-  ylimRV

  
  ###all points
  
  x <- glf 
  
  for (j in 1:length(yr.fig)) {
    
    y <- x [x$year == yr.fig [j],  ]
    
    
    plot(xlim.map, ylim.map, type = "n", xlab = "", ylab = "", 
         xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n")
    rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "white")
    plot(map_data[1],xlim=xlim.map, ylim=ylim.map, col= "grey75", border="grey40" , main="", add=TRUE)
    
    #points(y$X, y$Y, pch = 19, bg = "black", col = "black", cex = 0.5)
    
    scale = 0.02    #play with scale for size of points here
    index <- y$density == 0
    
    points(y$X [index], y$Y[index], pch = 4, cex = 0.3, lwd = 0.5)  #zero points
   # points(y$X[!index], y$Y[!index], pch = 1, cex = scale*sqrt(y$density[!index]))  #non-zero points adjusted to scale
    
    col_circles <- grDevices::adjustcolor("black", alpha = 0.50)
    points(y$X[!index], y$Y[!index], pch = 1, cex = scale*sqrt(y$density[!index]), col=col_circles)  #non-zero points adjusted to scale
    
    
    
    v <- c(0, 500, 2500, 5000, 15000, 30000)  #vector of the kg/km2 to identify in legend. keep same number of categories
    
    
    # Display survey year:  top right    for NS and RV, bottom right for MI
    
    

      text(0.98*(xlim.map[2]-xlim.map[1])+xlim.map[1],
           0.9 * (ylim.map[2]-ylim.map[1]) + ylim.map [1],
           yr.fig [j], cex = 2, pos = 2)

    
    
    ##add label with number of sets
    
    nsets <- paste0("n = ", nrow(y))

      text(0.98*(xlim.map[2]-xlim.map[1])+xlim.map[1],
           0.82 * (ylim.map[2]-ylim.map[1]) + ylim.map [1],
           nsets , cex =1, pos = 2)

      ##add probability of occurence 
      
      p <- round((sum(y$commercial.weight > 0) )/ length(y$commercial.weight), 2)
      
      p.occ = bquote(italic(P) == .(format(p, digits = 2)))##make label with italics for probability of occurrence
      #p.occ<-  paste0("P = ", p)
      
      
      text(0.98*(xlim.map[2]-xlim.map[1])+xlim.map[1],
           0.75 * (ylim.map[2]-ylim.map[1]) + ylim.map [1],
           p.occ , cex =1, pos = 2)
      
      
    box(col = "grey50")
    
  }
  

  
  ##legend for points
  
  leg_title <- expression(Kg/km^{2}) 
  plot.new() 
  
  leg.v <- format(v, big.mark = ",", scientific = FALSE)
  
  legend("center",
         legend = leg.v,   #use vector v for legend
         pch = c(4, 1, 1, 1, 1, 1),  #symbols to use for the numbers (0, 1, 5, 25, 100)
         pt.cex = c(0.5, scale*sqrt(v[2:6])),  #size of the symbols with adjustment with scale for non-zero values
         bg = "white",     #background color of the legend
         bty="n",
         ncol=1,
         #text.width = c(0.3, 0.3),
         x.intersp = 1.8,
        # y.intersp= c(1, 2.1, 2.1, 2.1, 2.1, 2.1) ,
        y.intersp= c(0.9, 1.9, 1.9, 1.9, 1.9, 1.9) ,
       #  y.intersp= c(0.5, 1.9, 1.9, 1.9, 1.9, 1.9) ,
         title = leg_title,
         cex=0.8)
  
  box(col = "grey50")
  
  ###predictions 
  
  ##color scheme for survey predictions
  cols <- colorRampPalette(c("blue4", "blue", "mediumturquoise", "yellow", "orange", "red", "firebrick4"))
 # cols <- colorRampPalette(c("lightsteelblue1", "skyblue1", "mediumturquoise", "yellow", "orange", "red", "firebrick4"))
  
  # cols <- colorRampPalette( c("aliceblue",
  #                             "lightskyblue",
  #                             "lemonchiffon",
  #                             "gold",
  #                             "darkorange",
  #                             "red",
  #                             "firebrick4"))
  
  
 # cols <- colorRampPalette(c("white", "lightyellow", "yellow", "orange", "darkorange", "red", "firebrick4"))
breaks <- c(0, 10,  50, 100, 500, 1000, 5000, 10000, 25000, 50000, 100000)
  
 # breaks <- c(0, 10,  50, 100, 500, 1000, 2500, 5000, 10000, 25000, 50000)
  
  x <- ALL
  
  ###prediction from integrated survey model
  
  for (j in 1:length(yr.fig)) {
    
    y <- x [x$year == yr.fig [j],  ]
    
    # Square out predictions: 
    zz <- y$pred_modALL 
    
    xgrid <- unique(y$X)
    ygrid <- unique(y$Y)
    dim(zz) <- c(length(xgrid), length(ygrid))
    
    
    # Plot map:
    plot(xlim.map, ylim.map, type = "n", xlab = "", ylab = "", 
         xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n")
    rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "white")
    
    image(xgrid, ygrid, zz, col = cols(length(breaks)-1), 
          
          xaxt = "n", yaxt = "n", breaks = breaks, 
          xlim = xlim.map, ylim = ylim.map, add = TRUE)
    
    ##block out areas of very low density
    # y$mask <- ifelse(y$pred_modALL > 10, NA, 1 )
    # 
    # zz <- y$mask
    # xgrid <- unique(y$X)
    # ygrid <- unique(y$Y)
    # dim(zz) <- c(length(xgrid), length(ygrid))
    # 
    # mask.col <- grDevices::adjustcolor("grey90", alpha = 0.50)
    # 
    # image(xgrid, ygrid, zz, col = mask.col,
    #       xaxt = "n", yaxt = "n",
    #       xlim = xlim.map, ylim = ylim.map, add = TRUE)
    # 
  
    
    plot(map_data[1],xlim=xlim.map, ylim=ylim.map, col= "grey75", border="grey15" , main="", add=TRUE)
    box(col = "grey50")
  }
  
  ##legend for points
  
  leg.breaks <- format(breaks, big.mark = ",", scientific = FALSE)
  
  leg_title<- expression(Kg/km^{2}) 
  plot.new() 
  legend("center", legend = leg.breaks, pch = 22, pt.bg = cols(length(breaks)), 
         bty="n", pt.cex = 2, pt.lwd = 0.5, cex = 0.8, title =  leg_title)
  
  box(col = "grey50")
 
  
   ###SE from integrated survey model
  
 breaks <- c(0, 10,  50, 100, 500, 1000, 5000, 10000, 25000, 50000, 100000)
  #breaks <- c(0, 10,  50, 100, 500, 1000, 2500, 5000, 10000, 25000, 50000)
  
  for (j in 1:length(yr.fig)) {
    
    y <- x [x$year == yr.fig [j],  ]
    
    # Square out predictions: 
    zz <-  y$pred_SE_modALL
    xgrid <- unique(y$X)
    ygrid <- unique(y$Y)
    dim(zz) <- c(length(xgrid), length(ygrid))
    
    
    # Plot map:
    plot(xlim.map, ylim.map, type = "n", xlab = "", ylab = "", 
         xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n")
    rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "white")
    
    image(xgrid, ygrid, zz, col = cols(length(breaks)-1), 
          
          xaxt = "n", yaxt = "n", breaks = breaks, 
          xlim = xlim.map, ylim = ylim.map, add = TRUE)
    
    # ##block out areas of very low density
    # y$mask <- ifelse(y$pred_modALL > 10, NA, 1 )
    # 
    # zz <- y$mask
    # xgrid <- unique(y$X)
    # ygrid <- unique(y$Y)
    # dim(zz) <- c(length(xgrid), length(ygrid))
    # 
    # image(xgrid, ygrid, zz, col = "grey90", 
    #       xaxt = "n", yaxt = "n", 
    #       xlim = xlim.map, ylim = ylim.map, add = TRUE)
    
    plot(map_data[1],xlim=xlim.map, ylim=ylim.map, col= "grey75", border="grey15" , main="", add=TRUE)
    box(col = "grey50")
  }
  
  ##legend for points
 
 leg.breaks <- format(breaks, big.mark = ",", scientific = FALSE) 
 
  leg_title<- expression(Kg/km^{2}) 
  plot.new() 
  legend("center", legend = leg.breaks, pch = 22, pt.bg = cols(length(breaks)), 
         bty="n", pt.cex = 2, pt.lwd = 0.5, cex = 0.8, title =  leg_title)
  
  box(col = "grey50")

  
  ##CV from integrated survey
  
  breaks <- c(0, 0.1,  0.2, 0.3, 0.4, 0.5 , 0.75, 1,2,5)
  
  for (j in 1:length(yr.fig)) {
    
    y <- x [x$year == yr.fig [j],  ]
    

    # Square out predictions: 
    zz <-  y$pred_CV_modALL
    xgrid <- unique(y$X)
    ygrid <- unique(y$Y)
    dim(zz) <- c(length(xgrid), length(ygrid))
    
    # Plot map:
    plot(xlim.map, ylim.map, type = "n", xlab = "", ylab = "", 
         xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n")
    rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "white")
    

    image(xgrid, ygrid, zz, col = cols(length(breaks)-1), 
          xaxt = "n", yaxt = "n", breaks = breaks, 
          xlim = xlim.map, ylim = ylim.map, add = TRUE)
    
    
    ##block out areas of very low density
    y$mask <- ifelse(y$pred_modALL > 10, NA, 1 )

    zz <- y$mask
    xgrid <- unique(y$X)
    ygrid <- unique(y$Y)
    dim(zz) <- c(length(xgrid), length(ygrid))

    image(xgrid, ygrid, zz, col = "grey90",
          xaxt = "n", yaxt = "n",
          xlim = xlim.map, ylim = ylim.map, add = TRUE)

    plot(map_data[1],xlim=xlim.map, ylim=ylim.map, col= "grey75", border="grey15" , main="", add=TRUE)
    
    box(col = "grey50")
    
    
    
  }
  
  ##legend for points
  
  leg_title<- "CV"
  plot.new() 
  #legend("center", legend = breaks, pch = 22, pt.bg = cols(length(breaks)),
   #      bty="n",pt.cex = 2, pt.lwd = 0.5, cex = 0.8, title =  leg_title)
  
  leg.cat <-c( breaks [c(1:(length(breaks))-1)], expression("">="5"))
  
  legend("center", legend = leg.cat, pch = 22, pt.bg = cols(length(breaks)),
         bty="n",pt.cex = 2, pt.lwd = 0.5, cex = 0.8, title =  leg_title)
  
  
  
  box(col = "grey50")
  
  dev.off()
}

