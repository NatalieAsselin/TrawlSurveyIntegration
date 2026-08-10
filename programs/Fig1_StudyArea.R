##Map of data used in spatiotemporal models

##Load packages from https://github.com/TobieSurette 

library(gulf.data)
library(gulf.spatial)

##Load gulf package (DFO package to access and manipulate data)

library(gulf)

library(TeachingDemos)

###set up paths
main.fp <- getwd()

path <- paste0(main.fp, "/programs/")
fig.path <- file.path(main.fp, "results/figures")
data.path <- file.path(main.fp, "data")
results.path <- file.path(main.fp, "results/tables")


##load  modelled dataset
glf <- read.csv(file = paste0(data.path, "/sGSL_SpatioTemporalModelDataset.csv"))

glf$survey <- ifelse (glf$survey %in% c("NS_Fall", "NS_Neph", "NS_Spring"), "NS", glf$survey)

##width and height in inches
width <- 18.2 / 2.54 # 

height <- 14.0 /2.54  

grDevices::pdf(file=paste0(fig.path , "/", "sGSL_CommercialBiomass_STM_Study Area",  ".pdf"),  height = height, width = width, pointsize=7)

par(mar=c(2,2,2,2)+0.1)##reduce the white area around the map
 
xlim = c(-66.3, -60.3) ###longitudinal extent - used in density maps
ylim = c(45.4, 49.0)###latitudinal extent - used in density maps

gulf.map(xlim=xlim, ylim=ylim,  land=FALSE, sea=FALSE, xtick=FALSE, ytick=FALSE, sea.col=NA, dem.sea=FALSE)

bathymetry(dem = FALSE, dem.col = c("lightsteelblue1", "skyblue3"),
            breaks = c( -100, 0, seq(-100, 0, by = 50)),
            contour = TRUE, contour.col = "skyblue4",
            contour.lwd = c(0.5,1), contour.cex = c(1, 1),
            levels = list(minor = seq(-100, 0, by = 50), major = seq(-100, 0, by = 50)))
   
## add LFAs   
   
   p <- fz.sf.polygons[fz.sf.polygons$species.code==2550, ]
   p <- p[p$region=="gulf" | (p$region =="quebec" & p$label == 22), ]


 plot(p, lwd=1, border="grey30",  col=NA, add=TRUE)

  ##plot land on top again to clean up edges
  par(new=T)

  gulf.map(xlim=xlim, ylim=ylim, land=TRUE, col= "grey75", border="grey40", sea=FALSE,  lwd=0.5, xtick=FALSE, ytick=FALSE,)
  
  
  cols <- grDevices::adjustcolor(c("orangered2", "darkgreen",  "darkblue"), alpha = 0.40)
  
  points(glf$lon [glf$survey=="NS"], glf$lat [glf$survey=="NS"],pch = 18, bg = cols[2], col = cols[2], lwd = 0.5, cex = 0.75)
  points(glf$lon [glf$survey=="MI"], glf$lat [glf$survey=="MI"], pch = 19, bg = cols[3], col = cols[3], lwd = 0.5, cex = 0.75)
  points(glf$lon [glf$survey=="RV"], glf$lat [glf$survey=="RV"] ,pch = 17, bg = cols[1], col = cols[1], lwd = 0.5, cex = 0.75)
  
  # Display survey year:  top right    
    shadowtext(-61, 47.5, "22", cex = 2, pos = 3, col="black", bg= "white")
    shadowtext(-61.5, 46.8, "24", cex = 2, pos = 3, col="black", bg= "white")
    shadowtext(-61.3, 46.4, "26B", cex = 2, pos = 3, col="black", bg= "white")
    shadowtext(-63, 45.8, "26A", cex = 2, pos = 3, col="black", bg= "white")
    shadowtext(-64.4, 46.4, "25", cex = 2, pos = 3, col="black", bg= "white")
    shadowtext(-64.3, 47.5, "23", cex = 2, pos = 3, col="black", bg= "white")
  
    ###NB-NS border
    b<-read.csv(paste0(data.path,"/NB_NS_Border.csv"))
    lines(b$Longitude, b$Latitude, col="grey15", lwd=1)
    

    
    ## add in province labels 
    
    text(-63.13, 46.36, "Prince Edward Island", cex =1.5)
    text(-65.5, 46.75, "New Brunswick", cex =1.5)
    text(-63.5, 45.55, "Nova Scotia", cex =1.5)
    text(-65.5, 48.5, "Québec", cex =1.5)
    text(-63, 47.25, "Gulf of St. Lawrence", cex =1.5, font=3, col="blue4")
    
    
    ### add in scale bar
    scale.bar(longitude=-61.8, latitude=45.45, length=100, units="km", convert = TRUE,
              number = 4, double = TRUE, thickness = 0.05, pos = 1 , margin = 0.25, cex = 0.6)
    

  legend <- c("September Research Vessel survey", "Northumberland Strait survey", "Magdalen Islands survey")
  
  cols2 <- c("orangered2", "darkgreen",  "darkblue")
  
  legend(x= "topright",
         legend = legend,
         col = cols2,
         bg =c("white"),
         cex = 1.5,
         pch=c(17, 18, 19),
         ncol=1)
  
  
  ### add an inset
  
  ##get a base map
  map_data <- rnaturalearth::ne_countries(
    scale = "medium",
    returnclass = "sf", continent = "north america")
  
  # calculate position of inset
  plotdim <- par("plt")
  xleft    = plotdim[2] *0.022
  xright   = plotdim[2]  *0.3
  ybottom  = plotdim[4]  *0.03
  ytop     = plotdim[4]  *0.31#
  
  xlim_map <- c(-131,-56)
  ylim_map <- c(22, 85)
  
  # set position for inset
  par(
    fig = c(xleft, xright, ybottom, ytop)
    , mar=c(1.5,1.5,1.5,1.5)+0.1##reduce the white area around the map
    , new=TRUE
  )

  plot(map_data[1], col= "grey80", border="grey15" , main="" , xlim= xlim_map, ylim=ylim_map, lwd=0.5, bg="white")
  
  par(
    fig = c(xleft, xright, ybottom, ytop)
  )

  ##center of study area
  
  xc <- (xlim[2]- xlim[1] )/2 + xlim[1]
  
  yc <- (ylim[2]- ylim[1] )/2 + ylim [1]
  
  points (-62.1, 46.7, col ="red", bg="red", pch=23, cex=1) ###the coordinates for xlim and ylim are in the units of the main map, not the inset, which makes no f'n sense
  
  
 dev.off()

 
