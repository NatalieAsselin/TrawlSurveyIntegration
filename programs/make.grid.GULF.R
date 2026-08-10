###Function to make grids for predictions in sdmTMB

make.grid.GULF <- function (res =1, xlim = c(280, 720) , ylim = c(5055, 5432), year,  mean.depth, sd.depth, mean.temp, sd.temp, temp.month){
  
   # Interpolation grid:
   if (length(res) == 1) res <- rep(res, 2)
   dx <- res[1]
   dy <- res[2]
   xgrid <- seq(xlim[1], xlim[2], by = dx)
   ygrid <- seq(ylim[1], ylim[2], by = dy)
  
  grid <- expand.grid(X = xgrid, Y  = ygrid, year = years)
   
  grid <- st_as_sf(grid, coords=c("X", "Y")) ##convert the grid to an sf object 
  grid <- st_set_crs(grid, 32620) #set crs to UTM 20N
 
  grid$X <- as.numeric( lapply(grid$geometry,`[[`, 1))
  grid$Y <- as.numeric(lapply(grid$geometry,`[[`, 2))
   
   ##add lat and long in decimal degrees
   x <- grid
   x$geometry <- x$geometry * 1000
   x <- st_set_crs(x, 32620) ##re-set crs to UTM 20N
   st_crs(x)
   x <- sf::st_transform(x, 4326)
   
   x$lon <- as.numeric( lapply(x$geometry,`[[`, 1))
   x$lat <- as.numeric(lapply(x$geometry,`[[`, 2))
   
  
   ##add LFA
  x$lfa <- NA 
  x$lfa <- gulf::fishing.zone(longitude=x$lon, latitude=x$lat, species = 2550) ##sGSL
  x$lfa <- ifelse(is.na(x$lfa), gulf.data::fishing.zone(longitude=x$lon, latitude=x$lat, species = 2550, region= "quebec"), x$lfa)
  
  
  ##put all of 23 together as one LFA
  x$lfa <- ifelse(x$lfa %in% c("23A", "23B", "23C", "23D"), "23", x$lfa)
  
  ##define study areas
  x$rv.stratum <- gulf::stratum(x$lon, x$lat, region = "gulf", survey="rv")
  x$rv.study.area <- ifelse(is.na(x$rv.stratum), FALSE, TRUE)
  
  x$ns.study.area <- ifelse(x$lfa %in% c("25", "26A") & (x$X  < 584), TRUE, FALSE) ##take out part of LFA 26A at western end - no samples
  
  ##bring in MI study area, from MCP analysis
  load(paste0(rdatapath, "MI_StudyArea",  ".rdata"))
  
 mi.mcp <- st_as_sf(mi.mcp, coords=c("X", "Y")) ##convert the study area to an sf object 
  
  mi.mcp <- sf::st_transform(mi.mcp, 4326)
  
  x <-  x %>%
    st_join(mi.mcp)
  
  x$mi.study.area <- ifelse(x$id %in% "1", TRUE, FALSE) 
  
   x <- as.data.frame(x)

   grid <- cbind(grid,x [c("lon","lat", "lfa", "rv.study.area", "ns.study.area", "mi.study.area")])
   
   grid$depth <- gulf::depth(grid$lon, grid$lat) 
   ##flip the sign of depth so that depth are positive
   grid$depth <- -1 *grid$depth
   
   grid$log.depth <- log(grid$depth)
   
 
   grid$temp.month <- temp.month
   
   grid <- get.bottom.temp(file=grid, years=years, month=temp.month)
   
   grid <- get.substrate.grid(file=grid, years=years)
   
   ##define prediction area
   grid$ix <- FALSE
   
   if (label == "commercial biomass_RV") { 
     
     grid$ix[which( (grid$depth >= 15)& (grid$depth <=max.depth) & !is.na(grid$depth))] <- TRUE 
     
   } else if (label == "commercial biomass_ALL"){ 
     
    grid$ix[which( (grid$depth >= 15)& (grid$depth <=max.depth) & !is.na(grid$depth))] <- TRUE 
   
    grid$ix[which( (grid$depth >= 5) & (grid$depth <=max.depth) & !is.na(grid$depth) & (grid$ns.study.area))] <- TRUE 
    
    grid$ix[which( (grid$depth >= 7) & (grid$depth <=max.depth) & !is.na(grid$depth) & (grid$mi.study.area))] <- TRUE 
    
   } else if (label == "commercial biomass_NS"){
     
     grid$ix[which( (grid$depth >= 5) & (grid$depth <=max.depth) & !is.na(grid$depth) & (grid$ns.study.area))] <- TRUE 
    
   } else if (label == "commercial biomass_MI"){
     
     grid$ix[which( (grid$depth >= 7) & (grid$depth <=max.depth) & !is.na(grid$depth) & (grid$mi.study.area))] <- TRUE 
  
}

   ##scale the depth and temp by mean and sd from integrated dataset
   
   grid$depth.s <- (grid$depth - mean.depth)/sd.depth
   grid$depth.s2 <- grid$depth.s ^ 2
   
   grid$bottom.temp.s <- (grid$bottom.temp - mean.temp)/sd.temp
   grid$bottom.temp.s2 <- grid$bottom.temp.s ^ 2
   
  grid$ix [which (grid$substrate == "glacial_deposits")] <- FALSE
   grid$ix [which(is.na(grid$bottom.temp))] <-FALSE ##take out any grid squares with FALSE temps. 

   return(grid)
}
