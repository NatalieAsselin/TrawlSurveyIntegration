##Script to determine a study area for magdalen islands based on buffering points

library(sp)
library(ptools)

###set up paths
main.fp <- getwd()

path <- paste0(main.fp, "/programs/")
data.path <- file.path(main.fp, "data")

rdatapath <- paste0(main.fp, "/programs/RData/")

load(paste0(data.path, "/MI_Lobster data.rdata"))

mi.sp <- setMI [, c("X", "Y")]


##multiply coordinates by 1000 to be in actual numbers (need this to work in grid script)

mi.sp$X <- mi.sp$X *1000
mi.sp$Y <- mi.sp$Y *1000

mi.sp$ID <- 1 ##animal ID, all the same in this case

mi.sp <- as.data.frame(mi.sp)

coordinates(mi.sp)<-mi.sp[c("X", "Y")]

proj4string(mi.sp) <- CRS( "EPSG:32620" )  #set crs to UTM 20N

mi.buffer <- ptools::buff_sp(mi.sp , 4000)##make a buffer big enough that points touch

mi.buffer$id <-  "1"

mi.mcp <- mi.buffer ##renaming to keep things consistent with other files
save(mi.mcp,   file = paste0(rdatapath, "MI_StudyArea",  ".rdata"))



