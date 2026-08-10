##Function to get bottom temperature in two steps
## 1- interpolate temperatures at locations
## 2- use nearest neighbour for NA values from step 1 (i.e., points on edge of spatial domain of temperature dataset

##modified from Tobie Surette gulf.data/R/temperature.R
#https://github.com/TobieSurette/gulf.data/blob/c3937441e0b40f424055da91396fc18f0598583f/R/temperature.R


temperatureNN <- function(x, year, month = "september", longitude, latitude, depth, radius.search = 0, polygon, polygon.longitude, polygon.latitude, ...){
  
  # Parse input arguments:
  if (!missing(x) & missing(year)) if (is.numeric(x)) year <-x
  if (missing(year)) stop("'year' must be specified.")
  if (length(year) != 1) stop("'year' must be a single value.")
  if (is.numeric(month)){
    month <- month[round(month) %in% 1:12]
    month <- month[!is.na(month)]
    if (length(month) > 0) month <- month.name[month]
  } 
  if (length(month) > 1) stop("'month' must be a single value.")
  month <- match.arg(tolower(month), c("june", "september"))
  
  
  # Locate data file:
  file <- locate(keywords = c("water", "temperature", month, year), package = "gulf.data")
  if (length(file) == 1){
    load(file)
  }else{
    stop("Unable to locate temperature data file.")
  }
  
  if (!missing(depth)){
    if (is.character(depth)){
      depth <- match.arg(depth, c("surface", "bottom"))
      
      # Surface temperatures:
      if (depth == "surface") depth <- 0:round(radius.search)
      
      # Bottom temperatures:
      if (depth == "bottom"){
        fun <- function(x){
          ix <- which(!is.na(x))
          if (length(ix) > 1) return(x[max(ix)]) else return(NA) 
        } 
        x <- apply(x, 1:2, fun)
        tmp <- dimnames(x)
        dim(x) <- c(dim(x), 1)
        dimnames(x) <- list(longitude = tmp[[1]], latitude = tmp[[2]], depth = NULL)
      }
    }
    
    # Look up numeric depths:
    if (is.numeric(depth)){
      depth <- as.character(round(depth))
      depth <- depth[which(as.character(depth) %in% dimnames(x)[[3]])]
      x <- x[,,depth, drop = FALSE]
    }
    
  }
  
  ##Function to get nearest but not take NA
  get_nearest <- function(x, y, z, x0, y0, max_dist = 5) {
    result <- numeric(length(x0))
    for (k in seq_along(x0)) {
      ix <- which.min(abs(x - x0[k]))
      iy <- which.min(abs(y - y0[k]))
      
      found <- FALSE
      for (r in 0:max_dist) { # Search box of increasing radius
        idxs <- expand.grid(dx = -r:r, dy = -r:r)
        for (j in seq_len(nrow(idxs))) {
          ix1 <- ix + idxs$dx[j]
          iy1 <- iy + idxs$dy[j]
          if (ix1 >= 1 && ix1 <= length(x) && iy1 >= 1 && iy1 <= length(y)) {
            val <- z[ix1, iy1]
            if (!is.na(val)) {
              result[k] <- val
              found <- TRUE
              break
            }
          }
        }
        if (found) break
      }
      if (!found) result[k] <- NA # All searched neighbourhoods are NA
    }
    return(result)
  }
  
  
  # Interpolate temperature data at specified coordinates:
  if (!missing(longitude) & !missing(latitude)){
    if (length(longitude) != length(latitude)) stop("'longitude' and 'latitude' must be the same length.")
    
    lons <- as.numeric(dimnames(x)[[1]])
    lats <- as.numeric(dimnames(x)[[2]])
    res <- matrix(NA, nrow = length(longitude), ncol = dim(x)[3])
    for (i in 1:dim(x)[3]){
      zz <- x[,,i]
    
      zz[is.na(zz)] <- -1000000  # Because function does not take NA values. ##from tobie
     
      res[,i] <- akima::bilinear(x = lons, y = lats, z = zz, x0 = longitude, y0 = latitude)$z ##original function, interpolates bottom temp based on nearby points
      res[res < -2] <- NA ##these are the values from the -1000000 zz values
      
      ix <- is.na(res)
      na.lons <- longitude[ix]
      na.lats <- latitude[ix]
      na.res <- matrix(NA, nrow = length(na.lons), ncol = dim(x)[3])
      
      zz[zz==-1000000] <- NA##put these back to na
      na.res[,i] <- get_nearest(x = lons, y = lats, z = zz, x0 = na.lons, y0 = na.lats) ##get nearest neighbour from original data for missing values
      res[ix] <- na.res

      }
   # res[res < -2] <- NA
    dimnames(res) <- list(NULL,  depth = dimnames(x)[[3]])
    x <- res
  }
  
  if (ncol(x) == 1) x <- x[, 1]
  
  return(x)   
}

