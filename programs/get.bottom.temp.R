##Function to add estimated bottom temperature to datapoints

get.bottom.temp <- function (file, years, month) {
  
  library(mgcv)
  y <- NULL
  
  for (i in 1:length(years)){
  
    for (j in 1:length(month)){ 
      
  a <- file [ file$year==years[i], ]
  
  x <- a [ a$temp.month == month [j],]
  
  x$bottom.temp <- temperatureNN( years[i], longitude=x$lon, latitude= x$lat, depth="bottom", month = month [j])
  
  y <- rbind (y, x)
    
    }
  }

  return(y)
}
