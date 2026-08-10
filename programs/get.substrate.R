##Function to add substrate to datapoints

##Based on data from 
##Loring and Nota, 1973, https://gcgeo.gc.ca/geonetwork/metadata/eng/8c269a91-d3a2-4f49-943d-6b2401c42cba, accessed April 10, 2025

get.substrate <- function (file, years) {

library(sf)
library(terra)
library(sp)
  
dat<- file
  
sub.path <- paste0(data.path, "/Substrat")
  
sub <- sf::read_sf(dsn = paste0(sub.path, "/Seafloor_SubstratBenthique_Loring"), layer="Seafloor_SubstratBenthique") 

proj.sub <- crs(sub) #identify projection of shapefile
  
  # Project our survey data coordinates in same projection as substrate data:
  
  s <- dat %>% dplyr::select(lon, lat, year, vessel.code, cruise.number, set.number) %>%
    sf::st_as_sf(crs = 4326, coords = c("lon", "lat")) %>%
    sf::st_transform(proj.sub)


 sel <-  sf::st_join(s,sub, left=TRUE, st_nearest_feature) 
 
  sel <- unique(sel[,c('year','vessel.code','cruise.number', "set.number", "DEPOT_EN")]) 

  names(sel) <- gsub("DEPOT_EN", "substrate",   names(sel))
  
  sel$substrate <-ifelse(sel$substrate %in% c("Sand"),
                         "sand",
                         sel$substrate)
  
  
  sel$substrate <-ifelse(sel$substrate %in% c("Pelite-sand",
                                              "Pelite" ),
                         "pelite",
                         sel$substrate)
  
  sel$substrate <-ifelse(sel$substrate %in% c("Gravel", 
                                              "Gravel-sand",
                                              "Sand-Gravier"),
                         "gravel",
                         sel$substrate)
  
  sel$substrate <-ifelse(sel$substrate %in% c("Glacial deposits" ),
                         "glacial_deposits",
                         sel$substrate)
  
  
  sel$substrate <-ifelse(is.na(sel$substrate),
                         "unknown",
                         sel$substrate)
  
  y <- merge(file, sel)
  return(y)
}
