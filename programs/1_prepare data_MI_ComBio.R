###File to prep data from Magdalen Islands trawl survey
#
##Data are from: https://ouvert.canada.ca/data/fr/dataset/4b595b87-ad3f-438c-8310-c3d545be1aee, accessed Feb 3, 2025

###set up paths
main.fp <- getwd()

path <- paste0(main.fp, "/programs/")
fig.path <- file.path(main.fp, "results/figures")
data.path <- file.path(main.fp, "data")
results.path <- file.path(main.fp, "results/tables")

##Load packages from https://github.com/TobieSurette 

library(gulf.data)
library(gulf.spatial)

##Load gulf package (DFO package to access and manipulate data)
library(gulf)

##Load other libraries
library(dplyr)
 
source(paste0(path, "get.substrate.R"))
source(paste0(path, "/BIO_functions.r"))
source(paste0(path, "/temperature_near.neigh.r"))
source(paste0(path, "get.bottom.temp.R"))


# Define survey years:
years <- c(2001:2022) 

# Load set data:

set <- read.csv(file = paste0(data.path, "/Data_Quebec/Information_stations.csv")) ##

set$year <- set$an 

set <- set[set$year %in% years, ]

set$survey <- "MI"

set$vessel.code <- set$bat

set <- base::sort(set, by= c("date", "heure.d"))

set$cruise.number <- 1

set <- set %>% 
    group_by(year) %>% 
      mutate(set.number = row_number())

##add in lat and long - midpoint of start and end
set$lat <- (set$lat.ddec + set$lat.fdec)/2
set$lon <- (set$long.ddec + set$long.fdec)/2

##put in start lat and long for sets without end lats and longs
set$lat <- ifelse(is.na(set$lat.fdec), set$lat.ddec, set$lat)
set$lon <- ifelse(is.na(set$long.fdec), set$long.ddec, set$lon)

##put in end lat and long for sets without end lats and longs
set$lat <- ifelse(is.na(set$lat.ddec), set$lat.fdec, set$lat)
set$lon <- ifelse(is.na(set$long.ddec), set$long.fdec, set$lon)

##make longitude -
set$lon <- -set$lon

##delete sets without positions as it seems they might not have been fished
set <- set [!is.na(set$lat), ]

# assign LFA:

set$lfa <- gulf.data::fishing.zone(longitude=set$lon, latitude=set$lat, species = 2550, region= "quebec")

mls <- read.csv(paste0(data.path, "/MLS_Quebec.csv"))

##Fix na values - put in 500 as no lobster ever that size. will work for min and max

mls$min.fem <- ifelse(is.na(mls$min.fem), 500, mls$min.fem) 
mls$max.fem <- ifelse(is.na(mls$max.fem), 500, mls$max.fem)
mls$max.size <- ifelse(is.na(mls$max.size), 500, mls$max.size)

mls$min.fem <- as.numeric(mls$min.fem)
mls$max.fem <- as.numeric(mls$max.fem)
mls$max.size <- as.numeric (mls$max.size)



##calculate average in lfas with more than one MLS per year (e.g. 26A)

tmp <- aggregate(mls[c("mls", "min.fem", "max.fem", "max.size")], by = mls[c("year", "lfa")], mean) 

##put in MLS 
ix <- gulf.utils::match(set[c("year", "lfa")], tmp[c("year", "lfa")])
set$mls <- tmp$mls[ix]
set$min.fem <- tmp$min.fem[ix]
set$max.fem <- tmp$max.fem[ix]
set$max.size <- tmp$max.size[ix]


set$depth <- round(gulf::depth(set$lon, set$lat)) 
##flip the sign of depth so that depth are positive
set$depth <- -1 *set$depth


##add in average swept area (surf) by boat for sets without swept.area values

tmp <- aggregate (set ["surf"], by = set["bat"], mean, na.rm=TRUE)

set$surf <- ifelse(is.na(set$surf) & set$bat == "C", tmp [1,2] , set$surf)
set$surf <- ifelse(is.na(set$surf) & set$bat == "L", tmp [2,2] , set$surf)
set$surf <- ifelse(is.na(set$surf) & set$bat == "M", tmp [3,2] , set$surf)

set$swept.area <- set$surf

# Log swept area:
set$log.swept.area <- log(set$swept.area)

# Calculate UTM coordinates:
tmp <- gulf.spatial::deg2km(set$lon, set$lat)
set$X <- tmp$x
set$Y  <- tmp$y

##add in bottom temp
set$temp.month <- "september"
month <- c( "september")
set <- get.bottom.temp(file=set, years=years, month=month)

##add in substrate

set <- get.substrate(file=set, years=years)

# Load in biological data:

len <- read.csv(file = paste0(data.path, "/Data_Quebec/Données_homard_lobster_survey.csv")) ##

len$year <- year(len$date)

len$size <- len$taille

##rounddown lengths to smaller mm

len$size <- floor(len$size)

## add mls

mls <- mls[mls$lfa==22,]

len <- merge(len, mls)

##key to identify individual sets
vars <- c("year", "st",  "tr")

by1 <- len$year
by2 <- len$st
by3 <- len$tr

### Identify berried females to remove from commercial

len$p.berried <- ifelse(len$sex %in% c(1,3),  0.5 / (1 + exp(-(-16.94+ 0.239 * as.numeric(len$size)))), 0) 
fem.berried <- aggregate (len$p.berried, by = list(by1,by2,by3), sum)
names (fem.berried) <- c(vars, "berried")
set <- merge (set, fem.berried, by=vars, all.x=TRUE)
set$berried <- ifelse(is.na(set$berried), 0, set$berried)


##commercial lobsters above mls, below max size and not berried
len$commercial.all <- ifelse(len$size >= len$mls & len$size < len$max.size , 1, 0)
len$commercial <- ifelse(len$commercial.all == 1, len$commercial.all - len$p.berried, len$commercial.all)


### commercial weight
len$weight <- NA
len$weight <- ifelse(len$sex==2, (0.0006 * len$size ^ 3.0782) / 1000, len$weight) ##males
len$weight <- ifelse(len$sex %in% c(1,3,9,12), (0.0013 * len$size ^ 2.8822) / 1000, len$weight) ##females (1,3) and all others as female weight is lower than male by formula
len$comm.weight <- len$commercial * len$weight
lob.comm.weight <- aggregate (len$comm.weight, by = list(by1,by2,by3), sum)
names (lob.comm.weight) <- c(vars, "commercial.weight")
set <- merge (set, lob.comm.weight, by=vars, all.x=TRUE)
set$commercial.weight <- ifelse(is.na(set$commercial.weight), 0, set$commercial.weight)


# Set year as numeric:
set$year <- as.numeric(as.character(set$year))

set <- set[set$depth >= 0,]

set <- set[set$depth <= 100,]

set <- set[set$substrate != "glacial_deposits",]

##make the folds for model validation
k_folds <-10
set$IDS = "I"
set = cv_SpaceTimeFolds(set,idCol = 'IDS',k_folds=k_folds)


setMI <- set [, c("year", "survey", "lfa", "vessel.code", "cruise.number", "set.number", "lon", "lat","X", "Y",  "swept.area"  , "log.swept.area",
                 "depth", 
                 "bottom.temp",
                 "substrate" ,
                 "fold_id",
                 "commercial.weight")  ] 


save(setMI, file = paste0(data.path, "/MI_Lobster data.rdata"))

