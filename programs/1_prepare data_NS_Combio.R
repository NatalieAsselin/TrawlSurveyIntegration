###File to prep data from Northumberland Strait trawl survey
#

graphics.off()

## Set up paths

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


source(paste0(path, "get.substrate.R"))
source(paste0(path, "/BIO_functions.r"))
source(paste0(path, "/temperature_near.neigh.r"))
source(paste0(path, "get.bottom.temp.R"))
 
###Put in password to access DFO internal database

years <- c(2001:2025)

cruise<-c(24,139,241,341,434,536,637,30,22,29,103,129,26,126,21,521,18,24,41,901, 140, 150,151,152,303, 402, 502, 
          103, 129, 
          329, 829, 953 
          )

# Load set data:
set <- gulf.data::read.gulf.set(year = years, species = 2550, password = password, survey = "ns", source = "ptran") ##have to be connected to VPN

set$cruise.number <- as.numeric(set$cruise.number)
set <- set[which(set$cruise.number %in% cruise), ]
set <- set[set$experiment != 3, ]
set$year <- year(set) 
a <- which(!is.na(set$year))
set <- set[a, ]

set$lat <- lat(set)
set$lon <- lon(set)


# Change gear variable:
set$gear.code <- set$gear

##corrections to gear codes from the small comparative

set$gear.code <- ifelse(set$gear.code == 9, 16, set$gear.code)
set$gear.code <- ifelse(set$gear.code == 13, 16, set$gear.code)


set$gear[set$gear.code == 16] <- "OT"
set$gear[set$gear.code == 19] <- "NT"
set$gear[set$gear.code == 14] <- "NephT"

# Identify LFA:
set$lfa <- gulf::fishing.zone(lon(set), lat(set), species = 2550) 

set$lfa <- ifelse(is.na(set$lfa), gulf.data::fishing.zone(longitude=set$lon, latitude=set$lat, species = 2550, region= "quebec"), set$lfa)

set$lfa <- ifelse((set$set.number == 31 & set$year == 2008), "25", set$lfa) ##mapped it out and it is very close to land inLFA25

set$lfa <- ifelse(set$lfa == "20B" , "20", set$lfa) 


#put in survey identifiers
set$survey <- NA
set$survey <- ifelse(set$gear.code %in% c(16,19) & set$lfa %in% c("25"), "NS_Fall", set$survey)
set$survey <- ifelse(set$gear.code %in% c(16,19) & set$lfa %in% c("20", "23A", "23B", "23C", "23D", "24", "26A", "26B"), "NS_Spring", set$survey)


##seperating the nephrops sets into spring and fall
set$survey <- ifelse(set$gear.code %in% c(14) & set$lfa %in% c("25"), "NS_Neph_Fall", set$survey)
set$survey <- ifelse(set$gear.code %in% c(14) & set$lfa %in% c("20", "23A", "23B", "23C", "23D", "24", "26A", "26B"), "NS_Neph_Spring", set$survey)


# Clean up station number:
set$station.number <- gsub(" ", "", set$station.number)
set$station.number <- gsub("^0+", "", set$station.number)

##bring in MLS file with female size restrictions
mls <- read.csv(paste0(data.path,"/mls_lobster_femSize.csv"))

#make fishing zones for 23 match lfa file
mls$lfa <- ifelse(mls$lfa == "23" & mls$sub.lfa=="A", "23A", mls$lfa)
mls$lfa <- ifelse(mls$lfa == "23" & mls$sub.lfa =="B", "23B", mls$lfa)
mls$lfa <- ifelse(mls$lfa== "23" & mls$sub.lfa =="C", "23C", mls$lfa)
mls$lfa <- ifelse(mls$lfa == "23" & mls$sub.lfa =="D", "23D", mls$lfa)

mls <- mls[c("lfa", "year", "mls", "min.fem", "max.fem", "max.size")]

##bring in Quebec file for stations in LFA 20
mls.pq <- read.csv(paste0(data.path,"/MLS_Quebec.csv"))

##merge the two files

mls <- rbind(mls, mls.pq)

##fix na values - put in 500 as no lobster over that size. will work for min and max

mls$min.fem <- ifelse(mls$min.fem=="n/a", NA, mls$min.fem) ##spelling
mls$max.fem <- ifelse(mls$max.fem=="n/a", NA, mls$max.fem)

mls$min.fem <- ifelse(is.na(mls$min.fem), 500, mls$min.fem) ##spelling
mls$max.fem <- ifelse(is.na(mls$max.fem), 500, mls$max.fem)
mls$max.size <- ifelse(is.na(mls$max.size), 500, mls$max.size)


mls$min.fem <- as.numeric(mls$min.fem)
mls$max.fem <- as.numeric(mls$max.fem)
mls$max.size <- as.numeric(mls$max.size)


##make a set$lfa_mls column to match with MLS file- mls prior to 2007 for 23 is one number, not 4

set$lfa<- ifelse(set$lfa %in% c("23A", "23B", "23C", "23D") & set$year < 2008, "23", set$lfa)

##calculate average in lfas with more than one MLS per year (e.g. 26A)

tmp <- aggregate(mls[c("mls", "min.fem", "max.fem", "max.size")], by = mls[c("year", "lfa")], mean) ##can use mean for all of them as fem size are throughout

##Add in MLS 
ix <- gulf.utils::match(set[c("year", "lfa")], tmp[c("year", "lfa")])
set$mls <- tmp$mls[ix]
set$min.fem <- tmp$min.fem[ix]
set$max.fem <- tmp$max.fem[ix]
set$max.size <- tmp$max.size[ix]

# # Fix depth variable from set card, to use for wingspread and calibration, as done in comparative:
 set$depth.recorded <- apply(set[c("depth.start", "depth.end")], 1, mean, na.rm = TRUE)
 set$depth.recorded[set$vessel.code == "O"] <- set$depth.recorded[set$vessel.code == "O"] + 3 ##correct for draft
 set$depth.recorded[set$vessel.code == "P"] <- set$depth.recorded[set$vessel.code == "P"] + 3 ##correct for draft
 ix <- which(is.na(set$depth.recorded))
 set$depth.recorded[ix] <- round(gulf.spatial::depth(lon(set[ix, ]), lat(set[ix, ]))) ##put in value from function if missin


set$depth <- round(gulf::depth(set$lon, set$lat)) ##gulf.package depth

set$depth <- -1 *set$depth

# Fix tow distance variable:
ix <- which(is.na(set$distance) & !is.na(set$longitude.end))
if(length(ix)>0){
## calculate distance from diagonal in naut. miles
  set$distance[ix] <- 0.539957 * diag(gulf.spatial::distance(set$longitude.start[ix], set$latitude.start[ix], set$longitude.end[ix], set$latitude.end[ix]))
}
set$distance[which(is.na(set$distance) & set$gear == "OT")] <- 0.625


##spot correction of one set over 10K long, 15 minutes, 2.5Kts
set$distance <- ifelse(set$cruise.number == 829 & set$set.number == 9, 0.625, set$distance)

##add in wingspread measurements for those that were measured directly
set$wingspread <- NA

# Calculate wing spread: ##for all sets - use recorded depth to be consistent with how formula originally calculated - Asselin et al. 2023 - comparative survey
set$wingspread[which(set$gear == "OT")] <- 0.093 * set$depth.recorded[which(set$gear == "OT")] + 6.269
set$wingspread[which(set$gear == "NT")] <- 0.070 * set$depth.recorded[which(set$gear == "NT")] + 7.294


##put in wingspread for Nephrops from Rondeau et al. 2014 (rock crab paper, 5.9 m)
set$wingspread[which(set$gear == "NephT")] <- 5.9


##put in actual winspread measurements from Notus sensors for when we have them 

 tmp <- read.csv(paste0(data.path,"/NSSurvey_WingspreadEstimates-2019to2023.csv"))
 tmp2 <- read.csv(paste0(data.path,"/NSSurvey_WingspreadEstimates-2024.csv"))
 tmp <- rbind(tmp, tmp2)
 
colnames(tmp) [6] <- "wingspread"
 vars <- c("year", "vessel.code", "cruise.number", "set.number")
 ix <- gulf.utils::match(tmp[vars], set[vars])
 tmp <- tmp[!is.na(ix), ]   
 ix <- ix[!is.na(ix)]
 set$wingspread[ix] <- tmp$wingspread 
 
   # Calculate swept area:
set$swept.area <- NA
set$swept.area <- 1000 * (set$distance / 0.539957) * set$wingspread ##calculate swept area in m2

##add in bottom temp
set$temp.month <- "september"
month <- c( "september")
set <- get.bottom.temp(file=set, years=years, month=month)

##add in substrate
 set <- get.substrate(file=set, years=years)
 
# Load length sampling data:
len <- gulf.data::read.gulf.len(year = years, species = 2550, password = password, survey = "ns")
ix <- gulf.utils::match(len[gulf.metadata::key(set)], set[gulf.metadata::key(set)])
len <- len[!is.na(ix), ]

# Log swept area:
set$log.swept.area <- log(set$swept.area)

# Calculate UTM coordinates:
tmp <- deg2km(lon(set), lat(set))
set$X <- tmp$x
set$Y  <- tmp$y

# Lobster size-frequencies:##freq translates length cards into a table 
f <- list()
f[[1]] <- gulf.data::freq(len[len$sex == 1, ])        # Total males. 
f[[2]] <- gulf.data::freq(len[len$sex == 2, ])        # Total females.
f[[3]] <- gulf.data::freq(len[len$sex %in% c(0,9), ]) # Unsexed.

fun <- function(x) return(names(x)[gsub("[0-9]", "", names(x)) == ""])
fvars <- sort(unique(unlist(lapply(f, fun))))
fvars <- fvars[order(as.numeric(fvars))]
vars <- c("year", "vessel.code", "cruise.number", "set.number")
for (i in 1:length(f)){
   f[[i]][, setdiff(fvars, fun(f[[i]]))] <- 0
   f[[i]] <- f[[i]][c(setdiff(names(f[[i]]), fvars), fvars)]
   
   tmp <- set
   import(tmp, fill = 0) <- f[[i]]
   f[[i]] <- tmp[c(vars, fvars)]
}
names(f) <- c("male", "female", "unsexed")

# Calculate berried females:
p.berried <- repvec(0.5 / (1 + exp(-(-16.94+ 0.239  * as.numeric(fvars)))), nrow = nrow(f$female))
f$berried <- f$female
f$berried[fvars] <- p.berried * f$female[fvars]

# Apply depth-size correction for lobster surveyed with OT
theta <- c('(Intercept)' = 2.9857551500, depth = 0.0346625680, size = -0.0657729990, size2 = 0.0003374518)  # Model parameters.
ix <- which(set$gear == "OT")
depth <- repvec(set$depth.recorded[ix], ncol = length(fvars)) ##use recorded depth to match methods from comparative survey
size  <- repvec(as.numeric(fvars), nrow = length(ix))
scale <- exp(theta[["(Intercept)"]] + theta[["depth"]] * depth + theta[["size"]] * size + theta[["size2"]] * size * size)
colnames(scale) <- fvars
for (j in fvars[as.numeric(fvars) > 90]) scale[, j] <- scale[, "90"]
for (i in 1:length(f)) f[[i]][ix, fvars] <- scale * f[[i]][ix, fvars]

##prep length restrictions
mls  <- repvec(set$mls, ncol = length(fvars)) 
min.fem <- repvec(set$min.fem , ncol = length(fvars)) 
max.fem <- repvec(set$max.fem , ncol = length(fvars))
max.size <- repvec(set$max.size , ncol = length(fvars))

lens <- repvec(as.numeric(fvars), nrow(set))


# Calculate total commercial lobster:
##males
ix <- lens >= mls & lens < max.size ## for males
tmp <- f$male[fvars]
set$commercial.males <- apply(ix * tmp, 1, sum) # Set non-legal sized catches to zero.
##females
ix <- lens >= mls & lens < max.size & (lens < min.fem | lens > max.fem) ## females above MLS and not in window
tmp <- f$female[fvars] - f$berried[fvars]
set$commercial.fem <- apply(ix  * tmp, 1, sum) # Set non-legal sized catches to zero.
set$commercial <- set$commercial.males + set$commercial.fem

# Calculate total commercial weight:
##males
ix <- lens >= mls & lens < max.size ## for males
weight.male   <- function(x) return((0.0006 * x ^ 3.0782) / 1000)
tmp <- f$male[fvars] * weight.male(lens) 
set$commercial.weight.males <- apply(ix * tmp, 1, sum) # Set non-legal sized catches to zero.

##females
ix <- lens >= mls & lens < max.size & (lens < min.fem | lens > max.fem) ## lengths of legal sized females no windows
weight.female <- function(x) return((0.0013 * x ^ 2.8822) / 1000)
tmp <- (f$female[fvars] - f$berried[fvars]) * weight.female(lens)
set$commercial.weight.fem <- apply(ix * tmp, 1, sum) # Set non-legal sized catches to zero.

set$commercial.weight <- set$commercial.weight.males + set$commercial.weight.fem

##remove sets not using in analysis before doing folds

set <- set[set$depth >= 0,]

##take out data points above 100 m - no lobster past 75 m. To check if this makes the models run better. 

set <- set[set$depth <= 100,]

##take out glacial deposits as no lobsters on glacial deposits - total of 88 sets. 

set <- set[set$substrate != "glacial_deposits",]

##add in of k folds for model validation
k_folds <-10

##make the folds for model validation
set$IDS = "I"
set = cv_SpaceTimeFolds(set,idCol = 'IDS',k_folds=k_folds)


##keep only columns needed for modelling

setNS <- set [, c("year", "survey", "lfa", "vessel.code", "cruise.number", "set.number", "lon", "lat","X", "Y",  "swept.area"  , "log.swept.area",
                  "depth", 
                  "bottom.temp",
                  "substrate" ,
                  "fold_id",
                "commercial.weight")  ]                              



save(setNS, file = paste0(data.path, "/NS_Lobster data.rdata"))


##Second dataset with only data in LFAs 25 and 26A for NS stand-alone model

set2 <- set2[set2$lfa %in% c("25", "26A"),]

##number of folds for cross validation, 90% of data to predict 10%

k_folds <- as.numeric(length(unique(set2$fold_id)))

##re-do folds since data removed

##make the folds for model validation 
set2$IDS = "I" 
set2 = cv_SpaceTimeFolds(set2,idCol = 'IDS',k_folds=k_folds)

setNS <- set2

save(setNS, file=paste0(data.path, "/NS_Lobster data_LFA25and26Aonly.rdata"))
