###File to prep data from September Research Vessel Survey 

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
 

source(paste0(path, "get.substrate.R"))
source(paste0(path, "/BIO_functions.r"))
source(paste0(path, "/temperature_near.neigh.r"))
source(paste0(path, "get.bottom.temp.R"))


# Define survey years:
years <- c(2001:2025)

# Load set data:

set <- read.card(year = years, card.type = "set", survey = "rv") 

set <- set[which(!(set$experiment %in% (c(3,9)))),] ### 3 is not a valid set, 9 is hydrography only

set <- set[which(!(set$stratum %in% (c(440,441,442)))),] ### take ou strata 440:442: not part of RV survey, very few sets

set$year <- year(set) 

##add in lat and long
set$lat <- lat(set)
set$lon <- lon(set)


# # Clean up station number:
set$station.number <- gsub(" ", "", set$station.number)
set$station.number <- gsub("^0+", "", set$station.number)

# Calculate LFA:
set$lfa <- gulf::fishing.zone(lon(set), lat(set), species = 2550) ##use gulf package to get sub-zones in 23
set$lfa <- ifelse(is.na(set$lfa), gulf.data::fishing.zone(longitude=set$lon, latitude=set$lat, species = 2550, region= "quebec"), set$lfa)
set$lfa <- ifelse(is.na(set$lfa), gulf.data::fishing.zone(longitude=set$lon, latitude=set$lat, species = 2550, region= "newfoundland"), set$lfa)


### fix up Québec LFAs to match mls file,, sub-lfas don't matter for MLS

set$lfa <- ifelse(set$lfa %in% c("17A", "17B"), "17", set$lfa)
set$lfa <- ifelse(set$lfa %in% c("19C"), "19", set$lfa)
set$lfa <- ifelse(set$lfa %in% c("20A", "20B"), "20", set$lfa)
set$lfa <- ifelse(set$lfa %in% c("21A"), "21", set$lfa)


##put 23 together prior to 2008 to match mls file
set$lfa<- ifelse(set$lfa %in% c("23A", "23B", "23C", "23D") & set$year < 2008, "23", set$lfa)

##put in 99 for areas outside LFAs 
set$lfa<- ifelse(is.na(set$lfa), 99, set$lfa) 

##bring in MLS file with female size restrictions - made from file on W drive
mls <- read.csv(paste0(data.path,"/mls_lobster_femSize.csv"))

#make fishing zones for 23 match lfa file
mls$lfa <- ifelse(mls$lfa == "23" & mls$sub.lfa=="A", "23A", mls$lfa)
mls$lfa <- ifelse(mls$lfa == "23" & mls$sub.lfa =="B", "23B", mls$lfa)
mls$lfa <- ifelse(mls$lfa== "23" & mls$sub.lfa =="C", "23C", mls$lfa)
mls$lfa <- ifelse(mls$lfa == "23" & mls$sub.lfa =="D", "23D", mls$lfa)

mls <- mls[, c("lfa", "year", "mls", "min.fem", "max.fem", "max.size")] ## max.size only in Québec but put in NAs to keep things consistent

##bring in Québec mls file

mls_PQ <- read.csv(paste0(data.path, "/MLS_Quebec.csv"))
mls_PQ <- mls_PQ[, c("lfa", "year", "mls", "min.fem", "max.fem", "max.size")] #female size restictions not in Québec but consistent file structure

mls <- rbind ( mls, mls_PQ)

##fix na values - put in 500 as no lobster ever that size. will work for min and max

mls$min.fem <- ifelse(mls$min.fem=="n/a", NA, mls$min.fem) ##spelling
mls$max.fem <- ifelse(mls$max.fem=="n/a", NA, mls$max.fem)

mls$min.fem <- ifelse(is.na(mls$min.fem), 500, mls$min.fem) ## replacing with 500 (larger than any lobster, so that math works)
mls$max.fem <- ifelse(is.na(mls$max.fem), 500, mls$max.fem)
mls$max.size <- ifelse(is.na(mls$max.size), 500, mls$max.size)

mls$min.fem <- as.numeric(mls$min.fem)
mls$max.fem <- as.numeric(mls$max.fem)
mls$max.size <- as.numeric (mls$max.size)


##calculate one mean mls by year to use when outside all LFAs (call it 99 to match above)
tmp <- aggregate(mls["mls"], by = mls[c("year")], mean)
tmp$lfa <- 99
tmp$min.fem <- 500##make these all 500 (i.e. no limit)
tmp$max.fem <- 500
tmp$max.size <- 500
tmp <- tmp[, c("lfa", "year", "mls", "min.fem", "max.fem", "max.size")]

mls <- rbind(mls, tmp)

##calculate average in lfas with more than one MLS per year (e.g. 26A)

tmp <- aggregate(mls[c("mls", "min.fem", "max.fem", "max.size")], by = mls[c("year", "lfa")], mean) ##can use mean for all of them as fem size are throughout##put in MLS - Tobie's way working

ix <- gulf.utils::match(set[c("year", "lfa")], tmp[c("year", "lfa")])
set$mls <- tmp$mls[ix]
set$min.fem <- tmp$min.fem[ix]
set$max.fem <- tmp$max.fem[ix]
set$max.size <- tmp$max.size[ix]

##add depth

set$depth <- round(gulf::depth(set$lon, set$lat)) ##gulf.package depth
set$depth <- -1 *set$depth


#Wingspread from (Hurlbut and Clay, 1990) - all data adjusted to this gear and wingspread later
set$wingspread <- 12.497 

# Calculate swept area in m2. all based on teleost, 1.76 naut.mile set (Benoît and Yin (2023)
set$swept.area <- 1000 * (1.75 / 0.539957) * set$wingspread 

# Log swept area:
set$log.swept.area <- log(set$swept.area)

# Calculate UTM coordinates:
tmp <- deg2km(lon(set), lat(set))
set$X <- tmp$x
set$Y  <- tmp$y

##add in bottom temp
set$temp.month <- "september"
month <- c( "september")
set <- get.bottom.temp(file=set, years=years, month=month)

##add in substrate

set <- get.substrate(file=set, years=years)

# Load length sampling data:

len <- read.card(year = years, card.type = "length", species = 2550, survey = "rv") 

##adjust length card to standard day tow for teleost, distance = true because correction factors include distance
len <- gulf::adjust(len, set, distance=TRUE, day.night = TRUE, vessel=TRUE, scale=TRUE)

ix <-gulf.utils:: match(len[gulf.metadata::key(set)], set[gulf.metadata::key(set)])
len <- len[!is.na(ix), ]

# Lobster size-frequencies:##freq translates length cards into a table 
f <- list()

vars <- c("year", "vessel.code",  "cruise.number", "set.number")
f[[1]] <- gulf::freq(len[len$sex == 1, ], by= vars   )     # Total males. 
f[[2]] <- gulf::freq(len[len$sex == 2, ], by=vars)        # Total females.
f[[3]] <- gulf::freq(len[len$sex %in% c(0,9), ], by=vars) # Unsexed.

fun <- function(x) return(names(x)[gsub("[0-9]", "", names(x)) == ""])
fvars <- sort(unique(unlist(lapply(f, fun))))
fvars <- fvars[order(as.numeric(fvars))]

for (i in 1:length(f)){
   f[[i]][, setdiff(fvars, fun(f[[i]]))] <- 0
   f[[i]] <- f[[i]][c(setdiff(names(f[[i]]), fvars), fvars)]
   
   tmp <- set
   import(tmp, fill = 0) <- f[[i]]
   f[[i]] <- tmp[c(vars, fvars)]
}
names(f) <- c("male", "female", "unsexed")

# Calculate berried females:
p.berried <- gulf.utils::repvec(0.5 / (1 + exp(-(-16.94+ 0.239  * as.numeric(fvars)))), nrow = nrow(f$female)) ##based on maturity curve and being berried 1/2 years
f$berried <- f$female
f$berried[fvars] <- p.berried * f$female[fvars]


##prep length restrictions
mls  <- repvec(set$mls, ncol = length(fvars)) 
min.fem <- repvec(set$min.fem , ncol = length(fvars)) 
max.fem <- repvec(set$max.fem , ncol = length(fvars))
max.size <- repvec(set$max.size , ncol = length(fvars))

lens <- gulf.utils::repvec(as.numeric(fvars), nrow(set))

##commercial lobsters
##males
ix <- lens >= mls & lens < max.size ## for males
tmp <- f$male[fvars] ##males
set$commercial.males <- apply(ix * tmp, 1, sum) # Set non-legal sized catches to zero.
##females
ix <- lens >= mls & lens < max.size & (lens < min.fem | lens > max.fem) ## lengths of legal sized females and excluding window-sized females
tmp <- f$female[fvars] - f$berried[fvars] ##females
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

set$survey <- "RV"

# Set year as numeric:
set$year <- as.numeric(as.character(set$year))

##remove sets not using in analysis before doing folds

set <- set[set$depth >= 0,]

set <- set[set$depth <= 100,]

set <- set[set$substrate != "glacial_deposits",]

##make the folds for model validation
k_folds <-10
set$IDS = "I"
set = cv_SpaceTimeFolds(set,idCol = 'IDS',k_folds=k_folds)


 setRV <- set [, c("year", "survey", "lfa", "vessel.code", "cruise.number", "set.number", "lon", "lat","X", "Y",  "swept.area"  , "log.swept.area",
                   "depth", 
                   "bottom.temp",
                   "substrate" ,
                   "fold_id",
                   "commercial.weight")  ]   


save(setRV, file = paste0(data.path, "/RV_Lobster data.rdata"))

