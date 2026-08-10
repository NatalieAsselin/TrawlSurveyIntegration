##Script to compare biomass prediction to landings and make figure

set.seed <- 787

library(Ecfun)

## Set paths
main.fp <- getwd()

fig.path <- file.path(main.fp, "results/figures")
data.path <- file.path(main.fp, "data")
results.path <- file.path(main.fp, "results/tables")

path <- paste0(main.fp, "/programs/")
rdatapath <- paste0(main.fp, "/programs/RData/")

label <- "commercial biomass_ALL"

model.lab <- "m13" ##to label maps

years <- c(2001:2024)

bio.est <- read.csv ( file = paste0(results.path, "/sGSL_", model.lab, "_CommercialBiomasssGSL_ 1000 tons_2001_2025_1Kres.csv"))

land <- read.csv( file = paste0(data.path, "/LandingsByRegion_1968_2024.csv"))

land <- land[,c("Year", "TOTAL")]

colnames(land) <- tolower(colnames(land)) 

a <- which(land$year %in% years)
land <- land[a,]

colnames (land) <- c("year", "landings_T")

land$landings_1000T <- land$landings_T/1000

land$year <- land$year - 1 ##bring landings back one year to match year of biomass estimate (i.e., predicting 2004 landings with 2003 biomass)

bio.est <- merge(bio.est, land)

  sim.res <- NULL
  sim.res$area <- "total"
  
 dn <- expand.grid(sim = seq(0, max(bio.est$upr *1.05),  by = 0.1)) ###data used in simulation of model parameter error
  
 sim.res<- cbind(sim.res, dn)
  
  ##set up inverse
  var.p <- as.data.frame(1/dn$sim)

  sims <-1000
  
  for (j in 1:sims){
  
  y <- bio.est
  
  b <- rnorm(nrow(bio.est), mean=bio.est$log_est, sd=bio.est$se)
  
  y$sim.log <- b
  
 y$sim <- (exp(y$sim.log)) / 1000000 ##back on response scale

 y$var <- 1/y$sim
 
 a  <- glm(landings_1000T~var,data=y,family=gaussian(link='inverse'), control=list(maxit=100)) ##from Adam ##added maxit to try more iterations
 
  d <-simulate.glm(a, newdata=var.p) 
  
  d$sim <- d$response ##flip back
   sim.res <- cbind(sim.res, as.data.frame(d$sim))

  }
  

mean <-apply (sim.res [2:sims+2], 
                1, mean) 
 
min <-apply (sim.res [2:sims+2],
               1, quantile, p=0.025) 

max <-apply (sim.res [2:sims+2], 
              1, quantile, p=0.975) 
 
results <- cbind(dn, mean,min, max)
  
write.csv(results, file = paste0(results.path, "/sGSL_", label, "_Landings_Biomass_Comp_SimResults",model.lab, ".csv"), row.names = FALSE)
  
save(results, bio.est, years, file = paste0(rdatapath, "sGSL_", label, "_Landings_Biomass_Comp_SimResults",model.lab,  ".rdata"))
  
 