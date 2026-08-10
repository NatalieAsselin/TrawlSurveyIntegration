##Script to get average values of complete dataset to use to scale depth and temperature


## Set paths
main.fp <- getwd()


fig.path <- file.path(main.fp, "results/figures")
data.path <- file.path(main.fp, "data")
results.path <- file.path(main.fp, "results/tables")

path <- paste0(main.fp, "/programs/")
rdatapath <- paste0(main.fp, "/programs/RData/")


# ##load  datasets
load(paste0(data.path, "/NS_Lobster data.rdata"))
load(paste0(data.path, "/RV_Lobster data.rdata"))
load(paste0(data.path, "/MI_Lobster data.rdata"))

xlim = c(280, 720)
ylim = c(5055, 5432)

# Define dataset:
data <- rbind(setNS, setRV, setMI)

mean.depth <- mean(data$depth)
sd.depth <- sd(data$depth)

mean.temp <- mean(data$bottom.temp)
sd.temp <- sd(data$bottom.temp)

save(mean.depth, sd.depth, mean.temp, sd.temp, xlim, ylim,   file = paste0(rdatapath, "sGSL_", "GeneralValues",  ".rdata"))
