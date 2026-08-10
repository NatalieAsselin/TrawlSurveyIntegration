###Script for figure comparing landings to biomass estimates
##Run Comparison to Landings.R first

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
load(file = paste0(rdatapath, "sGSL_", label, "_Landings_Biomass_Comp_SimResults",model.lab,  ".rdata"))

##width and height in inches
width <- 8.84/ 2.54 # 

height <- 8.84 /2.54  

grDevices::pdf(file = paste0(fig.path, "/sGSL_", model.lab, "_", label,"_Compare to landings"," _", min(years), "_", max(years),  ".pdf"),  height = height, width = width, pointsize=7)

cols <- colorRampPalette(c("blue", "mediumturquoise", "yellow", "orange", "red"))

bio.est$pt.color <- cols(diff(range(bio.est$year))+1)[bio.est$year-min(years)+1]


xlim <- c(0, max(bio.est$upr )*1.05)
ylim <- c(0, max(bio.est$landings_1000T)*1.10)


plot(xlim, ylim, type = "n", xlab = "", ylab = "", xaxs = "i", yaxs = "i")

grid()

polygon(c(results$sim ,rev(results$sim)),c(results$min ,rev(results$max)),col = "grey75", border = FALSE)

lines(results$sim, results$mean)
lines(results$sim, results$min, col="grey50")
lines(results$sim, results$max, col="grey50")

for (k in 1:nrow(bio.est)){
  lines( c(bio.est$lwr[k], bio.est$upr[k]), rep(bio.est$landings_1000T[k], 2), col = bio.est$pt.color[k],  lwd = 1)
  lines(rep(bio.est$lwr[k], 2),  c(bio.est$landings_1000T[k] - (0.025*ylim[2]), bio.est$landings_1000T[k] + (0.025*ylim[2])) , col = bio.est$pt.color[k],  lwd = 1)
  lines(rep(bio.est$upr[k], 2),  c(bio.est$landings_1000T[k] - (0.025*ylim[2]), bio.est$landings_1000T[k] + (0.025*ylim[2])) , col = bio.est$pt.color[k],  lwd = 1)
}

points( bio.est$est, bio.est$landings_1000T, pch=20,col = bio.est$pt.color , cex=1.5)

text(x=bio.est$est[bio.est$year %in% c(2020)], y=bio.est$landings_1000T[bio.est$year %in% c(2020)], labels=bio.est$year[bio.est$year %in% c(2020)], cex=0.65, pos=3) #

legend("bottomright", legend = (min(bio.est$year):max(bio.est$year)), lwd = 2, col = (cols(diff(range(bio.est$year))+1)), ncol=2)
mtext("Landings (1000s t)", 2, 2.5, cex = 1.25)
mtext("Biomass (1000s t)", 1, 2.5, cex = 1.25)

box()

dev.off()


