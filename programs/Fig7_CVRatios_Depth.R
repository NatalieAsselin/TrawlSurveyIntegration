library(mgcv)


## everything below is based on the main file path set here
main.fp <- getwd()

fig.path <- file.path(main.fp, "results/figures")
data.path <- file.path(main.fp, "data")
results.path <- file.path(main.fp, "results/tables")

path <- paste0(main.fp, "/programs/")
rdatapath <- paste0(main.fp, "/programs/RData/")

mod.label <- "m13"



##load in individual survey files
SRVS <- read.csv (paste0(results.path, "/sGSL_2001_2025_1Kres_RV_Predictions_", mod.label,"ALL.csv"))
SRVS$survey <- "SRVS"

SRVS <- SRVS[!is.na(SRVS$SE_ratio),]

NS <- read.csv (paste0(results.path, "/sGSL_2001_2025_1Kres_NS_Predictions_", mod.label,"ALL.csv"))
NS$survey <- "NS"

NS <- NS[!is.na(NS$SE_ratio),]

MI <- read.csv (paste0(results.path, "/sGSL_2001_2025_0.5Kres_MI_Predictions_", mod.label,"ALL.csv"))
MI$survey <- "MI"
MI <- MI[!is.na(MI$SE_ratio),]
 
mean(SRVS$pred_SE_modRV, na.rm=TRUE)
mean(SRVS$pred_SE_modALL, na.rm=TRUE)

mean(NS$pred_SE_modNS, na.rm=TRUE)
mean(NS$pred_SE_modALL, na.rm=TRUE)

mean(MI$pred_SE_modMI, na.rm=TRUE)
mean(MI$pred_SE_modALL, na.rm=TRUE)


plot(SRVS$depth, SRVS$CV_ratio, ylim=c(0,10))
plot(NS$depth, NS$SE_ratio)
plot(MI$depth, MI$SE_ratio)

###gamm

label<- c("SRVS", "NSS", "MIS")


### three panel figure of CV ratio for SRVS, NS and MI

width <- 8.4/2.54

height <- 10/2.54

xlim <- c(0,100)


#png(file = paste0(fig.path, "/sGSL_CVRatios", mod.label, ".png"), res = 500, units = "in", height = height, width = width)

grDevices::pdf(file = paste0(fig.path, "/sGSL_CVRatios", mod.label,  ".pdf"),  height = height, width = width, pointsize=7)


m <- kronecker(matrix(c(1:3), ncol = 1, byrow=TRUE), matrix(1, ncol = 5, nrow = 5)) ##larger squares for the maps

m <- rbind(0, cbind(0, m,  0), 0)
layout(m)
par(mar = c(1, 0, 0, 0))


for (i in 1:length(label) ){

data <- if (label[i] == "SRVS") {
  SRVS
}else if  (label[i] == "NSS"){
  NS
}else if  (label[i] == "MIS"){
  MI
}

data$log.cvratio <- log(data$CV_ratio)

## determine percentage of biomass in grid points with <10kg/km2

a <- data[data$pred_modALL < 10, ]
  
if (nrow(a)>0){ ##throws up error for MI as no rows predicted to have less than 10 kg/km2

by1 <- a$year

a2 <- aggregate (a$pred_modALL, by=list (by1), sum)

by1 <- data$year

colnames(a2) <- c("year", "Biomass_lessThan10kg_kg")

a3 <- aggregate (data$pred_modALL, by=list (by1), sum)

colnames(a3) <- c("year", "TotalBiomass_kg")

a4 <- merge(a2, a3)

a4$perc_low_dens <- (a4$Biomass_lessThan10kg_kg/a4$TotalBiomass_kg) *100

print(label[i])
print(max(a4$perc_low_dens))
}

###take out points where predicted biomass is <10kg/km2

#data <- data[data$pred_modALL >= 10, ]

#a <- gam(CV_ratio ~ s(depth), data = data, family = gaussian)

#a <- gam(log.cvratio ~ s(depth), data = data, family = gaussian)

a <- gam(log.cvratio ~ s(depth, k=3), data = data, family = gaussian) ##logging or the SRVS one is no good due to high values

dn <- expand.grid(depth = seq(round(min(data$depth),0), round(max(data$depth),0),  by = 0.01))
#dn <- expand.grid(depth = seq(round(min(data$depth),0), 60,  by = 0.01))
res2  = predict(a, newdata = dn, se.fit=TRUE, type="response")


#results <- data.frame(mean = exp(res2$fit),
 #                     lower.ci = (exp(res2$fit) - 1.96 * exp(res2$se.fit)),
  #                    upper.ci = (exp(res2$fit) + 1.96 * exp(res2$se.fit)))

results <- data.frame(mean = (res2$fit),
                      lower.ci = ((res2$fit) - 1.96 * (res2$se.fit)),
                      upper.ci = ((res2$fit) + 1.96 * (res2$se.fit)))


results <- cbind(dn, results)
#results$col <- ifelse(results$mean >= 1, "red", "blue")

results$col <- ifelse(results$mean >= 0, "red", "blue")

#xlim3 <- c(0, max(data$depth )*1.05)

#ylim3 <- c(0, ceiling(max(data$CV_ratio)))
#ylim3 <- c(0, max(results$mean)*1.05)
# ylim3 <- if (label[i] == "SRVS") {
#   c(0,7)
# }else if  (label[i] == "NSS"){
#   c(0,5)
# }else if  (label[i] == "MIS"){
#   c(0,2)
# }

#ylim3 <- c(-2,2)

ylim3 <- c(floor(min(data$log.cvratio)), ceiling(max(data$log.cvratio)))

###graph of cv ratio to depth

#png(file = paste0(fig.path, "/sGSL_", label[i], min(data$year), "_", max(data$year), "_CVRatio_", mod.label, ".png"), res = 500, units = "in", height = 8.5, width = 10)

plot(xlim, ylim3, type = "n", xlab = "", ylab = "", xaxs = "i", yaxs = "i", xaxt="n", axes = FALSE)
grid()


cols <- grDevices::adjustcolor(c("grey50"), alpha = 0.1)

points( data$depth, data$log.cvratio, pch=20,col = cols , cex=0.1)

axis(2, las=2) ##las=2 makes horizontal labels

#polygon(c(results$depth,rev(results$depth)),c(results$lower.ci ,rev(results$upper.ci)),col = "grey75", border = FALSE)

#lines(results$depth, results$lower.ci, col="grey50")
#lines(results$depth, results$upper.ci, col="grey50")
#lines(results$depth, results$mean, col="red", lwd=3)



points( results$depth, results$mean, pch=20,col = results$col , cex=0.4) ##plotting line as points to change colors based on position

#lines(c(0, xlim[2]),c(1,1), col="grey20", lwd=1.5, lty=2)

lines(c(0, xlim[2]),c(0,0), col="grey20", lwd=1.5, lty=2)

text(0.02*(xlim[2]-xlim[1])+xlim[1],
     0.90 * (ylim3[2]-ylim3[1]) + ylim3 [1],
     label[i], cex = 2, pos = 4)


#legend("bottomright", legend = (min(bio.est.a$year):max(bio.est.a$year)), lwd = 2, col = (cols(diff(range(bio.est.a$year))+1)), ncol=2)

if (i==2) {
mtext("Log CV ratio", 2, 2.5, cex = 1.25)
}


if (i==3) {
mtext("Depth (m)", 1, 2.8, cex = 1.25)
  
  ax.lab <- sequence(nvec=11, from=xlim[1], to=xlim[2], by=10)
  
  axis(1, labels=ax.lab,at=ax.lab)
}

box(col="grey50")


}

dev.off()


