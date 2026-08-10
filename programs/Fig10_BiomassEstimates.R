##Script to make figure of predictions


##Set paths
main.fp <- getwd()

fig.path <- file.path(main.fp, "results/figures")
data.path <- file.path(main.fp, "data")
results.path <- file.path(main.fp, "results/tables")

path <- paste0(main.fp, "/programs/")
rdatapath <- paste0(main.fp, "/programs/RData/")

mod.label <- "m13"


##load in individual index files


SRVS_index <- read.csv (paste0(results.path, "/sGSL_2001_2025_1Kres_RV_Index_", mod.label,"ALL.csv"))
vars <- colnames(SRVS_index)
colnames(SRVS_index) <- c(vars[1], "est_modSurvey", "lwr_modSurvey",  "upr_modSurvey" ,     "log_est_modSurvey",  "se_modSurvey",
                          vars[7:11])

NS_index <- read.csv (paste0(results.path, "/sGSL_2001_2025_1Kres_NS_Index_", mod.label,"ALL.csv"))
colnames(NS_index) <- c(vars[1], "est_modSurvey", "lwr_modSurvey",  "upr_modSurvey" ,     "log_est_modSurvey",  "se_modSurvey",
                        vars[7:11])
MI_index <- read.csv (paste0(results.path, "/sGSL_2001_2025_0.5Kres_MI_Index_", mod.label,"ALL.csv"))
colnames(MI_index) <- c(vars[1], "est_modSurvey", "lwr_modSurvey",  "upr_modSurvey" ,     "log_est_modSurvey",  "se_modSurvey",
                        vars[7:11])

All_index <- read.csv (paste0(results.path, "/sGSL_2001_2025_1Kres_ALL_Index_", mod.label,"ALL.csv"))

##4 panel vertical figure

width <- 8.84 / 2.54 # 8.84 cm for one column figure

height <- 20 /2.54  

grDevices::pdf(file = paste0(fig.path, "/sGSL_BiomassPredictions_",  mod.label,  ".pdf"),  height = height, width = width, pointsize=7)

m <- kronecker(matrix(c(1:4), ncol = 1, byrow=TRUE), matrix(1, ncol = 5, nrow = 5)) ##larger squares

m <- rbind(0, cbind(0, m, 0), 0)
layout(m)
par(mar = c(0.75, 0,0,0 ))
par(oma = c(0.0, 0.75, 0, 0))

xlim<-  c(2000, 2026)

surveys <- c("SRVS", "NSS", "MIS")

cols <- c("goldenrod3", "slateblue1" )

for (i in 1:length(surveys)) {
  
  surveys[i]
  
  if (surveys [i] == "SRVS"){
    d <- SRVS_index
    }else if (surveys [i] == "NSS") {
      d <- NS_index
    }else if (surveys [i] == "MIS") {
      d <- MI_index
    }
 
  ylim <- c(0, ((max(d$upr_modSurvey, d$upr_modALL) * 1.05)))
  
  plot( xlim,  ylim, type = "n", xlab = "", ylab = "", xaxs = "i", yaxs = "i", axes=FALSE)
  #grid()
  
  

  pol.col <- grDevices::adjustcolor(cols[1], alpha = 0.25)
  polygon(c(d$year,rev(d$year)),c(d$upr_modSurvey,rev(d$lwr_modSurvey)),col = pol.col, border = FALSE) 
  lines(d$year, d$est_modSurvey, lwd = 2, col=cols[1])
  
  
  pol.col <- grDevices::adjustcolor(cols[2], alpha = 0.25)
  polygon(c(d$year,rev(d$year)),c(d$upr_modALL ,rev(d$lwr_modALL )),col = pol.col, border = FALSE) 
  lines(d$year, d$est_modALL, lwd = 2, col=cols[2])

  text(0.05*(xlim[2]-xlim[1])+xlim[1],
       0.90 * (ylim[2]-ylim[1]) + ylim [1],
       surveys[i], cex = 2.5, pos = 4)
  
    axis(2, las=1, cex.axis=1.5)

  box(col = "grey50")
  
}
d <-  All_index

ylim <- c(0, ((max(d$upr_modALL) * 1.05)))

plot( xlim,  ylim, type = "n", xlab = "", ylab = "", xaxs = "i", yaxs = "i", axes=FALSE)

pol.col <- grDevices::adjustcolor(cols[2], alpha = 0.25)
polygon(c(d$year,rev(d$year)),c(d$upr_modALL ,rev(d$lwr_modALL )),col = pol.col, border = FALSE) 
lines(d$year, d$est_modALL, lwd = 2, col=cols[2])

## sGSL stock estimate
e <- read.csv (paste0(results.path, "/sGSL_", mod.label, "_CommercialBiomasssGSL_ 1000 tons_2001_2025_1Kres.csv"))

e <- e[e$area =="total",]


pol.col <- grDevices::adjustcolor("grey50", alpha = 0.50)
polygon(c(e$year,rev(e$year)),c(e$upr ,rev(e$lwr )),col = NULL, border = "grey50", lty=2, lwd=0.5) 
lines(e$year, e$est, lwd = 1, col="grey25", lty=2)



axis(2, las=1, cex.axis=1.5)
box(col = "grey50")
axis (1, cex.axis=1.5)

text(0.05*(xlim[2]-xlim[1])+xlim[1],
     0.90 * (ylim[2]-ylim[1]) + ylim [1],
     "ID", cex = 2.5, pos = 4)

mtext("Biomass (1000s t)", 2, -3, cex = 1.75, outer=TRUE)

mtext ("Year", 1, -1.2, cex = 1.75, outer=TRUE)

dev.off()

