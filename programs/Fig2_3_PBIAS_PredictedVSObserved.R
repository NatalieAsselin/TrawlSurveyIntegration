##script to make predicted vs observed figure and pbias over time figure

## Set paths
main.fp <- getwd()

fig.path <- file.path(main.fp, "results/figures")
data.path <- file.path(main.fp, "data")
results.path <- file.path(main.fp, "results/tables")

path <- paste0(main.fp, "/programs/")
rdatapath <- paste0(main.fp, "/programs/RData/")

library(mgcv)

##labels for file names  
label <- "commercial biomass"
mod.label<- "m13"

load( file = paste0(rdatapath, "sGSL_", label, "_MultiSurveyComparisonResults_",mod.label,  ".rdata"))
load(file = paste0(rdatapath, "sGSL_", label, "_PBIAS_By_Year_",mod.label,  ".rdata"))

##width and height in inches
width <- 18.2 / 2.54 # 18.2 cm for two column figure 7.165 inches

height <- 21.6 /2.54  # 

grDevices::pdf(file = paste0(fig.path, "/sGSL_ObserverdVsPredicted", mod.label, ".pdf"),  height = height, width = width, pointsize=7)

m <- kronecker(matrix(c(1:8), ncol = 2, byrow=TRUE), matrix(1, ncol = 5, nrow = 5)) ##larger squares for the maps

m <- rbind(0, cbind(0, m,  0), 0)
layout(m)
par(mar = c(2.5, 1.5, 0, 0))

surveys <- c("SRVS", "NSS", "MIS", "ID")


for (i in 1:length(surveys)) {
  s <- surveys[i]
  
  d <- if (s=="NSS") {
    NS_nd
  } else if(s=="SRVS")  {
    RV_nd
  } else if(s=="MIS")  {
    MI_nd
  } else if(s=="ID")  {
    ALL_nd
  }
  
  
  tmp <- d[d$commercial.weight <=50,]
  
  print ("percent sets <=50kg")
  
  print(s)
  
  print (nrow(tmp)/nrow(d)) *100 
  
  res <- d
  
  xlim <- c(0, plyr::round_any(max(res$commercial.weight), f=ceiling, accuracy=10) )
  
  ylim <- c(0,plyr::round_any(max(c(res$predicted_ALL, res$predicted_survey)), f=ceiling, accuracy=10) ) 
  
  
  for(j in 1:2) {
    
    
    if(j==1){
      
      if (surveys[i] %in% c("SRVS", "NSS", "MIS")){ ##survey specific model
        
        plot( xlim,  ylim, type = "n", xlab = "", ylab = "", xaxs = "i", yaxs = "i", axes=FALSE)
        
        grid()
        
        if (surveys[i] == "SRVS") {
          mtext("Single-survey", side=3, line=1, cex = 2.5)
        }
        
        cols <- grDevices::adjustcolor(c("slateblue4"), alpha = 0.20)
        
        points(res$commercial.weight, res$predicted_survey, pch=19,  col = cols[1], bg= cols[1], cex=0.5, lwd=0)
        
        x <- gam(predicted_survey ~ s(commercial.weight, k=3), data = res, family = gaussian)
        
        dn <- expand.grid(commercial.weight = seq(round(min(res$commercial.weight),0), round(max(res$commercial.weight),0),  by = 1))
        
        res2  = predict(x, newdata = dn, se.fit=TRUE, type="response")
        
        results <- data.frame(mean = (res2$fit),
                              lower.ci = ((res2$fit) - 1.96 * (res2$se.fit)),
                              upper.ci = ((res2$fit) + 1.96 * (res2$se.fit)))
        results <- cbind(dn, results)
        
        res_survey <- results
        
        colnames(res_survey) <- c("observed weight (kg)", "mean_survey", "lower.ci.survey", "upper.ci.survey")
        
        lines(results$commercial.weight, results$mean, lwd = 2, lty=1, col= "orange3")
        
        lines(c(0, ylim[2]), c(0, ylim[2]), lty= 2, lwd=2, col="grey25")
        
        
      }
      
      if (surveys[i] %in% c("ID")){ ##survey specific model
        
        plot( xlim,  ylim, type = "n", xlab = "", ylab = "", xaxs = "i", yaxs = "i", axes=FALSE)
        
        
      }
      
      if (surveys[i] %in% c("SRVS", "NSS", "MIS")){ ##survey specific model
        
        axis(1, las=1, cex.axis=1.5)
      }
      axis(2, las=1, cex.axis=1.5)
      
      
      text(0.02*(xlim[2]-xlim[1])+xlim[1],
           0.90 * (ylim[2]-ylim[1]) + ylim [1],
           s, cex = 2.5, pos = 4)
      
      box(col="grey5")
      
    }
    
    if(j==2){
      
      plot( xlim,  ylim, type = "n", xlab = "", ylab = "", xaxs = "i", yaxs = "i", axes=FALSE)
      grid()
      if (surveys[i] == "SRVS") {
        mtext("ID", side=3, line=1, cex = 2.5)
      }
      
      points(res$commercial.weight, res$predicted_ALL, pch=19,  col = cols[1], bg=cols[1], cex=0.5, lwd=0)
      
      
      x <- gam(predicted_ALL ~ s(commercial.weight, k=3), data = res, family = gaussian) 
      
      dn <- expand.grid(commercial.weight = seq(round(min(res$commercial.weight),0), round(max(res$commercial.weight),0),  by = 1))
      res2  = predict(x, newdata = dn, se.fit=TRUE, type="response")
      
      results <- data.frame(mean = (res2$fit),
                            lower.ci = ((res2$fit) - 1.96 * (res2$se.fit)),
                            upper.ci = ((res2$fit) + 1.96 * (res2$se.fit)))
      results <- cbind(dn, results)
      
      res_ALL <- results
      colnames(res_ALL) <- c("observed weight (kg)", "mean_ALL", "lower.ci.ALL", "upper.ci.ALL")
      
      res_both <- merge(res_survey, res_ALL, all=TRUE)
      
      write.csv(res_both, file = paste0(results.path, "/sGSL_", label, "_observed vs predicted_", s, "_", mod.label, ".csv"), row.names = FALSE) ##save results to look at details
      
      lines(results$commercial.weight, results$mean, lwd = 2, lty=1, col= "orange3")
      
      lines(c(0, ylim[2]), c(0, ylim[2]), lty= 2, lwd=2, col="grey25")
      
      axis(1, las=1, cex.axis=1.5)
      box(col="grey5")
    }
    
  }
  
}


mtext("Predicted catch (kg)", side=2, line=-4.5, cex = 2.5, outer=TRUE)
mtext("Observed catch (kg)", side = 1, line=-3, cex = 2.5, outer=TRUE)

dev.off()

##Graph of PBIAS by year


width <- 18.2 / 2.54 # 18.2 cm for two column figure 7.165 inches

height <- 18.2 / 2.54  

grDevices::pdf(file = paste0(fig.path, "/sGSL_PBIAS_Year_IntegratedvsSingle_lineGraph_FourPanel_", mod.label,  ".pdf"),  height = height, width = width, pointsize=7)

m <- kronecker(matrix(c(1:4), ncol = 2, byrow=TRUE), matrix(1, ncol = 5, nrow = 5)) ##larger squares for the maps

m <- rbind(0, cbind(0, m,  0), 0)
layout(m)

par(mar = c(2, 1.5, 0, 0))


surveys <- c("SRVS", "NSS", "MIS", "ID")


xlim <- c(2000, 2026)
ylim <- c(-40, 40)


#cols <- c("blue", "red", "green", "black") 

cols <- c( "blue3", "black")

res.gls <- NULL

for (i in 1:length(surveys)) {
  plot( xlim,  ylim, type = "n", xlab = "", ylab = "", xaxs = "i", yaxs = "i", axes=FALSE)
  
  
  s <- surveys[i]
  
  text(0.02*(xlim[2]-xlim[1])+xlim[1],
       0.90 * (ylim[2]-ylim[1]) + ylim [1],
       s, cex = 3, pos = 4)
  
  x <- res_year[res_year$survey==s,]
  
  if (s %in% c("SRVS", "NSS", "MIS")){
    
    val <-x$pBIAS_survey
    yr <-  x$year-2001
    xx <- nlme::gls(val ~ yr,   correlation = corAR1())
    
    yy <- as.data.frame(summary(xx)$tTable)
    
    yyy <- yy
    
    int <- round(yy[1,1],0)
    slope <- round(yy[2,1], 3)
    pval <- round(yy[2,4],3)
    
    symb  <- ifelse (slope > 0, " + ", " - ")
    
    eq_label <- paste0("y = ", round(int, 2), symb,
                       abs(slope), "x")
    
    pval_label <-bquote(italic(p) == .(format(pval, digits=3, nsmall=3)))
    
    text(0.02*(xlim[2]-xlim[1])+xlim[1],
         0.15 * (ylim[2]-ylim[1]) + ylim [1],
         "Single-survey", cex = 2, pos = 4, col=cols[1])
    
    text(0.02*(xlim[2]-xlim[1])+xlim[1],
         0.1 * (ylim[2]-ylim[1]) + ylim [1],
         eq_label, cex = 1.5, pos = 4, col=cols[1])
    
    text(0.02*(xlim[2]-xlim[1])+xlim[1],
         0.05 * (ylim[2]-ylim[1]) + ylim [1],
         pval_label , cex = 1.5, pos = 4,col=cols[1])
    
    
    yy <-yy[2,]
    yy$survey <- s
    yy$model <- "single"
    
    res.gls <- rbind(res.gls, yy)
    
    points(x$year, x$pBIAS_survey, col=cols[1],  pch=23, cex=1.5)
    
    
    if(yy$`p-value` < 0.05) {
      
      lines(x=xlim, y= c(yyy[1,1], yyy[1,1] + yyy[2,1]*(xlim[2]-2001)) ,  col=cols[1], lty=1, lwd=2)
      
    }
    
  }
  
  
  
  lwd <-  2
  lty <- 1
  
  
  val <-x$pBIAS_All
  yr <-  x$year-2001
  xx <- nlme::gls(val ~ yr,   correlation = corAR1())
  
  yy <- as.data.frame(summary(xx)$tTable)
  
  yyy <- yy
  
  int <- round(yy[1,1],0)
  slope <- round(yy[2,1], 3)
  pval <- round(yy[2,4],3)
  
  symb  <- ifelse (slope > 0, " + ", " - ")
  
  eq_label <- paste0("y = ", round(int, 2), symb,
                     abs(slope), "x")
  
  pval_label <-bquote(italic(p) == .(format(pval, digits=3, nsmall=3)))
  
  text(0.40*(xlim[2]-xlim[1])+xlim[1],
       0.15 * (ylim[2]-ylim[1]) + ylim [1],
       "ID", cex = 2, pos = 4, col=cols[2])
  
  text(0.40*(xlim[2]-xlim[1])+xlim[1],
       0.1 * (ylim[2]-ylim[1]) + ylim [1],
       eq_label, cex = 1.5, pos = 4, col=cols[2])
  
  text(0.40*(xlim[2]-xlim[1])+xlim[1],
       0.05 * (ylim[2]-ylim[1]) + ylim [1],
       pval_label , cex = 1.5, pos = 4, col=cols[2])
  
  
  yy <-yy[2,]
  yy$survey <- s
  yy$model <- "ID"
  
  res.gls <- rbind(res.gls, yy)
  
  
  if(yy$`p-value` < 0.05) {
    
    lines(x=xlim, y= c(yyy[1,1], yyy[1,1] + yyy[2,1]*(xlim[2]-2001)) , col=cols[2], lty=1, lwd=2)
    
  }
  
  
  points(x$year, x$pBIAS_All, col=cols[2],  bg= cols[2], pch=21, cex=1)
  
  
  
  lines(c(xlim), c(0.0, 0.0), lty= 1, col="grey50")
  lines(c(xlim), c(-10, -10), lty= 2, col="grey50")
  lines(c(xlim), c(10, 10), lty= 2, col="grey50")
  
  box(col="grey5")
  
  
  if (i %in% c(1,3)){  
    axis(2, las=1, col="grey5", cex.axis=1.5)
  }
  if (i %in% c(3,4)){ 
    axis(1, col="grey5", cex.axis=1.5)
  }
}

mtext("PBIAS", 2, -4.5, cex = 2.5, outer=TRUE)

mtext("Year", 1, -3.5, cex = 2.5, outer=TRUE)


dev.off()




