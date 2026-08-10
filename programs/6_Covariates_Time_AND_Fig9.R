##Script to make figure of covariates over time, and check for significant changes

library(Hmisc)

## Set paths
main.fp <- getwd()

fig.path <- file.path(main.fp, "results/figures")
data.path <- file.path(main.fp, "data")
results.path <- file.path(main.fp, "results/tables")

path <- paste0(main.fp, "/programs/")
rdatapath <- paste0(main.fp, "/programs/RData/")

mod.label <- "m13"


data <- read.csv (paste0(results.path, "/sGSL_2001_2025_1Kres_ALL_Predictions_", mod.label,"ALL.csv"))

data <- data[!is.na(data$pred_modALL),]

years <- c(min(data$year): max(data$year))


data2 <- data
data2$depth <- round(data2$depth, 0 )

by1 <- data2$year
by2 <- data2$depth

data.sum <- aggregate(data2$pred_modALL, by=list(by1,by2), sum)

colnames(data.sum) <- c("year", "depth", "pred_modALL")



##four panel vertical figure
width <- 8.84 / 2.54 # 8.84 cm for one column figure

height <- 15 /2.54  

grDevices::pdf(file = paste0(fig.path, "/sGSL_CovariatesOverTime_", mod.label, ".pdf"),  height = height, width = width, pointsize=7)

m <- kronecker(matrix(c(1:4), ncol = 1, byrow=TRUE), matrix(1, ncol = 5, nrow = 5)) ##larger squares

m <- rbind(0, cbind(0, m, 0), 0)
layout(m)
par(mar = c(1, 0,0,0 ))
par(oma = c(0.0, 0.75, 0, 0))

xlim<-  c(2000, 2026)

cols <- c("royalblue3", "grey25", "darkorange3")

pch <- c(21,22,23)

#Depth quantiles
res <- NULL

probs <- c(0.05, 0.5, 0.95)

for (i in min(years): max(years)) {
  
  print(i)

  data.yr <- data[data$year==i,]

  a <- NULL

  a$year <- i

 
b <- as.data.frame(t(wtd.quantile( data.yr$depth,  data.yr$pred_modALL, prob =probs)))

colnames(b)<- probs

 a <- cbind (a,b)
 
 res <- rbind(res, a)

}

res2 <- res ##make a copy to not mess up depth symbol everywhere

##make all depths negative for figure

res$`0.05` <- - (res$`0.05`)
res$`0.5` <- - (res$`0.5`)
res$`0.95`<- - (res$`0.95`)

##rename categories for gls 
res$shallow <- res$`0.05`
res$mid <- res$`0.5`
res$deep <-  res$`0.95`

ylim <- c(-53,0)

plot( xlim,  ylim, type = "n", xlab = "", ylab = "", xaxs = "i", yaxs = "i", axes=FALSE)

grid()


x <- nlme::gls(deep ~ year, data = res, correlation = corAR1()) 
dn <- expand.grid(year = seq(min(res$year), max(res$year),  by = 1))
top.line  = predict(x, newdata = dn, type="response")
xx <- as.data.frame(summary(x)$tTable)
x.pval <- xx[2,4]


y <- nlme::gls(mid ~ year, data = res, correlation = corAR1()) 
half.line  = predict(y, newdata = dn, type="response")
yy <- as.data.frame(summary(y)$tTable)
y.pval <- yy[2,4]


z <- nlme::gls(shallow ~ year, data = res, correlation = corAR1()) 
zz <- as.data.frame(summary(z)$tTable)
z.pval <- zz[2,4]


bottom.line  = predict(z, newdata = dn, type="response")
pol.col <- grDevices::adjustcolor("grey25", alpha = 0.25)
polygon(c(res$year,rev(res$year)),c(top.line ,rev(bottom.line)),col = pol.col, border = FALSE) ##showing range


points(res$year, res$`0.05`, pch=pch[1], col = cols[1], bg = cols [1],cex=1) 
if(z.pval <=0.05){
abline(z, lwd = 2, col = cols[1] , lty=2)
}
points(res$year, res$`0.5`, , pch=pch[2], col = cols[2] , bg = cols [2],cex=1) 
if(y.pval <=0.05){
abline( y, lwd = 2, col = cols[2] )
}
points(res$year, res$`0.95`, pch=pch[3], col = cols[3] , bg = cols [3],cex=1) 

if(x.pval <=0.05){
abline(x, lwd = 2, col = cols[3] )
}

axis(2, labels=c(0, 10, 20 , 30, 40, 50 ),at=c(0, -10, -20, -30, -40, -50), las=2)


mtext("Depth (m)", 2, 3, cex = 1.25)

box(col = "grey50")

leg_title <- "Quantile" 

legend("bottomright", legend = probs, pch = pch,  pt.bg = cols, col = cols, ncol=1, bty="0", bg="white", box.col="white")

# text(0.05*(xlim[2]-xlim[1])+xlim[1],
#      0.05 * (ylim[2]-ylim[1]) + ylim [1],
#      "A", cex = 2.5, pos = 4)

mtext("A", 4, 1, cex = 1.25, las=2)

box(col = "grey50")


##temp quantiles

res.temp <- NULL

probs <- c(0.05, 0.5, 0.95)

for (i in min(years): max(years)) {
  
  print(i)
  
  data.yr <- data[data$year==i,]
  
  a <- NULL
  
  a$year <- i
  
  
  b <- as.data.frame(t(wtd.quantile( data.yr$bottom.temp,  data.yr$pred_modALL, prob =probs)))
  
  colnames(b)<- probs
  
  a <- cbind (a,b)
  
  res.temp <- rbind(res.temp, a)
  
}

ylim <- c(0,20)
plot( xlim,  ylim, type = "n", xlab = "", ylab = "", xaxs = "i", yaxs = "i", add=TRUE, xaxt="n", axes=FALSE)

grid()

res.temp$warm <- res.temp$`0.95`
res.temp$mid <- res.temp$`0.5`
res.temp$cold <- res.temp$`0.05`

x <- nlme::gls(cold ~ year, data = res.temp, correlation = corAR1()) 

dn <- expand.grid(year = seq(min(res$year), max(res$year),  by = 1))
top.line  = predict(x, newdata = dn, type="response")
xx <- as.data.frame(summary(x)$tTable)
x.pval <- xx[2,4]


y <-nlme::gls(mid ~ year, data = res.temp, correlation = corAR1()) 
mid.line  = predict(y, newdata = dn, type="response")
yy <- as.data.frame(summary(y)$tTable)
y.pval <- yy[2,4]


z <-nlme::gls(warm ~ year, data = res.temp, correlation = corAR1()) 
zz <- as.data.frame(summary(z)$tTable)
z.pval <- zz[2,4]

dn <- expand.grid(year = seq(min(res$year), max(res$year),  by = 1))
bottom.line  = predict(z, newdata = dn, type="response")
pol.col <- grDevices::adjustcolor("grey25", alpha = 0.25)
polygon(c(res.temp$year,rev(res.temp$year)),c(top.line ,rev(bottom.line)),col = pol.col, border = FALSE) ##showing range

points(res.temp$year, res.temp$`0.95`, pch=pch[1], col = cols[1] , bg = cols [1],cex=1) 
if(z.pval <=0.05){
  abline(z, lwd = 2, col = cols[1] , lty=2)
}

points(res.temp$year, res.temp$`0.5`, pch=pch[2], col = cols[2] , bg = cols [2],cex=1)  

if(y.pval <=0.05){
  abline( y, lwd = 2, col = cols[2] )
}

points(res.temp$year, res.temp$`0.05`, pch=pch[3], col = cols[3] , bg = cols [3],cex=1) 

if(x.pval <=0.05){
  abline(x, lwd = 2, col = cols[3] )
}

ylab  = expression("Temperature ("*~degree*C*")") 

mtext(ylab, 2, 3, cex = 1.25)

axis(2, las=2)

leg_title <- "Quantile" 

legend("bottomright", legend = probs, pch = pch,  pt.bg = cols, col = cols, ncol=1, bty="o", bg="white", box.col="white")
# 
# text(0.05*(xlim[2]-xlim[1])+xlim[1],
#      0.05 * (ylim[2]-ylim[1]) + ylim [1],
#      "B", cex = 2.5, pos = 4)

mtext("B", 4, 1, cex = 1.25, las=2)

box(col = "grey50")

##proportion by substrate

by1 <- data$year
by2 <- data$substrate

a <- aggregate (data$pred_modALL, by = list (by1,by2), sum)

colnames (a) <- c("year", "substrate", "kg")

b <- aggregate (data$pred_modALL, by = list (by1), sum)

colnames (b) <- c("year",  "total_kg")

a <- merge(a,b)

a$"ratio" <- a$kg/a$total_kg

ylim <- c(0,1)

plot(xlim, ylim, type = "n", xlab = "", ylab = "", xaxs = "i", yaxs = "i", axes=FALSE)

grid()

by1 <- a$year
check <- aggregate(a$ratio, by=list(by1), sum)


p <- a[a$substrate=="pelite",]
s <- a[a$substrate=="sand",]
g <- a[a$substrate=="gravel",]

p$topline <- 1
p$bottomline <- s$ratio + g$ratio

s$topline <- p$bottomline
s$bottomline <-g$ratio

g$topline <- g$ratio
g$bottomline <- 0

substrates <- c("p", "s", "g")

cols <- c("#8b4513","#ddd6c2","#69747c" )


for (i in 1:length(substrates)){
  
  x <-  if (substrates[i] == "p"){
    p
  }else if (substrates[i] == "s"){
    s
  }else if (substrates[i] == "g"){
    g
  }
  
  
  pol.col <- grDevices::adjustcolor(cols[i], alpha = 0.5)
  polygon(c(x$year,rev(x$year)),c(x$topline,rev(x$bottomline)),col = pol.col, border = FALSE) ##showing range
  
}

axis(2, las=2)

mtext("Proportion", 2, 3, cex = 1.25)

box(col = "grey50")

leg.labels <- c("Mud", "Sand", "Gravel")

leg.col <- grDevices::adjustcolor(cols, alpha = 0.5)

legend("bottomright", legend = leg.labels , pch = 22,  pt.bg = leg.col, col = leg.col, ncol=1,  bty="o", bg="white", box.col="white")

# text(0.05*(xlim[2]-xlim[1])+xlim[1],
#      0.05 * (ylim[2]-ylim[1]) + ylim [1],
#      "C", cex = 2.5, pos = 4)

mtext("C", 4, 1, cex = 1.25, las=2)

box(col = "grey50")




cols <- c("royalblue3", "grey25", "darkorange3")

a <- data

dens <- c(100,500,1000)

c <- as.data.frame(years)

for (i in 1:length (dens) ) {

a$lobster.hab <- ifelse(a$pred_modALL >= dens[i], 1,0 ) ##checking what 500kg / km2 looks like

by1 <- a$year

b <- aggregate(a$lobster.hab, by=list(by1), sum)

colnames(b) <- c("year", paste0("km2 with ", dens[i], "_kg"))

c <- cbind(c,b [2] )

}


ylim <- c(0, max(c [, 2]) * 1.2)

plot(xlim, ylim, type = "n", xlab = "", ylab = "", xaxs = "i", yaxs = "i", axes=FALSE)

grid()

##lines
for (i in 1:length (dens) ) {

yval <- c [, i+1]  
yr <- c[,1]
xx <- nlme::gls(yval ~ yr,   correlation = corAR1()) 
area  = predict(xx, newdata = dn, type="response")
abline(xx, lwd = 1, col = cols[i])
#points(yr,yval, pch=pch[i], col = cols [i] , bg = cols [i], cex=1) 
}


for (i in 1:length (dens) ) {
  
  yval <- c [, i+1]  
  points(c$year,yval, pch=pch[i], col = cols [i] , bg = cols [i], cex=1) 
}

ylab <- bquote("Area (" * 1000 * "s km"^2 * ")")

ticks <- axTicks(2)/1000

axis(2,
     at = axTicks(2),
     labels = ticks,
     las = 2)

axis(1)

mtext(ylab, 2, 3, cex = 1.25)
mtext("Year", 1, 3, cex = 1.25)

leg_title<- expression(Kg/km^{2}) 

leg.dens <- format(dens, big.mark = ",", scientific = FALSE)

legend("bottomright", legend = leg.dens, pch = pch,  pt.bg = cols, col = cols, ncol=1, title = leg_title, bty="o", bg="white", box.col="white")

# text(0.05*(xlim[2]-xlim[1])+xlim[1],
#      0.05 * (ylim[2]-ylim[1]) + ylim [1],
#      "D", cex = 2.5, pos = 4)

mtext("D", 4, 1, cex = 1.25, las=2)

box(col = "grey50")

dev.off()







##check for significant trends for substrates
by1 <- data$year
by2 <- data$substrate

a <- aggregate (data$pred_modALL, by = list (by1,by2), sum)

colnames (a) <- c("year", "substrate", "kg")

b <- aggregate (data$pred_modALL, by = list (by1), sum)

colnames (b) <- c("year",  "total_kg")

a <- merge(a,b)

a$"ratio" <- a$kg/a$total_kg



plot(a$year [a$substrate =="gravel"],a$ratio [a$substrate =="gravel"])
#xx <- lm(a$ratio [a$substrate =="gravel"] ~ a$year [a$substrate =="gravel"])
val <-a$ratio [a$substrate =="gravel"]
yr <-  a$year [a$substrate =="gravel"]
xx <- nlme::gls(val ~ yr,   correlation = corAR1()) ##regression with ar1 temperal correlation ##not significant

gravel = predict(xx, newdata = dn, type="response")
abline(xx, lwd = 2, col = "grey40")##not significant


plot(a$year [a$substrate =="pelite"],a$ratio [a$substrate =="pelite"])

val <-a$ratio [a$substrate =="pelite"]
yr <-  a$year [a$substrate =="pelite"]
xx <- nlme::gls(val ~ yr,   correlation = corAR1()) ##regression with ar1 temperal correlation ##not significant
pelite = predict(xx, newdata = dn, type="response")
abline(xx, lwd = 2, col = "grey40")##not significant

plot(a$year [a$substrate =="sand"],a$ratio [a$substrate =="sand"])
val <-a$ratio [a$substrate =="sand"]
yr <-  a$year [a$substrate =="sand"]
xx <- nlme::gls(val ~ yr,   correlation = corAR1()) ##regression with ar1 temperal correlation ##not significant
sand = predict(xx, newdata = dn, type="response")
abline(xx, lwd = 2, col = "grey40")##not significant
