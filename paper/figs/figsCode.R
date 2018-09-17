#plots for the PMC paper
library(stats)
library(zoo)
par(mfrow=c(1,1))

### Plot 1 ###
coef_truth <- read.table("/home/sjones/Documents/PMCs/paper/figs/ss48008/ss48008_90.77.txt")
spec_truth <- coef_truth$V2

coefs <- read.table("/home/sjones/Documents/PMCs/paper/figs/calc_coefs_132.8_2.2.txt")
spec <- exp(-(coefs$V2*(436.43*(1e-12))*(100000*(1e6))))

plot(rollmean(coef_truth$V1, 51), rollmean(coef_truth$V2, 51), ylim=c(0.96, 1.005), 
     xlim=c(2800,3550), type="l", xlab=parse(text="Wavenumber (cm^{-1})"), ylab="Transmission")
lines(1/(coefs$V1/10000), spec, col=2, lwd=2)


### Overplot ss48008 measurements ###
coef1 <- read.table("/home/sjones/Documents/PMCs/paper/figs/ss48008/ss48008_90.77.txt")
coef2 <- read.table("/home/sjones/Documents/PMCs/paper/figs/ss48008/ss48008_85.47.txt")
coef3 <- read.table("/home/sjones/Documents/PMCs/paper/figs/ss48008/ss48008_80.13.txt")
coef4 <- read.table("/home/sjones/Documents/PMCs/paper/figs/ss48008/ss48008_74.76.txt")

png(filename="/home/sjones/Documents/PMCs/paper/figs/ss48008_all.png", width=600, height = 400)
plot(rollmean(coef1$V1, 51), rollmean(coef1$V2, 51), ylim=c(0.96, 1.005),
     xlim=c(2800,3550), type="l", xlab=parse(text="Wavenumber (cm^{-1})"), ylab="Transmission")
lines(rollmean(coef2$V1, 51), rollmean(coef2$V2, 51), col=2)
lines(rollmean(coef3$V1, 51), rollmean(coef3$V2, 51), col=3)
lines(rollmean(coef4$V1, 51), rollmean(coef4$V2, 51), col=4)
legend(2830, 0.98, c("90.77 km", "85.47 km", "80.13 km", "74.76 km"), col=c(1,2,3,4), lty=c(1,1,1), text.width = 100)
dev.off()

### Plot 2 ###

par(mfrow=c(3,1))
library(ggplot2)
require(gridExtra)

setwd("/home/sjones/Documents/PMCs")
data <- read.table("occOutputs.txt")
detections <- read.table("pmclist.txt")

names(data) <- c("Occultation", "Temp", "TempErr", "eps", "epsErr", "Dens", "DensErr")
names(detections) <- c("Occultation", "Altitude", "Ratio", "Latitude")

isEq <- which(detections$Occultation %in% data$Occultation)
data$Altitude <- detections$Altitude[isEq]
data$Latitude <- detections$Latitude[isEq]

dataDens <- subset(data, data$Dens>101.3)
dataTemp <- subset(dataDens, dataDens$Temp>120.1 & dataDens$Temp<159.9)
dataEps <- subset(dataTemp, dataTemp$eps>1.11 & dataTemp$eps<4.94)
dataEps2 <- subset(dataEps, dataEps$epsErr>0)
# Filter to higher than 79 km, like SOFIE, ice cannot form much lower
dataEps3 <- subset(dataEps2, Altitude > 79)

g1 <- ggplot(dataEps3, aes(Temp)) + geom_histogram(binwidth=3) + labs(x="Temperature (K)")
g2 <- ggplot(dataEps3, aes(eps)) + geom_histogram(binwidth=0.4) + labs(x="Axial Ratio")
g3 <- ggplot(dataEps3, aes(Dens)) + geom_histogram(binwidth=35) + labs(x="Density (/cm^3)")
g4 <- ggplot(dataEps3, aes(Altitude)) + geom_histogram(binwidth=5) + labs(x="Altitude (km)")

g <- arrangeGrob(g1, g2, g3, g4, nrow=4, ncol=1)
ggsave("/home/sjones/Documents/PMCs/paper/figs/allHistograms.png", g)

latAlt <- ggplot(dataEps2, aes(Latitude, Altitude)) + geom_point()