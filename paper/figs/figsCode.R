#plots for the PMC paper
library(stats)

### Plot 1 ###
coef_truth <- read.table("/home/sjones/Documents/PMCs/paper/figs/ss48008/ss48008_90.77.txt")
spec_truth <- coef_truth$V2

coefs <- read.table("/home/sjones/Documents/PMCs/paper/figs/calc_coefs_132.8_2.2.txt")
spec <- exp(-(coefs$V2*(436.43*(1e-12))*(100000*(1e6))))

plot(rollmean(coef_truth$V1, 51), rollmean(coef_truth$V2, 51), ylim=c(0.96, 1.005), 
     xlim=c(2800,3550), type="l", xlab=parse(text="Wavenumber (cm^{-1})"), ylab="Transmission")
lines(1/(coefs$V1/10000), spec, col=2, lwd=2)

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

dataDens <- subset(data, data$Dens>101.3)
dataTemp <- subset(dataDens, dataDens$Temp>120.1 & dataDens$Temp<159.9)
dataEps <- subset(dataTemp, dataTemp$eps>1.11 & dataTemp$eps<4.94)
dataEps2 <- subset(dataEps, dataEps$epsErr>0)

#hist(data$Temp)
g1 <- ggplot(dataEps2, aes(Temp)) + geom_histogram(binwidth=3) + labs(x="Temperature (K)")
g2 <- ggplot(dataEps2, aes(eps)) + geom_histogram(binwidth=0.4) + labs(x="Axial Ratio")
g3 <- ggplot(dataEps2, aes(Dens)) + geom_histogram(binwidth=35) + labs(x="Density (/cm^3)")

grid.arrange(g1, g2, g3, nrow=3)