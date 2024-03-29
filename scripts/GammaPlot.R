library(ggplot2)
library(reshape2)

rm(list=ls())

Gamma <- read.csv("Bin_2F_Intensity.csv",header=TRUE,row.names=1)
Gamma <- t(Gamma)
GammaP <- Gamma/rowSums(Gamma)

GammaMelt <- melt(GammaP)

colnames(GammaMelt) <- c('Sample','Strain','Freq')

GammaMelt$Strain <- as.factor(GammaMelt$Strain)

GammaMelt$Sample <- gsub('X','S', GammaMelt$Sample)

p <- ggplot(data=GammaMelt,aes(x=Sample,y=Freq,colour=Strain,group=Strain)) + geom_point()

p <- p + theme_bw() + ylab("Relative frequency")

p <- p + theme(axis.title = element_text(size=12, face="bold")) 

pdf("TimeSeries.pdf")
plot(p + geom_smooth(se=F))
dev.off()
