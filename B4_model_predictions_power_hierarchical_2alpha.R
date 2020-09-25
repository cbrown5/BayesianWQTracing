# Summary stats of model predictions from run of Bayes water quality model
#
# CJ Brown 17 May 2016
# This version uses a power function and source contributions are sampled from a heirarchical model.


rm(list = ls())

library(splines)
library(RColorBrewer)
library(ggplot2)
library(tidyr)
library(corrplot)
library(rjags)

setwd('/Users/s2973410/Databases/Fiji reefs analysis/wq_models/large_region_data_v3')
source('~/Code/Useful R Code/figlabel.R')

# runname <- 'large_region_hp_21Mar16_unscaled'
# runname <- 'large_region_2Mar16_exp_v3'
 # runname <- 'large_region_hp_1Apr16_unscaled'
  # runname <- 'large_region_hp_3Apr16_unscaled'
 # runname <- 'large_region_hp_28Apr16'
# runname <- 'large_region_hp_17May16_2alpha'
runname <- 'large_region_hp_24May16_2alpha'



# dathighres <- read.csv('WQ_dat_highres_grouped.csv', header = T)
# save(list = 'dathighres', file = 'WQ_dat_highres_grouped.RData')
load('WQ_dat_highres_grouped.RData') #use load to save some time

icellhr <- read.csv('sampled_cell_ids.csv')$x

#Validation data
ikeep <- which(!dathighres$cells %in% icellhr)
datv <- dathighres[ikeep,]

#Check validation data
# datv2 <-dathighres[which(dathighres$cells %in% icellhr),]
# library(sp)
# datsp <- datv[,c(30:32)]
# datsp2 <- datv2[,c(30:32)]
# coordinates(datsp) <- ~x+y
# coordinates(datsp2) <- ~x+y
# plot(datsp)
# plot(datsp2, add = T, col ='red')

rm(dathighres)

#
#  Some key catchments
#
icatch <- c(8, 15, 7, 12, 23, 10, 26, 27)
# catchnams <- c('Savusavu Bay catchments','Dreketi river',
# 'Bua Bay catchments', 'Dawara river' ,'Labasa', 'Wainunu Bay catchments')
 catchnams <- c('','Dreketi river',
 'Bua Bay catchments', '' ,'Labasa', '', 'Nasauu river','Rukuruku Bay')


#
# Load last run
#
load(paste(runname,'RData', sep ='.'))

#
# Data summaries
#

smc <- summary(mcout)
smc

#Which rows of summary are which params?
dimnam <- attr(smc$stat, 'dimnames')[[1]]

ibeta <- grep('lnbeta', dimnam)
ialpha <- grep('alpha', dimnam)
ilnb <- grep('lnb', dimnam)[1]

a_est <- -smc$statistics[ialpha,1]
b_est <- exp(smc$statistics[ilnb,1])
betaest <- exp(smc$statistics[ibeta,1])
betaCI <- exp(smc$quantiles[ibeta,c(1,5)])
sdest <- sqrt(1/smc$statistics[grep('tau', dimnam),1][1])

#
# Plot traces for alpha and ascale
#
dev.new()
par(mfrow =c(2,3))
plot(as.numeric(mcout[[1]][,ialpha[1]]),type = 'l', xlab = '', ylab ='', main = 'alpha')
plot(as.numeric(mcout[[1]][,ialpha[2]]),type = 'l', xlab = '', ylab ='', main = 'alpha')

plot(as.numeric(mcout[[1]][,ilnb]),type = 'l', xlab = '', ylab ='', main = 'ln(w)')
plot(as.numeric(mcout[[1]][,ialpha[1]]), as.numeric(mcout[[1]][, ilnb]),type = 'p', xlab = '', ylab ='', main = 'correlation')
plot(as.numeric(mcout[[1]][,ialpha[2]]), as.numeric(mcout[[1]][, ilnb]),type = 'p', xlab = '', ylab ='', main = 'correlation')

#
# Plot betas vs yields
#
pchcols <- grey(0.2, 0.4)
figlab <- LETTERS[1:4]
xoff <- 0.05
yoff <- 0.08
letcex <- 0.6
letcex1 <- 1.5
#Offsets for labels
xofflab <- rep(-0.03, length(catchnams))
xofflab[8] <- -0.1
yofflab <- c(-0.04, 0, 0.1, 0.02, 0, 0, 0, -0.15)

dev.new(width = 9, height = 4)
par(mfrow = c(1,2), mar = c(5,5,2,2))
plot(bmean, betaest/b_est, col = pchcols, ylim = c(0,4), xlim = c(0,1.5), pch = 16, xlab = 'Relative influence on turbidity', ylab = 'Relative influence on turbidity (posterior)')
arrows(bmean, betaCI[,1]/b_est, bmean, betaCI[,2]/b_est, len =0, col = pchcols)
abline(0, 1, col = 'grey80')
text(bmean[icatch]+xofflab, (betaest/b_est)[icatch]+yofflab, catchnams, cex = letcex, pos = 4)
# text(bmean, betaest/b_est, 1:28)
figlabel(label = figlab[1], font =2, xoffset=xoff, yoffset = yoff, cex =letcex1)


plot(log(bmean), log(betaest/b_est), ylim = c(-14, 2), xlim = c(-10, 0), pch = 16, col = pchcols, xlab = 'Relative influence on turbidity (log)', ylab = 'Relative influence on turbidity (posterior, log)')

arrows(log(bmean), log(betaCI[,1]/b_est), log(bmean), log(betaCI[,2]/b_est), len =0, col = pchcols)
abline(0,1, col = 'grey80')
figlabel(label = figlab[2], font =2, xoffset=xoff, yoffset = yoff, cex =letcex1)


#Plot for presentation

dev.new()
par(mar =c(5,6,4,2))
plot(log10(bmean), log10(betaest/b_est), ylim = c(-6, 2), xlim = c(-3, 0),
    pch = 16, col = pchcols, cex = 1.5, xlab = "", ylab = "", cex.axis = 2, las = 1, xaxt = 'n', yaxt = 'n')
#arrows(log10(bmean), log10(betaCI[,1]/b_est), log10(bmean), log10(betaCI[,2]/b_est), len =0, col = pchcols,
    #lwd = 2)
abline(0,1, col = 'blue', lwd = 2)
axis(1, cex = 1.5, at = -3:0, labels = 10^(-3:0), cex.axis = 2)
axis(2, cex = 1.5, at = seq(-6, 2, by = 2), labels = 10^seq(-6, 2, by = 2), cex.axis = 2, las = 1)



dev.new()
par(mar =c(5,6,4,2))
plot(b1$mnyld_tons, as.numeric(betaest/b_est), ylim = c(0.000001, 150),
    xlim = c(100, 65000), pch = 16, col = pchcols,
    xlab = '', ylab = '', log = 'xy', las = 1, cex.axis = 2, cex = 1.5)
arrows(b1$mnyld_tons, betaCI[,1]/b_est, b1$mnyld_tons, betaCI[,2]/b_est, len =0, col = pchcols, lwd = 2)

# text(bmean, betaest, 1:28)

#
# Plot correlation between parameters
#

cor(mcout[[1]])
corrplot(cor(mcout[[1]]), diag = F, type = 'lower')

mcall <- rbind(mcout[[1]], mcout[[2]], mcout[[3]])
corrplot(cor(mcall), diag = F, type = 'lower')

#
# Convergence
#
gelman.diag(mcout, autoburnin = F)

#
# Predictions
#

bmat <- matrix(rep(betaest, length(y)), ncol = nriv, byrow = T)
amat <- matrix(rep(a_est[icoast],length(y)), ncol = nriv, byrow = T)

ymat <- bmat *(xdat^ amat)

ypred <- rowSums(ymat)*exp(sdest/2)
yresid <- yscale - ypred


cor(ypred, yscale)
cor(log(ypred), lny)

plot(log(ypred), lny, xlab = 'predictions', ylab = 'observed'); abline(0,1)

plot(ypred, yscale, xlab = 'predictions', ylab = 'observed', xlim = c(0,1)); abline(0,1)

plot(ypred, yresid, xlab = 'predictions', ylab = 'residuals'); abline(0,0)


#
# Predictive deviance
#
lnypred <- log(ypred)
Gm <- sum((lnypred - lny)^2)
Pm <- mean(lnypred ^2) - (mean(lnypred)^2)
var(lny - lnypred)
var(lnypred)

Gm + Pm

#
# Calculate DIC
#

predictmodel <- function(betavalf, yf, nrivf, xmapf, a_estf, icoastf){
	# Create predictions df
	bmapf <- matrix(rep(betavalf, length(yf)), ncol = nrivf, byrow = T)
	amapf <- matrix(rep(a_estf[icoastf], length(yf)), ncol = nrivf, byrow = T)
	log(rowSums(bmapf * (xmapf ^ amapf)))
}

predictmodel_mean <- function(betavalf, yf, nrivf, xmapf, a_estf, icoastf, sdest){
	# Create predictions df
	bmapf <- matrix(rep(betavalf, length(yf)), ncol = nrivf, byrow = T)
	amapf <- matrix(rep(a_estf[icoastf], length(yf)), ncol = nrivf, byrow = T)
	log(rowSums(bmapf * (xmapf ^ amapf))*exp(sdest/2))
}


nchain <- length(mcout)
nstep <- nrow(mcout[[1]])
dimnam <- attr(mcout[[1]], 'dimnames')[[2]]
ibeta <- grep('lnbeta',dimnam)
ialpha <- grep('alpha',dimnam)
itaubeta <- grep('taubeta',dimnam)
itau <- grep('tau',dimnam)
itau <- itau[itau!=itaubeta]


# Posterior deviance for each step
loglik <- matrix(NA, nrow = nstep, ncol = nchain)
for (ichain in 1:nchain){
	xchain <- as.matrix(mcout[[ichain]])
	for (i in 1:nstep){
		lnypred <- predictmodel(exp(xchain[i,ibeta]), y, nriv, xdat, -xchain[i, ialpha], icoast)
		loglik[i, ichain] <-  sum(-2*dnorm(lny, mean = lnypred, sd = sqrt(1/xchain[i,itau]), log = T))
	}
}

# Posterior averaged deviance
lnypred <- predictmodel_mean(betaest, y, nriv, xdat, a_est, icoast, sdest)
loglik_bar <- sum(2*dnorm(lny, mean = lnypred, sd = sdest, log = T))

# Effective number of params
pD <- mean(loglik) + loglik_bar

# DIC
mean(loglik) + pD


#
# Predict to cells not used in fitting
#
yobshr <- as.numeric(datv$wq)
preddat <- as.numeric(unlist(datv[,icol]))
preddat <- matrix(preddat, ncol = nriv, byrow = F)*xscale

betaest2 <- betaest
# betaest2[1] <- 1

bmathr <- matrix(rep(betaest2, length(yobshr)), ncol = nriv, byrow = T)
amathr <- matrix(rep(a_est[icoast], length(yobshr)), ncol = nriv, byrow = T)
ymathr <- bmathr *(preddat^ amathr)

ypredhr <- rowSums(ymathr)*exp(sdest/2)

yntuhr <- ypredhr + (ymin-cadd)

yresidhr <- yobshr - yntuhr

plot(yntuhr, yresidhr); abline(0,0)
plot(log(yntuhr), log(yobshr)); abline(0,1)
cor(yobshr, yntuhr)
cor(log(yntuhr), log(yobshr))

# Save df of high res data
dfsav <- data.frame(predntu = yntuhr, obsntu = yobshr)
write.csv(dfsav, paste0(runname, 'datpred_highres.csv'))


#Root mean-squared error - error in NTUs
#This equals sd if there is no bias
sqrt(sum((yntuhr - yobshr)^2)/length(yntuhr))

rm(preddat)
rm(bmathr)
rm(ymathr)


#
# Create plot of predictions versus distances
#
nx <- 100
xdist <- seq(1, 10, length.out = nx)

y1 <- betaest[10] * (xdist ^ a_est[icoast[10]])
y2 <- betaest[19] * (xdist ^ a_est[icoast[10]])

par(las = 1, mar = c(6,6,2,2))
plot(xdist, y1, ylim = c(0, 4), type = 'l', col = 'tomato3', lwd = 5, yaxt = 'n', xlab = 'Distance from river mouth (km)', cex.axis =2.53, tck = 0.02, cex.lab = 2.8, ylab = 'Contribution to turbidity', bty = 'l')
lines(xdist, y2, lwd = 5, col = 'steelblue')

lines(xdist, y1+y2, lwd = 5, col = 'black')

#
# Save summary data
#

save(list = ls(), file = paste(runname, '_stats','.RData', sep =''))

#
#Plot posterior for beta
#
library(ggplot2)
library(tidyr)

nrows <- length((mcout[[1]][,1]))
nchains <- length(mcout)

ybeta <- data.frame(matrix(NA, nrow = nrows * nchains, ncol = nriv-1))
for (i in 1:nriv){
	for (j in 1:nchains){
		ybeta[(j*nrows - nrows+1) : (j*nrows),i] <-
		 exp(as.numeric(mcout[[j]][,ibeta[i]]))
	}
}


colnames(ybeta) <- paste('Source', 1:nriv, sep ='_')
ybetas <- gather(ybeta, betavar, Beta)

blines <- data.frame(betavar = unique(ybetas$betavar), bpost = betaest[1:nriv], bprior = bmean[1:nriv]*b_est)

ybetas$sourcenum <- as.numeric(unlist(lapply(strsplit(ybetas$betavar, '_'), function(x) x[2])))

ybetas <- ybetas[order(ybetas$sourcenum),]

ggplot(ybetas, aes(x=Beta)) +
geom_density(fill = 'grey', alpha = 0.2, color = 'grey') +
ylab("Density") +
facet_wrap( ~ betavar,nrow = 5, scales = 'free_y') +
theme_bw() +
geom_vline(data = blines, aes(xintercept = bpost), color = 'grey20') +
geom_vline(data = blines, aes(xintercept = bprior), color = 'grey20', lty = 2) +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
