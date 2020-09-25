# GLM for relating water quality to distance to river mouths.
#
# CJ Brown
# 16 Feb 2017

rm(list = ls())

library(raster)
library(RColorBrewer)
library(splines)
library(mgcv)

 setwd('/Users/s2973410/Databases/Fiji reefs analysis/wq_models/large_region_data_v3')

# ***************** #
# DATA PREPARATION
# ***************** #

#
# Load data
#

#WQ data
dat <- read.csv('WQ_dat_sampled_grouped.csv', header = T)
#river data
bfull <- read.csv('RiversPoints_with_groups.csv', header = T)
#River summarised by groups
b1 <- read.csv('River_groups_sampled.csv', header = T)

icellhr <- read.csv('sampled_cell_ids.csv')

nriv <- nrow(b1)

refriv <- which.max(b1$mnyld_tons)

#Scaling for distances - distances as loaded are in metres
xscale <- 1/1000 #turn into km

#Make sediment loadings numeric (JAGS compile fails if it's a matrix or dataframe)
byld <- as.numeric(b1$mnyld_tons)

#Rescale sediment contributions for using in Bayesian model. Easier to guess initial values then.
bmean <- byld/byld[refriv]

#
# Select data inputs for distance and turbidity
#

ywq <- as.numeric(dat$wq)
icol <- grep(pattern = 'groupID', names(dat))
xdat <- as.numeric(unlist(dat[,icol]))
xdat <- matrix(xdat, ncol = nriv, byrow = F)*xscale

N <- nrow(dat)
print(N)

#
# Rescale y to be in (0, 1) and add a small constant
#

ymin <- min(ywq)
cadd <- 0.001 #constant so there are no zeros
dat$yscale <- (ywq-min(ywq)+cadd)/max(ywq-min(ywq)+cadd)

dat$mdist <- apply(xdat, 1, min)

dat$lon <- (dat$x-min(dat$x))/max(dat$x-min(dat$x))
dat$lat <- (dat$y-min(dat$y))/max(dat$y-min(dat$y))

#
# Binomial GAM with log link
#

m1 <- gam(yscale ~ s(mdist) + s(lat) + s(lon), family = binomial(link = log), data = dat)
summary(m1)

par(mfrow = c(1,3))
plot(m1)
plot(dat$yscale, predict(m1, type = 'response'))
abline(0,1)
cor(dat$yscale, predict(m1, type = 'response'))

newdat <- data.frame(mdist = seq(min(dat$mdist), max(dat$mdist), length.out = 100),
    lat = 0.5, lon = 0.5)
newdat$p <- predict(m1, newdata = newdat, type = 'response')

plot(dat$mdist, dat$yscale)
lines(newdat$mdist, newdat$p, col = 'red')


#
# Lognormal GAM
#
dat$lny <- log(dat$yscale)

m2a <- gam(yscale ~ s(mdist, k = 20) + s(lat) + s(lon), family = gaussian(link = log), data = dat)
m2 <- gam(yscale ~ s(mdist, k = 8) + s(lat, lon, k = 35), family = gaussian(link = log), data = dat)

summary(m2)
gam.check(m2)

#EMFig2
dev.new(width = 8, height = 4)
par(mfrow = c(1,2))
plot(m2, xlab = "Distance to nearest river (km)", ylab ="Effect on log(turbidity)", select = 1, main = "A")
plot(m2, xlab = "x-coordinate", ylab ="y-coordinate", select = 2, scheme = 2, main = "B", contour.col = "black", n2 = 50, too.far = 0.03, labcex = 1.1, vfont = c("sans serif", "bold"))


plot(dat$yscale, predict(m2, type = 'response'))
abline(0,1)
cor(dat$yscale, predict(m2, type = 'response'))

newdat <- data.frame(mdist = seq(min(dat$mdist), max(dat$mdist), length.out = 100),
    lat = 0, lon = 0)
newdat$p <- predict(m2, newdata = newdat, type = 'response')

plot(dat$mdist, dat$yscale)
lines(newdat$mdist, newdat$p, col = 'red')

# ---------------
# Map turbidity across region
# ---------------
rwq <- raster('Turbidity_not_sampled.grd')
wq2 <- (rwq[]-min(ywq)+cadd)/max(ywq-min(ywq)+cadd)
rwq2 <- raster(rwq)
rwq2[] <- wq2
plot(rwq2)

rdists_agg <- stack('river_distances/agg_river_dist_highres.grd')
rmdist <- stackApply(rdists_agg, rep(1, nlayers(rdists_agg)), min)
rmdist2 <- rmdist/1000
plot(rmdist2)

xycell <- xyFromCell(rmdist2, 1:ncell(rmdist2))
xcell <- (xycell[,1]-min(dat$x))/max(dat$x-min(dat$x))
ycell <- (xycell[,2]-min(dat$y))/max(dat$y-min(dat$y))
#rx <- raster(rmdist2); rx[] <- xcell; plot(rx)

datreg <- data.frame(mdist = rmdist2[], lon = xcell, lat = ycell)

datreg$pred <- as.numeric(predict(m2, newdata = datreg, type = 'response'))
rpred <- raster(rmdist2)
rpred[] <- datreg$pred

rresid <- rpred - rwq2
par(mfrow = c(1,2)); plot(rpred); plot(rwq2)
plot(rpred - rwq2)

#
# Save model
#

llminmax <- list(lonmin = min(dat$x), lonmax = max(dat$x),
    latmin = min(dat$y), latmax = max(dat$y))

mlist <- list(m1 = m1, m2 = m2, llminmax = llminmax, rpred = rpred, rresid = rresid)

save(mlist, file = "~/Code/Fiji reefs analysis/Ridge2ReefFish/data-raw/gam_turbidity_modelv2.RData")
