# Bayesian model for relating water quality to distance to river mouths. 
#
# CJ Brown with Ben Stewart-Koster
# 17 May 2016
# 
# This version uses a power function and source contributions are sampled from a heirarchical model. 
#
# Allows for two alpha values
#
# Model specification: 
#	
# y[i] ~ sum_j(z[i,]) * exp(dnorm(0,tau))
# z[i,j] = beta[j] * d[i,j] ^ -alpha
# beta[j] ~ theta * sediment[j] * exp(dnorm(0,tau_beta))
#
#y[i] ~ sum_j((theta * sediment[j] * exp(dnorm(0,tau_beta))) * d[i,j] ^ -alpha) * exp(dnorm(0,tau))
#
#
# Greek letters are estimated. English letters are fixed data.
#
# y[i] is observations of turbidity at a point in the ocean
# tau is precision for normal error on turbidity observations
# z[i,j] is latent state representing the influence of river j on ocean point i
# beta[j] is sediment yield from river j rescaled to turbidity units. 
# d[i,j] is distance from river j to ocean site i
# alpha is power decay of sediment yield with distance. 
# theta is scaling paramater from sediment[j] at a river to its turbdity
# tau_beta is precision for scaling model. 
#
# In practice we use logging liberally so we can sample from normal distributions


rm(list = ls())

library(rjags)
library(raster)
library(RColorBrewer)

 setwd('/Users/s2973410/Databases/Fiji reefs analysis/wq_models/large_region_data_v3')
# setwd('/Users/s2973410/Databases/Fiji reefs analysis/wq_models/north_data')

#Name to save run under
# runname <- 'large_region_hp_1Apr16_unscaled'
# runname <- 'large_region_hp_3Apr16_unscaled'
runname <- 'large_region_hp_24May16_2alpha'

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

nriv <- nrow(b1)

refriv <- which.max(b1$mnyld_tons)

#Scaling for distances - distances as loaded are in metres
xscale <- 1/1000 #turn into km

#Make sediment loadings numeric (JAGS compile fails if it's a matrix or dataframe)
byld <- as.numeric(b1$mnyld_tons)

#Rescale sediment contributions for using in Bayesian model. Easier to guess initial values then. 
bmean <- byld/byld[refriv]

#
# Specify coasts
#
icoast <- rep(1, nriv)
northcoast <- c(27, 28, 14, 13, 16, 15, 17, 18, 24, 23, 19, 20, 21, 22, 25, 26, 7)
icoast[northcoast] <- 2


#
# Select data inputs for distance and turbidity
#

y <- as.numeric(dat$wq)
icol <- grep(pattern = 'groupID', names(dat))
xdat <- as.numeric(unlist(dat[,icol]))
xdat <- matrix(xdat, ncol = nriv, byrow = F)*xscale

N <- length(y)
print(N)

#
# Rescale y to be in (0, 1) and add a small constant
#

# ymin <- min(y)
# ymax <- max(y-ymin)
# cadd <- 0.0001 #constant so there are no zeros
# yscale <- ((y - ymin)/ymax) + cadd

ymin <- min(y)
cadd <- 0.001 #constant so there are no zeros
yscale <- y-(ymin-cadd)

# lny <- log(yscale) #log response, so errors are normal
lny <- log(yscale)

#
# Estimate starting values for variance and 
#
aguess <- -1.38

bguess <- rep(NA, nriv)
for (iriv in 1:nriv){
	imin <- which(xdat[,iriv] == min(xdat[,iriv]))[1]
	bguess[iriv] <- yscale[imin]/((xdat[imin,iriv]^aguess))
	
}
mod1 <- lm(bguess ~ bmean + 0)
# plot(bmean, bguess); abline(mod1);  summary(mod1)
lnbinit <- as.numeric(log(coef(mod1)))

taubetainit <- 1/var(log(bguess))

#
# Specify initial values
# 
# Initial values derived from fitting model to just the northern region first 

   # inits <- list(alpha = rep(1.38,2), tau = 2, lnb = lnbinit, taubeta = taubetainit)

#Inits for 3 chains
    inits <- list(
   list(alpha = rep(1.38,2), tau = 2, lnb = lnbinit*2, taubeta = taubetainit*2),
  list(alpha = rep(1.38,2), tau = 2, lnb = lnbinit, taubeta = taubetainit) ,
    list(alpha = rep(1.38,2), tau = 2, lnb = lnbinit*0.5, taubeta = taubetainit*0.5)
   )
 

# ***************** #
# SPECIFY BAYESIAN MODEL
# ***************** #

cat("
model {
	#loop through turbidity obs at reef sites
	for (i in 1:N) {
		
		# Contributions of all rivers to turbidity at a reef site
		# 
		# z[i,j] = beta * d ^ -alpha
		#
		for (j in 1:nriv){ 
			z[i,j] <- exp(lnbeta[j])*(d[i,j]^(-alpha[icoast[j]]))
			}
		
		# Likelihood of turbidity observation at a reef site
		#
		# y = sum_j(z[i,]) * exp(error)
		#
		lny[i] ~ dnorm(log(sum(z[i,])), tau)
	}
		
	#
	#Priors
	#
	
	tau ~ dgamma(0.001, 0.001) #precision on turbidity obs
	alpha[1] ~ dgamma(0.001, 0.001) #Prior for power, must be >0
	alpha[2] ~ dgamma(0.001, 0.001) #Prior for power, must be >0
	
	# Heirarchical model for sediment to turbidity scaling param
	# Also wraps up uncertainty in actual sediment loading
	# beta = b * sediment * exp(error)
	#
	# Loop through each river
	for (jriv in 1:nriv){
		lnbeta[jriv] ~ dnorm(lnb + log(sedcontrib[jriv]), taubeta)
		}
		
	#
	# Prior on heirarchical terms
	#
	# Scaling parameter
	lnb ~ dnorm(0, 0.001)
	#Error in scaling
	taubeta ~ dgamma(0.001, 0.001)
	
}
", 
file = 'wqmod_lnorm.bug')

# ***************** #
# RUNNING BAYESIAN MODEL
# ***************** #

# Specify model

ptime <- proc.time()

jags <- jags.model('wqmod_lnorm.bug', 
data = list('lny' = lny,
'N' = N, 
'nriv' = nriv, 
'd' = xdat, 
'sedcontrib' = bmean, 
'icoast' = icoast
), inits = inits,
n.chains = 3)
print(proc.time() - ptime)

#Burn in
update(jags, 5000)
print(proc.time() - ptime)

#Extract samples
niter <- 18000
nthin <- 15
nsamp <- round(niter/nthin)
system.time(mcout <- coda.samples(jags, variable.names=c("tau", "alpha", "lnbeta", "lnb", "taubeta"), n.iter=niter, thin = nthin))

print(proc.time() - ptime)

#
# Save the run
#

save(list = ls(), file = paste(runname,'.RData', sep =''))

# load(paste(runname,'.RData', sep =''))
#
# Quick check summary statistics
#

gelman.diag(mcout)

# load(paste(runname,'.RData', sep =''))

library(corrplot)

dev.new()
corrplot(cor(mcout[[1]]), diag = F, type = 'lower')
cor(mcout[[1]])[1:3,1:3]

hist(mcout[[1]][,1])

(smc <- summary(mcout))
dimnam <- attr(smc$stat, 'dimnames')[[1]]

ibeta <- grep('lnbeta',dimnam)
#Which params are betas? 
a_est <- -smc$statistics[1:2,1]
betaest <- exp(smc$statistics[ibeta,1])
betaCI <-  exp(smc$quantiles[ibeta,c(1,5)])
best <- exp(smc$statistics[grep('lnb', dimnam),1][1])
sdest <- sqrt(1/smc$statistics[grep('tau', dimnam),1][1])

dev.new()
plot(bmean, betaest, xlab = 'prior yield', ylab = 'Predicted yield', ylim = c(0, max(betaCI[,2])))
arrows(bmean, betaCI[,1] , bmean, betaCI[,2] , len = 0)
abline(0,best)
text(bmean, betaest, b1$groupID, pos = 4)


bmat <- matrix(rep(betaest, length(y)), ncol = nriv, byrow = T)
amat <- matrix(rep(a_est[icoast],length(y)), ncol = nriv, byrow = T)
ymat <- bmat *(xdat^ amat)

i <- 15
dev.new(); plot(xdat[,i],ymat[,i])

yvals <- rowSums(ymat)*exp(sdest/2)

cor(yvals, yscale)
cor(log(yvals), lny)

dev.new()
plot(log(yvals), lny); abline(0,1)

plot(yvals, exp(lny), xlim = c(0,1)); abline(0,1)

plot(yvals, exp(lny) - yvals); abline(0,0)
sum((yscale - yvals)^2)/N

#
# Predictive deviance
#
lnypred <- log(yvals)
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
dimnam <- attr(xchain, 'dimnames')[[2]]
ibeta <- grep('lnbeta',dimnam)
ialpha <- grep('alpha',dimnam) 
itaubeta <- grep('taubeta',dimnam) 
itau <- grep('tau',dimnam)
itau <- itau[itau!=itaubeta]


# Posterior deviance for each step
loglik <- matrix(NA, nrow = nstep, ncol = nchain)
nstep <- nrow(xchain)
for (ichain in 1:nchain){
	xchain <- as.matrix(mcout[[ichain]])
	for (i in 1:nstep){
		lnypred <- predictmodel(exp(xchain[i,ibeta]), y, nriv, xdat, -xchain[i, ialpha], icoast)
		loglik[i, ichain] <-  sum(-2*dnorm(lny, mean = lnypred, sd = sqrt(1/xchain[i,itau]), log = T))
	}
}

# Posterior averaged deviance
lnypred <- predictmodel_mean(betaest, y, nriv, xdat, a_est, icoast, sdest)
# lnypred2 <- predictmodel(exp(smc$quantiles[ibeta,3]), y, nriv, xdat,-smc$quantiles[ialpha,3], icoast)
loglik_bar <- sum(2*dnorm(lny, mean = lnypred, sd = sdest, log = T))

# Effective number of params
pD <- mean(loglik) + loglik_bar 

# DIC
mean(loglik) + pD



#
#Plot posterior for beta
#
library(ggplot2)
library(tidyr)

nrows <- length((mcout[[1]][,1]))

ybeta <- data.frame(matrix(NA, nrow = nrows, ncol = nriv-1))

for (i in 1:(nriv-1)){
	ybeta[,i] <- exp(as.numeric(mcout[[1]][,ibeta[i]]))
	}

colnames(ybeta) <- paste('Source', 2:nriv, sep ='_')
ybetas <- gather(ybeta, betavar, Beta)

blines <- data.frame(betavar = unique(ybetas$betavar), bpost = betaest[2:nriv], bprior = bmean[2:nriv])


ggplot(ybetas, aes(x=Beta)) + 
geom_density(fill = 'purple', alpha = 0.2) + 
ylab("Density") + 
facet_wrap( ~ betavar,nrow = 5, scales = 'free') + 
theme_bw() + 
geom_vline(data = blines, aes(xintercept = bpost), color = 'grey20') + 
geom_vline(data = blines, aes(xintercept = bprior), color = 'grey20', lty = 2) + 
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


#
#Plot traces for beta
#
colnms <- attr(mcout[[1]], "dimnames")[[2]]
ibeta <- grep('beta', colnms)

dev.new(width = 10, height = 6)
par(mfrow = c(5,6), mar = c(1,2,1,1))
for (i in 2:nriv){
	ybeta <- exp(as.numeric(mcout[[1]][,ibeta[i-1]]))
	plot(ybeta,type = 'l', xlab = '', ylab ='', ylim = c(0, max(ybeta)))
	
	abline(h = bmean[i], col ='red')
	abline(h = betaest[i], col = 'green')
	if(i==nriv){ 
		legend('topright', legend = c('prior mean','posterior mean'), col = c('red','green'), lty =1, bg = 'white')
		}
	}

#
# Plot traces for other params
#
dev.new(width = 10, height = 6)
par(mfrow = c(1,3), mar = c(1,2,1,1))
alphasamp <-as.numeric(mcout[[1]][,1])
plot(alphasamp,type = 'l', xlab = '', ylab ='')
bsamp <-as.numeric(mcout[[1]][,2])
plot(bsamp,type = 'l', xlab = '', ylab ='')	
plot(alphasamp, bsamp)
cor(alphasamp, bsamp)	