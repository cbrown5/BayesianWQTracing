# Functions for running test simulations of bayesian data analysis
#CJ Brown 17 Oct 2015



# ********************** #
# Data creation functions
# ********************** #

calcBeta <- function(S, b, nb, sd_beta){
	# Calculates wq at pixels
	beta_est <- log(exp(b*S)*exp(rnorm(nb, mean = 0, sd = sd_beta)))
	return(beta_est)
	}


calcwq <- function(N, a, b, x){
	# Calculates wq at pixels
	brep <- matrix(rep(b, N), ncol = length(b), byrow = T)
	ymean <- rowSums(brep * (x^a))
	return(ymean)
	}


genseascape <- function(nc, nr, pos, ext, dscale=1, dmin = 1){
	#Function to generate a seascape
	nx <- nrow(pos)
	r <- raster(extent(ext), ncol = nc, nrow = nr)
	N <- ncell(r)
	r[] <- NA
	dists <- NULL
	for (i in 1:nx){
		rtemp <- r
		icell <- cellFromXY(rtemp, pos[i,])
		rtemp[icell] <- 1
		dtemp <- distance(rtemp)
		dtemp[dtemp==0] <-0.01
		dists <- c(dists, list(dtemp[]))
		}
	xmat <- matrix(unlist(dists), ncol = nx)
	xmat <- (dscale * xmat/max(xmat))+dmin
	
	#Return data
	return(list(r = r, xmat = xmat))
	}

simdata <- function(N, ymean, sd, cadd = 0.0001){
	# Simulates data based on a base layer
	yobs <- exp(rnorm(N, mean = log(ymean), sd = sd))
	 # ymin <- min(yobs)
	 # ymax <- max(yobs-ymin)
	 # ((yobs - ymin)/ymax) + cadd
	 yobs
	}



# ********************** #
# Modelling functions
# ********************** #

runmod <- function(modelname, N, nx, x, lny, initparams, S, nchains = 1, nadapt = 1000, nupdate = 5000, niter = 20000, nthin = 10, inits = NULL){
	# Randomly generates a dataset and runs JAGS
	inits <- list(alpha = -initparams$a, tau = initparams$tau, lnb = initparams$lnb, taubeta = initparams$taubeta)
	
		jags <- jags.model(modelname, 
	data = list('lny' = lny,
	'N' = N, 
	'nriv' = nx, 
	'd' = x, 
	'sedcontrib' = S
	), inits = inits,
	n.chains = nchains, 
	n.adapt = nadapt)


	#Burn in
	update(jags, nupdate)

	#Extract samples
	mcs <- coda.samples(jags, variable.names=c("tau", "alpha", "lnbeta", "lnb", "taubeta"), n.iter=niter, thin = nthin)
	
	return(mcs)
}

#
# Extract estimates for prediction
#
modextract <- function(smc){
	dimnam <- attr(smc$stat, 'dimnames')[[1]]
	ibeta <- grep('lnbeta', dimnam)
	ialpha <- grep('alpha', dimnam)
	ilnb <- grep('lnb', dimnam)[1]

	a_est <- -smc$statistics[ialpha,1]
	b_est <- exp(smc$statistics[ilnb,1])
	betaest <- exp(smc$statistics[ibeta,1])
	sdest <- sqrt(1/smc$statistics[grep('tau', dimnam),1][1])
	
	a_sdest <- smc$statistics[ialpha,2]
	b_sdest <- exp(smc$statistics[ilnb,2])
	beta_sdest <- exp(smc$statistics[ibeta,2])
	
	
	list(a_est = a_est, b_est = b_est, betaest = betaest, sdest = sdest, a_sdest= a_sdest, b_sdest= b_sdest, beta_sdest= beta_sdest)
	
	}

#
# Predict from model
#

predictmodel <- function(betaval, N, nriv, x, a_est, sdest){
	bmap <- matrix(rep(betaval, N), ncol = nriv, byrow = T)
	y <- bmap * (x ^ a_est)
	y <-  rowSums(y, na.rm = T) *exp(sdest/2) 
	return(y)
}

rmse <- function(ypred, ytrue){
	sqrt(sum((ypred - ytrue)^2)/length(ypred))
	}
	
calcbiasvar <- function(ahat, atrue, asd){
	bias <- (ahat - atrue) / atrue
	variance <- asd/abs(ahat)
	c(bias, variance)
	}	
	
