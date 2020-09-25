# Runs test simulations for varying location of rivers and their contributions
# CJ Brown 21 Mar 2016


rm(list = ls())
library(rjags)
library(raster)
setwd('/Users/s2973410/Code/Fiji reefs analysis/wq_bayes_testsim')

source('A_testsim_functions.R')


#
# Seacape params
#
dscale <- 100 #scale distances
dmin <- 2 #min distance
ext <- c(0,100, 0, 400)
nc <- 10
nr <- 30
N <- nc * nr
nx <- 3

spacing <- round((1:nx) * ext[4]/(nx+1))
pos <- matrix(c(rep(0,3),spacing), ncol = 2, byrow = F)

niter <- 20 #number of random trials

#
# WQ params
#
modnam <- 'wq_simtesting.bug'
a <- -2
b <- 1
sd_beta <- 0.1
sd_data_mult <- 30

#
#Params to vary in sims
#
n_S <- 3
Sunq <- cbind(c(1, 2, 4),c(rep(1, n_S)))
Sbase <- rep(1, nx)
n_pos <- 3
posrep <- c(100, 150, 175)

#Create column of var indices
paramsrep <- expand.grid(iS = 1:n_S, ipos = 1:n_pos)
nreps <- nrow(paramsrep)
#
# Generate seascape
#


# Plot some simulated data to visualise
dev.new(width = 10, height = 6)
par(mfrow = c(5, n_S), mar = c(1,3,1,1))

for (irep in 1:nreps){
	thispos <- pos
	thispos[1,2] <- posrep[[paramsrep[irep,'ipos']]]
	
	xout <- genseascape(nc = nc, nr = nr, pos = thispos, ext = ext, dscale = dscale, dmin)	
	Suse <- Sbase
	Suse[1] <- Sunq[[paramsrep[irep,'iS']]]
	beta_est <- calcBeta(Suse, b, nx, sd_beta)
	ymean <- calcwq(N, a, beta_est, xout$xmat)
	yobs <- simdata(N, ymean = ymean, sd = sd_data_mult*mean(ymean))
	
	r <- xout$r
	r[] <- yobs
	plot(r, asp = NA)
	
}


#
# Run model, varying position and bint
#
ptime <- proc.time()

asave <- NULL
bsave <- NULL
betasave <- NULL
rmserun <- matrix(nrow = niter, ncol = nreps)

for (irep in 1:nreps){
	print(irep)
	
	thispos <- pos
	thispos[1,2] <- posrep[[paramsrep[irep,'ipos']]]
	xout <- genseascape(nc = nc, nr = nr, pos = thispos, ext = ext, dscale = dscale, dmin)	
	Suse <- Sbase
	Suse[1] <- Sunq[[paramsrep[irep,'iS']]]
	
	betas <- calcBeta(Suse, b, nx, sd_beta)

	aout <- matrix(NA, ncol = 2, nrow = niter)
	bout <- matrix(NA, ncol = 2, nrow = niter)
	betaout <- matrix(NA, ncol = 2*nx, nrow = niter)

	for (iter in 1:niter){
		set.seed(iter)
		print(paste('rep', iter))
		
		ymean <- calcwq(N, a, betas, xout$x)
		yobs <- simdata(N, ymean, sd = sd_data_mult*mean(ymean))
		lny <- log(yobs)
		initparams <- list(a = a, tau = 1/((sd_data_mult*mean(ymean))^2), lnb = log(b), taubeta = 1/(sd_beta ^2))
		
		set.seed(iter)
		modout <- runmod(modnam, N, nx,lny = lny, x = xout$xmat, initparams, Suse, nchains = 1, niter = 10000, nthin = 5)
		
		# Save bias and variance
		smc <- summary(modout)
		
		paramest <- modextract(smc)
		
		ypred <- predictmodel(paramest$betaest, N, nx, xout$xmat, paramest$a_est, paramest$sdest)
		
		rmserun[iter, irep] <- rmse(ypred, ymean)
		
		aout[iter,] <- calcbiasvar(paramest$a_est, a, paramest$a_sdest)
		bout[iter,] <- calcbiasvar(paramest$b_est, b, paramest$b_sdest)
		betaout[iter,] <- calcbiasvar(paramest$betaest, betas, paramest$beta_sdest)
		# crosscorr(modout)
		
		#Clean up
		rm(smc)
		rm(modout)
		
	}
		
	#Save results
	asave <- c(asave, list(aout))
	bsave <- c(bsave, list(bout))	
	betasave <- c(betasave, list(betaout))	
	
	#Clean up
	rm(aout)
	rm(bout)
	rm(betaout)
	rm(ymean)
	rm(betas)
	rm(thispos)
	rm(Suse)
	
	}

proc.time() - ptime

# Save results 

save(list = ls(), file = '/Users/s2973410/Databases/Fiji reefs analysis/wq_models_simtesting/varybint_pos.RData')

load('/Users/s2973410/Databases/Fiji reefs analysis/wq_models_simtesting/varybint_pos.RData')


# source('/Users/s2973410/Code/Fiji reefs analysis/wq_bayes_model_v4/B3_Turbidity_model_power_hierarchical.R')

#
# Summary stats
#
 datout <- paramsrep

datout$a_bias <- unlist(lapply(asave, function(x, i) mean(x[,i]), 1))
 datout$a_cv <- unlist(lapply(asave, function(x, i) mean(x[,i]), 2))

 datout$b_bias <- unlist(lapply(bsave, function(x, i) mean(x[,i]), 1))
 datout$b_cv <- unlist(lapply(bsave, function(x, i) mean(x[,i]), 2))

 datout$beta_bias1 <- unlist(lapply(betasave, function(x, i) mean(x[,i]), 1))
  datout$beta_bias2 <- unlist(lapply(betasave, function(x, i) mean(x[,i]), 2))
 datout$beta_bias3 <- unlist(lapply(betasave, function(x, i) mean(x[,i]), 3))

 datout$beta_cv1 <- unlist(lapply(betasave, function(x, i) mean(x[,i]), 4))
 datout$beta_cv2 <- unlist(lapply(betasave, function(x, i) mean(x[,i]), 5))
 datout$beta_cv3 <- unlist(lapply(betasave, function(x, i) mean(x[,i]), 6))


#RMSE

rmsedat <- colMeans(rmserun)
paramsrep$rmsedat <- rmsedat

# #
# # Plots
# #
 library(ggplot2)

mytheme <-   theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

 ggplot(datout, aes(x = iS, y = a_bias, colour = factor(ipos))) + geom_line(size = 1.5) + mytheme

 ggplot(datout, aes(x = iS, y = a_cv, colour = factor(ipos))) + geom_line(size = 1.5) + mytheme# + ylim(0,0.1)

 ggplot(datout, aes(x = iS, y = b_bias, colour = factor(ipos))) + geom_line(size = 1.5) + mytheme 
 ggplot(datout, aes(x = iS, y = b_cv, colour = factor(ipos))) + geom_line(size = 1.5) + mytheme 

 ggplot(datout, aes(x = iS, y = beta_bias1, colour = factor(ipos))) + geom_line(size = 1.5) + mytheme 
 ggplot(datout, aes(x = iS, y = beta_cv1, colour = factor(ipos))) + geom_line(size = 1.5) + mytheme 

 ggplot(datout, aes(x = iS, y = beta_bias2, colour = factor(ipos))) + geom_line(size = 1.5) + mytheme 
 ggplot(datout, aes(x = iS, y = beta_cv2, colour = factor(ipos))) + geom_line(size = 1.5) + mytheme 

 ggplot(datout, aes(x = iS, y = beta_bias3, colour = factor(ipos))) + geom_line(size = 1.5) + mytheme 
 ggplot(datout, aes(x = iS, y = beta_cv3, colour = factor(ipos))) + geom_line(size = 1.5) + mytheme 

 ggplot(datout, aes(x = iS, y = rmsedat, colour = factor(ipos))) + geom_line(size = 1.5) + mytheme 








