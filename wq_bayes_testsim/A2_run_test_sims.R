# Runs test simulations for varying sd and a
# CJ Brown 21 Mar 2016


rm(list = ls())
library(rjags)
library(raster)
setwd('/Users/s2973410/Code/Fiji reefs analysis/wq_bayes_testsim')

source('A_testsim_functions.R')


#
# Seacape params
#
dscale <- 400 #scale distances
dmin <- 2 #min distance
ext <- c(0,100, 0, 100)
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
Sint <- c(1, 1)

modnam <- 'wq_simtesting.bug'

#Params to vary in sims
S <- rep(1, nx)
sd_beta <- 0.1
b <- 1
n_a <- 3
auni <- seq(-2, -0.5, length.out = n_a)
n_sd <- 3
sdunimult <- c(0.5, 1, 2) #multiples of the mean
paramsrep <- expand.grid(a = auni, sd = sdunimult)
nreps <- nrow(paramsrep)

#
# Generate seascape
#
xout <- genseascape(nc = nc, nr = nr, pos = pos, ext = ext, dscale = dscale, dmin)

# Plot some simulated data to visualise
dev.new(width = 10, height = 6)
par(mfrow = c(3, n_a), mar = c(1,3,1,1))

for (irep in 1:nreps){
	beta_est <- calcBeta(S, b, nx, sd_beta)
	ymean <- calcwq(N, paramsrep[irep,'a'], beta_est, xout$xmat)
	yobs <- simdata(N, ymean = ymean, sd = mean(ymean)*paramsrep[irep,'sd'])
	# with(xout, plot(xmat[,1], yobs))
	r <- xout$r
	r[] <- yobs
	plot(r, asp=NA)
}


#
# Run model, varying a and sdlog
#
ptime <- proc.time()
asave <- NULL
bsave <- NULL
betasave <- NULL
rmserun <- matrix(nrow = niter, ncol = nreps)

for (irep in 1:nreps){
	print(irep)
	a <- paramsrep[irep, 'a']
	sdmult <- paramsrep[irep, 'sd']
	
	betas <- calcBeta(S, b, nx, sd_beta)

	aout <- matrix(NA, ncol = 2, nrow = niter)
	bout <- matrix(NA, ncol = 2, nrow = niter)
	betaout <- matrix(NA, ncol = 2*nx, nrow = niter)

	for (iter in 1:niter){
			
		ymean <- calcwq(N, a, betas, xout$x)
		yobs <- simdata(N, ymean, sd = sdmult*mean(ymean))
		lny <- log(yobs)
		initparams <- list(a = a, tau = 1/((sdmult*mean(ymean))^2), lnb = log(b), taubeta = 1/(sd_beta ^2))
		
		set.seed(iter)
		print(paste('rep', iter))
		modout <- runmod(modnam, N, nx,lny = lny, x = xout$xmat, initparams, S, nchains = 1, niter = 10000, nthin = 5)
	
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
	rm(a)
	rm(betas)
		
	}

proc.time() - ptime

# Save results 

save(list = ls(), file = '/Users/s2973410/Databases/Fiji reefs analysis/wq_models_simtesting/varya_sd.RData')

load('/Users/s2973410/Databases/Fiji reefs analysis/wq_models_simtesting/varya_sd.RData')

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

 ggplot(datout, aes(x = sd, y = a_bias, colour = factor(a))) + geom_line(size = 1.5) + mytheme

 ggplot(datout, aes(x = sd, y = a_cv, colour = factor(a))) + geom_line(size = 1.5) + mytheme# + ylim(0,0.1)

 ggplot(datout, aes(x = sd, y = b_bias, colour = factor(a))) + geom_line(size = 1.5) + mytheme 
 ggplot(datout, aes(x = sd, y = b_cv, colour = factor(a))) + geom_line(size = 1.5) + mytheme 

 ggplot(datout, aes(x = sd, y = beta_bias1, colour = factor(a))) + geom_line(size = 1.5) + mytheme 
 ggplot(datout, aes(x = sd, y = beta_cv1, colour = factor(a))) + geom_line(size = 1.5) + mytheme 

 ggplot(datout, aes(x = sd, y = beta_bias2, colour = factor(a))) + geom_line(size = 1.5) + mytheme 
 ggplot(datout, aes(x = sd, y = beta_cv2, colour = factor(a))) + geom_line(size = 1.5) + mytheme 

 ggplot(datout, aes(x = sd, y = beta_bias3, colour = factor(a))) + geom_line(size = 1.5) + mytheme 
 ggplot(datout, aes(x = sd, y = beta_cv3, colour = factor(a))) + geom_line(size = 1.5) + mytheme 

 ggplot(datout, aes(x = sd, y = rmsedat, colour = factor(a))) + geom_line(size = 1.5) + mytheme 


# source('/Users/s2973410/Code/Fiji reefs analysis/wq_bayes_testsim/A3_run_test_sims_loc.R')



