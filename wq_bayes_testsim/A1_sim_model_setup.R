# Set-up bayesian model of water quality
#CJ Brown 17 Oct 2015
#
# Bayes - model 2
#

setwd('/Users/s2973410/Code/Fiji reefs analysis/wq_bayes_testsim')


#
# base model
#

cat("
model {
	#loop through turbidity obs at reef sites
	for (i in 1:N) {
		
		# Contributions of all rivers to turbidity at a reef site
		# 
		# z[i,j] = beta * d ^ -alpha
		#
		for (j in 1:nriv){ 
			z[i,j] <- exp(lnbeta[j])*(d[i,j]^(-alpha))
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
	alpha ~ dgamma(0.001, 0.001) #Prior for power, must be >0
	
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
file = 'wq_simtesting.bug')



