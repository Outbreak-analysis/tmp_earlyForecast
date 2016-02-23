library(plyr)
library(adaptivetau) 

lappend <- function (lst, ...){
	lst <- c(lst, list(...))
	return(lst)
}

trans.SEmInR <- function(nE,nI){
	###
	### Generate the list of transitions for SEmInR 
	###
	
	# infection:
	z <- list(c(S=-1, E1=1))
	
	# transition through E_k compartments
	for(i in 2:nE){
		tmp <- c(-1,1)
		names(tmp) <- c(paste0("E",i-1),paste0("E",i))
		z <- lappend(z,tmp)
	}
	
	# infectiousness triggered:
	tmp <- c(-1,1)
	names(tmp) <-  c(paste0("E",nE),"I1")
	z <- lappend(z,tmp)
	
	# transition through I_k compartments
	for(i in 2:nI){
		tmp <- c(-1,1)
		names(tmp) <- c(paste0("I",i-1),paste0("I",i))
		z <- lappend(z,tmp)
	}
	
	# Recovery
	tmp <- c(-1,1)
	names(tmp) <-  c(paste0("I",nI),"R")
	z <- lappend(z,tmp)
	
	return(z)
}

trans.rate.SEmInR <- function(x, params, t){
	###
	### Generate transition rates for SEmInR 
	###
	###   * WARNING * : same order as transitions!
	
	# total population
	N <- sum(x)
	
	# prevalence infectious
	idx.I <- which(grepl(pattern = "I",x = names(x)))
	Itot <- sum(x[idx.I])
	
	# infection rate
	inf.rate <- params$beta*x["S"]*Itot/N*(x["S"]>=1)
	
	# transition through E_k
	trans.E <- numeric(params$nE)
	for(i in 1:params$nE){
		trans.E[i] <- params$sigma*params$nE * x[paste0("E",i)]*(x[paste0("E",i)]>=1)
	}
	
	# transition through I_k
	trans.I <- numeric(params$nI)
	for(i in 1:params$nI){
		trans.I[i] <- params$sigma*params$nI * x[paste0("I",i)]*(x[paste0("I",i)]>=1)
	}
	
	# Gather all transition rates that match
	# transition events defined elsewhere:
	return(c(inf.rate,trans.E,trans.I))
}


simul.SEmInR <- function(horizon.years,# horizon of the simulation 
						 DOL.days,     # duration of latency in DAYS
						 DOI.days,     # duration of infectiousness in DAYS
						 R0,
						 pop.size,     # population size (Population size needs to be very large for convergence stochastic -> deterministic)
						 nE,           # Number of "exposed" compartment for Erlang distribution
						 nI,           # Number of "infectious" compartment for Erlang distribution
						 I.init,       # Initial number of infectious individuals
						 n.MC = 1,
						 time.bucket = 1/365,     # Aggregation of all events in a time bucket unit
						 remove.fizzles = FALSE,
						 thres.fizz = 0.1,        # Threshold (fraction of max incidence) to identify fizzles
						 do.adaptivetau = TRUE,
						 epsilon = 0.05,
						 seed = 1234,
						 save.to.Rdata.file = TRUE 
						 ) {
	
	### Simulate several epidemics with 
	### a SEmInR stochastic model.
	### Returns a data frame with all simulations.
	
	set.seed(seed)
	
	# SEIR parameters
	gamma <- 365/DOL.days
	sigma <- 365/DOI.days
	
	# R0(SEmInR) is a complicated expression (see Feng et al, Epidemiological models with nonexponentially distributed disease stages and applications to disease control. Bulletin of Mathematical Biology, 69(5):1511–1536, 2007)
	# but when the DOI and DOL 
	# very small compared to average life of host (e.g. 80 years for humans)
	# then good approx with:
	# R0(SEmInR) ~ DOI * beta   (beta:contact rate)
	beta <- R0*gamma
	
	params <- list(beta=beta, 
				   sigma=sigma,
				   gamma=gamma, 
				   nE=nE, nI=nI)
	
	# Initial values
	I0 <- I.init   
	x0 <- c(pop.size-I0,
			rep(0,nE), 
			I0,
			rep(0,nI-1),
			0)
	names(x0) <- c("S",paste0("E",1:nE),paste0("I",1:nI),"R")
	
	# Run all Monte Carlo iterations of simulations 
	for(mc in 1:n.MC){
		message(paste("MC",mc,"/",n.MC))
		if (!do.adaptivetau){
			res.ATAU <- ssa.exact(init.values = x0,
								  transitions = trans.SEmInR(nE,nI), 
								  rateFunc = trans.rate.SEmInR, 
								  params = params,
								  tf=horizon.years)
		}
		if(do.adaptivetau){
			res.ATAU <- ssa.adaptivetau(init.values = x0,
										transitions = trans.SEmInR(nE,nI), 
										rateFunc = trans.rate.SEmInR, 
										params = params,
										tl.params = list(epsilon=epsilon),  
										tf=horizon.years)
		}
		# Store in data frame:
		idx.E <- grepl("E",colnames(res.ATAU))
		idx.I <- grepl("I",colnames(res.ATAU))
		idx.IE <- as.logical(idx.E+idx.I)
		idx.S <- which(grepl("S",colnames(res.ATAU)))
		tmp <- data.frame(t = res.ATAU[,1], 
						  inc = c(I0,-diff(res.ATAU[,idx.S])),
						  prev = apply(res.ATAU[,idx.IE],1,sum),
						  prevI = apply(res.ATAU[,idx.I],1,sum),
						  mc = rep(mc,nrow(res.ATAU)))
		
		# Aggregate in time buckets
		tmp$tb <- ceiling(tmp$t/time.bucket)
		tmp$tb[tmp$tb==0] <- 1
		
		# Store all  aggregated results in dataframe:
		if(mc==1) all.sim <- tmp
		if(mc>1) all.sim <- rbind(all.sim, tmp)
	}
	
	# Remove fizzles:
	if(remove.fizzles){
		 
		sim.nofizz = ddply(all.sim,.variables = c("mc"),
						   summarize,
						   i.max = max(prev))
		imax.all <- max(sim.nofizz$i.max)
		
		mc.nofizz <- which(sim.nofizz$i.max>thres.fizz*imax.all)
		p.fizz <- 1-length(mc.nofizz)/n.MC
		
		all.sim.nofizz <- all.sim[all.sim$mc %in% mc.nofizz,]
		all.sim <- all.sim.nofizz
	}
	
	# Incidence only at time buckets
	inc.tb <- ddply(all.sim, c("tb","mc"), summarize, inc=sum(inc))
	
	# All params saved here:
	param.synthetic.sim <- list(DOL.days = DOL.days,
								DOI.days = DOI.days,
								R0 = R0,
								nE = nE,
								nI = nI,
								I0 = I0,
								pop.size = pop.size)
	
	if(save.to.Rdata.file){
		fname <- paste0("SEmInR_",
						"DOL_",DOL.days,"_",
						"DOI_",DOI.days,"_",
						"R_",R0,"_",
						"nE_",nE,"_",
						"nI_",nI,"_",
						"popsize_",pop.size/1000,"k",
						".RData")
		save(list=c("all.sim","inc.tb","param.synthetic.sim"),
			 file = fname)	
		message(paste("objects saved in",fname))
	}
	
	# Objects returned:
	return(list(inc = inc.tb, 
				param = param.synthetic.sim)
				)

}