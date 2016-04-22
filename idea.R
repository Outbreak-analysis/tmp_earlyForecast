### 
### IMPLEMENT 'IDEA' MODEL BY FISMAN ET AL.
### Fisman DN, Hauck TS, Tuite AR, Greer AL (2013) 
### An IDEA for Short Term Outbreak Projection: Nearcasting Using the Basic Reproduction Number. PLoS ONE 8(12): e83622. doi:10.1371/journal.pone.0083622
###

idea.inc <- function(t,R0,d,GI,I0=1,CI=0.95) {
	### CALCULATE INCIDENCE WITH IDEA FORMULA
	# GI is the generation interval
	# t is calendar time
	# tau is time in generation interval units
	# This implements the stochastic version
	# of IDEA (assumes Poisson distrib)
	
	tau <- t/GI
	lambda <- I0 * (R0/(1+d)^tau)^tau
	a <- (1-CI)/2
	qlo <- qpois(p=a,lambda = lambda)
	qhi <- qpois(p=1-a,lambda = lambda)
	median <- qpois(p=0.5,lambda = lambda) 
	return( list(mean=lambda,
				 median = median,
				 qlo=qlo, 
				 qhi=qhi) )
}

idea.fit <- function(inc, GI, I0, 
					 ignore.data = 0,
					 do.plot = FALSE) {
	### FIT 'IDEA' MODEL
	
	t0 <- 1:length(inc)
	
	# Ignore the first data points (for example if they are too noisy)
	if(ignore.data>0){
		t <- t0[ignore.data:length(t0)]
		inc <- inc[ignore.data:length(t0)]
	}
	# Work with log incidence
	loginc <- log(inc)
	loginc[is.infinite(loginc)] <- NA
	
	# Do a simple polynomial regression
	# on log incidence:
	fit <- lm(formula = loginc ~ poly(t,degree = 2,raw = TRUE))
	a <- -unname(fit$coefficients[3])
	b <- unname(fit$coefficients[2])
	c <- unname(fit$coefficients[1])
	
	if (a<0){
		# If a<0 then super-exponential growth.
		# that's not in the spirit of the model,
		# so force d to 0
		fit <- lm(formula = loginc ~ t)
		b <- unname(fit$coefficients[2])
		c <- unname(fit$coefficients[1])
		a <- 0
		warning("In IDEA fit, discount forced to 0 because its previous fitted value was negative")
	}
	
	# Confidence intervals:
	ci <- confint(object = fit, level = 0.95)
	c.lo <- ci[1,1] 
	c.hi <- ci[1,2]
	b.lo <- ci[2,1] 
	b.hi <- ci[2,2]
	a.lo <- a.hi <- 0;
	if (a>0){
		a.lo <- -ci[3,1]
		a.hi <- -ci[3,2]
	}
	
	I0 <- exp(c)
	d <- exp(a*GI^2)-1
	R0 <- exp(b*GI)
	
	I0.lo <- exp(c.lo)
	d.lo <- exp(a.lo*GI^2)-1
	R0.lo <- exp(b.lo*GI)
	
	I0.hi <- exp(c.hi)
	d.hi <- exp(a.hi*GI^2)-1
	R0.hi <- exp(b.hi*GI)
	
	if(do.plot){
		par(mfrow=c(1,2))
		plot(t,inc,typ="s")
		inc.fit <- c + b*t -a*t^2
		inc.fit.lo <- c.lo + b.lo*t -a.lo*t^2
		inc.fit.hi <- c.hi + b.hi*t -a.hi*t^2
		plot(t, loginc,
			 main = "Internal fit",
			 pch=16, 
			 ylim=range(inc.fit,inc.fit.lo,inc.fit.hi))
		lines(t,loginc,typ="s")
		lines(t,inc.fit, col="blue",lwd=2)
		lines(t,inc.fit.lo, col="blue")
		lines(t,inc.fit.hi, col="blue")
	}
	return(c(I0=I0, R0=R0, d=d,
			 I0.lo=I0.lo, R0.lo=R0.lo, d.lo=d.lo,
			 I0.hi=I0.hi, R0.hi=R0.hi, d.hi=d.hi))
}

idea.fit.stoch <- function(inc, GI, I0, 
						   ignore.data = 0,
						   do.plot = FALSE){
	# FIT STOCHASTIC IDEA
	
	t0 <- 1:length(inc)
	# Ignore the first data points (for example if they are too noisy)
	if(ignore.data>0){
		t <- t0[ignore.data:length(t0)]
		inc <- inc[ignore.data:length(t0)]
	}
	
	tau <- t/GI
	
	# Calculate likelihood
	# given incidence data
	loglik <- function(x,I0,tau) {
		# x[1]: R0
		# x[2]: d
		lambda <- I0 * (x[1]/(1+x[2])^tau)^tau
		res <- dpois(x=inc, lambda = lambda,log = TRUE)
		return(-res)
	}
	
	# Maximize loglikelihood:
	# with some constraints on parameters:
	lower <- c(0,0)
	upper <- c(20,0.99)
	startval <- c(2,0.01)
	estim <- nlminb(start = startval,
					objective = loglik,
					lower = lower,
					upper = upper, 
					I0=I0,tau=tau)
	R0.fit <- estim$par[1]
	d.fit <- estim$par[2]
	
	return(c(R0=R0.fit, d=d.fit, I0 = I0))
}

idea.forecast <- function(data, 
						  horiz.forecast, # how many time units ahead to forecast
						  GI, 
						  CI, # confidence interval of the forecasts
						  ignore.data= 0, 
						  stoch = TRUE, # Use stochastic(T) or deterministic(F) version of IDEA
						  do.plot = FALSE,
						  data.full = NULL){
	
	n.data <- length(data)
	
	t <- c(1:(n.data+horiz.forecast))
	
	# Fit to data
	if(!stoch){
		x <- idea.fit(inc = data, 
					  GI = GI, 
					  I0 = I0,
					  ignore.data = ignore.data,
					  do.plot = do.plot)
	}
	if(stoch){
		x <- idea.fit.stoch(inc = data, 
							GI = GI,
							I0 = I0,
							ignore.data = ignore.data,
							do.plot = do.plot)
	}
	
	# Calculate incidence from IDEA formula
	# and fitted parameters:
	tmp <- idea.inc(t = t, 
					R0 = x["R0"], 
					d = x["d"], 
					I0 = x["I0"],
					CI = CI,
					GI = GI)
	inc.idea <- tmp[["mean"]]
	
	# Keep only the forecasted points:
	fcast.rng <- (n.data+1):(n.data+horiz.forecast)
	inc.fcast <- inc.idea[fcast.rng]
	inc.fcast.lo <- tmp[["qlo"]][fcast.rng]
	inc.fcast.hi <- tmp[["qhi"]][fcast.rng]
	

	
	if(do.plot){
		par(mfrow=c(1,2))
		
		plot(t, (inc.idea), 
			 main = "IDEA forecast",
			 typ="l", 
			 col="blue", lwd=2,
			 ylim = c(0,(max(inc.idea,data))) )
		points(t[fcast.rng],(inc.fcast),col="blue", lwd=3)
		lines(t[fcast.rng],inc.fcast.lo,col="blue",lty=2)
		lines(t[fcast.rng],inc.fcast.hi,col="blue",lty=2)
		points(1:length(data),(data),pch=16,cex=1)
		lines(1:length(data),(data),typ="s")
		abline(v=ignore.data, lty=2)
		abline(v=last.date, lty=2)
		if(!is.null(data.full)){
			points(x = t[fcast.rng],
				   y = data.full$inc[fcast.rng],
				   pch=15, col="gray")
		}
		grid()
		
		plot(t, log(inc.idea), 
			 main = "IDEA forecast (log scale)",
			 typ="l", 
			 col="blue", lwd=2,
			 ylim = c(0,log(max(inc.idea,data))))
		points(t[fcast.rng],log(inc.fcast),col="blue", lwd=3)
		lines(t[fcast.rng],log(inc.fcast.lo),col="blue",lty=2)
		lines(t[fcast.rng],log(inc.fcast.hi),col="blue",lty=2)
		
		points(1:length(data),log(data),pch=16,cex=1)
		lines(1:length(data),log(data),typ="s")
		abline(v=ignore.data, lty=2)
		abline(v=last.date, lty=2)
		if(!is.null(data.full)){
			points(x = t[fcast.rng],
				   y = log(data.full$inc[fcast.rng]),
				   pch=15, col="gray")
		}
		grid()
	}
	inc.f.m <- c(data, inc.fcast)
	inc.f.hi <- c(data, inc.fcast.hi)
	inc.f.lo <- c(data, inc.fcast.lo)
	
	return(list(inc.f.m  = inc.f.m,
				inc.f.hi = inc.f.hi,
				inc.f.lo = inc.f.lo,
				R0 = x["R0"],
				d = x["d"]) )
}




# - - - - - -  - - - - - - - - - - - - - -  - - - 
# - - - - - -  - - - - - - - - - - - - - -  - - - 

do.test <- FALSE 

if(do.test){
	load("./data/SEmInR_sim.Rdata")
	data.full <- subset(inc.tb, mc==2)
	t <- 1:nrow(data.full)
	
	ignore.data <- 5
	last.date <- 12
	horiz.forecast <- 13
	GI <- 4
	
	data  <- data.full$inc[1:last.date]
	
	x <- idea.forecast(data, 
					   horiz.forecast,
					   GI, 
					   stoch = F,
					   CI = 0.99999,
					   ignore.data, 
					   do.plot = TRUE,
					   data.full = data.full)
	
	zz = x[["inc.fcast.hi"]] - x[["inc.fcast.lo"]]
	
}