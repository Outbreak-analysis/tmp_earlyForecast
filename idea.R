### 
### IMPLEMENT 'IDEA' MODEL BY FISMAN ET AL.
###

idea.inc <- function(t,R0,d,GI,I0=1) {
	### CALCULATE INCIDENCE WITH IDEA FORMULA
	
	# GI is the generation interval
	# t is calendar time
	# tau is time in generation interval units
	tau <- t/GI
	
	return(I0 * (R0/(1+d)^tau)^tau )
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


idea.forecast <- function(data, 
						  horiz.forecast,
						  GI, 
						  ignore.data= 0, 
						  do.plot = FALSE){
	
	t <- c(1:horiz.forecast)
	
	stopifnot(horiz.forecast>length(data))
	
	x <- idea.fit(inc=data, 
				  GI=GI, 
				  ignore.data = ignore.data,
				  do.plot = do.plot)
	
	inc.fcast <- idea.inc(t = t, 
						  R0 = x["R0"], 
						  d = x["d"], 
						  I0 = x["I0"],
						  GI = GI)
	
	if(do.plot){
		par(mfrow=c(1,2))
		plot(t,log(inc.fcast), typ="o", col="blue", lwd=2,
			 ylim = c(0,log(max(inc.fcast,data))))
		points(1:length(data),log(data),pch=16,cex=1)
		lines(1:length(data),log(data),typ="s")
		abline(v=ignore.data, lty=2)
		abline(v=last.date, lty=2)
		
		plot(t,(inc.fcast), typ="o", col="blue", lwd=2,
			 ylim = c(0,(max(inc.fcast,data))))
		points(1:length(data),(data),pch=16,cex=1)
		lines(1:length(data),(data),typ="s")
		abline(v=ignore.data, lty=2)
		abline(v=last.date, lty=2)
	}
	
	return(inc.fcast)
}


# - - - - - -  - - - - - - - - - - - - - -  - - - 

load("./data/SEmInR_sim.Rdata")
data.full <- subset(inc.tb, mc==2)
t <- 1:nrow(data.full)

ignore.data <- 6
last.date <- 15
horiz.forecast <- last.date + 7

data  <- data.full$inc[1:last.date]

GI <- 4

idea.forecast(data, 
			  horiz.forecast,
			  GI, 
			  ignore.data, 
			  do.plot = TRUE)



# 
# x <- idea.fit(inc=inc, GI=GI, 
# 			  ignore.data = ignore.data,
# 			  do.plot = T)
# 
# print(x)
# inc.fcast <- idea.inc(t = t, 
# 					  R0 = x["R0"], 
# 					  d = x["d"], 
# 					  I0 = x["I0"],
# 					  GI = GI)
# 
# inc.fcast.lo <- idea.inc(t = t, 
# 						 R0 = x["R0.lo"], 
# 						 d = x["d.lo"], 
# 						 I0 = x["I0.lo"],
# 						 GI = GI)
# 
# inc.fcast.hi <- idea.inc(t = t, 
# 						 R0 = x["R0.hi"], 
# 						 d = x["d.hi"],
# 						 I0 = x["I0.hi"],
# 						 GI = GI)
# 
# par(mfrow=c(1,1))
# plot(t,log(data.full$inc), typ="s",
# 	 ylim = c(0,log(max(inc.fcast,data.full$inc))))
# points(t,log(data.full$inc),pch=16,cex=0.6)
# lines(t,log(inc.fcast),col="blue",lwd=2)
# lines(t,log(inc.fcast.lo),col="blue", lty=3)
# lines(t,log(inc.fcast.hi),col="blue", lty=3)
# abline(v=ignore.data, lty=2)
# abline(v=last.date, lty=2)
# 
