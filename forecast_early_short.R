library(R0)
library(EpiEstim)


plot.GI <- function(g, GI.dist){
	n <- length(g)
	m <- sum(c(1:n)*g)
	v <- sum((c(1:n)^2)*g) - m^2
	plot(x=1:n, y=g, 
		 type="o",lwd=4,
		 pch=5,
		 xlab="time since infection",ylab="pdf",
		 main=paste0("GI distribution - ",GI.dist))
	abline(v=m, lwd=2, lty=2)
	abline(v=m+sqrt(v)*c(-1,1), lwd=1, lty=2)
	grid()
}

plot.exp.phase <- function(dat){
	# Check exponential phase is more or less present:
	loginc <- log(dat$inc+1)
	tt <-dat$t 
	m <- lm(loginc~tt)
	r2 <- summary(m)$r.squared
	plot(tt,loginc,
		 pch=16,type="o",
		 xlab="time",
		 main=paste0("Check exp phase: R2=",round(r2,3)))
	text(x = tt,y=loginc, labels = dat$inc, pos = 3)
	abline(a = m$coefficients[1], b=m$coefficients[2],lty=2)
	
	if( r2<0.5 || r2 > 0.91) 
		warning("Exponential phase either not satisfied, or too good (reaching peak incidence?)")
	
}

plot.fcast <- function(model,
					   dat,
					   inc.f.m,
					   inc.f.lo,
					   inc.f.hi,
					   Restim,
					   dat.full=NULL,
					   log=FALSE){
	
	tt <-dat$t 
	nf <- length(inc.f.m)
	ntt <- length(tt)
	f.rng <- (ntt+1):nf
	inc <- dat$inc
	ylab <- "incidence"
	
	if(log){
		inc.f.m <- log(inc.f.m)
		inc.f.lo <- log(inc.f.lo)
		inc.f.hi<- log(inc.f.hi)
		inc <- log(inc)
		ylab <- "LOG incidence"
	}
	
	title <- paste(model,"forecast time series\n",Restim)
	
	plot(x= c(tt,(length(tt)+1):nf),
		 y=inc.f.m,
		 cex=2, lwd=2,
		 main = title,
		 xlab="time", ylab=ylab,
		 ylim=c(0,max(inc.f.hi)))
	
	polygon(x=c(f.rng,rev(f.rng)),
			y = c(inc.f.lo[f.rng],rev(inc.f.hi[f.rng])),
			border = NA,
			col = rgb(0,0,0,0.2))
	points(inc,pch=16)
	
	if(!is.null(dat.full)){
		inc.full <- dat.full$inc[f.rng]
		if(log) inc.full <- log(dat.full$inc[f.rng])
		
		points(x=dat.full$t[f.rng], 
			   y=inc.full,
			   pch=3, col="red",lwd=3)
		
		legend(x = "topleft",bty = "n",
			   pch=c(1,3), col=c("black","red"), 
			   pt.cex = 2,
			   legend = c("Forecast","Target"))
	}
	abline(v=tt[ntt],lty=2)
}

plot.fcast.vs.actual <- function(dat,
								 dat.full,
								 inc.f.m,
								 inc.f.lo,
								 inc.f.hi,
								 log=FALSE){
	
	tt <-dat$t 
	nf <- length(inc.f.m)
	ntt <- length(tt)
	f.rng <- (ntt+1):nf
	inc <- dat$inc
	inc.full <- dat.full$inc[f.rng]
	title <- "Incidence Forecast vs Actual"
	if(log){
		inc.f.m <- log(inc.f.m)
		inc.f.lo <- log(inc.f.lo)
		inc.f.hi<- log(inc.f.hi)
		inc <- log(inc)
		inc.full <- log(inc.full)
		title <- "LOG Incidence Forecast vs Actual"
	}
	
	plot(x=inc.full,
		 y=inc.f.m[f.rng],
		 xlab="Actual incidence", 
		 ylab="Forecast incidence",
		 main = title,
		 ylim=range(inc.f.lo[f.rng],inc.f.hi[f.rng]),
		 cex=2,lwd=3)
	
	segments(x0=inc.full, x1=inc.full,
			 y0=inc.f.lo[f.rng], y1=inc.f.hi[f.rng],
			 lwd=3)
	grid()
	abline(0,1,col="red",lty=2,lwd=2)
}

plot.all <- function(model,
					 Restim,
					 dat,
					 dat.full,
					 inc.f.m,
					 inc.f.lo,
					 inc.f.hi){
	
	plot.fcast(model,
			   Restim,
			   dat = dat,
			   inc.f.m = inc.f.m,
			   inc.f.lo = inc.f.lo,
			   inc.f.hi = inc.f.hi,
			   dat.full = dat.full,
			   log=FALSE)
	plot.fcast.vs.actual(dat,
						 dat.full,
						 inc.f.m,
						 inc.f.lo,
						 inc.f.hi,
						 log=FALSE)
	plot.fcast(model,
			   Restim,
			   dat = dat,
			   inc.f.m = inc.f.m,
			   inc.f.lo = inc.f.lo,
			   inc.f.hi = inc.f.hi,
			   dat.full = dat.full,
			   log=TRUE)
	plot.fcast.vs.actual(dat,
						 dat.full,
						 inc.f.m,
						 inc.f.lo,
						 inc.f.hi,
						 log=TRUE)
}


translate.model <- function(x){
	res <- NA
	
	if(x=="WalLip") res <- "EG"
	if(x=="WhiPag") res <- "ML"
	if(x=="SeqBay") res <- "SB"
	if(x=="CoriParam") res <- "ParametricSI"
	if(x=="CoriNonParam") res <- "NonParametricSI"
	if(x=="CoriUncertain") res <- "UncertainSI"
	
	if(is.na(res)) stop(paste("Model name unknown:",x))
	return(res)
}

fcast.inc.early.short <- function(prms, do.plot=FALSE){
	### FORECAST INCIDENCE WITH A CHOICE OF MODELS
	### - EARLY IN THE OUTBREAK (MUST BE EXPONENTIAL PHASE)
	### - FORCAST FOR A SHORT HORIZON ONLY
	
	# Unpack parameters depending on model 
	model <- prms[["model"]]
	pname <- c("dat","dat.full", "horiz.forecast","GI.dist","GI.val","GI.truncate")
	pname.cori <- c("cori.window","cori.mean.prior","cori.std.prior","Std.Mean.SI","Min.Mean.SI",
					"Max.Mean.SI","Std.SI", "Std.Std.SI", "Min.Std.SI", "Max.Std.SI", 
					"n.coriUnc.meanstdv", "n.coriUnc.postR")
	if(grepl("Cori",model)) pname <- c(pname,pname.cori)
	for(p in pname) assign(x=p,value=prms[[p]])
	
	#  ==== Estimate R ====
	
	model2 <- translate.model(model)
	
	# R0 package
	if(model %in% c("WalLip","WhiPag","SeqBay")){
		# Create generation time 
		GI <- generation.time(type = GI.dist, 
							  val = GI.val,
							  truncate = GI.truncate,
							  step = 1 )
		# Estimate R
		
		inc <- dat$inc
		if (model=="SeqBay"){
			# There is a conceptual problem with this model
			# when incidence data = 0, because:
			# I[t+1] ~ Poisson(*I[t])
			# so when I[t]=0, probability(I[t+1]>0)=0
			# (bc the poisson intensity is 0)
			#
			# So, try to go round this issue by
			# artificially nudging incidence away from 0 (e.g. to 1)
			if(any(inc==0)) {
				warning("Incidence was changed for SeqBay method")
				inc[inc==0] <- 1
			}
		}
		R <- estimate.R(epid = inc ,
						GT = GI, 
						methods = model2)
	}
	# EpiEstim package
	if(grepl(pattern = "Cori",x = model)){
		### Reference: Cori 2013 American Journal of Epidemiology
		
		# We just want to estimate R at the 
		# last incidence point available:
		N <- length(dat$inc)
		
		# deal with default values
		if(is.null(cori.mean.prior)) cori.mean.prior <- 5
		if(is.null(cori.std.prior)) cori.std.prior <- 5
		
		R <- EstimateR(I = dat$inc,
					   T.Start = (N-cori.window+1),
					   T.End = N,
					   method = model2,
					   SI.Distr = GI.val,
					   Mean.SI = GI.val[1],
					   Std.SI = GI.val[2],
					   Mean.Prior = cori.mean.prior,
					   Std.Prior = cori.std.prior,
					   plot = FALSE,
					   # Param below for "CoriUncertain" model only:
					   Std.Mean.SI= Std.Mean.SI,
					   Min.Mean.SI=Min.Mean.SI, 
					   Max.Mean.SI=Max.Mean.SI, 
					   Std.Std.SI=Std.Std.SI, 
					   Min.Std.SI=Min.Std.SI, 
					   Max.Std.SI=Max.Std.SI, 
					   n1=n.coriUnc.meanstdv, 
					   n2=n.coriUnc.postR
		)
	}
	
	### Retrieve estimates
	
	if(model %in% c("WalLip","WhiPag")){
		R.m <- R$estimates[[model2]]$R
		R.lo <- R$estimates[[model2]]$conf.int[1]
		R.hi <- R$estimates[[model2]]$conf.int[2]
	}
	if(model %in% c("SeqBay")){
		R.m <- R$estimates[[model2]]$R
		# Take last estimate:
		nr <- length(R.m)
		R.m <- R.m[nr]
		R.lo <- R$estimates[[model2]]$conf.int[nr,1]
		R.hi <- R$estimates[[model2]]$conf.int[nr,2]
	}
	if(grepl("Cori",model)){
		R.m <- R$R$`Mean(R)`
		R.lo <- R$R$`Quantile.0.025(R)`
		R.hi <- R$R$`Quantile.0.975(R)`
	}
	
	
	# ==== Forecast ====
	
	inc.f.m <- inc.f.lo <- inc.f.hi <- dat$inc
	
	# -- Retrieve GI distribution:
	
	if(model %in% c("WalLip","WhiPag","SeqBay")){
		#  Remove day 0 of GI distribution
		if(GI$time[1]==0) g <- GI$GT[-1]  
	}
	if(grepl("Cori",model)){
		if(model=="CoriNonParam") {
			g <- GI.val
			GI.dist <- "empirical"
		}
		if(model=="CoriParam") {
			g <- R$SIDistr
			GI.dist <- "gamma"
		}
		if(model=="CoriUncertain") 
			stop("Model CoriUncertain implementation not finished!")
		
		#  Remove day 0 of GI distribution
		if(g[1,1]==0) g <- g[-1,] 
		g <- g[,2]
	}
	
	n.g <- length(g)
	n.i <- length(inc.f.m)
	
	if(model %in% c("WalLip","WhiPag") || 
	   grepl("Cori",model)){
		# Apply renewal equation at each time step 
		# until forecast horizon:
		for(i in 1:horiz.forecast){
			n.g2 <- min(n.g,n.i)
			inc.rev.m <- rev(inc.f.m)[1:n.g2]
			inc.rev.lo <- rev(inc.f.lo)[1:n.g2]
			inc.rev.hi <- rev(inc.f.hi)[1:n.g2]
			
			# renewal equation:
			inc.f.m <- c(inc.f.m, R.m*sum(g[1:n.g2]*inc.rev.m[1:n.g2]))
			inc.f.lo <- c(inc.f.lo, R.lo*sum(g[1:n.g2]*inc.rev.lo[1:n.g2]))
			inc.f.hi <- c(inc.f.hi, R.hi*sum(g[1:n.g2]*inc.rev.hi[1:n.g2]))
		}
	}
	
	if(model %in% c("SeqBay")){
		# Apply the formula I(t+h) = exp(h*gamma(R-1))*I(t)
		# from Bettencourt 2008 PLoS ONE
		
		tmp.m <- exp( (1:horiz.forecast)*(R.m-1.0)/GI$mean )
		tmp.lo <- exp( (1:horiz.forecast)*(R.lo-1.0)/GI$mean )
		tmp.hi <- exp( (1:horiz.forecast)*(R.hi-1.0)/GI$mean )
		
		inc.f.m <- c(inc.f.m, inc.f.m[n.i] * tmp.m)
		inc.f.lo <- c(inc.f.lo, inc.f.lo[n.i] * tmp.lo)
		inc.f.hi <- c(inc.f.hi, inc.f.hi[n.i] * tmp.hi)
	}
	
	Restim <- paste0(round(R.m,2)," (",
					 round(R.lo,2),";",
					 round(R.hi,2),")")
	
	# ==== Plots ====
	if (do.plot){
		par(mfrow=c(3,2))
		plot.exp.phase(dat)
		plot.GI(g,GI.dist)
		plot.all(model,Restim,dat, dat.full, inc.f.m, inc.f.lo, inc.f.hi)
	}
	
	# Target data (if exists, for testing purpose)
	target.dat <- NULL
	if(!is.null(dat.full)) target.dat <- dat.full$inc[(length(dat$t)+1):length(inc.f.m)]
	
		
	return(list(R=R,
				inc.f.m = inc.f.m,
				inc.f.lo = inc.f.lo,
				inc.f.hi = inc.f.hi,
				target.dat = target.dat
	)
	)
}

