
dolog <- TRUE

compare.fcast.early.2 <- function(fcast, dolog){
	### Compare forecasts
	
	n.model <- length(fcast)
	
	tgt <- fcast[[1]]$target.dat
	n.tgt <- length(tgt)
	n <- length(fcast[[1]]$inc.f.m)
	frng <- (n-n.tgt+1):n
	n.obs <- length(fcast[[1]]$obs.dat)
	
	obs.dat <- fcast[[1]]$obs.dat
	tgt.dat <- fcast[[1]]$target.dat
	
	ymax <- 0
	for (i in 1:length(fcast)) 
		ymax <- max(ymax, fcast[[i]]$inc.f.hi)
	
	if(dolog){
		obs.dat <- log10(obs.dat)
		tgt.dat <- log10(tgt.dat)
		ymax <- log10(ymax)
	}
	
	plot(0, cex=0,
		 xlim = c(1,length(fcast[[i]]$inc.f.m)),
		 ylim = c(0,ymax),
		 las = 1,
		 xlab = "time",
		 ylab = "")
	grid(lty=2,col="lightgrey")
	
	nudge <- 0.15*c(0:(length(fcast)-1))
	nudge <- nudge-mean(nudge)
	
	# Plot observed data:
	points(x = 1:n.obs, y = obs.dat,typ="s")
	points(x = 1:n.obs, y = obs.dat,typ="p", pch=16)
	# Plot future data if available:
	if(!is.null(tgt.dat)){
		points(x = frng, y = tgt.dat, 
			   cex = 3,
			   col="lightgrey", pch=15)
	}
	
	lwd.fcast <- 4
	
	for (i in 1:n.model) {
		
		f.m <- fcast[[i]]$inc.f.m[frng]
		f.hi <- fcast[[i]]$inc.f.hi[frng]
		f.lo <- fcast[[i]]$inc.f.lo[frng]
		
		if(dolog){
			f.m <- log10(f.m)
			f.hi <- log10(f.hi)
			f.lo <- log10(f.lo)
		}
		
		points(x = frng+nudge[i],
			   y = f.m,
			   col = i, 
			   lwd = lwd.fcast, cex = 1.3)
		segments(x0=frng+nudge[i], x1=frng+nudge[i],
				 y0 = f.lo, y1 = f.hi, col=i, lwd=lwd.fcast/2)
	}
	legend(x="topleft",legend = names(fcast),col=c(1:n.model),lwd=lwd.fcast,pch = 1)
}


# compare.fcast.early.2(fcast, dolog = T)



compare.fcast.early <- function(fcast){
	### Compare forecasts
	
	tgt <- fcast[[1]]$target.dat
	n.tgt <- length(tgt)
	n <- length(fcast[[1]]$inc.f.m)
	frng <- (n-n.tgt+1):n
	
	cex <- 2
	
	nudge <- 0.015*c(0:length(fcast))
	nudge <- nudge-mean(nudge)
	dolog <- TRUE
	if(dolog) tgt <- log(tgt)
	
	yrng <- range(unlist(fcast[[1:3]]))
	
	rng <- list()
	for(i in 1:length(fcast)) 
		rng[[i]]<- range(fcast[[i]][["inc.f.lo"]][frng],fcast[[i]][["inc.f.hi"]][frng])
	yrng <- range(unlist(rng))
	if(dolog) yrng<-log(yrng)
	
	
	for(i in 1:length(fcast)){
		
		f.m <- fcast[[i]]$inc.f.m[frng]
		f.lo <- fcast[[i]]$inc.f.lo[frng]
		f.hi <- fcast[[i]]$inc.f.hi[frng]
		
		if(dolog){
			f.m <- log(f.m)
			f.lo <- log(f.lo)
			f.hi <- log(f.hi)
		}
		
		pch <- 14+i
		
		if(i==1){
			plot(x=tgt+nudge[i],
				 y=f.m,
				 xlab="actual incidence",
				 ylab="forecast incidence",
				 ylim=yrng,
				 xlim=c(min(tgt),1.05*max(tgt)),
				 pch=pch, col=i, typ="o",cex=cex)
			grid()
			abline(0,1,lwd=9,lty=1,col="gray")
		}
		if(i>1){
			lines(x=tgt+nudge[i],y=f.m,pch=pch, col=i, typ="o",cex=cex)
		}
		segments(x0=tgt+nudge[i],x1=tgt+nudge[i],y0=f.lo,y1=f.hi,col=i,lwd=3)
		text(x=tgt[length(tgt)],y=f.m[length(f.m)],labels = names(fcast)[i],col=i,pos=4)
	}
	irng <- c(1:length(fcast))
	legend(x="topleft",col = irng, pch=14+irng, legend = names(fcast))
}

dist.target <- function(fcast){
	
	M <- length(fcast) # Number of models
	
	ME <- vector(length = M) # Mean Error
	MAE <- vector(length = M) # Mean Absobule Error
	MQE <- vector(length = M)   # Mean Quantile Error
	
	for(i in 1:M){
		x <- fcast[[i]]
		tg <- x$target.dat
		n.tg <- length(tg)
		n.f <- length(x$inc.f.m)  
		
		# Mean and quantiles of forecasts
		frng <- (n.f-n.tg+1):n.f
		n <- length(frng) # number of actual forecasts
		fm <- x$inc.f.m[frng]
		flo <- x$inc.f.lo[frng]
		fhi <- x$inc.f.hi[frng]
		ME[i] <- sum(fm/tg-1)/n
		MAE[i] <- sum(abs(fm/tg-1))/n
		MQE[i] <- sum(max(flo/tg-1,0))/n + sum(max(1-fhi/tg,0))/n 
	}
	df <- data.frame(model=names(fcast),ME,MAE,MQE)
	return(df)
}


get.model.prm <- function(dat,
						  dat.full,
						  horiz.forecast ,  
						  GI.mean,GI.stdv,
						  GI.dist,
						  cori.window){
	
	PRM <- list(Cori = list(model = "CoriParam",  
							dat = dat,
							dat.full = dat.full,
							horiz.forecast = horiz.forecast,  
							GI.val = c(GI.mean,GI.stdv),
							cori.window = cori.window),
				
				WalLip = list(model = "WalLip",
							  dat = dat,
							  dat.full = dat.full,
							  horiz.forecast = horiz.forecast,  
							  GI.dist = "gamma",  # gamma, lognormal, weibull
							  GI.val = c(GI.mean,GI.stdv)),
				
				WhiPag = list(model = "WhiPag",
							  dat = dat,
							  dat.full = dat.full,
							  horiz.forecast = horiz.forecast,  
							  GI.dist = "gamma",  # gamma, lognormal, weibull
							  GI.val = c(GI.mean,GI.stdv)),
				
				SeqBay = list(model = "SeqBay",
							  dat = dat,
							  dat.full = dat.full,
							  horiz.forecast = horiz.forecast,  
							  GI.dist = "gamma",  # gamma, lognormal, weibull
							  GI.val = c(GI.mean,GI.stdv))
				,
				
				IDEA = list(model = "IDEA",
							dat = dat,
							dat.full = dat.full,
							horiz.forecast = horiz.forecast,  
							GI.val = GI.mean)
	)
	return(PRM)
}






fcast.wrap <- function(mc, datafilename, 
					   trunc,
					   horiz.forecast,
					   GI.mean, GI.stdv,
					   GI.dist = "gamma",
					   cori.window = 3,
					   do.plot = FALSE){
	
	# Read incidence data:
	x <- read.incidence(filename = datafilename,
						  objname = "inc.tb",
						  type = "simulated",
						  find.epi.start.window = horiz.forecast + 3,
						  find.epi.start.thresrate = 0.5,
						  truncate.date = trunc,
						  mc.choose = mc)
	
	if(is.na(x)) return(NULL)
	
	dat <- x[["dat"]]
	dat.full <- x[["dat.full"]]
	
	# Set parameters for every models:
	PRM <- get.model.prm(dat,
						 dat.full,
						 horiz.forecast ,  
						 GI.mean, 
						 GI.stdv,
						 GI.dist = GI.dist,
						 cori.window = cori.window)
	# Forecast:
	fcast <- try(lapply(PRM,
						fcast.inc.early.short,
						do.plot=do.plot),
				 silent = TRUE)
	
	# Merge all results in one data frame:	
	df.tmp <- NULL
	if(class(fcast)!="try-error"){
		df.tmp <- dist.target(fcast)
		df.tmp$mc <- mc
	}
	return(df.tmp)
}