
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
}

dist.target <- function(fcast){
	M <- length(fcast)
	
	b <- vector(length = M)
	s <- vector(length = M)
	
	for(i in 1:M){
		x <- fcast[[i]]
		tg <- x$target.dat
		n.tg <- length(tg)
		n.f <- length(x$inc.f.m)
		fm <- x$inc.f.m[(n.f-n.tg+1):n.f]
		b[i] <- sum(fm/tg-1)
		s[i] <- sum((fm/tg-1)^2)
	}
	df <- data.frame(model=names(fcast),b,s)
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
	)
	return(PRM)
}