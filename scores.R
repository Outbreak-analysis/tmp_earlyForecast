scores.fcast <- function(fcast, rel.err = FALSE){
	
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
		if(rel.err){
			tg[tg==0] <- 1 
			ME[i] <- sum(fm/tg-1)/n
			MAE[i] <- sum(abs(fm/tg-1))/n
			MQE[i] <- sum(max(flo/tg-1,0))/n + sum(max(1-fhi/tg,0))/n 
		}
		if(!rel.err){
			ME[i] <- sum(fm-tg)/n
			MAE[i] <- sum(abs(fm-tg))/n
			MQE[i] <- sum(max(flo-tg,0))/n + sum(max(tg-fhi,0))/n 
		}
	}
	df <- data.frame(model=names(fcast),ME,MAE,MQE)
	return(df)
}


calc.scores <- function(res.parallel, rel.err) {
	
	df.tmp <- list()
	for(k in 1:length(res.parallel)){
		if(class(res.parallel[[k]])!="try-error"){
			df.tmp[[k]] <- scores.fcast(fcast = res.parallel[[k]],
										rel.err = TRUE)
			df.tmp[[k]]$mc <- k
		}
	}
	return(do.call('rbind',df.tmp))
}



merge.sum.scores <- function(sc.tmp,CI){
	# Merge all backtesting scenarios:
	sc <- do.call('rbind',sc.tmp)
	
	# Summary statistics:
	scsum <- ddply(sc,c('model','source','modelsyndata',
						'R0','GI.mean','nMCdata'),
				   summarize, 
				   ME.med    = median(ME),
				   ME.lo     = quantile(ME,probs = 0.5+CI/2),
				   ME.hi     = quantile(ME,probs = 0.5-CI/2),
				   MAE.med   = median(MAE),
				   MAE.lo    = quantile(MAE,probs = 0.5+CI/2),
				   MAE.hi    = quantile(MAE,probs = 0.5-CI/2))
	
	return(scsum)
}

plot.scores <-function(scsum){
	
	pdf(file = 'scores-summary.pdf', width=25,height =15)
	
	### OVERALL RANKING
	CI <- 0.8
	z <- ddply(scsum,c('model'),summarize,
			   ME.med2 = median(ME.med), 
			   ME.med.lo = quantile(ME.med, probs = 0.5-CI/2), 
			   ME.med.hi = quantile(ME.med, probs = 0.5+CI/2), 
			   MAE.med2 = median(MAE.med),
			   MAE.med.lo = quantile(MAE.med, probs = 0.5-CI/2), 
			   MAE.med.hi = quantile(MAE.med, probs = 0.5+CI/2))
	g <- ggplot(z)
	g <- g + geom_point(aes(x=MAE.med2,y=ME.med2,shape=model,colour=model),size=8)
	g <- g + geom_segment(aes(x=MAE.med.lo,xend=MAE.med.hi,y=ME.med2,yend=ME.med2,
							  shape=model,colour=model),size=6,alpha=0.3)
	g <- g + geom_segment(aes(x=MAE.med2,xend=MAE.med2,y=ME.med.lo,yend=ME.med.hi,
							  shape=model,colour=model),size=6,alpha=0.3)
	g <- g + geom_hline(yintercept=0,linetype=2,size=2)
	g <- g + scale_x_log10()
	g <- g + ggtitle('Overall Summary - Median of Median')
	plot(g)
	
	### DETAILS BY SCENARIO:
	
	scsum$source2 <- paste(scsum$source,"\nMC",scsum$nMCdata)
	
	# cosmetics
	size.pt   <- 3
	size.seg  <- 2
	alpha.seg <- 0.4
	
	g <- ggplot(scsum)+geom_point(aes(x=MAE.med,y=ME.med,
									  shape = model, colour=model),
								  size = size.pt)
	g <- g + geom_segment(aes(x=MAE.lo,xend=MAE.hi,y=ME.med,yend=ME.med, colour=model),
						  size = size.seg, alpha=alpha.seg)
	g <- g + geom_segment(aes(y=ME.lo,yend=ME.hi,x=MAE.med,xend=MAE.med, colour=model),
						  size = size.seg, alpha=alpha.seg)
	g <- g + facet_wrap(~source2)
	g <- g + geom_hline(yintercept=0,linetype = 2)
	g <- g + scale_x_log10()
	g <- g + theme(strip.text.x = element_text(size = 8))
	plot(g)
	g <- g + facet_wrap(~source2,scales = 'free')
	plot(g)
	
	g.R0 <-  ggplot(scsum)+geom_point(aes(x=factor(R0), y=MAE.med, 
										  colour=model,
										  shape=model),size=size.pt)
	g.R0 <- g.R0 + facet_wrap(~modelsyndata+GI.mean)
	g.R0 <- g.R0 + scale_y_log10()
	g.R0 <- g.R0 + ggtitle('MAE w.r.t. R0 (GI mean faceted)')
	plot(g.R0)
	
	g.GI <-  ggplot(scsum)+geom_point(aes(x=(GI.mean), y=MAE.med, 
										  colour=model,
										  shape=model),size=size.pt)
	g.GI <- g.GI + facet_wrap(~modelsyndata+R0)
	g.GI <- g.GI + scale_y_log10()
	g.GI <- g.GI + ggtitle('MAE w.r.t. GI (R0 faceted)')
	plot(g.GI)
	
	# ME
	g.R0 <-  ggplot(scsum)+geom_point(aes(x=factor(R0), y=ME.med, 
										  colour=model,
										  shape=model),size=size.pt)
	g.R0 <- g.R0 + facet_wrap(~modelsyndata+GI.mean, scales='free')
	g.R0 <- g.R0 + ggtitle('ME w.r.t. R0 (GI mean faceted)')
	g.R0 <- g.R0 + geom_hline(yintercept=0,linetype=2)
	plot(g.R0)
	
	g.GI <-  ggplot(scsum)+geom_point(aes(x = GI.mean, 
										  y = ME.med, 
										  colour = model,
										  shape = model),size=size.pt)
	g.GI <- g.GI + facet_wrap(~modelsyndata+R0, scales='free')
	g.GI <- g.GI + ggtitle('ME w.r.t. GI (R0 faceted)')
	g.GI <- g.GI + geom_hline(yintercept=0,linetype=2)
	plot(g.GI)
	
	
	dev.off()
}



