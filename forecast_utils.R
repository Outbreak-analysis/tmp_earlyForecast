
model.colour <- function(model.name){
	col <- "grey"
	if(model.name == "Cori") col <- "black"
	if(model.name == "WalLip") col <- "blue"
	if(model.name == "WhiPag") col <- "green"
	if(model.name == "SeqBay") col <- "red"
	if(model.name == "IDEA") col <- "orange"
	return(col)
}

compare.fcast.early.shiny <- function(fcast, dolog){
	### Compare forecasts for Shiny GUI
	
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
		ymax <- max(ymax, fcast[[i]]$inc.f.hi, tgt.dat)
	
	ymin <- ifelse(dolog,1,0)
	plot(0, cex=0,
		 xlim = c(1,length(fcast[[i]]$inc.f.m)),
		 ylim = c(ymin,ymax),
		 las = 1,
		 log = ifelse(dolog,"y",""),
		 xlab = "time",
		 ylab = "")
	grid(lty=2,col="lightgrey")
	
	nudge <- 0.15*c(0:(length(fcast)-1))
	nudge <- nudge-mean(nudge)
	
	# Plot observed data:
	points(x = 1:n.obs, y = obs.dat, typ="s")
	points(x = 1:n.obs, y = obs.dat, typ="p", pch=16)
	# Plot future data if available:
	if(!is.null(tgt.dat)){
		points(x = frng, y = tgt.dat, 
			   cex = 3,
			   col = "lightgrey",
			   pch = 15)
	}
	lwd.fcast <- 4
	
	for (i in 1:n.model) {
		f.m <- fcast[[i]]$inc.f.m[frng]
		f.hi <- fcast[[i]]$inc.f.hi[frng]
		f.lo <- fcast[[i]]$inc.f.lo[frng]
		
		points(x = frng+nudge[i],
			   y = f.m,
			   col = model.colour(names(fcast[i])), 
			   lwd = lwd.fcast, 
			   cex = 1.3)
		segments(x0 = frng+nudge[i], 
				 x1=frng+nudge[i],
				 y0 = f.lo, y1 = f.hi, 
				 col = model.colour(names(fcast[i])), 
				 lwd = lwd.fcast/2)
	}
	legend(x="topleft",
		   legend = names(fcast),
		   col = sapply(names(fcast),FUN = model.colour),
		   lwd = lwd.fcast, 
		   pch = 1)
}

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

dist.target <- function(fcast, rel.err = FALSE){
	
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


create.model.prm <- function(dat,
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

fcast.wrap <- function(mc, inc.tb, 
					   trunc.date,
					   trunc.generation,
					   horiz.forecast,
					   GI.mean, GI.stdv,
					   GI.dist = "gamma",
					   cori.window = 3,
					   do.plot = FALSE){
	### FORECASTING FUNCTION:
	### CALLS EVERY MODELS
	
	# Read incidence data:
	x <- read.incidence.obj(inc.tb = inc.tb,
							type = "simulated",
							find.epi.start.window = horiz.forecast + 3,
							find.epi.start.thresrate = 0.5,
							truncate.date = trunc.date,
							truncate.generation = trunc.generation,
							mc.choose = mc)
	if(is.na(x)) return(NULL)
	
	dat <- x[["dat"]]
	dat.full <- x[["dat.full"]]
	t.epi.start <- x[["tstart"]]
	t.epi.trunc <- x[["ttrunc"]]
	
	# Set parameters for every models:
	PRM <- create.model.prm(dat,
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
	return(list(df = df.tmp, 
				t.epi.start = t.epi.start,
				t.epi.trunc = t.epi.trunc))
}


fcast.wrap.mc <- function(m, dat.all,ttrunc,horiz.forecast,GI.mean,GI.stdv){
	
	### WRAP FOR PARALLEL CODE IN THE MONTE CARLO LOOP
	
	# Extract one realization:
	dat.full <- subset(dat.all, mc == m)
	# Find practical start of epidemic:
	tstart <- find.epi.start(dat.full$inc, w=4, thres.rate = 0.6)
	# shift times such that t=tstart --> t=1
	dat.no.trunc     <- subset(dat.full, tb >= tstart)
	dat.no.trunc$tb  <- dat.no.trunc$tb - tstart + 1
	# Truncates:
	dat.chopped <- subset(dat.no.trunc, tb <= ttrunc)
	
	# Forecast scores:
	res <- fcast.wrap.new(mc       = m, 
						  dat.full = dat.no.trunc, 
						  dat      = dat.chopped,
						  horiz.forecast = horiz.forecast,
						  GI.mean  = as.numeric(GI.mean), 
						  GI.stdv  = as.numeric(GI.stdv),
						  GI.dist  = "gamma",
						  cori.window = 3,
						  rel.err = T,
						  do.plot = FALSE)
	return(res)
}




fcast.wrap.new <- function(mc, dat.full, dat,
						   horiz.forecast,
						   GI.mean, GI.stdv,
						   GI.dist = "gamma",
						   cori.window = 3,
						   rel.err = FALSE,
						   do.plot = FALSE){
	### FORECASTING FUNCTION:
	### CALLS EVERY MODELS
	
	# Set parameters for every models:
	PRM <- create.model.prm(dat,
							dat.full,
							horiz.forecast ,  
							GI.mean, 
							GI.stdv,
							GI.dist = GI.dist,
							cori.window = cori.window)
	# Forecast:
	fcast <- try(lapply(PRM,
						fcast.inc.early.short,
						do.plot = do.plot),
				 silent = TRUE)
	
	# Merge all results in one data frame:	
	df.tmp <- NULL
	if(class(fcast)!="try-error"){
		df.tmp <- dist.target(fcast,rel.err)
		df.tmp$mc <- mc
	}
	return(df.tmp)
}

get.synthetic.data.db <- function(db.path,
								  country = NULL,
								  disease = NULL,
								  synthetic = NULL,
								  source.keys = NULL,
								  eventtype = NULL,
								  eventtype2 = NULL,
								  social.struct = NULL,
								  save.to.file = TRUE){
	# Get the incidence curves from data base
	# (there are potentially several Monte-Carlo realizations)
	
	# Load synthetic data:
	dat0 <- read.database(db.path,
						  country,
						  disease,
						  synthetic,
						  source.keys,
						  eventtype,
						  eventtype2,
						  social.struct)
	# Reformat (keeps only what's used later on)
	inc.tb <- convert.for.backtest(dat0)
	return(inc.tb)
}


get.synthetic.epi.param <- function(source){
	# Extract parameter information 
	# located in column named 'source'
	prm.name <- get.prm.names.from.source(source.string = source)
	x <- vector()
	for(i in 1:length(prm.name)){
		x[i] <- get.prm.value.from.source(source, prm.name = prm.name[i])[1]
	}
	names(x) <- prm.name
	return(x)
}


get.GI <- function(trueparam){
	### Retrieve GI from model parameters
	###
	seir <- sum(grepl('DOL.days',names(trueparam)))
	if(seir>0) {
		DOL.true <- trueparam['DOL.days']
		DOI.true <- trueparam['DOI.days']
		nI.true <- trueparam['nI']
		nE.true <- trueparam['nI']
		GI.mean <- (DOL.true+DOI.true/2/nI.true*(nI.true+1))
		GI.var <- GI.mean/mean(nI.true,nE.true) # <-- not sure about that!!!
	}
	if(seir==0) {
		GI.mean <- trueparam['GImean']
		GI.var <- trueparam['GIvar']
	}
	res <- as.numeric(c(GI.mean,GI.var))
	return(c(GI.mean=res[1], GI.var=res[2]))
}

get.trunc.time <- function(file,trueparam){
	
	# Parameters for the backtesting:
	prm <- read.csv(file,header = FALSE)
	
	# Truncation date (synthetic beyond
	# this time is supposed unknown)
	trunc.type <- as.character(prm[prm[,1]=="trunc.type",2])
	trunc.date <- as.numeric(as.character(prm[prm[,1]=="trunc.date",2]))
	trunc.gen  <- as.numeric(as.character(prm[prm[,1]=="trunc.gen",2]))
	if(trunc.type=="date") trunc.gen <- NULL
	if(trunc.type=="generation") trunc.date <- NULL
	if(trunc.type!="generation" & trunc.type!="date") trunc.gen <- trunc.date <- NULL
	
	# # Generation interval
	# seir <- sum(grepl('DOL.days',names(trueparam)))
	# if(seir>0) {
	# 	DOL.true <- trueparam['DOL.days']
	# 	DOI.true <- trueparam['DOI.days']
	# 	nI.true <- trueparam['nI']
	# 	GI.mean <- (DOL.true+DOI.true/2/nI.true*(nI.true+1))
	# }
	# if(seir==0)	 GI.mean <- trueparam['GImean']
	GI.mean <- get.GI(trueparam)['GI.mean']
	# truncated time:
	if(is.null(trunc.date)) ttr <- as.numeric(GI.mean*trunc.gen)
	if(is.null(trunc.gen))  ttr <- as.numeric(trunc.date)
	return(ttr)
}

plot.scores <-function(scsum){
	
	pdf(file = 'scores-summary.pdf', width=25,height =15)
	g <- ggplot(scsum)+geom_point(aes(x=MAE.med,y=ME.med,
									  shape = model, colour=model),
								  size = 1.5)
	g <- g + geom_segment(aes(x=MAE.lo,xend=MAE.hi,y=ME.med,yend=ME.med, colour=model),alpha=0.5)
	g <- g + geom_segment(aes(y=ME.lo,yend=ME.hi,x=MAE.med,xend=MAE.med, colour=model),alpha=0.5)
	g <- g + facet_wrap(~source)
	g <- g + geom_hline(yintercept=0,linetype = 2)
	g <- g + scale_x_log10()
	plot(g)
	g <- g + facet_wrap(~source,scales = 'free')
	plot(g)
	
	g.R0 <-  ggplot(scsum)+geom_point(aes(x=factor(R0), y=MAE.med, 
										  colour=model,
										  shape=model))
	g.R0 <- g.R0 + facet_wrap(~modelsyndata+GI.mean)
	g.R0 <- g.R0 + scale_y_log10()
	g.R0 <- g.R0 + ggtitle('MAE w.r.t. R0 (GI mean faceted)')
	plot(g.R0)
	
	g.GI <-  ggplot(scsum)+geom_point(aes(x=(GI.mean), y=MAE.med, 
										  colour=model,
										  shape=model))
	g.GI <- g.GI + facet_wrap(~modelsyndata+R0)
	g.GI <- g.GI + scale_y_log10()
	g.GI <- g.GI + ggtitle('MAE w.r.t. GI (R0 faceted)')
	plot(g.GI)
	
	# ME
	g.R0 <-  ggplot(scsum)+geom_point(aes(x=factor(R0), y=ME.med, 
										  colour=model,
										  shape=model))
	g.R0 <- g.R0 + facet_wrap(~modelsyndata+GI.mean, scales='free')
	g.R0 <- g.R0 + ggtitle('ME w.r.t. R0 (GI mean faceted)')
	g.R0 <- g.R0 + geom_hline(yintercept=0,linetype=2)
	plot(g.R0)
	
	g.GI <-  ggplot(scsum)+geom_point(aes(x = GI.mean, 
										  y = ME.med, 
										  colour = model,
										  shape = model))
	g.GI <- g.GI + facet_wrap(~modelsyndata+R0, scales='free')
	g.GI <- g.GI + ggtitle('ME w.r.t. GI (R0 faceted)')
	g.GI <- g.GI + geom_hline(yintercept=0,linetype=2)
	plot(g.GI)
	
	
	dev.off()
}



#
#
#
#
#
#
#
#

plot.data.old <- function(flist, prm.bcktest.file, backtest){
	### PLOT FULL DATA SETS 
	### AS WELL AS THE WINDOWS OF DATA USED
	
	prm <- read.csv(prm.bcktest.file,header = FALSE)
	
	trunc.type <- as.character(prm[prm[,1]=="trunc.type",2])
	trunc.date <- as.numeric(as.character(prm[prm[,1]=="trunc.date",2]))
	trunc.gen <- as.numeric(as.character(prm[prm[,1]=="trunc.gen",2]))
	if(trunc.type=="date") trunc <- trunc.date
	if(trunc.type=="generation") trunc <- trunc.gen
	
	horiz.forecast <- as.numeric(as.character(prm[prm[,1]=="horiz.forecast",2]))
	
	df <- list()
	ts <- list()
	for(i in 1:length(flist)){
		load(flist[[i]])	
		df[[i]] <- inc.tb
		ds <- gsub(".RData","",flist[[i]])
		ds <- gsub("./data/SEmInR_","",ds)
		df[[i]]$dataset <- ds
		if(!is.na(backtest[[i]])[1]) ts[[i]] <- backtest[[i]]$df.t.bounds
	}
	D <- do.call("rbind",df)
	D.ts <- do.call("rbind",ts)
	D.ts$dataset <- gsub(".RData","",D.ts$datafile)
	D.ts$dataset <- gsub("./data/SEmInR_","",D.ts$dataset )
	
	Davg <- ddply(D,c("tb","dataset"),summarize,
				  inc.m = mean(inc),
				  inc.lo = quantile(inc,probs = 0.05),
				  inc.hi = quantile(inc,probs = 0.95))
	
	Davg1 <- Davg[Davg$inc.m > 0.1, ]
	
	D.ts.avg <- ddply(D.ts,"dataset",summarize,
					  tstart.m = mean(tstart),
					  ttrunc.m = mean(ttrunc)
	)
	
	pdf("plot_data.pdf",width=24,height = 16)
	
	g <- ggplot(Davg1,aes(x=tb,y=inc.m)) + geom_pointrange(aes(ymax=inc.hi,ymin=inc.lo),size=0.2) + geom_line()
	g <- g + geom_vline(data = D.ts.avg,aes(xintercept=tstart.m), colour="red")
	g <- g + geom_vline(data = D.ts.avg,aes(xintercept=ttrunc.m), colour="red")
	g <- g + scale_y_log10()
	g <- g + facet_wrap(~dataset, scales = "free")
	plot(g)
	
	Davg2 <- subset(Davg1,tb <= trunc + horiz.forecast*2)
	g <- ggplot(Davg2) + geom_line(aes(x=tb,y=inc.m))
	g <- g + facet_wrap(~dataset, scales = "free")+scale_y_log10()
	g <- g + geom_text(aes(x=tb,y=inc.m,label=round(inc.m,1)),size=3)
	g <- g + geom_point(data = D.ts,aes(x=tstart,y=0.1),colour="red",shape=2)
	plot(g)
	dev.off()
}

plot.data <- function(db.path,
					  source.keys.vec,
					  country = NULL,
					  disease = NULL,
					  synthetic = NULL,
					  eventtype = NULL,
					  eventtype2 = NULL,
					  social.struct = NULL){
	
	# Load synthetic data:
	dat <- list()
	for(i in 1:length(source.keys.vec)){
		dat[[i]] <- read.database(db.path = db.path,
								  country = country,
								  disease = disease,
								  synthetic = synthetic,
								  source.keys = source.keys.vec[i],
								  eventtype =  eventtype,
								  eventtype2 =  eventtype2,
								  social.struct = social.struct)	
	}
	dat0 <- do.call("rbind", dat)
	
	inc.tb <- convert.for.backtest(dat0)
	df <- subset(inc.tb, mc <9)
	g <- ggplot(df) + geom_line(aes(x=tb, y=inc, colour=factor(mc)))
	g <- g + facet_wrap(~source,scales="free")
	g <- g + scale_y_log10()
	plot(g)
}


backtest.fcast <- function(RData.file,
						   save.to.file = TRUE) {
	### BACKTEST FORECASTING MODELS 
	### FROM PRE-SIMULATED SYNTHETIC DATA
	
	# Load synthetic data:
	load(RData.file)
	nE.true <- param.synthetic.sim[["nE"]]
	nI.true <- param.synthetic.sim[["nI"]]
	DOL.true <- param.synthetic.sim[["DOL.days"]]
	DOI.true <- param.synthetic.sim[["DOI.days"]]
	
	# Parameters for the backtesting:
	prm <- read.csv("prm_multi_bcktest.csv",header = FALSE)
	
	# Truncation date (synthetic beyond
	# this time is supposed unknown)
	trunc.type <- as.character(prm[prm[,1]=="trunc.type",2])
	trunc.date <- as.numeric(as.character(prm[prm[,1]=="trunc.date",2]))
	trunc.gen <- as.numeric(as.character(prm[prm[,1]=="trunc.gen",2]))
	
	if(trunc.type=="date") trunc.gen <- NULL
	if(trunc.type=="generation") trunc.date <- NULL
	if(trunc.type!="generation" & trunc.type!="date") trunc.gen <- trunc.date <- NULL
	
	# Forecast horizon (time units after last know data)
	horiz.forecast <- as.numeric(as.character(prm[prm[,1]=="horiz.forecast",2]))
	# Maximum of Monte Carlo realizations backtested:
	n.MC.max <- as.numeric(as.character(prm[prm[,1]=="n.MC.max",2]))
	
	# Generation interval
	bias <- as.numeric(as.character(prm[prm[,1]=="GI.bias",2]))
	GI.mean <- bias * (DOL.true+DOI.true/2)  # approx!
	GI.stdv <- GI.mean/sqrt((mean(nE.true,nI.true))) # approx!
	
	# This loop performs the forecast 
	# on every synthetic data set.
	# Each forecast is evaluated with 
	# specified metrics.
	# All forecasts are merged in one data frame.
	#
	n.cores <- detectCores()
	sfInit(parallel = (n.cores>1), 
		   cpu = n.cores)
	sfLibrary(R0)
	if(!use.DC.version.of.EpiEstim) sfLibrary(EpiEstim)
	
	idx.apply <- unique(inc.tb$mc)
	message(paste("Synthetic data contains",length(idx.apply),"MC iterations"))
	# Reduce backtesting to specified MC realizations:
	if(n.MC.max>0) {
		idx.apply <- idx.apply[1:n.MC.max]
		message(paste("but not more than",length(idx.apply),"are used."))
		inc.tb <- subset(inc.tb, mc %in% idx.apply)
	}
	
	inc.tb$datafile <- RData.file
	# approximate generation interval length:
	inc.tb$GI <- DOL.true + DOI.true/2
	
	### Parallel execution:
	
	# Object 'all.sim' can be very large
	# when many synthetic data were generated
	# but this object is not used (for now).
	# So, removed from memory before being exported
	# because it can crash the memory
	rm("all.sim")
	
	sfExportAll()
	res0 <- sfSapply(idx.apply, 
					 simplify = FALSE,
					 fcast.wrap, 
					 inc.tb = inc.tb,
					 trunc.date = trunc.date,
					 trunc.generation = trunc.gen,
					 horiz.forecast = horiz.forecast ,
					 GI.mean = GI.mean,
					 GI.stdv = GI.stdv ,
					 GI.dist = "gamma" ,
					 cori.window = 3,
					 do.plot = FALSE)
	sfStop()
	
	res <- list()
	tstart <- list()
	ttrunc <- list()
	for(i in 1:length(res0)) {
		res[[i]] <- res0[[i]]$df
		tstart[[i]] <- res0[[i]]$t.epi.start
		ttrunc[[i]] <- res0[[i]]$t.epi.trunc
	}
	
	# If all df are NULL results, then something
	# went wrong with this data set:
	if(length(res)==0) {
		warning(paste("---> WARNING: backtesting problems with",RData.file))
		return(NA)
	}
	
	# Remove results that gave NULL:
	nullres <- unlist(lapply(res,is.null))
	for(i in 1:length(res)){
		if(nullres[i]) res[[i]] <- NULL # <-- assigning NULL to a list element _removes_ it
		warning(paste("---> WARNING: backtesting problems with",RData.file))
	}
	
	df <- do.call("rbind", res)
	
	df.t.bounds <- data.frame(tstart = unlist(tstart), 
							  ttrunc = unlist(ttrunc),
							  datafile = RData.file)
	
	# If all NULL results, then something
	# went wrong with this data set:
	if(is.null(df)) return(NA)
	
	df.m <- NULL
	
	if(!is.null(df)){
		# Specify the (modified) measures to be plotted:
		df$b <- sign(df$ME)*(abs(df$ME))^(1/4)
		df$s <- df$MAE + 1*df$MQE
		
		# Summarize forecast performance across
		# all synthetic data sets:
		df.m <- ddply(df,c("model"),summarize, 
					  b.m=mean(b, na.rm = TRUE), 
					  s.m=mean(s, na.rm = TRUE),
					  b.md=median(b, na.rm = TRUE), 
					  s.md=median(s, na.rm = TRUE),
					  b.lo=quantile(b,probs = 0.1, na.rm = TRUE),
					  b.hi=quantile(b,probs = 0.9, na.rm = TRUE),
					  s.lo=quantile(s,probs = 0.1, na.rm = TRUE),
					  s.hi=quantile(s,probs = 0.9, na.rm = TRUE)
		)
	}
	return(list(stat.errors = df.m, 
				param.synthetic.sim = param.synthetic.sim,
				bias = bias,
				n.mc = length(unique(df$mc)),
				df.t.bounds = df.t.bounds)
	)
}




backtest.fcast.db <- function(db.path,
							  country = NULL,
							  disease = NULL,
							  synthetic = NULL,
							  source.keys = NULL,
							  eventtype = NULL,
							  eventtype2 = NULL,
							  social.struct = NULL,
							  save.to.file = TRUE) {
	### BACKTEST FORECASTING MODELS 
	### FROM PRE-SIMULATED SYNTHETIC DATA
	### READING DATA FROM DATABASE
	
	# Load synthetic data:
	dat0 <- read.database(db.path,
						  country,
						  disease,
						  synthetic,
						  source.keys,
						  eventtype,
						  eventtype2,
						  social.struct)
	
	inc.tb <- convert.for.backtest(dat0)
	
	# Extract parameter information 
	# located in column named 'source'
	prm.name.raw <- get.prm.names.from.source(source.string = inc.tb$source)
	prm.name <- unique(prm.name.raw)
	param.synthetic.sim <- vector()
	for(i in 1:length(prm.name)){
		param.synthetic.sim[i] <- get.prm.value.from.source(inc.tb$source, 
															prm.name = prm.name[i])[1]
	}
	nE.true  <- get.prm.value.from.source(inc.tb$source, prm.name = "nE")[1]
	nI.true  <- get.prm.value.from.source(inc.tb$source, prm.name = "nI")[1]
	DOL.true <- get.prm.value.from.source(inc.tb$source, prm.name = "DOL.days")[1]
	DOI.true <- get.prm.value.from.source(inc.tb$source, prm.name = "DOI.days")[1]
	GImean   <- get.prm.value.from.source(inc.tb$source, prm.name = "GImean")[1]
	GIvar    <- get.prm.value.from.source(inc.tb$source, prm.name = "GIvar")[1]
	
	# Parameters for the backtesting:
	prm <- read.csv("prm_multi_bcktest.csv",header = FALSE)
	
	# Truncation date (synthetic beyond
	# this time is supposed unknown)
	trunc.type <- as.character(prm[prm[,1]=="trunc.type",2])
	trunc.date <- as.numeric(as.character(prm[prm[,1]=="trunc.date",2]))
	trunc.gen  <- as.numeric(as.character(prm[prm[,1]=="trunc.gen",2]))
	
	if(trunc.type=="date") trunc.gen <- NULL
	if(trunc.type=="generation") trunc.date <- NULL
	if(trunc.type!="generation" & trunc.type!="date") trunc.gen <- trunc.date <- NULL
	
	# Forecast horizon (time units after last know data)
	horiz.forecast <- as.numeric(as.character(prm[prm[,1]=="horiz.forecast",2]))
	# Maximum of Monte Carlo realizations backtested:
	n.MC.max <- as.numeric(as.character(prm[prm[,1]=="n.MC.max",2]))
	
	# Generation interval
	bias <- as.numeric(as.character(prm[prm[,1]=="GI.bias",2]))
	if(!is.na(DOL.true)){
		GI.mean <- bias * (DOL.true+DOI.true/2)  # approx!
		GI.mean.true <- (DOL.true+DOI.true/2)  # approx!
		GI.stdv <- GI.mean/sqrt((mean(nE.true,nI.true))) # approx!
	}
	if(is.na(DOL.true)){
		GI.mean <- bias * GImean
		GI.mean.true <- GImean
		GI.stdv <- sqrt(GIvar)
	}
	
	# This loop performs the forecast 
	# on every synthetic data set.
	# Each forecast is evaluated with 
	# specified metrics.
	# All forecasts are merged in one data frame.
	#
	n.cores <- detectCores()
	sfInit(parallel = (n.cores>1), 
		   cpu = n.cores)
	sfLibrary(R0)
	if(!use.DC.version.of.EpiEstim) sfLibrary(EpiEstim)
	
	idx.apply <- unique(inc.tb$mc)
	message(paste("Synthetic data contains",length(idx.apply),"MC iterations"))
	# Reduce backtesting to specified MC realizations:
	if(n.MC.max>0 & n.MC.max<length(idx.apply)) {
		idx.apply <- idx.apply[1:n.MC.max]
		message(paste("but not more than",length(idx.apply),"are used."))
		inc.tb <- subset(inc.tb, mc %in% idx.apply)
	}
	
	inc.tb$datafile <- source.keys
	
	# approximate generation interval length:
	inc.tb$GI <- GI.mean.true
	
	### Parallel execution:
	
	# Object 'all.sim' can be very large
	# when many synthetic data were generated
	# but this object is not used (for now).
	# So, removed from memory before being exported
	# because it can crash the memory
	rm("all.sim")
	
	sfExportAll()
	res0 <- sfSapply(idx.apply, 
					 simplify = FALSE,
					 fcast.wrap, 
					 inc.tb = inc.tb,
					 trunc.date = trunc.date,
					 trunc.generation = trunc.gen,
					 horiz.forecast = horiz.forecast ,
					 GI.mean = GI.mean,
					 GI.stdv = GI.stdv ,
					 GI.dist = "gamma" ,
					 cori.window = 3,
					 do.plot = FALSE)
	sfStop()
	
	res <- list()
	tstart <- list()
	ttrunc <- list()
	for(i in 1:length(res0)) {
		res[[i]] <- res0[[i]]$df
		tstart[[i]] <- res0[[i]]$t.epi.start
		ttrunc[[i]] <- res0[[i]]$t.epi.trunc
	}
	
	# If all df are NULL results, then something
	# went wrong with this data set:
	if(length(res)==0) {
		warning(paste("---> WARNING: backtesting problems with",source.keys,": Empty simulation!"))
		return(NA)
	}
	
	# Remove results that gave NULL:
	nullres <- unlist(lapply(res,is.null))
	for(i in 1:length(res)){
		if(nullres[i]) res[[i]] <- NULL # <-- assigning NULL to a list element _removes_ it
		warning(paste("---> WARNING: backtesting problems with",source.keys,": NULL simulation!"))
	}
	
	df <- do.call("rbind", res)
	
	df.t.bounds <- data.frame(tstart = unlist(tstart), 
							  ttrunc = unlist(ttrunc),
							  datafile = source.keys)
	
	# If all NULL results, then something
	# went wrong with this data set:
	if(is.null(df)) return(NA)
	
	df.m <- NULL
	
	if(!is.null(df)){
		# Specify the (modified) measures to be plotted:
		df$b <- sign(df$ME)*(abs(df$ME))^(1/4)
		df$s <- df$MAE + 1*df$MQE
		
		# Summarize forecast performance across
		# all synthetic data sets:
		df.m <- ddply(df,c("model"),summarize, 
					  b.m=mean(b, na.rm = TRUE), 
					  s.m=mean(s, na.rm = TRUE),
					  b.md=median(b, na.rm = TRUE), 
					  s.md=median(s, na.rm = TRUE),
					  b.lo=quantile(b,probs = 0.1, na.rm = TRUE),
					  b.hi=quantile(b,probs = 0.9, na.rm = TRUE),
					  s.lo=quantile(s,probs = 0.1, na.rm = TRUE),
					  s.hi=quantile(s,probs = 0.9, na.rm = TRUE)
		)
	}
	return(list(stat.errors = df.m, 
				param.synthetic.sim = param.synthetic.sim,
				bias = bias,
				n.mc = length(unique(df$mc)),
				df.t.bounds = df.t.bounds)
	)
}


make.title <- function(x) {
	prm <- x[["param.synthetic.sim"]]
	bias <- x[["bias"]]
	n.mc <- x[["n.mc"]]
	title <- paste(names(prm),prm,sep="=",collapse = " ; ")
	title <- paste0(title, "; GIbias=",bias,"; MC=",n.mc)
	return(title)
}

make.title.db <- function(x) {
	prm <- as.character(x$df.t.bounds$datafile[1])
	bias <- x[["bias"]]
	n.mc <- x[["n.mc"]]
	title <- paste0(prm, "; GIbias=",bias,"; MC=",n.mc)
	return(title)
}

plot.backtest <- function(x) {
	### PLOT THE RESULTS OF THE BACKTESTS
	
	df.m <- x[["stat.errors"]]
	# 	prm <- x[["param.synthetic.sim"]]
	# 	bias <- x[["bias"]]
	# 	n.mc <- x[["n.mc"]]
	# 	title <- paste(names(prm),prm,sep="=",collapse = " ; ")
	# 	title <- paste(title, "; GI bias =",bias,"; n.MC=",n.mc)
	# 	
	title <- make.title(x)
	g <- ggplot(df.m)
	g <- g + geom_segment(data = df.m, 
						  aes(x=log(s.lo),xend=log(s.hi),
						  	y=b.md,yend=b.md,
						  	colour=model, shape=model),
						  size=1)
	g <- g + ggtitle(title)
	
	g <- g + geom_segment(data = df.m, 
						  aes(x=log(s.md),xend=log(s.md),
						  	y=b.lo,yend=b.hi,
						  	colour=model, shape=model),
						  size=1)
	
	g <- g + geom_point(data = df.m, 
						aes(x=log(s.md),y=(b.md),
							colour=model, shape=model),
						size=4)
	
	g <- g + geom_point(data = df.m, 
						aes(x=log(s.m),y=(b.m),
							colour=model, shape=model),
						size=6)
	
	g <- g + geom_hline(yintercept=0,size=2,alpha=0.5,linetype=2)
	
	return(g)
}

plot.backtest.db <- function(x) {
	### PLOT THE RESULTS OF THE BACKTESTS
	
	df.m <- x[["stat.errors"]]
	# 	prm <- x[["param.synthetic.sim"]]
	# 	bias <- x[["bias"]]
	# 	n.mc <- x[["n.mc"]]
	# 	title <- paste(names(prm),prm,sep="=",collapse = " ; ")
	# 	title <- paste(title, "; GI bias =",bias,"; n.MC=",n.mc)
	# 	
	title <- make.title.db(x)
	g <- ggplot(df.m)
	g <- g + geom_segment(data = df.m, 
						  aes(x=log(s.lo),xend=log(s.hi),
						  	y=b.md,yend=b.md,
						  	colour=model, shape=model),
						  size=1)
	g <- g + ggtitle(title)
	
	g <- g + geom_segment(data = df.m, 
						  aes(x=log(s.md),xend=log(s.md),
						  	y=b.lo,yend=b.hi,
						  	colour=model, shape=model),
						  size=1)
	
	g <- g + geom_point(data = df.m, 
						aes(x=log(s.md),y=(b.md),
							colour=model, shape=model),
						size=4)
	
	g <- g + geom_point(data = df.m, 
						aes(x=log(s.m),y=(b.m),
							colour=model, shape=model),
						size=6)
	
	g <- g + geom_hline(yintercept=0,size=2,alpha=0.5,linetype=2)
	
	return(g)
}

plot.backtest.all <- function(x) {
	### PLOT BACKTESTS FROM ALL DATASETS (using facet_wrap)
	
	# Merge all data frame into a single one:
	df <- list()
	for(i in 1:length(x)){
		if( is.list(x[[i]]) ){
			df[[i]] <- x[[i]]$stat.errors
			df[[i]]$dataset <- make.title.db(x[[i]])
			df[[i]]$R0 <- x[[i]]$param.synthetic.sim$R0
			df[[i]]$DOLI <- x[[i]]$param.synthetic.sim$DOL.days+x[[i]]$param.synthetic.sim$DOI.days
		}
	}
	D <- do.call("rbind",df)
	mean.ok <- ( sum(is.infinite(D$b.m)) + sum(is.infinite(D$s.m)) ) ==0
	
	# --- Plots ---
	pdf(paste0("plot_backtest_all_B",x[[1]]$bias,".pdf"),
		width=28, height=20)
	
	g <- ggplot(D) 
	if(mean.ok) g <- g + geom_point(aes(x=s.m,xend=s.m,y=b.m,yend=b.m,
										colour=model,shape=model),size=1)
	g <- g + geom_point(aes(x=s.md,y=b.md,colour=model,shape=model),size=4)
	g <- g + geom_segment(aes(x=s.lo,xend=s.hi,y=b.md,yend=b.md,colour=model))
	g <- g + geom_segment(aes(x=s.md,xend=s.md,y=b.lo,yend=b.hi,colour=model))
	g <- g + geom_hline(yintercept = 0, linetype=2, colour="black") 
	g <- g + scale_x_log10()
	g <- g + scale_colour_brewer(palette = "Set1")
	g <- g + theme(text = element_text(size = 28),
				   strip.text.x = element_text(size = 10))
	g <- g + xlab("MAQE") + ylab("Bias")
	
	g.ds <- g + facet_wrap(~dataset)
	plot(g.ds)
	
	g2 <- ggplot(D) + theme(text = element_text(size=28))
	
	g.R0.b <- g2 + geom_pointrange(aes(x=factor(R0),
									   y=b.md,
									   ymin=b.lo,
									   ymax=b.hi,
									   colour=model,shape=model),size=1,
								   position =position_dodge(width = 0.3)) 
	g.R0.b <- g.R0.b + facet_grid(~DOLI) + geom_hline(yintercept = 0)+ ylab("Bias")
	plot(g.R0.b)
	
	g.R0.s <- g2 + geom_pointrange(aes(x=factor(R0),
									   y=s.md,
									   ymin=s.lo,
									   ymax=s.hi,
									   colour=model,shape=model),size=1,
								   position =position_dodge(width = 0.3)) 
	g.R0.s <- g.R0.s + facet_grid(~DOLI)
	g.R0.s <- g.R0.s + scale_y_log10()+ ylab("MAQE")
	plot(g.R0.s)
	
	g.DOLI.s <- g2 + geom_pointrange(aes(x=factor(DOLI),
										 y=s.md,
										 ymin=s.lo,
										 ymax=s.hi,
										 colour=model,shape=model),size=1,
									 position =position_dodge(width = 0.3)) 
	g.DOLI.s <- g.DOLI.s + facet_grid(~R0) 
	g.DOLI.s <- g.DOLI.s + scale_y_log10()+ ylab("MAQE")
	plot(g.DOLI.s)
	dev.off()
}

plot.backtest.all.db <- function(x) {
	### PLOT BACKTESTS FROM ALL DATASETS (using facet_wrap)
	
	# Merge all data frame into a single one:
	df <- list()
	for(i in 1:length(x)){
		if( is.list(x[[i]]) ){
			df[[i]] <- x[[i]]$stat.errors
			synthetic_i <- make.title.db(x[[i]])
			df[[i]]$dataset <- synthetic_i
			df[[i]]$R0 <- get.prm.value.from.source(synthetic_i,"R0")  
			df[[i]]$DOLI <- get.prm.value.from.source(synthetic_i,"DOL.days")+get.prm.value.from.source(synthetic_i,"DOI.days")
		}
	}
	D <- do.call("rbind",df)
	
	# --- Plots ---
	
	pdf(paste0("plot_backtest_all_B",x[[1]]$bias,".pdf"),
		width=28, height=20)
	
	g <- ggplot(D) 
	g <- g + geom_point(aes(x=s.md,y=b.md,colour=model,shape=model),size=4)
	g <- g + geom_segment(aes(x=s.lo,xend=s.hi,y=b.md,yend=b.md,colour=model))
	g <- g + geom_segment(aes(x=s.md,xend=s.md,y=b.lo,yend=b.hi,colour=model))
	g <- g + geom_hline(yintercept = 0, linetype=2, colour="black") 
	g <- g + scale_x_log10()
	g <- g + scale_colour_brewer(palette = "Set1")
	g <- g + theme(text = element_text(size = 28),
				   strip.text.x = element_text(size = 10))
	g <- g + xlab("MAQE") + ylab("Bias")
	
	g.ds <- g + facet_wrap(~dataset)
	plot(g.ds)
	
	g2 <- ggplot(D) + theme(text = element_text(size=28))
	
	g.R0.b <- g2 + geom_pointrange(aes(x=factor(R0),
									   y=b.md,
									   ymin=b.lo,
									   ymax=b.hi,
									   colour=model,shape=model),size=1,
								   position =position_dodge(width = 0.3)) 
	g.R0.b <- g.R0.b + facet_grid(~DOLI) + geom_hline(yintercept = 0)+ ylab("Bias")
	plot(g.R0.b)
	
	g.R0.s <- g2 + geom_pointrange(aes(x=factor(R0),
									   y=s.md,
									   ymin=s.lo,
									   ymax=s.hi,
									   colour=model,shape=model),size=1,
								   position =position_dodge(width = 0.3)) 
	g.R0.s <- g.R0.s + facet_grid(~DOLI)
	g.R0.s <- g.R0.s + scale_y_log10()+ ylab("MAQE")
	plot(g.R0.s)
	
	g.DOLI.s <- g2 + geom_pointrange(aes(x=factor(DOLI),
										 y=s.md,
										 ymin=s.lo,
										 ymax=s.hi,
										 colour=model,shape=model),size=1,
									 position =position_dodge(width = 0.3)) 
	g.DOLI.s <- g.DOLI.s + facet_grid(~R0) 
	g.DOLI.s <- g.DOLI.s + scale_y_log10()+ ylab("MAQE")
	plot(g.DOLI.s)
	dev.off()
}



