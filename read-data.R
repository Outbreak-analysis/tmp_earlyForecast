
find.epi.start <- function(inc,w, thres.rate, doplot = FALSE) {
	
	### Identify the start of an epidemic
	### Try to discard "noisy/fizzly" begining and
	### locate the start of significant outbreak.
	### Use simple linear regression on 'slices'
	### of incidence data
	
	tt <- 1:length(inc)
	K <- vector()
	
	if(doplot){
		plot(tt,inc,typ="s",lwd=3)
		abline(v=t1,lty=2)
	}
	
	for(start in 2:(length(inc)-w-1)){
		# define the 'slice':
		rng <- start:(start+w)
		ts <- tt[rng]-tt[start-1]
		incs <- inc[rng] - inc[start]
		# linear regression:
		m <- lm(incs~ts+0)
		K[start] <- m$coefficients
		if(doplot){
			yy = K[start] * (tt-tt[start]) + inc[start]
			lines(tt[rng],yy[rng],col="red",typ="o")
			text(x=tt[start],y=yy[start],labels = round(K[start],2),pos = 1,col="black",cex = 0.9)
		}
	}
	# Retrieve slopes that are larger than threshold
	idx <- which(K>thres.rate)
	d <- diff(idx)
	# Only the first 2 consecutive slopes
	# above the threshold qualify:
	idx2 <- min(which(d==1))
	TSTART <- idx[idx2]
	
	# If start is when incidence=0
	# then move forward until it's >0
	if(inc[TSTART]==0 & !is.na(TSTART)){
		idx.nozero <- min(which(inc[(TSTART+1):length(inc)]>0))
		TSTART <- TSTART + idx.nozero
	}
	
	if(doplot) abline(v=TSTART,lwd=6,lty=3,col="red")
	
	return(TSTART)
}



read.incidence <- function(filename, # RData file
						   objname, # object storing incidence
						   type, # simulated or real
						   find.epi.start.window = NULL,
						   find.epi.start.thresrate = NULL,
						   truncate.date = NULL,
						   mc.choose = 1
){
	load(filename)
	tmp <- get(objname)
	if(type=="simulated") tmp <- subset(get(objname),mc==mc.choose)
	
	# find the start of the significant growth of the epidemic
	# (ignores the fizzles at the start)
	if(!is.null(find.epi.start.window)){
		tstart <- find.epi.start(inc = tmp$inc,
								 w = find.epi.start.window,
								 thres.rate = find.epi.start.thresrate,
								 doplot = F)
		if(is.na(tstart)) {
			warning(paste("Cannot find start of epidemic growth",filename,"MC:",mc.choose))
			return(NA)
		}
		tmp <- subset(tmp, tb>=tstart)
		tmp$tb <- tmp$tb - tstart+1
	}
	dat.full <- data.frame(t=tmp$tb, inc=tmp$inc)
	dat <- dat.full
	# Truncate
	if(!is.null(truncate.date)) dat <- dat.full[1:truncate.date,]	
	return(list(dat=dat, dat.full=dat.full))
}



read.incidence2 <- function(inc.tb, 
						   type, # simulated or real
						   find.epi.start.window = NULL,
						   find.epi.start.thresrate = NULL,
						   truncate.date = NULL,
						   mc.choose = 1
){
	if(type=="simulated") tmp <- subset(inc.tb,mc==mc.choose)
	
	if(nrow(tmp)==0 | is.null(tmp)) {
		warning(paste("No epidemic found for MC:",mc.choose))
		return(NA)
	}
	
	# find the start of the significant growth of the epidemic
	# (ignores the fizzles at the start)
	if(!is.null(find.epi.start.window)){
		tstart <- find.epi.start(inc = tmp$inc,
								 w = find.epi.start.window,
								 thres.rate = find.epi.start.thresrate,
								 doplot = F)
		if(is.na(tstart)) {
			warning(paste("Cannot find start of epidemic growth MC:",mc.choose))
			return(NA)
		}
		tmp <- subset(tmp, tb>=tstart)
		tmp$tb <- tmp$tb - tstart+1
	}
	dat.full <- data.frame(t=tmp$tb, inc=tmp$inc)
	dat <- dat.full
	# Truncate
	if(!is.null(truncate.date)) dat <- dat.full[1:truncate.date,]	
	return(list(dat=dat, dat.full=dat.full, tstart=tstart))
}


read.incidence3 <- function(inc.tb, 
							type, # simulated or real
							find.epi.start.window = NULL,
							find.epi.start.thresrate = NULL,
							truncate.date = NULL,
							truncate.generation = NULL,
							mc.choose = 1
){
	if(type=="simulated") tmp <- subset(inc.tb,mc==mc.choose)
	
	if(nrow(tmp)==0 | is.null(tmp)) {
		msg <- paste("No epidemic found for MC:",mc.choose,". MC available:")
		msg <- c(msg,paste(unique(inc.tb$mc),collapse = ";"))
		warning(msg)
		return(NA)
	}
	
	# find the start of the significant growth of the epidemic
	# (ignores the fizzles at the start)
	if(!is.null(find.epi.start.window)){
		tstart <- find.epi.start(inc = tmp$inc,
								 w = find.epi.start.window,
								 thres.rate = find.epi.start.thresrate,
								 doplot = F)
		if(is.na(tstart)) {
			warning(paste("Cannot find start of epidemic growth MC:",mc.choose))
			return(NA)
		}
		# Rebase starting time
		# (everything before is ignored)
		tmp <- subset(tmp, tb>=tstart)
		tmp$tb <- tmp$tb - tstart+1
	}
	dat.full <- data.frame(t=tmp$tb, inc=tmp$inc)
	dat <- dat.full
	
	### Truncate
	
	# If truncation date specified:
	if(!is.null(truncate.date)) {
		td <- truncate.date
		dat <- subset(dat.full, t <= td)
	}
	# If truncation specified in terms of generations number:
	if(!is.null(truncate.generation)) {
		td <- round(truncate.generation*inc.tb$GI[1],0)
		dat <- subset(dat.full, t <= td)
	}
	
	return(list(dat = dat, 
				dat.full = dat.full, 
				# These dates are meant
				# to be used on the _original_ data
				# to visualize what was kept:
				tstart = tstart,
				ttrunc = tstart + td))
}


