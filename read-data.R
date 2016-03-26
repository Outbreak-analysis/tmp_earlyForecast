library(plyr)

# Assume sister folder contains database
# and database utils scripts:
source("../Datsid/read_db.R")

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
			text(x=tt[start],y=yy[start],
				 labels = round(K[start],2),
				 pos = 1,col="black",cex = 0.9)
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



read.incidence.file <- function(filename, # RData file
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
			warning(paste("Cannot find start of epidemic growth",
						  filename,"MC:",mc.choose))
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

read.incidence.obj <- function(inc.tb, 
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



# datevec <- zz2$date
# typevec <- zz2$mc

conv.date.to.duration <- function(datevec, typevec){
	df <- data.frame(date=datevec, type=typevec)	
	df2 <- ddply(df,"type",summarize, m=min(date))
	f <- function(i){return(as.numeric(df2[df2[,1]==i,2]))}
	df$t <- as.numeric(df$date) - sapply(df$type, f) 
	df$t
	return(df$t)
}

convert.for.backtest <- function(x){
	### Converts a query in a simple format
	### used in backtesting codes.
	stopifnot(unique(x$eventtype)=="incidence")
	
	dd <- as.Date(x$reportdate)
	t <- conv.date.to.duration(datevec = dd, 
							   typevec = x$synthetic)
	return(data.frame(date = dd,
					  tb   = as.numeric(t),
					  inc  = as.numeric(x$count),
					  mc   = as.numeric(x$synthetic),
					  source = as.character(x$source))
	)
}

read.database <- function(db.path,
						  country = NULL,
						  disease = NULL,
						  synthetic = NULL,
						  source.keys = NULL,
						  eventtype = NULL,
						  eventtype2 = NULL,
						  social.struct = NULL) {
	### Retrieve all fields of a given query
	
	# First (base) query:
	x <- get.epi.ts(db.path,country,disease,synthetic)
	
	# Then, filter what is asked:
	if(!is.null(source.keys))   x <- x[grepl(pattern = source.keys,x = x$source),]
	if(!is.null(eventtype))     x <- x[x$eventtype %in% eventtype,]
	if(!is.null(eventtype2))    x <- x[x$eventtype2 %in% eventtype2,]
	if(!is.null(social.struct)) x <- x[x$socialstruct %in% social.struct,]
	
	return(x)
}



get.prm.value.from.source <- function(source.string,prm.name) {
	### Retrieve the value of a parameter used to 
	### generate synthetic data.
	### For synthetic data, the 'source' field 
	### has (should have) a string describing these
	### values. 
	### For example "blabla_XYZ_1.23;blabla" indicates
	### the parameter named "XYZ" was "1.23"
	### (the undescore and semi-colon is the expected separator)
	
	# find the position of the parameter name in the string:
	pstart <- as.numeric(gregexpr(prm.name,source.string,fixed=TRUE))
	# discard what's before:
	x2 <- substring(text=source.string, first = pstart)
	# Retrieve positions of all underscores:
	uds <- gregexpr("_",x2,fixed=TRUE)[[1]]
	# Retrieve positions of all semicolon:
	sc <- gregexpr(";",x2,fixed=TRUE)[[1]]
	# Retrieve the value between the next 2 underscores (following param name):
	val <- as.numeric(substr(x2,start = uds[1]+1, stop= sc[1]-1))
	return(val)
}

get.prm.names.from.source <- function(source.string){
	uds <- as.numeric(gregexpr("_",source.string,fixed = T)[[1]])[-1]
	sc <- as.numeric(gregexpr(";",source.string,fixed = T)[[1]])
	prm.name <- substring(text = source.string, first = sc+1, last = uds-1)
	# prm.name <- substr(x = source.string, start = sc+1, stop = uds-1)
	return(prm.name)
}

# DEBUG - - - - 
# db.path <- "../Datsid/a.db"
# country <-  'synthetic'  #  'synthetic' 'LIBERIA'
# disease <-  'synthetic'  #'synthetic' NULL
# synthetic <- NULL
# source.keys <-  "BACKTEST_2" # "SET 2" NULL
# eventtype <- "incidence"  #c("incidence","deaths")
# eventtype2 <- NULL # "confirmed"
# social.struct <- NULL   #"HCW"
# 
# zz <- read.database(db.path,
# 					country,
# 					disease,
# 					synthetic,
# 					source.keys,
# 					eventtype,
# 					eventtype2,
# 					social.struct)
# 
# zz2 <- convert.for.backtest(zz)
# 
# yy <- get.prm.value.from.source(source.string = zz2$source,
# 								prm.name = "nE")
# 
# 
# ss <- get.list.sources(db.path)
# - - - - - - - 

