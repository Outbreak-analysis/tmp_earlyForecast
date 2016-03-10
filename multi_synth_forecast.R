###
###   FORECAST MODELS BACKTESTING
###

library(ggplot2);theme_set(theme_bw())
library(gridExtra)
library(plyr)
library(snowfall)
library(parallel)

t1 <- as.numeric(Sys.time())

# - - - - - - - - - - - - - - - - - - - - - - -
# There is a problem in 'OverallInfectivity' 
# function when data set is short.
use.DC.version.of.EpiEstim <- TRUE  
if(!use.DC.version.of.EpiEstim) library(EpiEstim)
if(use.DC.version.of.EpiEstim) source("EstimationR.R")
# - - - - - - - - - - - - - - - - - - - - - - -

source("read-data.R")
source("forecast_early_short.R")
source("forecast_utils.R")


plot.data <- function(flist, prm.bcktest.file, backtest){
	
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
					fcast.wrap3, 
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

make.title <- function(x) {
	prm <- x[["param.synthetic.sim"]]
	bias <- x[["bias"]]
	n.mc <- x[["n.mc"]]
	title <- paste(names(prm),prm,sep="=",collapse = " ; ")
	title <- paste0(title, "; GIbias=",bias,"; MC=",n.mc)
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



plot.backtest.all <- function(x) {
	### PLOT BACKTESTS FROM ALL DATASETS (using facet_wrap)
	
	# Merge all data frame into a single one:
	df <- list()
	for(i in 1:length(flist)){
		if( is.list(x[[i]]) ){
			df[[i]] <- x[[i]]$stat.errors
			df[[i]]$dataset <- make.title(x[[i]])
			df[[i]]$R0 <- x[[i]]$param.synthetic.sim$R0
			df[[i]]$DOLI <- x[[i]]$param.synthetic.sim$DOL.days+x[[i]]$param.synthetic.sim$DOI.days
		}
	}
	D <- do.call("rbind",df)
	mean.ok <- ( sum(is.infinite(D$b.m)) + sum(is.infinite(D$s.m)) ) ==0
	# Plots:
	pdf("plot_backtest_all.pdf",width=28,height=20)
	g <- ggplot(D) 
	if(mean.ok) g <- g + geom_point(aes(x=s.m,xend=s.m,y=b.m,yend=b.m,colour=model,shape=model),size=1)
	g <- g + geom_point(aes(x=s.md,y=b.md,colour=model,shape=model),size=4)
	g <- g + geom_segment(aes(x=s.lo,xend=s.hi,y=b.md,yend=b.md,colour=model))
	g <- g + geom_segment(aes(x=s.md,xend=s.md,y=b.lo,yend=b.hi,colour=model))
	g <- g + geom_hline(yintercept = 0, linetype=2, colour="black") 
	g <- g + scale_x_log10()
	g <- g + scale_colour_brewer(palette = "Set1")
	
	g.ds <- g + facet_wrap(~dataset)
	plot(g.ds)
	
	g.R0.b <- ggplot(D) + geom_pointrange(aes(x=factor(R0),
											  y=b.md,
											  ymin=b.lo,
											  ymax=b.hi,
											  colour=model,shape=model),size=0.7,
									 position =position_dodge(width = 0.3)) 
	g.R0.b <- g.R0.b + facet_grid(~DOLI) + geom_hline(yintercept = 0)
	plot(g.R0.b)
	
	g.R0.s <- ggplot(D) + geom_pointrange(aes(x=factor(R0),
											  y=s.md,
											  ymin=s.lo,
											  ymax=s.hi,
											  colour=model,shape=model),size=0.7,
										  position =position_dodge(width = 0.3)) 
	g.R0.s <- g.R0.s + facet_grid(~DOLI)
	g.R0.s <- g.R0.s + scale_y_log10()
	plot(g.R0.s)
	
	
	g.DOLI.s <- ggplot(D) + geom_pointrange(aes(x=factor(DOLI),
											  y=s.md,
											  ymin=s.lo,
											  ymax=s.hi,
											  colour=model,shape=model),size=0.7,
										  position =position_dodge(width = 0.3)) 
	g.DOLI.s <- g.DOLI.s + facet_grid(~R0) 
	g.DOLI.s <- g.DOLI.s + scale_y_log10()
	plot(g.DOLI.s)
	
	
	dev.off()
}



### --- Run the backtesting ---


# Read all data available:
cmd <- "ls ./data/*.RData"
flist <- system(command = cmd, intern = TRUE)

# Backtest every data sets:
x <- list()
for(i in 1:length(flist)){
	message(paste("data sets:",i,"/",length(flist),flist[i]))
	x[[i]] <- backtest.fcast(RData.file = flist[i])
}
save.image("backtest.RData")

pdf("plot_backtest.pdf",width=15,height = 15)

g <- list()
for(i in 1:length(flist)){
	msg.ok <- 	paste("plotting results from data sets:",i,"/",length(flist),flist[i],": OK")
	msg.fail <-	paste("plotting results from data sets:",i,"/",length(flist),flist[i],": Failed!")
	
	if(is.na(x[[i]])) message(msg.fail)
	
	if(!is.na(x[[i]])){
		g[[i]] <- plot.backtest(x[[i]])
		try.plot <- try(plot(g[[i]]), silent = TRUE)
		if(class(try.plot)!="try-error") message(msg.ok)
		if(class(try.plot)=="try-error") message(msg.fail)
	}
}
dev.off()

try(plot.backtest.all(x), silent=TRUE)

plot.data(flist = flist, 
		  prm.bcktest.file = "prm_multi_bcktest.csv",
		  backtest = x)

# ==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
t2 <- as.numeric(Sys.time())
message(paste("Finished in",round((t2-t1)/60,2),"min"))
