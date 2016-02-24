library(ggplot2);theme_set(theme_bw())
library(plyr)
library(snowfall)
library(parallel)

use.DC.version.of.EpiEstim <- TRUE  # There is a problem in 'OverallInfectivity' function when data set is short

if(!use.DC.version.of.EpiEstim) library(EpiEstim)
if(use.DC.version.of.EpiEstim) source("EstimationR.R")

source("read-data.R")
source("forecast_early_short.R")
source("forecast_utils.R")


backtest.fcast <- function(RData.file,
						   n.MC = NULL,
						   save.to.file = TRUE) {
	load(RData.file)
	nE.true <- param.synthetic.sim[["nE"]]
	nI.true <- param.synthetic.sim[["nI"]]
	DOL.true <- param.synthetic.sim[["DOL.days"]]
	DOI.true <- param.synthetic.sim[["DOI.days"]]
	
	prm <- read.csv("prm_multi_bcktest.csv",header = FALSE)
	
	
	# Truncation date (synthetic beyond
	# this time is supposed unknown)
	trunc <- prm[prm[,1]=="trunc",2]
	
	# Forecast horizon (time units after last know data)
	horiz.forecast <- prm[prm[,1]=="horiz.forecast",2]
	
	# Generation interval
	bias <- prm[prm[,1]=="GI.bias",2]
	GI.mean <- bias * (DOL.true+DOI.true/2)  # approx!
	GI.stdv <- GI.mean/sqrt((mean(nE.true,nI.true))) # approx!
	
	# This loop performs the forecast 
	# on every synthetic data set.
	# Each forecast is evaluated with 
	# specified metrics.
	# All forecasts are merged in one data frame.
	#
	
	n.cores <- detectCores()
	sfInit(parallel = (n.cores>1), cpu = n.cores)
	sfLibrary(R0)
	if(!use.DC.version.of.EpiEstim) sfLibrary(EpiEstim)
	
	idx.apply <- unique(inc.tb$mc)
	if(!is.null(n.MC)) idx.apply <- idx.apply[1:n.MC]
	
	### Parallel execution:
	sfExportAll()
	res <- sfSapply(idx.apply, 
					simplify = FALSE,
					fcast.wrap, 
					datafilename = RData.file,
					trunc = trunc,
					horiz.forecast = horiz.forecast ,
					GI.mean = GI.mean,
					GI.stdv = GI.stdv ,
					GI.dist = "gamma" ,
					cori.window = 3,
					do.plot = FALSE)
	sfStop()
	df <- do.call("rbind", res)
	
	# Specify the (modified) measures to be plotted:
	df$b <- sign(df$ME)*(abs(df$ME))^(1/4)
	df$s <- df$MAE + 1*df$MQE
	
	# Summarize forecast performance across
	# all synthetic data sets:
	df.m <- ddply(df,c("model"),summarize, 
				  b.m=mean(b), 
				  s.m=mean(s),
				  b.md=median(b), 
				  s.md=median(s),
				  b.lo=quantile(b,probs = 0.1),
				  b.hi=quantile(b,probs = 0.9),
				  s.lo=quantile(s,probs = 0.1),
				  s.hi=quantile(s,probs = 0.9)
	)
	return(list(stat.errors = df.m, 
				param.synthetic.sim = param.synthetic.sim,
				bias = bias,
				n.mc = length(idx.apply))
	)
}


plot.backtest <- function(x) {
	
	df.m <- x[["stat.errors"]]
	prm <- x[["param.synthetic.sim"]]
	bias <- x[["bias"]]
	n.mc <- x[["n.mc"]]
	
	g <- ggplot(df.m)
	g <- g + geom_segment(data = df.m, 
						  aes(x=log(s.lo),xend=log(s.hi),
						  	y=b.md,yend=b.md,
						  	colour=model, shape=model),
						  size=1)
	title <- paste(names(prm),prm,sep="=",collapse = " ; ")
	title <- paste(title, "; GI bias =",bias,"; n.MC=",n.mc)
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
	
	plot(g)
}



t1 <- as.numeric(Sys.time())


cmd <- "ls ./data/*.RData"
flist <- system(command = cmd, intern = TRUE)

pdf("plot_backtest.pdf",width=15,height = 15)

for(i in 1:length(flist)){
	message(paste("data sets:",i,"/",length(flist),flist[i]))
	x <- backtest.fcast(RData.file = flist[i], n.MC=10)
	plot.backtest(x)
}
dev.off()


t2 <- as.numeric(Sys.time())
message(paste("Finished in",round((t2-t1)/60,2),"min"))
