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
use.DC.version.of.EpiEstim <- TRUE  # There is a problem in 'OverallInfectivity' function when data set is short
if(!use.DC.version.of.EpiEstim) library(EpiEstim)
if(use.DC.version.of.EpiEstim) source("EstimationR.R")
# - - - - - - - - - - - - - - - - - - - - - - -

source("read-data.R")
source("forecast_early_short.R")
source("forecast_utils.R")



plot.data <- function(flist, prm.bcktest.file){
	
	prm <- read.csv(prm.bcktest.file,header = FALSE)
	trunc <- prm[prm[,1]=="trunc",2]
	horiz.forecast <- prm[prm[,1]=="horiz.forecast",2]
	
	df <- list()
	for(i in 1:length(flist)){
		load(flist[[i]])	
		df[[i]] <- inc.tb
		ds <- gsub(".RData","",flist[[i]])
		ds <- gsub("./data/SEmInR_","",ds)
		df[[i]]$dataset <- ds
	}
	D <- do.call("rbind",df)
	Davg <- ddply(D,c("tb","dataset"),summarize,inc.m=mean(inc))
	
	pdf("plot_data.pdf",width=24,height = 16)
	g <- ggplot(Davg) + geom_line(aes(x=tb,y=inc.m))
	g <- g + facet_wrap(~dataset, scales = "free_y")
	plot(g)
	Davg2 <- subset(Davg,tb <= trunc + horiz.forecast*2)
	g <- ggplot(Davg2) + geom_line(aes(x=tb,y=inc.m))
	g <- g + facet_wrap(~dataset, scales = "free_y")+scale_y_log10()
	g <- g + geom_text(aes(x=tb,y=inc.m,label=round(inc.m,1)),size=3)
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
	trunc <- prm[prm[,1]=="trunc",2]
	# Forecast horizon (time units after last know data)
	horiz.forecast <- prm[prm[,1]=="horiz.forecast",2]
	# Maximum of Monte Carlo realizations backtested:
	n.MC.max <- prm[prm[,1]=="n.MC.max",2]
	
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
	# Reduce backtesting to specified MC realizations:
	if(n.MC.max>0) idx.apply <- idx.apply[1:n.MC.max]
	
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
		df[[i]] <- x[[i]]$stat.errors
		df[[i]]$dataset <- make.title(x[[i]])
	}
	D <- do.call("rbind",df)
	# Plots:
	pdf("plot_backtest_all.pdf",width=24,height=16)
	g <- ggplot(D) 
	g <- g + geom_point(aes(x=s.m,xend=s.m,y=b.m,yend=b.m,colour=model,shape=model),size=4)
	g <- g + geom_segment(aes(x=s.lo,xend=s.hi,y=b.m,yend=b.m,colour=model))
	g <- g + geom_segment(aes(x=s.m,xend=s.m,y=b.lo,yend=b.hi,colour=model))
	g <- g + geom_hline(yintercept = 0, linetype=2, colour="black") 
	g <- g + scale_x_log10()
	g <- g + scale_colour_brewer(palette = "Set1")
	g <- g + facet_wrap(~dataset)
	plot(g)
	dev.off()
}



### --- Run the backtesting ---


# Read all data available:
cmd <- "ls ./data/*.RData"
flist <- system(command = cmd, intern = TRUE)

plot.data(flist = flist, prm.bcktest.file = "prm_multi_bcktest.csv")


# Backtest every data sets:
g <- list()
x <- list()
pdf("plot_backtest.pdf",width=15,height = 15)
for(i in 1:length(flist)){
	message(paste("data sets:",i,"/",length(flist),flist[i]))
	x[[i]] <- backtest.fcast(RData.file = flist[i])
	g[[i]] <- plot.backtest(x[[i]])
	plot(g[[i]])
}
dev.off()

plot.backtest.all(x)

# ==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
t2 <- as.numeric(Sys.time())
message(paste("Finished in",round((t2-t1)/60,2),"min"))
