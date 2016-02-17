library(ggplot2);theme_set(theme_bw())
library(plyr)
library(snowfall)

source("read-data.R")
source("forecast_early_short.R")
source("forecast_utils.R")

save.to.file <- TRUE

t1 <- as.numeric(Sys.time())

load("./data/SEmInR_sim.Rdata")
nE.true <- param.synthetic.sim[["nE"]]
nI.true <- param.synthetic.sim[["nI"]]
DOL.true <- param.synthetic.sim[["DOL.days"]]
DOI.true <- param.synthetic.sim[["DOI.days"]]

### Unpack parameters used for backtesting:
prm <- read.csv("prm_multi_bcktest.csv",header = FALSE)

### Large data set definition

# Number of synthetic data set 
# that will be forecasted:
n.mc <- prm[prm[,1]=="n.mc",2]
n.cores <- prm[prm[,1]=="n.cores",2]

# Truncation date (synthetic beyond
# this time is supposed unknown)
trunc <- prm[prm[,1]=="trunc",2]

# Forecast horizon (time units after last know data)
horiz.forecast <- prm[prm[,1]=="horiz.forecast",2]

# Generation interval
bias <- prm[prm[,1]=="GI.bias",2]
GI.mean <- bias * (DOL.true+DOI.true/2)  # approx!
GI.stdv<- 1

# This loop performs the forecast 
# on every synthetic data set.
# Each forecast is evaluated with 
# specified metrics.
# All forecasts are merged in one data frame.
#
if (save.to.file) pdf("plots.pdf")

sfInit(parallel = (n.cores>1), cpu = n.cores)
sfLibrary(R0)
sfLibrary(EpiEstim)

idx.apply <- 1:n.mc

### Parallel execution:
sfExportAll()
res <- sfSapply(idx.apply, 
				simplify = FALSE,
				fcast.wrap, 
				datafilename = "./data/SEmInR_sim.Rdata",
				trunc = trunc,
				horiz.forecast = horiz.forecast ,
				GI.mean = GI.mean,
				GI.stdv = GI.stdv ,
				GI.dist = "gamma" ,
				cori.window = 3,
				do.plot = TRUE)
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

g <- ggplot(df.m)
g <- g + geom_segment(data = df.m, 
					  aes(x=log(s.lo),xend=log(s.hi),
					  	y=b.md,yend=b.md,
					  	colour=model, shape=model),
					  size=1)

g <- g + ggtitle(paste("GI bias =",bias))

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

g <- g + geom_point(data = df, aes(x=log(s),y=b,
						colour=model, shape=model),
					size=2, alpha=0.2)
plot(g)

if (save.to.file) dev.off()

t2 <- as.numeric(Sys.time())
message(paste("Finished in",round((t2-t1)/60,2),"min"))
