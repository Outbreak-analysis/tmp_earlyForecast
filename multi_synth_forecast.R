library(ggplot2);theme_set(theme_bw())
library(plyr)
source("read-data.R")
source("forecast_early_short.R")
source("forecast_utils.R")

save.to.file <- FALSE


### Large data set definition

n.mc <- 9
trunc <- 10

### Model choice and associated parameters:
# Models are:
# WalLip  WhiPag  SeqBay 
# CoriParam CoriNonParam CoriUncertain

horiz.forecast <- 7

GI.mean <- 2.3
GI.stdv<- 1


for(mc in 1:n.mc){
	# Read incidence data:
	dat <- read.incidence(filename = "./data/SEmInR_sim.Rdata",
						  objname = "inc.tb",
						  type = "simulated",
						  truncate.date = trunc,
						  mc.choose = mc)
	dat.full <- read.incidence(filename = "./data/SEmInR_sim.Rdata",
							   objname = "inc.tb",
							   type = "simulated",
							   truncate.date = NULL,
							   mc.choose = mc)
	# Set model parameters:
	PRM <- get.model.prm(dat,
						 dat.full,
						 horiz.forecast ,  
						 GI.mean,GI.stdv,
						 GI.dist="gamma",
						 cori.window=3)
	# Forecast:
	fcast <- try(lapply(PRM,fcast.inc.early.short,do.plot=FALSE),silent = TRUE)
	
	if(class(fcast)!="try-error"){
		print(paste("forecast",mc,"/",n.mc,"done."))
		
		df.tmp <- dist.target(fcast)
		df.tmp$mc <- mc
		
		if (mc==1) df <- df.tmp
		if (mc>1) df <- rbind(df,df.tmp)
	}
}


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

g <- ggplot(df)
g <- g + geom_point(aes(x=log(s),y=b,
						colour=model, shape=model),
					size=3, alpha=0.5)
g <- g + geom_point(data = df.m, 
					aes(x=log(s.m),y=(b.m),
						colour=model, shape=model),
					size=6)

g <- g + geom_segment(data = df.m, 
					  aes(x=log(s.lo),xend=log(s.hi),
					  	y=b.m,yend=b.m,
					  	colour=model, shape=model),
					  size=1)

g <- g + geom_segment(data = df.m, 
					  aes(x=log(s.md),xend=log(s.md),
					  	y=b.lo,yend=b.hi,
					  	colour=model, shape=model),
					  size=1)


plot(g)
