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




### - - - - - - - - - - - - - - -
### --- Run the backtesting ---
### - - - - - - - - - - - - - - -

# Models used to generate synthetic data:
syn.models <- list("SEmInR", "RESuDe")
# Identify the source names of  synthetic data
db.path <- "../Datsid/bcktest.db"
bcktest <- get.list.sources(db.path = db.path)
idx <- lapply(syn.models, grepl,x = bcktest )
idx <- rowSums(matrix(unlist(idx),ncol=length(syn.models)))
bcktest <- bcktest[as.logical(idx)]

# Backtest every data sets:
file.prm.bcktst <- 'prm_multi_bcktest.csv'
fpb <- read.csv(file.prm.bcktst,header = F)
horiz.forecast <- as.numeric(as.character(fpb[as.character(fpb[,1])=="horiz.forecast",2]))
GI.bias <- as.numeric(as.character(fpb[as.character(fpb[,1])=="GI.bias",2]))
n.bcktest <- length(bcktest)

sc.tmp <- list()

for(i in 1:n.bcktest){
	# Retrieve all synthetic epidemics from a model parameter set:
	dat.all <- get.synthetic.data.db(db.path     = db.path,
									 source.keys = bcktest[i],
									 eventtype   = 'incidence')
	source <- dat.all$source[1]
	# The 'true' parameters that generated these epidemics:
	trueparam <- get.synthetic.epi.param(source = source)
	GI.mean <- get.GI(trueparam)['GI.mean']
	GI.stdv <- sqrt(get.GI(trueparam)['GI.var'])
	
	print(paste(i,"/",n.bcktest,":",source))
	
	# Find time (after start date) 
	# to truncate full data (to make forecasts):
	ttrunc <- ceiling(get.trunc.time(file = file.prm.bcktst,
									 trueparam = trueparam))
	
	mcvec <- unique(dat.all$mc)
	x <- list()
	dat.chopped <- list()
	for (m in mcvec[1:2]){
		print(paste('m =',m))
		# Extract one realization:
		dat.full <- subset(dat.all, mc == m)
		# Find practical start of epidemic:
		tstart <- find.epi.start(dat.full$inc, w=4, thres.rate = 0.5)
		# shift times such that t=tstart --> t=1
		dat.no.trunc <- subset(dat.full, tb >= tstart)
		dat.no.trunc$tb <- dat.no.trunc$tb - tstart + 1
		# Truncates:
		dat.chopped[[m]] <- subset(dat.no.trunc, tb <= ttrunc)
		
		par(mfrow=c(1,1))
		try(plot(dat.chopped[[m]]$tb,dat.chopped[[m]]$inc, typ='s', main=source),
			silent = TRUE)
		
		# Forecast scores:
		x[[m]] <- fcast.wrap.new(mc = m, 
								 dat.full = dat.no.trunc, 
								 dat = dat.chopped[[m]],
								 horiz.forecast = horiz.forecast,
								 GI.mean = as.numeric(GI.mean), 
								 GI.stdv = as.numeric(GI.stdv),
								 GI.dist = "gamma",
								 cori.window = 3,
								 rel.err = T,
								 do.plot = T)
	}
	sc.tmp[[i]] <- do.call('rbind',x)
	sc.tmp[[i]]$source <- source
}
sc <- do.call('rbind',sc.tmp)
scsum <- ddply(sc,c('model','source'),summarize, ME.med=median(ME),MAE.med = median(MAE))

pdf(file='scores-summary.pdf', width=15,height =15)
g <- ggplot(scsum)+geom_point(aes(x=MAE.med,y=ME.med,
								  shape = model, colour=model))+facet_wrap(~source)
g <- g + scale_x_log10()
plot(g)
dev.off()


# ==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
t2 <- as.numeric(Sys.time())
message(paste("Finished in",round((t2-t1)/60,2),"min"))
