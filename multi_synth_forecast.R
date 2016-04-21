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
n.bcktest <- length(bcktest)

# Backtest every data sets:
file.prm.bcktst <- 'prm_multi_bcktest.csv'
fpb <- read.csv(file.prm.bcktst,header = F)
read_prm <- function(x){
	as.numeric(as.character(fpb[as.character(fpb[,1])==x,2]))
}
horiz.forecast <- read_prm('horiz.forecast') 
GI.bias        <- read_prm('GI.bias') 
n.MC.max       <- read_prm('n.MC.max')
CI             <- read_prm('CI') 
multicores     <- read_prm('parallel')

sc.tmp <- list()

for(i in 1:n.bcktest){
	# Retrieve all synthetic epidemics from a model parameter set:
	dat.all <- get.synthetic.data.db(db.path     = db.path,
									 source.keys = bcktest[i],
									 eventtype   = 'incidence')
	mcvec <- unique(dat.all$mc)
	source <- dat.all$source[1]
	print(paste(i,"/",n.bcktest,":",source))
	# The 'true' parameters that generated these epidemics:
	trueparam  <- get.synthetic.epi.param(source = source)
	GI.mean    <- get.GI(trueparam)['GI.mean']
	GI.stdv    <- sqrt(get.GI(trueparam)['GI.var'])
	
	# Find time (after start date) 
	# to truncate full data (to make forecasts):
	ttrunc     <- ceiling(get.trunc.time(file = file.prm.bcktst,
									 trueparam = trueparam))
	# Parallel execution for a given scenario
	# across all MC realizations:
	n.cores <- detectCores()
	if(!multicores) n.cores <- 1
	sfInit(parallel = (n.cores>1), 
		   cpu = n.cores)
	sfLibrary(R0)
	if(!use.DC.version.of.EpiEstim) sfLibrary(EpiEstim)
	
	idx.apply <- mcvec
	message(paste("Synthetic data contains",length(idx.apply),"MC iterations"))
	# Reduce backtesting to 
	# specified MC realizations:
	if(n.MC.max>0) {
		idx.apply <- idx.apply[1:n.MC.max]
		message(paste("but not more than",length(idx.apply),"are used."))
	}
	sfExportAll()
	res.parallel <- sfSapply(idx.apply, 
							 simplify = FALSE,
							 fcast.wrap.mc,
							 dat.all  = dat.all,
							 ttrunc   = ttrunc,
							 horiz.forecast = horiz.forecast,
							 GI.mean  = GI.mean,
							 GI.stdv  = GI.stdv,
							 do.plot = (n.cores==1)
	)
	sfStop()
	# Store each scenario in a list
	sc.tmp[[i]] <- do.call('rbind',res.parallel)
	sc.tmp[[i]]$modelsyndata <- substr(x = source, start = 1, stop=6)
	sc.tmp[[i]]$source <- source
	sc.tmp[[i]]$R0 <- trueparam['R0']
	sc.tmp[[i]]$GI.mean <- GI.mean
}
# Merge everything:
sc <- do.call('rbind',sc.tmp)

# Summary statistics:
scsum <- ddply(sc,c('model','source','modelsyndata','R0','GI.mean'),summarize, 
			   ME.med    = median(ME),
			   ME.lo     = quantile(ME,probs = 0.5+CI/2),
			   ME.hi     = quantile(ME,probs = 0.5-CI/2),
			   MAE.med   = median(MAE),
			   MAE.lo    = quantile(MAE,probs = 0.5+CI/2),
			   MAE.hi    = quantile(MAE,probs = 0.5-CI/2))

# Plot:
plot.scores(scsum)


# ==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
t2 <- as.numeric(Sys.time())
message(paste("Finished in",round((t2-t1)/60,2),"min"))
