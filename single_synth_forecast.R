###
###   COMPARE VARIOUS FORECASTING MODEL
###   ON A SINGLE SYNTHETIC DATA SET
###


source("read-data.R")
source("forecast_early_short.R")
source("forecast_utils.R")

save.to.file <- FALSE
do.all.plot <- TRUE 

### Read incidence data:

# Data set number
# (several have been pre-simulated):
mc <- 20

# Synthetic data beyond this date
# are assumed unknown:
trunc <- 19

# How far beyond the last known date
# should the forecast be performed:
horiz.forecast <- 7

# Load pre-simulated synthetic data:
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

# GI values inputed in models:
GI.mean <- 2.3
GI.stdv<- 1

### Model choice and associated parameters:
# Models are:
# WalLip  WhiPag  SeqBay 
# CoriParam CoriNonParam CoriUncertain
PRM <- get.model.prm(dat,
					 dat.full,
					 horiz.forecast ,  
					 GI.mean,GI.stdv,
					 GI.dist,
					 cori.window = 3)

### Forecast

if (save.to.file) pdf("forecast_early.pdf", width=15,height=10)

fcast <- lapply(PRM,
				fcast.inc.early.short,
				do.plot=do.all.plot)

par(mfrow=c(1,1))
compare.fcast.early(fcast)

if (save.to.file) dev.off()

### Note:
### Method "SeqBay" ("SB" in R0 pckg) needs a fix when incidence is zero



