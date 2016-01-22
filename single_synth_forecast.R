

source("read-data.R")
source("forecast_early_short.R")
source("forecast_utils.R")

save.to.file <- FALSE
do.all.plot <- TRUE 

### Read incidence data:

trunc <- 19
mc <- 20

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

### Model choice and associated parameters:
# Models are:
# WalLip  WhiPag  SeqBay 
# CoriParam CoriNonParam CoriUncertain

horiz.forecast <- 12
GI.mean <- 2.3
GI.stdv<- 1

PRM <- list(Cori = list(model = "CoriParam",  
						dat = dat,
						dat.full = dat.full,
						horiz.forecast = horiz.forecast,  
						GI.val = c(GI.mean,GI.stdv),
						cori.window = 3),
			
			WalLip = list(model = "WalLip",
						  dat = dat,
						  dat.full = dat.full,
						  horiz.forecast = horiz.forecast,  
						  GI.dist = "gamma",  # gamma, lognormal, weibull
						  GI.val = c(GI.mean,GI.stdv)),
			
			WhiPag = list(model = "WhiPag",
						  dat = dat,
						  dat.full = dat.full,
						  horiz.forecast = horiz.forecast,  
						  GI.dist = "gamma",  # gamma, lognormal, weibull
						  GI.val = c(GI.mean,GI.stdv)),
			
			SeqBay = list(model = "SeqBay",
						  dat = dat,
						  dat.full = dat.full,
						  horiz.forecast = horiz.forecast,  
						  GI.dist = "gamma",  # gamma, lognormal, weibull
						  GI.val = c(GI.mean,GI.stdv))
)

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



