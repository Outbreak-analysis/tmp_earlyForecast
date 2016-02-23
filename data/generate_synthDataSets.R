#### 
####    GENERATE MULTIPLE SYNTHETIC DATA SETS
####    USING A SEmInR STOCHASTIC MODEL
####

library(snowfall)
library(parallel)
library(ggplot2);theme_set(theme_bw())
source("SEmInR_Gillespie_FCT.R")

pdf.options(width=12)

wrap.sim <- function(prm,prmfxd) {
	
	# unpack fixed parameters:
	horizon.years <- prmfxd[["horizon.years"]]
	pop.size <- prmfxd[["pop.size"]]
	I.init <- prmfxd[["I.init"]]
	n.MC <- prmfxd[["n.MC"]]
	remove.fizzles <- prmfxd[["remove.fizzles"]]
	
	# unpack variable parameters:
	DOL.days <- prm[["DOL.days"]]
	DOI.days <- prm[["DOI.days"]]
	R0 <- prm[["R0"]]
	nE <- prm[["nE"]]
	nI <- prm[["nI"]]
	
	sim <- simul.SEmInR(horizon.years=horizon.years ,
						DOL.days=DOL.days,
						DOI.days=DOI.days,
						R0=R0 ,
						pop.size=pop.size,
						nE=nE,
						nI=nI,
						I.init=I.init,
						n.MC=n.MC,
						remove.fizzles=remove.fizzles,
						save.to.Rdata.file = TRUE)
	return(sim)
}

prmfxd <- list(horizon.years = 1.3,
			   pop.size = 1E4,
			   I.init = 2,
			   n.MC = 100,
			   remove.fizzles = TRUE)

# Define the various model parameters (data sets):
prm <- list()

prm[[1]] <- list(DOL.days = 2,
				 DOI.days = 2,
				 R0 = 3.5,
				 nE = 3,
				 nI = 3)

prm[[2]] <- list(DOL.days = 7,
				 DOI.days = 5,
				 R0 = 3.5,
				 nE = 3,
				 nI = 3)

prm[[3]] <- list(DOL.days = 2,
				 DOI.days = 2,
				 R0 = 1.6,
				 nE = 3,
				 nI = 3)

prm[[4]] <- list(DOL.days = 7,
				 DOI.days = 5,
				 R0 = 1.6,
				 nE = 3,
				 nI = 3)

# Run all data sets 
ncpus <- detectCores()
sfInit(parallel = TRUE, cpu = ncpus)
sfLibrary(adaptivetau)
sfLibrary(plyr)
sfExportAll()
SIM <- sfSapply(prm, wrap.sim, prmfxd=prmfxd, simplify = FALSE)
sfStop()

# Plots some simulations:
for(i in 1:length(prm)){
	mc.chosen <- c(1,2,3)
	title <- paste(names(SIM[[i]][["param"]]),SIM[[i]][["param"]],sep="=",collapse = ";")
	title <- paste("SET",i,":",title)
	df <- subset(SIM[[i]][["inc"]], mc %in% mc.chosen)
	g <- ggplot(df) + geom_step(aes(x=tb,y=inc,colour=factor(mc)),size=1)
	g <- g + ggtitle(title)
	plot(g)
}





