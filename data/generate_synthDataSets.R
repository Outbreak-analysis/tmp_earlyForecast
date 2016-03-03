#### 
####    GENERATE MULTIPLE SYNTHETIC DATA SETS
####    USING A SEmInR STOCHASTIC MODEL
####

library(snowfall)
library(parallel)
library(ggplot2);theme_set(theme_bw())
source("SEmInR_Gillespie_FCT.R")

args <- commandArgs(trailingOnly = TRUE)
n.MC <- as.numeric(args[1])

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
			   n.MC = n.MC,
			   remove.fizzles = TRUE)

# Define the various model parameters (data sets):
prm <- list()


Dvec <- c(2, 4, 8)
R0vec <- c(1.5, 3, 6, 12)
nI <- nE <- 5

cnt <- 1
for(d in Dvec){
	for(r in R0vec){
		prm[[cnt]] <- list(DOL.days = d,
						 DOI.days = d,
						 R0 = r,
						 nE = nE,
						 nI = nI)
		cnt <- cnt + 1
	}
}

message(paste("   ===> Simulating",length(prm),"x",prmfxd[["n.MC"]],"parameter sets...\n"))

t1 <- as.numeric(Sys.time())
# Run all data sets 
sfInit(parallel = TRUE, cpu = detectCores())
sfLibrary(adaptivetau)
sfLibrary(plyr)
sfExportAll()
SIM <- sfSapply(prm, wrap.sim, prmfxd=prmfxd, simplify = FALSE)
sfStop()

message("... done.")

# Plots some simulations:
df <- data.frame()
mc.chosen <- 1:4

for(i in 1:length(prm)){
	title <- paste(names(SIM[[i]][["param"]]),SIM[[i]][["param"]],sep="=",collapse = ";")
	title <- paste0("SET ",i,": ",title)
	tmp <- SIM[[i]][["inc"]]
	tmp <- subset(tmp, mc %in% mc.chosen)
	tmp$title <-factor(title)
	tmp$title2 <- i
	df <- rbind(df, tmp)
}
pdf("plot_data.pdf",width=22,height = 15)
g <- ggplot(df) + geom_step(aes(x=tb,y=inc,colour=factor(mc)),size=1)
g <- g + facet_wrap(~title) + scale_y_log10()
plot(g)
dev.off()



# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

t2 <- as.numeric(Sys.time())
message(paste("Completed in",round((t2-t1)/60,2),"minutes"))





