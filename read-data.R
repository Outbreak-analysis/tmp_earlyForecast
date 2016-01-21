
read.incidence <- function(filename, # RData file
						   objname, # object storing incidence
						   type, # simulated or real
						   truncate.date = NULL,
						   mc.choose = 1
){
	load(filename)
	tmp <- get(objname)
	if(type=="simulated") tmp <- subset(get(objname),mc==mc.choose)
	
	dat.full <- data.frame(t=tmp$tb, inc=tmp$inc)
	dat <- dat.full
	# Truncate
	if(!is.null(truncate.date)) dat <- dat.full[1:truncate.date,]	
	return(dat)
}
	
# 	# Read incidence
# 	load(paste0(data.dir,"SEmInR_sim.Rdata"))
# mc.choose <- 4   # choose the iteration
# tmp <- subset(inc.tb,mc==mc.choose)
# dat.full <- data.frame(t=tmp$tb, inc=tmp$inc)
# 
# # Truncate
# dat <- dat.full[1:9,]