library(R0)

# loads example dataset
data(Germany.1918)
dat.full <- Germany.1918

# truncate full data:
dat <- dat.full[1:30]

# plot
par(mfrow=c(1,2))
plot(dat.full,typ="s")
lines(dat,typ="s",lwd=6)
plot(log(dat.full),typ="o")
points(log(dat),pch=16,col="black",cex=1.3)

# create generation time : gamma distribution
# with mean 3.5 time units and standard deviation 1 time unit
GT <- generation.time("gamma", c(2.5,1))


# applies methods EG, ML, SB, TD to the dataset
method.R.est <- c("EG","ML","SB","TD")
res.R <- estimate.R(dat, 
					GT=GT, 
					methods=method.R.est,
					force.prior=5)


# Display estimates from all methods

par(mfrow=c(2,2))

for(m in method.R.est){
	print(m)
	x <- res.R$estimates[[m]]
	
	ci.lo <- x$conf.int[1]
	if(class(x$conf.int)=="data.frame") ci.lo <- x$conf.int[,1]
	ci.hi <- x$conf.int[2]
	if(class(x$conf.int)=="data.frame") ci.hi <- x$conf.int[,2]
	
	plot(x$R,typ="o", 
		 pch = 16,
		 main=m, ylab="R estimate")	
	grid()
	lines(ci.lo)
	lines(ci.hi)
	
}



