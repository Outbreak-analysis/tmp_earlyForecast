library(EpiEstim)


## load data on pandemic flu in a school in 2009
data("Flu2009")
plot(Flu2009$Incidence,type="s")

N <- length(Flu2009$Incidence)
first.date <- 5
window.size <- 7

t.start <- first.date:(N-window.size+1)
t.end <- (first.date+window.size-1):N

# N2 <- 15
# t.start <- (N2-window.size+1)
# t.end <- N2


## estimate the instantaneous reproduction number (method "NonParametricSI")
# the second plot produced shows, at each each day, 
# the estimate of the instantaneous reproduction number over the 7-day window finishing on that day.
x <- EstimateR(Flu2009$Incidence, 
			   T.Start=t.start, 
			   T.End=t.end, 
			   method= "NonParametricSI", 
			   SI.Distr=Flu2009$SI.Distr, 
			   plot=TRUE, 
			   leg.pos=xy.coords(1,3))

Rbar <- x$R$`Mean(R)`
x$R$`Quantile.0.025(R)`
x$R$`Quantile.0.975(R)`

data.frame(x$R$T.End,Rbar)

x$SIDistr

x2 <- EstimateR(Flu2009$Incidence, 
				T.Start=t.start, 
				T.End=t.end, 
				method= "ParametricSI", 
				Mean.SI=4.6, Std.SI=1.5, 
				plot=TRUE)

x2$SIDistr

x3 <- EstimateR(Flu2009$Incidence, 
				T.Start=t.start, 
				T.End=t.end, 
				method= "UncertainSI", 
				Mean.SI=2.6, Std.Mean.SI=1, 
				Min.Mean.SI=1, Max.Mean.SI=4.2, 
				Std.SI=1.5, Std.Std.SI=0.5, 
				Min.Std.SI=0.5, Max.Std.SI=2.5, 
				n1=100, n2=100, 
				plot=TRUE)

x3$SIDistr


