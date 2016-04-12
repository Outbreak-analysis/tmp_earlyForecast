xx <- RESuDe.generate.data(pop_size = 1e5, 
						   I.init = 2,
						   R0 = 6, 
						   alpha = 0, 
						   kappa = 0.00, 
						   GI_span = 20, 
						   GI_mean = 3, 
						   GI_var = 1.5,
						   horizon = 200,
						   seed=1234)

par(mfrow=c(1,2))
plot(xx$time, xx$I, typ='s',log='y')
plot(xx$time, xx$S, typ='s',log='')

xx$I
xx$S
