library(shiny)

source("read-data.R")
source("forecast_early_short.R")
source("forecast_utils.R")

shinyServer(function(input, output) {
	
	output$forecast.plot <- renderPlot(
		{
			# Read all inputs from UI:
			data.type <- input$data.type
			mc <- input$synth.dataset
			
			trunc <- input$trunc
			first.date <- input$first.date
			
			models <- input$model
			horiz.forecast <- input$fcast.horiz
			GI.dist <- input$GI.dist
			GI.mean <- input$GI.mean
			GI.stdv <- input$GI.stdv
			cori.window <- input$cori.window
			dolog <- input$dolog
			
			# Load data:
			if(data.type=="synthetic"){
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
				
				dat <- dat[first.date:nrow(dat),]
				dat.full <- dat.full[first.date:nrow(dat.full),]
				
			}
			
			if(data.type=="real"){
				d0 <- read.csv(file = "./data/2009_Flu_Mexico.csv",header = F)
				names(d0) <- c("t","inc")
				dat <- d0[first.date:trunc,]
				dat.full <- d0[first.date:nrow(d0),]
			}
			# Set-up model parameters
			# (wether used or not)
			PRM <- get.model.prm(dat,
								 dat.full,
								 horiz.forecast ,  
								 GI.mean,GI.stdv,
								 GI.dist,
								 cori.window = cori.window)
			
			PRM.chosen <- PRM[models]
			
			### Forecast
			fcast <- lapply(PRM.chosen,
							fcast.inc.early.short,
							do.plot= FALSE)
			# plot forecast
			compare.fcast.early.2(fcast, dolog = dolog)
		},
		height=700, 
		width=700)
	
	
})

