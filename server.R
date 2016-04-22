library(shiny)
source("EstimationR.R")

source("read-data.R")
source("forecast_early_short.R")
source("forecast_utils.R")

full.filename <- function(filename,datatype){
	if(datatype=="synthetic") filetype <- ".RData"
	if(datatype=="real") filetype <- ".csv"
	return(paste0("./data/SEmInR_",filename,filetype))
}

shinyServer(function(input, output) {
	
	output$forecast.plot <- renderPlot(
		{
			# Read all inputs from UI:
			data.type <- input$data.type
			data.file <- input$data.file
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
				
				load(full.filename(data.file,data.type))
				x <- read.incidence.obj(inc.tb = inc.tb,
									 type = "simulated",
									 find.epi.start.window = horiz.forecast + 3,
									 find.epi.start.thresrate = 0.5,
									 truncate.date = trunc,
									 truncate.generation = NULL,
									 mc.choose = mc)
				
				if(is.na(x)) return(NULL)
				
				dat <- x[["dat"]]
				dat.full <- x[["dat.full"]]
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
			PRM <- create.model.prm(dat,
								 dat.full,
								 horiz.forecast ,  
								 GI.mean,GI.stdv,
								 GI.dist,
								 cori.window = cori.window)
			
			PRM.chosen <- PRM[models]
			
			### Forecast
			fcast <- lapply(PRM.chosen,
							fcast_incidence,
							do.plot= FALSE)
			# plot forecast
			compare.fcast.early.shiny(fcast, dolog = dolog)
		},
		height=700, 
		width=700)
	
	
})

