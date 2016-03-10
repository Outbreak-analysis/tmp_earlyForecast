library(shiny)



model.names <- c("Cori","WalLip","WhiPag","SeqBay","IDEA")


get.file.list <- function(filetype){
	cmd <- paste0("ls ./data/*.",filetype)
	flist <- system(command = cmd, intern = TRUE)
	return(flist)
}

trim.filename <- function(x){
	y <- gsub("./data/","",x)
	y <- gsub(".RData","",y)
	y <- gsub(".csv","",y)
	y <- gsub("SEmInR_","",y)
	return(y)
}



# Define UI for application 
shinyUI(fluidPage(
	
	# Application title
	titlePanel("Epidemic Short-Term Early Forecast"),
	
	# Sidebar with a slider input for the number of bins
	sidebarLayout(
		sidebarPanel(width = 3,
					 
					 selectInput("data.type",
					 			"Data type: ",
					 			choices = c("synthetic","real"),
					 			selected= "synthetic")
					 ,
					 selectInput("data.file",
					 			"Data set: ",
					 			choices = trim.filename(get.file.list("RData")) )
					 ,
					 numericInput("synth.dataset",
					 			 "Synthetic data set number :",
					 			 min = 1,
					 			 step = 1,
					 			 value = 1)
					 ,
					 numericInput("first.date",
					 			 "Data truncation first date :",
					 			 min = 1,
					 			 step = 1,
					 			 value = 1)
					 ,
					 numericInput("trunc",
					 			 "Data truncation last date :",
					 			 min = 1,
					 			 step = 1,
					 			 value = 10)
					 ,
					 numericInput("fcast.horiz",
					 			 "Forecast horizon :",
					 			 min = 1,
					 			 step = 1,
					 			 value = 7) 
					 ,
					 checkboxInput("dolog",
					 			  "Log scale")
					 ,
					 checkboxGroupInput("model", 
					 				   "Model:", 
					 				   choices = model.names,
					 				   selected = "Cori")
					 ,
					 numericInput("cori.window",
					 			 "Cori window size:",
					 			 min = 1,
					 			 step = 1,
					 			 value = 3)
					 ,
					 
					 selectInput("GI.dist",
					 			"GI distribution :",
					 			choices = c("gamma","lognormal"),
					 			selected = "gamma") 
					 ,
					 numericInput("GI.mean",
					 			 "GI mean :",
					 			 min = 1,
					 			 step = 0.1,
					 			 value = 3.5)
					 ,
					 numericInput("GI.stdv",
					 			 "GI std dev:",
					 			 min = 1,
					 			 step = 0.1,
					 			 value = 1.1)
					 
		) # END SIDEBAR PANEL
		
		,
		
		# MAIN PANEL: 
		# Show a plot of forecasts
		mainPanel(
			plotOutput("forecast.plot")
		)
	)
))