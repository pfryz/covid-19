library(shiny)
library(rvest)
library(tidyverse)
library(forecast)
library(not)



wbs.sdll.cpt <- function(x, sigma = stats::mad(diff(x)/sqrt(2)), universal = TRUE, M = NULL, th.const = NULL, th.const.min.mult = 0.3, lambda = 0.9) {

	

	n <- length(x)

    if (n <= 1) {

        no.of.cpt <- 0

        cpt <- integer(0)

    }

    else {

		if (sigma == 0) stop("Noise level estimated at zero; therefore no change-points to estimate.")

		if (universal) {

        	u <- universal.M.th.v3(n, lambda)

        	th.const <- u$th.const

        	M <- u$M

    	}

    	else if (is.null(M) || is.null(th.const)) stop("If universal is FALSE, then M and th.const must be specified.")

    	th.const.min <- th.const * th.const.min.mult

    	th <- th.const * sqrt(2 * log(n)) * sigma

    	th.min <- th.const.min * sqrt(2 * log(n)) * sigma



 		rc <- t(wbs.K.int(x, M))

 		if (max(abs(rc[,4])) < th) {

    	    no.of.cpt <- 0

        	cpt <- integer(0)



 		}

		else {

			indices <- which(abs(rc[,4]) > th.min)

			if (length(indices) == 1) {

				cpt <- rc[indices, 3]

				no.of.cpt <- 1

			}

			else {

				rc.sel <- rc[indices,,drop=F]

				ord <- order(abs(rc.sel[,4]), decreasing=T)

				z <- abs(rc.sel[ord,4])

				z.l <- length(z)

				dif <- -diff(log(z))

				dif.ord <- order(dif, decreasing=T)

				j <- 1

				while ((j < z.l) & (z[dif.ord[j]+1] > th)) j <- j+1

				if (j < z.l) no.of.cpt <- dif.ord[j] else no.of.cpt <- z.l

				cpt <- sort((rc.sel[ord,3])[1:no.of.cpt])			

			}

		} 

    }

    est <- mean.from.cpt(x, cpt)

	list(est=est, no.of.cpt=no.of.cpt, cpt=cpt)

}





wbs.sdll.cpt.rep <- function(x, sigma = stats::mad(diff(x)/sqrt(2)), universal = TRUE, M = NULL, th.const = NULL, th.const.min.mult = 0.3, lambda = 0.9, repeats = 9) {



	res <- vector("list", repeats)

	

	cpt.combined <- integer(0)

	

	nos.of.cpts <- rep(0, repeats)

	

	for (i in 1:repeats) {

		

		res[[i]] <- wbs.sdll.cpt(x, sigma, universal, M, th.const, th.const.min.mult, lambda)

		cpt.combined <- c(cpt.combined, res[[i]]$cpt)

		nos.of.cpts[i] <- res[[i]]$no.of.cpt				

		

	}



	med.no.of.cpt <- median(nos.of.cpts)

	

	med.index <- which.min(abs(nos.of.cpts - med.no.of.cpt))

	

	med.run <- res[[med.index]]

	

	list(med.run = med.run, cpt.combined = sort(cpt.combined))



}





universal.M.th.v3 <- function(n, lambda = 0.9) {

		

	mat.90 <- matrix(0, 24, 3)

	mat.90[,1] <- c(10, 50, 100, 150, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1500, 2000, 2500, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000)

	mat.90[,2] <- c(1.420, 1.310, 1.280, 1.270, 1.250, 1.220, 1.205, 1.205, 1.200, 1.200, 1.200, 1.185, 1.185, 1.170, 1.170, 1.160, 1.150, 1.150, 1.150, 1.150, 1.145, 1.145, 1.135, 1.135)

	mat.90[,3] <- rep(100, 24)

	

	mat.95 <- matrix(0, 24, 3)

	mat.95[,1] <- c(10, 50, 100, 150, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1500, 2000, 2500, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000)

	mat.95[,2] <- c(1.550, 1.370, 1.340, 1.320, 1.300, 1.290, 1.265, 1.265, 1.247, 1.247, 1.247, 1.225, 1.225, 1.220, 1.210, 1.190, 1.190, 1.190, 1.190, 1.190, 1.190, 1.180, 1.170, 1.170)

	mat.95[,3] <- rep(100, 24)



	if (lambda == 0.9) A <- mat.90 else A <- mat.95



	d <- dim(A)

	if (n < A[1,1]) {

		th <- A[1,2]

		M <- A[1,3]

	}

	else if (n > A[d[1],1]) {

		th <- A[d[1],2]

		M <- A[d[1],3]

	}

	else {

		ind <- order(abs(n - A[,1]))[1:2]

		s <- min(ind)

		e <- max(ind)

		th <- A[s,2] * (A[e,1] - n)/(A[e,1] - A[s,1]) + A[e,2] * (n - A[s,1])/(A[e,1] - A[s,1])

		M <- A[s,3] * (A[e,1] - n)/(A[e,1] - A[s,1]) + A[e,3] * (n - A[s,1])/(A[e,1] - A[s,1])

	}



	list(th.const=th, M=M)

}





wbs.K.int <- function(x, M) {

	

	n <- length(x)

	if (n == 1) return(matrix(NA, 4, 0))

	else {

		cpt <- t(random.cusums(x, M)$max.val)

		return(cbind(cpt, wbs.K.int(x[1:cpt[3]], M), wbs.K.int(x[(cpt[3]+1):n], M) + c(rep(cpt[3], 3), 0)            ))

	}

	

}





mean.from.cpt <- function(x, cpt) {



	n <- length(x)

	len.cpt <- length(cpt)

	if (len.cpt) cpt <- sort(cpt)

	beg <- endd <- rep(0, len.cpt+1)

	beg[1] <- 1

	endd[len.cpt+1] <- n

	if (len.cpt) {

		beg[2:(len.cpt+1)] <- cpt+1

		endd[1:len.cpt] <- cpt

	}

	means <- rep(0, len.cpt+1)

	for (i in 1:(len.cpt+1)) means[i] <- mean(x[beg[i]:endd[i]])

	rep(means, endd-beg+1)

}





random.cusums <- function(x, M) {



	y <- c(0, cumsum(x))



	n <- length(x)

	

	M <- min(M, (n-1)*n/2)

		

	res <- matrix(0, M, 4)

	

	if (n==2) ind <- matrix(c(1, 2), 2, 1)

	else if (M == (n-1)*n/2) {

		ind <- matrix(0, 2, M)

		ind[1,] <- rep(1:(n-1), (n-1):1)

		ind[2,] <- 2:(M+1) - rep(cumsum(c(0, (n-2):1)), (n-1):1)

	}

	else {

		ind <- ind2 <- matrix(floor(runif(2*M) * (n-1)), nrow=2)

		ind2[1,] <- apply(ind, 2, min)

		ind2[2,] <- apply(ind, 2, max)

		ind <- ind2 + c(1, 2)

	}



	res[,1:2] <- t(ind)

	res[,3:4] <- t(apply(ind, 2, max.cusum, y))



	max.ind <- which.max(abs(res[,4]))



	max.val <- res[max.ind,,drop=F]



	list(res=res, max.val=max.val, M.eff=M)



}





max.cusum <- function(ind, y) {

	

		z <- y[(ind[1]+1):(ind[2]+1)] - y[ind[1]]

		m <- ind[2]-ind[1]+1

		ip <- sqrt(((m-1):1) / m / (1:(m-1))) * z[1:(m-1)] - sqrt((1:(m-1)) / m / ((m-1):1)) * (z[m] - z[1:(m-1)])

		ip.max <- which.max(abs(ip))

		

		c(ip.max + ind[1] - 1, ip[ip.max])



}



clipped <- function(x, minn, maxx) {
	
	min(max(x, minn), maxx)
		
}

read_data_wiki <- function() {
	
	cv.page <- "https://en.wikipedia.org/wiki/2020_coronavirus_pandemic_in_the_United_Kingdom"

	i <- 1
	
	repeat {
		
		xp <- paste('//*[@id="mw-content-text"]/div/table[', as.character(i), ']', sep="")
		read_html(cv.page) %>% html_node(xpath=xp) %>% html_table(fill=TRUE) -> dd
		if (dim(dd)[2] >= 19) break
		i <- i+1
		
	}

	gsub(",", "", dd[[13]]) -> cases_str
	n <- length(cases_str)
	cases_int <- as.numeric(cases_str[2:(n-2)])

	gsub(",", "", dd[[18]]) -> tested_str
	tested_int <- c(rep(0, 6), diff(as.numeric(tested_str[7:(n-2)])))


	tested_actual <- tested_int[7:(n-3)]
	cases_actual <- cases_int[7:(n-3)]
	
	list(tested_actual=tested_actual, cases_actual=cases_actual)
	
}


read_data_emma <- function() {
	
	f <- read_csv("https://raw.githubusercontent.com/emmadoughty/Daily_COVID-19/master/Data/COVID19_by_day.csv", col_types = cols())
	cases_int <- f %>% pull(2)
	n <- length(cases_int)
	cases_actual <- cases_int[35:n]
	tested_int <- f %>% pull(4)
	tested_actual <- tested_int[35:n]

	list(tested_actual=tested_actual, cases_actual=cases_actual)
	
}



fcast_cases <- function() {
	
#	uses libraries: rvest, tidyverse, forecast, not
#	also uses: https://github.com/pfryz/wild-binary-segmentation-2.0/blob/master/WBS2_SDLL_github_v3.R

	d_wiki <- read_data_wiki()
	d_emma <- read_data_emma()
	
	if (length(d_emma$tested_actual) > length(d_wiki$tested_actual)) {
		tested_actual <- d_emma$tested_actual
		cases_actual <- d_emma$cases_actual
	} else {
		tested_actual <- d_wiki$tested_actual
		cases_actual <- d_wiki$cases_actual
	}

	n <- length(tested_actual)
	prop_actual <- cases_actual / tested_actual
	
	tested_fit_lc <- round(mean.from.cpt(tested_actual, wbs.sdll.cpt(sqrt(tested_actual))$cpt))
	tested_fcast_lc <- clipped(tested_fit_lc[n], 0, Inf)

	tested_fit_fcast_tv <- forecast(tested_actual, 1)
	tested_fit_tv <- round(as.numeric(tested_fit_fcast_tv$fitted))
	tested_fcast_tv <- clipped(round(as.numeric(tested_fit_fcast_tv$mean)), 0, Inf)
	
	tested_fit_pl <- predict(not(tested_actual, contrast="pcwsLinContMean"))
	tested_fcast_pl <- clipped(round(2 * tested_fit_pl[n] - tested_fit_pl[n-1]), 0, Inf)


	prop_fit_lc <- wbs.sdll.cpt(prop_actual)$est
	prop_fcast_lc <- clipped(prop_fit_lc[n], 0, 1)

	prop_fit_fcast_tv <- forecast(prop_actual, 1)
	prop_fit_tv <- as.numeric(prop_fit_fcast_tv$fitted)
	prop_fcast_tv <- clipped(as.numeric(prop_fit_fcast_tv$mean), 0, 1)

	prop_fit_pl <- predict(not(prop_actual, contrast="pcwsLinContMean"))
	prop_fcast_pl <- clipped(2 * prop_fit_pl[n] - prop_fit_pl[n-1], 0, 1)
	

	cases_fcast_lc <- round(prop_fcast_lc * tested_fcast_lc)
	cases_fcast_tv <- round(prop_fcast_tv * tested_fcast_tv)
	cases_fcast_pl <- round(prop_fcast_pl * tested_fcast_pl)
	
	list(tested_actual=tested_actual, tested_fit_lc=tested_fit_lc, tested_fcast_lc=tested_fcast_lc, tested_fit_tv=tested_fit_tv, tested_fcast_tv=tested_fcast_tv,
	 tested_fit_pl=tested_fit_pl, tested_fcast_pl=tested_fcast_pl, prop_actual=prop_actual, prop_fit_lc=prop_fit_lc, prop_fcast_lc=prop_fcast_lc,
	 prop_fit_tv=prop_fit_tv, prop_fcast_tv=prop_fcast_tv, prop_fit_pl=prop_fit_pl, prop_fcast_pl=prop_fcast_pl, cases_fcast_lc=cases_fcast_lc, cases_fcast_tv=cases_fcast_tv, cases_fcast_pl=cases_fcast_pl)	
	
}



ui <- fluidPage(

  titlePanel("Trends and next day forecasts for the number of detected Covid-19 cases in the UK"),

  sidebarLayout(

    sidebarPanel(

radioButtons("radio", h3("Trend estimates (see the bottom of the page for references)"),
                        choices = list("piecewise linear" = 1, "smoothly varying" = 2, "piecewise constant" = 3),
                                       ,selected = 1)
                                     




    ),

    mainPanel(

	h3(textOutput("f_tests")),
	h3(textOutput("f_prop")),
	h3(textOutput("f_cases")),
      plotOutput(outputId = "ts_plot"),
      plotOutput(outputId = "ts_plot_1"),
      			h4("black: actual figures", align="center", style = "color:black"),
			h4("brown: statistical trend estimates", align="center", style = "color:brown"),
			h6("References:"),
			h6("[piecewise linear trend] R. Baranowski, Y. Chen and P. Fryzlewicz (2019), ", tags$a(href="https://rss.onlinelibrary.wiley.com/doi/full/10.1111/rssb.12322", "Narrowest‐over‐threshold detection of multiple change points and change‐point‐like features"), ", JRSSB, 81, 649-672"),
			h6("[smoothly varying trend] R package ", tags$a(href="https://CRAN.R-project.org/package=forecast", "forecast")),
			h6("[piecewise constant trend] P. Fryzlewicz (2020), ", tags$a(href="https://link.springer.com/article/10.1007/s42952-020-00060-x", "Detecting possibly frequent change-points: Wild Binary Segmentation 2 and steepest-drop model selection"), ", JKSS, to appear"),
			h6("[data sources]", tags$a(href="https://en.wikipedia.org/wiki/2020_coronavirus_pandemic_in_the_United_Kingdom", "https://en.wikipedia.org/wiki/2020_coronavirus_pandemic_in_the_United_Kingdom"), "and", tags$a(href="https://github.com/emmadoughty/Daily_COVID-19/blob/master/Data/COVID19_by_day.csv", "https://github.com/emmadoughty/Daily_COVID-19/blob/master/Data/COVID19_by_day.csv")),
			h6("[this app]", tags$a(href="https://github.com/pfryz/covid-19", "https://github.com/pfryz/covid-19")),
			h6("[author]", tags$a(href="http://stats.lse.ac.uk/fryzlewicz/", "Piotr Fryzlewicz"))



    )
  )
)


server <- function(input, output) {

	dd <- fcast_cases()

	
	output$f_tests <- renderText({
		
		if (input$radio == 1) pred_tests <- dd$tested_fcast_pl else if (input$radio == 2) pred_tests <- dd$tested_fcast_tv else pred_tests <- dd$tested_fcast_lc
		
		paste("Next day's predicted no. of tests:", pred_tests)
		
		
	})

	output$f_prop <- renderText({

		if (input$radio == 1) pred_prop <- dd$prop_fcast_pl else if (input$radio == 2) pred_prop <- dd$prop_fcast_tv else pred_prop <- dd$prop_fcast_lc
				
		paste("Next day's predicted proportion of +ve test results:", round(pred_prop, 2))
		
		
	})
	
	
	output$f_cases <- renderText({

		if (input$radio == 1) pred_cases <- dd$cases_fcast_pl else if (input$radio == 2) pred_cases <- dd$cases_fcast_tv else pred_cases <- dd$cases_fcast_lc
				
		paste("Next day's predicted no. of detected cases:", pred_cases)
		
		
	})
	
	
	output$ts_plot <- renderPlot({

    		ts.plot(dd$tested_actual, main="Daily number of tests performed, starting from 28/02/2020", ylab="", xlab="Day number")
   if (input$radio == 1) 		lines(dd$tested_fit_pl, col="brown", lwd=2)
   if (input$radio == 2)	lines(dd$tested_fit_tv, col="brown", lwd=2)
   if (input$radio == 3)	lines(dd$tested_fit_lc, col="brown", lwd=2)


    })


	output$ts_plot_1 <- renderPlot({

  		ts.plot(dd$prop_actual, main="Daily number of confirmed cases as a proportion of daily tests, starting from 28/02/2020", ylab="", xlab="Day number")

   if (input$radio == 1) 		lines(dd$prop_fit_pl, col="brown", lwd=2)
   if (input$radio == 2)	lines(dd$prop_fit_tv, col="brown", lwd=2)
   if (input$radio == 3)	lines(dd$prop_fit_lc, col="brown", lwd=2)



    })

	
	
}






shinyApp(ui = ui, server = server)
