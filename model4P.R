model4P <- function(x, y, W.coef = 0, Plot=FALSE, Title="", pCol = "grey50", lCol = "grey25"){

	x <- as.numeric(x)
	y <- as.numeric(y)

	if(any(is.na(x) | is.na(y))){
		NAs <- which(is.na(x) | is.na(y))
		x <- x[-NAs]
		y <- y[-NAs]
		}

	# Fonction logistique 4P
		nlm.fit.4P <- function(bottom, top, xmid, scal, X){
			Y.fit <-bottom+(top-bottom)/(1+10^((xmid-X)*scal))
			return(Y.fit)
			}

	# Fonction sce (somme carré résidus) avec pondérations
		sce.4P <- function(param, X, yobs, Weights, W.coef) {
			bottom <- param[1]
			top <- param[2]
			xmid <- param[3]
			scal <- param[4]
			ytheo <- nlm.fit.4P(bottom, top, xmid, scal, X)
			residus <- yobs - ytheo
			Weights <- (1/(residus^2))^(W.coef)
			return(sum(Weights*(yobs - ytheo)^2))
			}


	# initialisation des valeurs
		w.coef = W.coef
		bottom.ini = min(y, na.rm=T); min(y)
		top.ini = max(y, na.rm=T); max(y)
		xmid.ini = (max(x)+min(x))/2; xmid.ini

		min.y <- min(y, na.rm=T)		# *0.995
						# if (min(y)<0) min.y <- min(y)*1.005
	
		max.y <- max(y, na.rm=T)		# *1.005
						# if (max(y)<0) max.y <- max(y)*0.995
				
		z <- (y-min.y)/(max.y - min.y)
		z[z==0] <- 0.05
		z[z==1] <- 0.95
		scal.ini = coef(lm(x~log(z/(1-z))))[2]

		scal.ini <- as.numeric(scal.ini)											# ; scal.ini
		e = 1e-5
		# weights <- (z*(1-z))^(W.coef)												# ; weights
		# weights <- weights/sum(weights, na.rm=T)
		weights <- rep(1, length(y))		

		init <- c(bottom=bottom.ini, top=top.ini, xmid=xmid.ini, scal=scal.ini)
		best<-nlm(f = sce.4P, p = init, X = x, yobs = y, Weights=weights, W.coef=w.coef)						# calcul des paramètres

	# Récupération des paramètres
		best.bottom <- best$estimate[1]
		best.top <- best$estimate[2]
		best.xmid<- best$estimate[3]
		best.scal <- best$estimate[4]

	# Estimation des valeurs
		newX <- seq(min(x)*0.75, max(x)*1.2, length=500)						
		Y.best <- nlm.fit.4P(best.bottom, best.top, best.xmid, best.scal, newX)

	# Score de régression
		y.fit <- nlm.fit.4P(best.bottom, best.top, best.xmid, best.scal, x)
		lm.test <- lm(y.fit~y)
		rob.test <- lm.rob(y.fit, y)
		rlm.test <- rlm(y.fit~y)
															# coef(summary(lm.test)); coef(summary(rob.test)); coef(summary(rlm.test))
		t <- coef(summary(rlm.test))[2,3]
		p <- (1 - pt(abs(t), length(y)-2))*2; p

		first <- nlm.fit.4P(best.bottom, best.top, best.xmid, best.scal, x[1])			# ; first
		last <- nlm.fit.4P(best.bottom, best.top, best.xmid, best.scal, x[length(y)])		# ; last

		delta <- last - first											#; delta
		FC <- ifelse(delta>=0, 10^delta, -1/10^delta)							# ; FC

	# Représentations graphiques
		if (Plot){
			leg.pos = "bottomright"
			plot(y~x, pch = 20, cex = 2.5, col = pCol, xlim = range(min(newX), max(newX)), ylim = range(min(y, Y.best),max(y, Y.best)))	# , ylim = range(0, 1))
			lines(Y.best ~ newX, col = lCol, lwd = 3)
			title (main = Title, sub = "Weighted 4P logistic regr.")
			legend1 <- paste("FC =", round(FC,3))
			legend2 <- paste("p-fitting =", signif(p,3))
			if(FC<0) leg.pos = "topright"
			legend(leg.pos, legend = c(legend1, legend2), bty="n")
			}

	# return(list(signif(first, 3), signif(last, 3), round(best.scal, 3), signif(FC, 3), signif(p, 4)))
	return(list(model = best, fitted.values = y.fit))
}

